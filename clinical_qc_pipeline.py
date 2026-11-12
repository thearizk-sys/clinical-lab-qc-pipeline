"""Multi-site clinical laboratory quality-control and streaming ETL pipeline."""
import csv
import json
import math
import time
from collections import Counter
from datetime import datetime
from functools import wraps
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIELDS = ["patient_id", "sample_date", "site", "analyte", "value", "unit"]


def timed(func):
    """Report elapsed wall-clock time for a pipeline consumer."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:
            wrapper.last_duration = time.perf_counter() - start
            print(f"{func.__name__}: {wrapper.last_duration:.3f} seconds")
    wrapper.last_duration = None
    return wrapper


def load_config(path: Path | str) -> dict:
    """Read JSON configuration; missing or malformed files stop execution."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_value(text: str) -> dict:
    """Return an exact numeric value or a named missing, invalid, or censored status."""
    token = text.strip().lower()
    result = {"value": None, "value_status": "missing", "bound": None}
    if token in {"", "na", "n/a", "nd", "?"}:
        return result
    censored = token.startswith("<")
    number = token[1:] if censored else token
    try:
        value = float(number.replace(",", "."))
    except ValueError:
        result["value_status"] = "invalid"
        return result
    if not math.isfinite(value):
        result["value_status"] = "nonfinite"
    elif censored:
        result.update(value_status="censored", bound=value)
    else:
        result.update(value_status="numeric", value=value)
    return result


def iso_date(text: str) -> str:
    """Parse supported date formats into canonical ISO 8601 (YYYY-MM-DD)."""
    for pattern in ("%Y-%m-%d", "%d/%m/%Y", "%d-%b-%Y"):
        try:
            return datetime.strptime(text.strip(), pattern).date().isoformat()
        except ValueError:
            continue
    raise ValueError(f"Unsupported or invalid date: {text!r}")


# =====================================================================
# 1. Lazy Streaming Generator Stages (Constant-Memory ETL)
# =====================================================================

def stream_rows(path: Path, tally: Counter):
    """Yield structured six-field dictionaries; count and skip malformed CSV lines."""
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header != FIELDS:
            raise ValueError(f"Unexpected CSV header: {header!r}")
        for row in reader:
            tally["rows_read"] += 1
            if len(row) != 6:
                reason = "too_few_fields" if len(row) < 6 else "too_many_fields"
                tally[reason] += 1
                continue
            record = {key: text.strip() for key, text in zip(FIELDS, row)}
            record["patient_id"] = record["patient_id"].upper()
            tally["structured"] += 1
            yield record


def parse_records(records, tally: Counter):
    """Yield copied records with parsed values while retaining raw value text."""
    for record in records:
        out = record.copy()
        out["raw_value"] = record["value"]
        out.update(parse_value(record["value"]))
        tally["value_" + out["value_status"]] += 1
        yield out


def convert_records(records, config: dict, tally: Counter):
    """Yield records converted to canonical units using JSON schema lookup."""
    for record in records:
        out = record.copy()
        name = record["analyte"].strip().lower()
        name = config["aliases"].get(name, name)
        unit = record["unit"].strip().lower()
        out["analyte"] = name
        if name not in config["analytes"]:
            out["unit_status"] = "unknown_analyte"
        elif unit not in config["analytes"][name]["factors"]:
            out["unit_status"] = "unknown_unit"
        else:
            entry = config["analytes"][name]
            factor = entry["factors"][unit]
            if factor <= 0 or not math.isfinite(factor):
                raise ValueError(f"Invalid conversion factor for {name}, {unit}")
            out["unit_status"] = "converted"
            out["unit"] = entry["canonical_unit"]
            if out["value_status"] == "numeric":
                out["value"] *= factor
            elif out["value_status"] == "censored":
                out["bound"] *= factor
        tally[out["unit_status"]] += 1
        yield out


def filter_usable(records, tally: Counter):
    """Keep converted numeric/censored records; exclude impossible concentrations."""
    for record in records:
        if record["unit_status"] != "converted":
            tally["excluded_unit"] += 1
            continue
        status = record["value_status"]
        if status not in {"numeric", "censored"}:
            tally["excluded_" + status] += 1
            continue
        value = record["value"] if status == "numeric" else record["bound"]
        if value < 0 or (status == "censored" and value == 0):
            tally["negative_or_impossible_bound"] += 1
            continue
        yield record


def clean_dates(records, tally: Counter):
    """Yield copied records with ISO dates; count and skip invalid dates."""
    for record in records:
        try:
            date = iso_date(record["sample_date"])
        except ValueError:
            tally["invalid_date"] += 1
            continue
        out = record.copy()
        out["sample_date"] = date
        yield out


def label_ranges(records, ranges: dict, tally: Counter):
    """Classify numeric measurements against clinical reference limits."""
    for record in records:
        out = record.copy()
        if record["analyte"] not in ranges:
            tally["unknown_range"] += 1
            continue
        limits = ranges[record["analyte"]]
        expected = limits["unit"].replace("µ", "u").lower()
        if record["unit"].replace("µ", "u").lower() != expected:
            raise ValueError("Reference-range unit does not match cleaned unit")
        if record["value_status"] == "censored":
            status = "censored"
        elif record["value"] < limits["low"]:
            status = "low"
        elif record["value"] > limits["high"]:
            status = "high"
        else:
            status = "normal"
        out["range_status"] = status
        tally["range_" + status] += 1
        tally["output"] += 1
        yield out


def deduplicate(records, tally: Counter) -> list[dict]:
    """Retain the first usable record per (patient_id, sample_date, analyte)."""
    seen = {}
    results = []
    for record in records:
        key = (record["patient_id"], record["sample_date"], record["analyte"])
        signature = (record["value_status"], record["value"], record["bound"])
        if key in seen:
            old = seen[key]
            same = signature[0] == old[0]
            if same and signature[0] == "numeric":
                same = math.isclose(signature[1], old[1], rel_tol=1e-9)
            elif same:
                same = math.isclose(signature[2], old[2], rel_tol=1e-9)
            reason = "duplicate_same" if same else "duplicate_conflict"
            tally[reason] += 1
            continue
        seen[key] = signature
        results.append(record)
    return results


@timed
def consume_stream(records) -> Counter:
    """Consume a lazy generator stream in constant additional memory."""
    return Counter(record["range_status"] for record in records)


# =====================================================================
# 2. Vectorised Statistical Summary & Multi-Panel QC Visualisation
# =====================================================================

def summarise_and_plot(clean_records: list[dict], ranges: dict, config: dict, out_dir: Path):
    """Generate per-analyte stats, site out-of-range rates, outlier flags, and QC figure."""
    all_clean = pd.DataFrame(clean_records)
    df = all_clean[all_clean["value_status"] == "numeric"].copy()

    groups = df.groupby("analyte")["value"]
    summary = pd.DataFrame({
        "count": groups.count(),
        "mean": groups.mean(),
        "std": groups.std(ddof=1),
        "median": groups.median(),
        "min": groups.min(),
        "max": groups.max(),
    })

    df["out_of_range"] = df["range_status"] != "normal"
    site_groups = df.groupby(["site", "analyte"])
    rates = pd.DataFrame({
        "count": site_groups.size(),
        "out_of_range_rate": site_groups["out_of_range"].mean(),
    })

    # Flag |z| > 3 per analyte for manual review without silent deletion
    df["outlier"] = False
    for analyte in sorted(df["analyte"].unique()):
        mask = df["analyte"] == analyte
        values = df.loc[mask, "value"].to_numpy()
        sd = values.std(ddof=0)
        if sd > 0:
            df.loc[mask, "outlier"] = np.abs((values - values.mean()) / sd) > 3

    out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_dir / "summary.csv")
    rates.to_csv(out_dir / "site_rates.csv")
    df.to_csv(out_dir / "cleaned_exact.csv", index=False)

    analytes = sorted(df["analyte"].unique())
    sites = sorted(df["site"].unique())
    fig, axes = plt.subplots(
        len(analytes), 2, figsize=(12, 3 * len(analytes)), squeeze=False
    )
    for row, analyte in enumerate(analytes):
        subset = df[df["analyte"] == analyte]
        limits = ranges[analyte]
        unit = config["analytes"][analyte]["canonical_unit"]
        left, right = axes[row]

        display_max = 2 * limits["high"]
        visible = subset[subset["value"] <= display_max]
        outside = len(subset) - len(visible)

        left.hist(
            visible["value"],
            bins=np.linspace(0, display_max, 26),
            color="steelblue",
            edgecolor="white",
        )
        left.axvline(limits["low"], color="darkorange", linestyle="--", label="lower limit")
        left.axvline(limits["high"], color="crimson", linestyle="--", label="upper limit")
        left.set_title(f"{analyte}: zoom 0-{display_max:g} (total n={len(subset)})")
        left.set_xlabel(f"{analyte} ({unit})")
        left.set_ylabel("Record count")
        left.legend(fontsize=8)
        left.text(
            0.98, 0.97,
            f"Above window: {outside}; full max: {subset['value'].max():.2f}",
            transform=left.transAxes, ha="right", va="top", fontsize=8,
        )

        counts = []
        percentages = []
        for site in sites:
            site_values = subset[subset["site"] == site]
            counts.append(len(site_values))
            percentages.append(100 * site_values["out_of_range"].mean())

        right.bar(sites, percentages, color="teal")
        right.set_ylim(0, 110)
        right.set_xlabel("Hospital site")
        right.set_ylabel("Outside reference range (%)")
        right.set_title(f"{analyte}: site comparison")
        for index, (percentage, count) in enumerate(zip(percentages, counts)):
            if count:
                right.text(index, percentage + 2, f"n={count}", ha="center")

    fig.suptitle(
        "Clinical laboratory QC: exact measurements only; first usable repeat kept",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "qc.png", dpi=150)
    plt.close(fig)
    return summary, rates


def main():
    """Run the full QC pipeline on raw clinical lab data."""
    data_dir = Path("data")
    config = load_config("units.json")
    ranges = load_config(data_dir / "reference_ranges.json")

    tally = Counter()
    stream = stream_rows(data_dir / "labs_raw.csv", tally)
    stream = parse_records(stream, tally)
    stream = convert_records(stream, config, tally)
    stream = filter_usable(stream, tally)
    stream = clean_dates(stream, tally)
    deduped = deduplicate(stream, tally)
    clean_records = list(label_ranges(deduped, ranges, tally))

    print("QC Provenance Tally:", dict(tally))
    summary, rates = summarise_and_plot(clean_records, ranges, config, Path("outputs"))
    print("\nPer-Analyte Summary:\n", summary.round(3))
    print("\nOut-of-Range Rates:\n", rates.round(3))


if __name__ == "__main__":
    main()
