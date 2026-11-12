"""Clinical reference-range classification, statistical QC summaries, and plotting."""
import math
import time
from collections import Counter
from functools import wraps
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def timed(func):
    """Decorator to report elapsed execution time for a stream consumer."""
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


def deduplicate(records, tally: Counter) -> list[dict]:
    """Keep the first usable record per (patient_id, sample_date, analyte) in input order."""
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


def label_ranges(records, ranges: dict, tally: Counter):
    """Lazily classify numeric measurements as low, normal, or high against clinical limits."""
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


@timed
def consume_ranges(records) -> Counter:
    """Consume a range-labelled generator stream in constant memory."""
    return Counter(record["range_status"] for record in records)


def summarize_clean_data(clean_records: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute vectorised per-analyte stats, site out-of-range rates, and z-score outlier flags."""
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

    df["outlier"] = False
    for analyte in sorted(df["analyte"].unique()):
        mask = df["analyte"] == analyte
        values = df.loc[mask, "value"].to_numpy()
        sd = values.std(ddof=0)
        if sd > 0:
            df.loc[mask, "outlier"] = np.abs((values - values.mean()) / sd) > 3

    return df, summary, rates


def plot_qc(df: pd.DataFrame, ranges: dict, config: dict, output_path: str | Path = "qc.png") -> None:
    """Generate a multi-panel QC figure showing analyte distributions and site out-of-range rates."""
    analytes = sorted(df["analyte"].unique())
    sites = sorted(df["site"].unique())
    fig, axes = plt.subplots(len(analytes), 2, figsize=(12, 3 * len(analytes)), squeeze=False)

    for row, analyte in enumerate(analytes):
        subset = df[df["analyte"] == analyte]
        limits = ranges[analyte]
        unit = config["analytes"][analyte]["canonical_unit"]
        left, right = axes[row]

        display_max = 2 * limits["high"]
        visible = subset[subset["value"] <= display_max]
        outside = len(subset) - len(visible)
        left.hist(visible["value"], bins=np.linspace(0, display_max, 26), color="steelblue", edgecolor="white")
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

        counts = [len(subset[subset["site"] == s]) for s in sites]
        percentages = [100 * subset[subset["site"] == s]["out_of_range"].mean() for s in sites]
        right.bar(sites, percentages, color="teal")
        right.set_ylim(0, 110)
        right.set_xlabel("Hospital site")
        right.set_ylabel("Outside reference range (%)")
        right.set_title(f"{analyte}: site comparison")
        for idx, (pct, cnt) in enumerate(zip(percentages, counts)):
            if cnt:
                right.text(idx, pct + 2, f"n={cnt}", ha="center")

    fig.suptitle("Clinical laboratory QC: exact measurements only; first usable repeat kept", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
