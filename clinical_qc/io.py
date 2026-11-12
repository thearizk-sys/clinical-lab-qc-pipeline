"""Input/output, schema-driven unit conversion, and lazy stream parsing."""
import csv
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path

FIELDS = ["patient_id", "sample_date", "site", "analyte", "value", "unit"]


def load_config(path: str | Path) -> dict:
    """Read JSON configuration; missing or malformed files halt the pipeline."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_value(text: str) -> dict:
    """Parse raw value strings into numeric, censored, missing, or invalid states."""
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
    """Normalise multiple date formats into ISO 8601 (YYYY-MM-DD)."""
    for pattern in ("%Y-%m-%d", "%d/%m/%Y", "%d-%b-%Y"):
        try:
            return datetime.strptime(text.strip(), pattern).date().isoformat()
        except ValueError:
            continue
    raise ValueError(f"Unsupported or invalid date: {text!r}")


def read_rows(path: str | Path, tally: Counter):
    """Lazily yield structured records while tracking malformed CSV lines."""
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        if next(reader, None) != FIELDS:
            raise ValueError("Unexpected CSV header")
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
    """Lazily parse measurement values and record status provenance."""
    for record in records:
        out = record.copy()
        out["raw_value"] = record["value"]
        out.update(parse_value(record["value"]))
        tally["value_" + out["value_status"]] += 1
        yield out


def convert_records(records, config: dict, tally: Counter):
    """Lazily standardise analyte aliases and convert values to canonical units."""
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
    """Filter out unconverted units, non-numeric entries, and negative bounds."""
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
    """Lazily convert sample dates to ISO format and drop invalid dates."""
    for record in records:
        try:
            date = iso_date(record["sample_date"])
        except ValueError:
            tally["invalid_date"] += 1
            continue
        out = record.copy()
        out["sample_date"] = date
        yield out
