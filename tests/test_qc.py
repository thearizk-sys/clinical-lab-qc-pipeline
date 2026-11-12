"""Unit tests for clinical laboratory parsing, unit conversion, and streaming QC."""
from collections import Counter

import pytest

from clinical_qc.analysis import consume_ranges, label_ranges
from clinical_qc.io import convert_records, iso_date, parse_value


def test_parse_comma_decimal():
    assert parse_value(" 5,4 ") == {"value": 5.4, "value_status": "numeric", "bound": None}


def test_parse_censored_value():
    assert parse_value("<0.1") == {"value": None, "value_status": "censored", "bound": 0.1}


def test_parse_missing_and_invalid_tokens():
    assert parse_value("NA")["value_status"] == "missing"
    assert parse_value("twelve")["value_status"] == "invalid"
    assert parse_value("nan")["value_status"] == "nonfinite"


def test_iso_date_formats():
    assert iso_date("2025-05-26") == "2025-05-26"
    assert iso_date("07/03/2025") == "2025-03-07"
    assert iso_date("05-May-2025") == "2025-05-05"
    with pytest.raises(ValueError):
        iso_date("not-a-date")


def test_unit_conversion_and_alias():
    config = {
        "aliases": {"hemoglobin": "haemoglobin"},
        "analytes": {
            "haemoglobin": {
                "canonical_unit": "g/L",
                "factors": {"g/l": 1.0, "g/dl": 10.0},
            }
        },
    }
    records = [{"analyte": "Hemoglobin", "unit": "g/dL", "value": 14.2, "value_status": "numeric", "bound": None}]
    tally = Counter()
    converted = list(convert_records(records, config, tally))
    assert converted[0]["analyte"] == "haemoglobin"
    assert converted[0]["unit"] == "g/L"
    assert converted[0]["value"] == pytest.approx(142.0)
    assert tally["converted"] == 1


def test_streaming_range_labelling():
    ranges = {"glucose": {"unit": "mmol/L", "low": 3.9, "high": 5.5}}
    records = [
        {"analyte": "glucose", "unit": "mmol/L", "value": 3.2, "value_status": "numeric"},
        {"analyte": "glucose", "unit": "mmol/L", "value": 4.8, "value_status": "numeric"},
        {"analyte": "glucose", "unit": "mmol/L", "value": 7.1, "value_status": "numeric"},
        {"analyte": "glucose", "unit": "mmol/L", "value": None, "value_status": "censored"},
    ]
    tally = Counter()
    counts = consume_ranges(label_ranges(records, ranges, tally))
    assert counts == Counter({"low": 1, "normal": 1, "high": 1, "censored": 1})
