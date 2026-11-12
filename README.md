# Clinical Laboratory Quality-Control & Streaming ETL Pipeline

A reproducible, memory-efficient Python data pipeline designed to ingest, validate, normalise, summarise, and visualise multi-site clinical biomarker records.

## Key Features
* **Data Provenance & Structural Validation:** Tracks exact rejection counts (`Counter`) across malformed CSV lines (`too_few_fields`, `too_many_fields`), missing tokens, non-finite numbers, invalid dates, and left-censored (`< LOD`) values without halting execution.
* **Schema-Driven Unit Normalisation:** Reads from `units.json` to resolve analyte aliases (e.g., `hemoglobin` -> `haemoglobin`, `k+` -> `potassium`) and convert disparate units (`mg/dL`, `g/dL`, `mEq/L`) into canonical clinical standards.
* **ISO-8601 Date Standardisation & Deduplication:** Normalises mixed date formats (`YYYY-MM-DD`, `DD/MM/YYYY`, `DD-Mon-YYYY`) and applies deterministic deduplication per `(patient_id, sample_date, analyte)` while distinguishing identical vs. conflicting repeats.
* **Constant-Memory Streaming Mode:** Implements a lazy Python generator chain profiled via a custom `@timed` (`time.perf_counter`) decorator, capable of processing **2,000,000+ rows** without loading the dataset into RAM.
* **Vectorised QC Summary & Visualisation:** Uses NumPy and Pandas to compute per-analyte summary statistics, site-level out-of-range rates, and z-score outlier flags (`|z| > 3`), generating a multi-panel Matplotlib QC report comparing analyte distributions against clinical reference ranges.

## Tech Stack
* **Language:** Python 3.10+
* **Libraries:** `csv`, `json`, `math`, `time`, `datetime`, `pathlib`, `numpy`, `pandas`, `matplotlib`

## Project Structure
```text
├── clinical_qc_pipeline.py   # End-to-end ETL, streaming generators, stats & plotting
├── units.json                # Canonical unit conversion factors and analyte aliases
└── README.md                 # Project documentation
```

## Usage
1. Install dependencies:
   ```bash
   pip install numpy pandas matplotlib
   ```
2. Place your raw CSV (`patient_id,sample_date,site,analyte,value,unit`) and `reference_ranges.json` inside a `data/` directory.
3. Run the pipeline:
   ```bash
   python clinical_qc_pipeline.py
   ```
