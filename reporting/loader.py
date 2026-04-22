"""Load generator artifacts for reporting."""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
GEN_OUTPUT = ROOT / "generator" / "output"
REPORT_OUTPUT = Path(__file__).resolve().parent / "output"

CSV_PATH = GEN_OUTPUT / "synthetic_dataset.csv"
PARAMS_PATH = GEN_OUTPUT / "parameters_used.json"


def load_dataset():
    return pd.read_csv(CSV_PATH)


def load_parameters():
    with open(PARAMS_PATH) as f:
        return json.load(f)


def ensure_output_dir():
    REPORT_OUTPUT.mkdir(parents=True, exist_ok=True)
    return REPORT_OUTPUT
