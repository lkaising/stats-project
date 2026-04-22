"""Load the synthetic dataset and generator parameters for analysis."""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
GEN_OUTPUT = ROOT / "generator" / "output"
ANALYSIS_OUTPUT = Path(__file__).resolve().parent / "output"

CSV_PATH = GEN_OUTPUT / "synthetic_dataset.csv"
PARAMS_PATH = GEN_OUTPUT / "parameters_used.json"

SITE_ORDER = ["dorsal_hand_L", "dorsal_hand_R", "antecubital_L", "antecubital_R"]
CONDITION_ORDER = ["850", "940", "1050", "850_940", "940_1050"]


def load_dataset():
    df = pd.read_csv(CSV_PATH)
    df["site"] = pd.Categorical(df["site"], categories=SITE_ORDER, ordered=True)
    df["condition"] = pd.Categorical(
        df["condition"].astype(str), categories=CONDITION_ORDER, ordered=True
    )
    return df


def load_parameters():
    if not PARAMS_PATH.exists():
        return {}
    with open(PARAMS_PATH) as f:
        return json.load(f)


def ensure_output_dir():
    ANALYSIS_OUTPUT.mkdir(parents=True, exist_ok=True)
    return ANALYSIS_OUTPUT
