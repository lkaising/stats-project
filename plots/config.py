"""Locked constants, styling, and small formatting helpers for ``plots/``."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent

ANALYSIS_OUTPUT_DIR = REPO_ROOT / "analysis" / "output"
GENERATOR_OUTPUT_DIR = REPO_ROOT / "generator" / "output"
PLOTS_OUTPUT_DIR = REPO_ROOT / "plots" / "output"

CONDITION_ORDER = ["850", "850_940", "940", "940_1050", "1050"]
SITE_ORDER = ["dorsal_hand_L", "dorsal_hand_R", "antecubital_L", "antecubital_R"]

CONDITION_DISPLAY_LABELS = {
    "850": "850",
    "940": "940",
    "1050": "1050",
    "850_940": "850/940",
    "940_1050": "940/1050",
}

CONDITION_PALETTE = {
    "850": "#1f77b4",
    "850_940": "#2ca02c",
    "940": "#ff7f0e",
    "940_1050": "#9467bd",
    "1050": "#d62728",
}

SITE_PALETTE = {
    "dorsal_hand_L": "#1f77b4",
    "dorsal_hand_R": "#ff7f0e",
    "antecubital_L": "#2ca02c",
    "antecubital_R": "#d62728",
}

FIGSIZE = (7.5, 5.0)
FIGSIZE_WIDE = (10.0, 4.5)
DPI = 150


def fmt_p(p: float) -> str:
    if p is None:
        return "p=n/a"
    if p < 0.001:
        return "p<0.001"
    return f"p={p:.3f}"


def fmt_num(x, digits: int = 3) -> str:
    if x is None:
        return "n/a"
    return f"{float(x):.{digits}f}"


def ensure_output_dir() -> Path:
    PLOTS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return PLOTS_OUTPUT_DIR
