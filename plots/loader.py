"""Artifact loading for the ``plots/`` layer.

Centralizes all reads from ``analysis/output/`` and the raw dataset in
``generator/output/``. Every loader performs an existence preflight and a
minimal schema guard so downstream rendering code can assume a valid payload.
"""

import json

import pandas as pd

from config import (
    ANALYSIS_OUTPUT_DIR,
    CONDITION_ORDER,
    GENERATOR_OUTPUT_DIR,
)


class MissingArtifactError(FileNotFoundError):
    pass


class ArtifactSchemaError(ValueError):
    pass


def _read_json(filename: str, required: bool = True):
    path = ANALYSIS_OUTPUT_DIR / filename
    if not path.exists():
        if required:
            raise MissingArtifactError(
                f"Required analysis artifact missing: {path}"
            )
        return None
    return json.loads(path.read_text())


def _require_keys(payload: dict, keys, source: str) -> None:
    missing = [k for k in keys if k not in payload]
    if missing:
        raise ArtifactSchemaError(
            f"{source} missing expected keys: {missing}"
        )


def _assert_conditions_covered(conditions, source: str) -> None:
    missing = [c for c in CONDITION_ORDER if c not in conditions]
    if missing:
        raise ArtifactSchemaError(
            f"{source}: conditions {missing} absent; cannot reindex to canonical order"
        )


def load_a1() -> dict:
    data = _read_json("a1_results.json", required=True)
    _require_keys(
        data,
        ["condition_means", "path", "omnibus", "effect_size",
         "assumption_path", "posthoc", "ranking_desc_by_mean"],
        "a1_results.json",
    )
    conds = [row["condition"] for row in data["condition_means"]]
    _assert_conditions_covered(conds, "a1_results.json::condition_means")
    return data


def load_a3() -> dict:
    data = _read_json("a3_results.json", required=True)
    _require_keys(
        data,
        ["per_condition_cv", "bootstrap_cis", "a3_mode",
         "unstable_cells", "usable_site_count", "bootstrap"],
        "a3_results.json",
    )
    conds = [row["condition"] for row in data["bootstrap_cis"]]
    _assert_conditions_covered(conds, "a3_results.json::bootstrap_cis")
    return data


def load_michelson() -> dict:
    data = _read_json("michelson_results.json", required=True)
    _require_keys(
        data,
        ["condition_means", "path", "omnibus", "effect_size",
         "assumption_path", "ranking_comparison_vs_weber"],
        "michelson_results.json",
    )
    conds = [row["condition"] for row in data["condition_means"]]
    _assert_conditions_covered(conds, "michelson_results.json::condition_means")
    return data


def load_a2(required: bool = False):
    data = _read_json("a2_results.json", required=required)
    if data is None:
        return None
    _require_keys(
        data,
        ["interaction_primary", "interaction_denominator_df",
         "interaction_primary_assessable", "interaction_primary_significant",
         "interaction_interpretation"],
        "a2_results.json",
    )
    return data


def load_assumption_checks():
    return _read_json("assumption_checks.json", required=False)


def load_raw_dataset() -> pd.DataFrame:
    path = GENERATOR_OUTPUT_DIR / "synthetic_dataset.csv"
    if not path.exists():
        raise MissingArtifactError(
            f"Required raw dataset missing: {path}"
        )
    df = pd.read_csv(path)
    needed = {"site", "condition", "weber", "michelson"}
    missing = needed - set(df.columns)
    if missing:
        raise ArtifactSchemaError(
            f"synthetic_dataset.csv missing columns: {sorted(missing)}"
        )
    return df


def condition_means_df(condition_means: list) -> pd.DataFrame:
    df = pd.DataFrame(condition_means).set_index("condition")
    return df.reindex(CONDITION_ORDER)


def bootstrap_cis_df(bootstrap_cis: list) -> pd.DataFrame:
    df = pd.DataFrame(bootstrap_cis).set_index("condition")
    return df.reindex(CONDITION_ORDER)
