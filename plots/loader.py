"""Artifact loading for the ``plots/`` layer.

Centralizes reads from ``analysis/output/``. Every loader performs an
existence preflight and a minimal schema guard so downstream rendering code
can assume a valid payload.
"""

import json
from collections.abc import Iterable
from typing import Any, Literal, TypeAlias, overload

import pandas as pd

from config import (
    ANALYSIS_OUTPUT_DIR,
    CONDITION_ORDER,
    SITE_ORDER,
)


class MissingArtifactError(FileNotFoundError):
    pass


class ArtifactSchemaError(ValueError):
    pass


JsonDict: TypeAlias = dict[str, Any]


@overload
def _read_json(filename: str, required: Literal[True] = True) -> JsonDict:
    ...


@overload
def _read_json(filename: str, required: Literal[False] = False) -> JsonDict | None:
    ...


def _read_json(filename: str, required: bool = True) -> JsonDict | None:
    path = ANALYSIS_OUTPUT_DIR / filename
    if not path.exists():
        if required:
            raise MissingArtifactError(
                f"Required analysis artifact missing: {path}"
            )
        return None
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ArtifactSchemaError(
            f"{filename} must contain a top-level JSON object"
        )
    return payload


def _require_keys(payload: JsonDict, keys: Iterable[str], source: str) -> None:
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


def _assert_sites_covered(sites, source: str) -> None:
    missing = [s for s in SITE_ORDER if s not in sites]
    if missing:
        raise ArtifactSchemaError(
            f"{source}: sites {missing} absent; cannot reindex to canonical order"
        )


def load_a1() -> JsonDict:
    data = _read_json("a1_results.json", required=True)
    _require_keys(
        data,
        ["condition_means", "path", "omnibus", "effect_size",
         "assumption_path", "posthoc", "ranking_desc_by_mean",
         "site_condition_means"],
        "a1_results.json",
    )
    conds = [row["condition"] for row in data["condition_means"]]
    _assert_conditions_covered(conds, "a1_results.json::condition_means")
    site_condition_means_df(data["site_condition_means"])
    return data


def load_a3() -> JsonDict:
    data = _read_json("a3_results.json", required=True)
    _require_keys(
        data,
        ["per_condition_cv", "bootstrap_cis", "cv_matrix", "a3_mode",
         "unstable_cells", "usable_site_count", "bootstrap"],
        "a3_results.json",
    )
    point_conds = [row["condition"] for row in data["per_condition_cv"]]
    _assert_conditions_covered(point_conds, "a3_results.json::per_condition_cv")
    conds = [row["condition"] for row in data["bootstrap_cis"]]
    _assert_conditions_covered(conds, "a3_results.json::bootstrap_cis")
    cv_matrix_df(data["cv_matrix"])
    return data


def load_michelson() -> JsonDict:
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


def load_a2(required: bool = False) -> JsonDict | None:
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


def load_assumption_checks() -> JsonDict | None:
    return _read_json("assumption_checks.json", required=False)


def condition_means_df(
    condition_means: list,
    condition_order: list[str] | None = None,
) -> pd.DataFrame:
    df = pd.DataFrame(condition_means).set_index("condition")
    return df.reindex(condition_order or CONDITION_ORDER)


def site_condition_means_df(
    site_condition_means: list,
    condition_order: list[str] | None = None,
) -> pd.DataFrame:
    source = "a1_results.json::site_condition_means"
    df = pd.DataFrame(site_condition_means)
    needed = {"site", "condition", "mean"}
    missing = needed - set(df.columns)
    if missing:
        raise ArtifactSchemaError(
            f"{source} missing columns: {sorted(missing)}"
        )
    _assert_conditions_covered(df["condition"].astype(str).tolist(), source)
    _assert_sites_covered(df["site"].astype(str).tolist(), source)
    matrix = df.pivot(index="site", columns="condition", values="mean")
    matrix = matrix.reindex(
        index=SITE_ORDER,
        columns=condition_order or CONDITION_ORDER,
    )
    if matrix.isna().any().any():
        missing_count = int(matrix.isna().sum().sum())
        raise ArtifactSchemaError(
            f"{source} has {missing_count} missing site-condition cells"
        )
    return matrix


def a2_cell_means_df(
    cell_means: list,
    condition_order: list[str] | None = None,
) -> pd.DataFrame:
    source = "a2_results.json::cell_means"
    df = pd.DataFrame(cell_means)
    needed = {"site", "condition", "mean"}
    missing = needed - set(df.columns)
    if missing:
        raise ArtifactSchemaError(
            f"{source} missing columns: {sorted(missing)}"
        )
    _assert_conditions_covered(df["condition"].astype(str).tolist(), source)
    _assert_sites_covered(df["site"].astype(str).tolist(), source)
    matrix = df.pivot(index="site", columns="condition", values="mean")
    matrix = matrix.reindex(
        index=SITE_ORDER,
        columns=condition_order or CONDITION_ORDER,
    )
    expected_shape = (len(SITE_ORDER), len(condition_order or CONDITION_ORDER))
    if matrix.shape != expected_shape:
        raise ArtifactSchemaError(
            f"{source} has shape {matrix.shape}; expected {expected_shape}"
        )
    if matrix.isna().any().any():
        missing_count = int(matrix.isna().sum().sum())
        raise ArtifactSchemaError(
            f"{source} has {missing_count} missing site-condition cells"
        )
    return matrix


def bootstrap_cis_df(bootstrap_cis: list) -> pd.DataFrame:
    df = pd.DataFrame(bootstrap_cis).set_index("condition")
    return df.reindex(CONDITION_ORDER)


def cv_matrix_df(cv_matrix: JsonDict) -> pd.DataFrame:
    source = "a3_results.json::cv_matrix"
    if not isinstance(cv_matrix, dict):
        raise ArtifactSchemaError(
            f"{source} must be a mapping of condition -> site -> value"
        )

    _assert_conditions_covered(cv_matrix.keys(), source)
    for cond in CONDITION_ORDER:
        site_map = cv_matrix[cond]
        if not isinstance(site_map, dict):
            raise ArtifactSchemaError(
                f"{source}[{cond!r}] must be a mapping of site -> value"
            )
        _assert_sites_covered(site_map.keys(), f"{source}[{cond!r}]")

    df = pd.DataFrame(cv_matrix)
    _assert_sites_covered(df.index, source)
    return df.reindex(index=SITE_ORDER, columns=CONDITION_ORDER)
