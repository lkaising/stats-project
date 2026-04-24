"""A3 — Condition-level repeatability (CV with bootstrap CIs)."""

from numbers import Real
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import (
    CONDITION_DISPLAY_LABELS,
    CONDITION_ORDER,
    CONDITION_PALETTE,
    DPI,
    FIGSIZE,
    SITE_ORDER,
    fmt_num,
)
from loader import bootstrap_cis_df, cv_matrix_df


_TOL = 1e-9


def _as_float(value: object, label: str) -> float:
    if not isinstance(value, Real):
        raise TypeError(f"{label} must be real-valued, got {type(value).__name__}")
    return float(value)


def _fmt_trimmed(x: float, digits: int = 3) -> str:
    text = f"{float(x):.{digits}f}"
    if text.startswith("0."):
        return text[1:]
    if text.startswith("-0."):
        return "-" + text[2:]
    return text


def _fmt_p_compact(p: float) -> str:
    if p is None:
        return "p=n/a"
    if p < 0.001:
        return "p<.001"
    return f"p={_fmt_trimmed(p)}"


def _stat_line(a3: dict) -> str:
    mode = a3["a3_mode"]
    if mode == "inferential":
        fr = a3["friedman"]
        return (
            "Friedman test on per-site CVs: "
            f"χ²({fr['dof']:g})={fmt_num(fr['Q'], 2)}, "
            f"{_fmt_p_compact(fr['p'])}, Kendall's W={fmt_num(fr['kendalls_w'], 3)}"
        )
    reason = a3.get("fallback_reason") or "unspecified"
    return f"Descriptive (fallback: {reason})"


def _provenance_line(a3: dict) -> str:
    ci_level = int(a3["bootstrap"]["ci"] * 100)
    return (
        "Colored points = mean CV across sites; "
        f"error bars = bootstrap {ci_level}% CI; "
        "gray lines = site CV profiles. "
        "Lower CV = greater repeatability."
    )


def _metadata_line(a3: dict) -> str:
    bs = a3["bootstrap"]
    bits = [
        f"bootstrap iters={bs['iterations']}",
        f"CI={int(bs['ci'] * 100)}%",
        f"seed={bs['seed']}",
        f"usable sites={a3['usable_site_count']}",
    ]
    n_unstable = len(a3.get("unstable_cells") or [])
    if n_unstable:
        bits.append(f"{n_unstable} unstable cells excluded")
    return "  ·  ".join(bits)


def _per_condition_cv_df(per_condition_cv: list) -> pd.DataFrame:
    df = pd.DataFrame(per_condition_cv).set_index("condition")
    return df.reindex(CONDITION_ORDER)


def _condition_labels() -> list[str]:
    missing = [c for c in CONDITION_ORDER if c not in CONDITION_DISPLAY_LABELS]
    if missing:
        raise ValueError(f"Missing A3 condition display labels: {missing}")
    return [CONDITION_DISPLAY_LABELS[c] for c in CONDITION_ORDER]


def _require_columns(df: pd.DataFrame, columns: list[str], source: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{source} missing required column(s): {missing}")


def _require_complete_inferential_inputs(
    point_df: pd.DataFrame,
    ci_df: pd.DataFrame,
    cv_matrix: pd.DataFrame,
) -> None:
    expected_index = list(CONDITION_ORDER)
    if list(point_df.index) != expected_index:
        raise ValueError(
            "A3 per_condition_cv is not aligned to CONDITION_ORDER: "
            f"{list(point_df.index)}"
        )
    if list(ci_df.index) != expected_index:
        raise ValueError(
            "A3 bootstrap_cis is not aligned to CONDITION_ORDER: "
            f"{list(ci_df.index)}"
        )
    if list(cv_matrix.columns) != expected_index:
        raise ValueError(
            "A3 cv_matrix columns are not aligned to CONDITION_ORDER: "
            f"{list(cv_matrix.columns)}"
        )
    if list(cv_matrix.index) != list(SITE_ORDER):
        raise ValueError(
            "A3 cv_matrix rows are not aligned to SITE_ORDER: "
            f"{list(cv_matrix.index)}"
        )

    summary_values = point_df["cv_point_estimate"].to_numpy(dtype=float)
    ci_values = ci_df[["cv_point_estimate", "ci_low", "ci_high"]].to_numpy(dtype=float)
    site_values = cv_matrix.to_numpy(dtype=float)
    if not np.isfinite(summary_values).all():
        raise ValueError("A3 inferential summary CV point estimates must be finite")
    if not np.isfinite(ci_values).all():
        raise ValueError("A3 inferential bootstrap CI values must be finite")
    if not np.isfinite(site_values).all():
        raise ValueError("A3 inferential cv_matrix values must be finite")


def _point_ci_guard(point_df: pd.DataFrame, ci_df: pd.DataFrame) -> None:
    for cond in CONDITION_ORDER:
        summary_point = _as_float(
            point_df.loc[cond, "cv_point_estimate"],
            f"{cond} per_condition_cv CV point estimate",
        )
        ci_point = _as_float(
            ci_df.loc[cond, "cv_point_estimate"],
            f"{cond} bootstrap_cis CV point estimate",
        )
        ci_low = _as_float(ci_df.loc[cond, "ci_low"], f"{cond} bootstrap CI low")
        ci_high = _as_float(ci_df.loc[cond, "ci_high"], f"{cond} bootstrap CI high")

        values = np.array([summary_point, ci_point, ci_low, ci_high], dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"A3 non-finite point/CI value for {cond}")
        if abs(summary_point - ci_point) > _TOL:
            raise ValueError(
                f"A3 point estimate mismatch for {cond}: "
                f"per_condition_cv={summary_point} bootstrap_cis={ci_point}"
            )
        if not ci_low <= summary_point <= ci_high:
            raise ValueError(
                f"A3 CI bounds do not contain point estimate for {cond}: "
                f"ci_low={ci_low} point={summary_point} ci_high={ci_high}"
            )


def _numerical_guard(cv_matrix: pd.DataFrame, point_df: pd.DataFrame) -> None:
    for cond in CONDITION_ORDER:
        stable_values = cv_matrix[cond].dropna().to_numpy(dtype=float)
        if not stable_values.size:
            raise ValueError(f"A3 cv_matrix has no finite values for {cond}")
        recomputed = float(np.mean(stable_values))
        summary_point = _as_float(
            point_df.loc[cond, "cv_point_estimate"],
            f"{cond} CV point estimate",
        )
        if not np.isfinite(recomputed) or not np.isfinite(summary_point):
            raise ValueError(f"A3 non-finite numerical guard value for {cond}")
        if abs(recomputed - summary_point) > _TOL:
            raise ValueError(
                f"A3 spot check failed for {cond}: "
                f"summary_point={summary_point} mean_from_cv_matrix={recomputed}"
            )


def _a3_ylim(
    point_df: pd.DataFrame,
    ci_df: pd.DataFrame,
    cv_matrix: pd.DataFrame,
) -> tuple[float, float]:
    values = np.concatenate(
        [
            point_df["cv_point_estimate"].to_numpy(dtype=float),
            ci_df["ci_high"].to_numpy(dtype=float),
            cv_matrix.to_numpy(dtype=float).ravel(),
        ]
    )
    finite_values = values[np.isfinite(values)]
    if not finite_values.size:
        raise ValueError("A3 has no finite plotted values for y-axis scaling")
    upper = float(np.max(finite_values))
    return 0.0, upper * 1.12


def render_a3_condition_cv(a3: dict, out_path: Path) -> None:
    xlabels = _condition_labels()
    ci_df = bootstrap_cis_df(a3["bootstrap_cis"]).reindex(CONDITION_ORDER)
    point_df = _per_condition_cv_df(a3["per_condition_cv"]).reindex(CONDITION_ORDER)
    cv_matrix = cv_matrix_df(a3["cv_matrix"]).reindex(
        index=SITE_ORDER,
        columns=CONDITION_ORDER,
    )
    _require_columns(
        point_df,
        ["cv_point_estimate"],
        "a3_results.json::per_condition_cv",
    )
    _require_columns(
        ci_df,
        ["cv_point_estimate", "ci_low", "ci_high"],
        "a3_results.json::bootstrap_cis",
    )
    if a3["a3_mode"] == "inferential":
        _require_complete_inferential_inputs(point_df, ci_df, cv_matrix)
    _point_ci_guard(point_df, ci_df)
    _numerical_guard(cv_matrix, point_df)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(CONDITION_ORDER))
    cv = point_df["cv_point_estimate"].to_numpy(dtype=float)
    lo = cv - ci_df["ci_low"].to_numpy(dtype=float)
    hi = ci_df["ci_high"].to_numpy(dtype=float) - cv
    colors = [CONDITION_PALETTE[c] for c in CONDITION_ORDER]

    for site in SITE_ORDER:
        vals = cv_matrix.loc[site].to_numpy(dtype=float)
        mask = ~np.isnan(vals)
        if not mask.any():
            continue
        ax.plot(
            x[mask],
            vals[mask],
            color="#8a8a8a",
            linewidth=1.0,
            alpha=0.45,
            zorder=1,
        )
        ax.scatter(
            x[mask],
            vals[mask],
            color="#8a8a8a",
            s=24,
            alpha=0.6,
            linewidth=0,
            zorder=2,
        )

    ax.errorbar(
        x, cv, yerr=[lo, hi], fmt="none",
        ecolor="#333333", elinewidth=1.6, capsize=4, zorder=3,
    )
    ax.scatter(x, cv, c=colors, s=95, zorder=4, edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Coefficient of variation (CV)")
    ax.set_title("A3 · Repeatability by condition (Weber CV)")
    ax.set_ylim(*_a3_ylim(point_df, ci_df, cv_matrix))
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.08)

    fig.text(0.5, 0.070, _stat_line(a3), ha="center", fontsize=9)
    fig.text(0.5, 0.040, _provenance_line(a3), ha="center", fontsize=8, color="#555555")
    fig.text(0.5, 0.014, _metadata_line(a3), ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
