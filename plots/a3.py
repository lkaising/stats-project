"""A3 — Condition-level repeatability (CV with bootstrap CIs)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import (
    CONDITION_ORDER,
    CONDITION_PALETTE,
    DPI,
    FIGSIZE,
    SITE_ORDER,
    SITE_PALETTE,
    fmt_num,
    fmt_p,
)
from loader import bootstrap_cis_df, cv_matrix_df


_TOL = 1e-9
_SITE_X_OFFSETS = np.linspace(-0.18, 0.18, len(SITE_ORDER))


def _stat_line(a3: dict) -> str:
    mode = a3["a3_mode"]
    if mode == "inferential":
        fr = a3["friedman"]
        return (
            f"Inferential (Friedman)  ·  Q({fr['dof']:g}) = "
            f"{fmt_num(fr['Q'], 2)},  {fmt_p(fr['p'])}  ·  "
            f"Kendall's W = {fmt_num(fr['kendalls_w'], 3)}"
        )
    reason = a3.get("fallback_reason") or "unspecified"
    return f"Descriptive (fallback: {reason})"


def _provenance_line(a3: dict) -> str:
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


def _numerical_guard(cv_matrix: pd.DataFrame, point_df: pd.DataFrame) -> None:
    for cond in CONDITION_ORDER:
        stable_values = cv_matrix[cond].dropna().to_numpy(dtype=float)
        recomputed = float(np.mean(stable_values)) if stable_values.size else float("nan")
        summary_point = float(point_df.loc[cond, "cv_point_estimate"])
        if np.isnan(recomputed) and np.isnan(summary_point):
            continue
        if abs(recomputed - summary_point) > _TOL:
            raise ValueError(
                f"A3 spot check failed for {cond}: "
                f"summary_point={summary_point} mean_from_cv_matrix={recomputed}"
            )


def render_a3_condition_cv(a3: dict, out_path: Path) -> None:
    ci_df = bootstrap_cis_df(a3["bootstrap_cis"])
    point_df = _per_condition_cv_df(a3["per_condition_cv"])
    cv_matrix = cv_matrix_df(a3["cv_matrix"])
    _numerical_guard(cv_matrix, point_df)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(CONDITION_ORDER))
    cv = point_df["cv_point_estimate"].to_numpy(dtype=float)
    lo = cv - ci_df["ci_low"].to_numpy(dtype=float)
    hi = ci_df["ci_high"].to_numpy(dtype=float) - cv
    colors = [CONDITION_PALETTE[c] for c in CONDITION_ORDER]

    ax.errorbar(
        x, cv, yerr=[lo, hi], fmt="none",
        ecolor="#444444", elinewidth=1.5, capsize=4, zorder=1,
    )

    legend_handles = []
    for offset, site in zip(_SITE_X_OFFSETS, SITE_ORDER):
        vals = cv_matrix.loc[site].to_numpy(dtype=float)
        mask = ~np.isnan(vals)
        if not mask.any():
            continue
        handle = ax.scatter(
            x[mask] + offset,
            vals[mask],
            color=SITE_PALETTE[site],
            s=42,
            alpha=0.75,
            zorder=2,
            edgecolor="white",
            linewidth=0.5,
            label=site,
        )
        legend_handles.append(handle)

    ax.scatter(x, cv, c=colors, s=90, zorder=3, edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Coefficient of variation (CV)")
    ax.set_title("A3  ·  Condition-level repeatability (point estimate ± bootstrap CI)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.08)
    ax.set_ylim(bottom=0)
    if legend_handles:
        ax.legend(title="Site", loc="best", fontsize=8, title_fontsize=9)

    fig.text(0.5, 0.035, _stat_line(a3), ha="center", fontsize=9)
    fig.text(0.5, 0.005, _provenance_line(a3), ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
