"""A1 — Condition main effect on Weber contrast.

Produces two polished figures from ``a1_results.json``:

* ``a1_condition_means_ci.png`` — site-level mean ± 95% CI per condition.
* ``a1_posthoc_matrix.png`` — compact pairwise Holm-corrected matrix.
"""

from numbers import Real
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from config import (
    CONDITION_ORDER,
    CONDITION_PALETTE,
    DPI,
    FIGSIZE,
    fmt_num,
    fmt_p,
)
from loader import condition_means_df


_TOL = 1e-9


def _as_float(value: object, label: str) -> float:
    if not isinstance(value, Real):
        raise TypeError(f"{label} must be real-valued, got {type(value).__name__}")
    return float(value)


def _omnibus_line(a1: dict) -> str:
    omni = a1["omnibus"]
    es = a1["effect_size"]
    path = a1["path"]
    if path == "rm_anova":
        gg = omni.get("sphericity_correction_applied", False)
        corr = " (Greenhouse-Geisser corrected)" if gg else ""
        stat = (
            f"RM-ANOVA  F({omni['df1']:g}, {omni['df2']:g}) = "
            f"{fmt_num(omni['F'], 2)},  {fmt_p(omni['p_reported'])}{corr}"
        )
    else:
        stat = (
            f"Friedman  Q({omni['dof']:g}) = {fmt_num(omni['Q'], 2)},  "
            f"{fmt_p(omni['p'])}"
        )
    if es.get("value") is not None:
        label = "η²ₚ" if es["type"] == "partial_eta_sq" else "W"
        stat += f";  {label} = {fmt_num(es['value'], 3)}"
    return stat


def _provenance_line(a1: dict) -> str:
    ap = a1["assumption_path"]
    bits = [f"n_sites={a1['n_sites']}", f"α={a1['omnibus']['alpha']}"]
    if ap.get("fallback_used"):
        bits.append("Friedman fallback")
    if ap.get("sphericity_correction_applied"):
        bits.append("GG-corrected")
    return "  ·  ".join(bits)


def _numerical_guard(a1: dict) -> None:
    df = condition_means_df(a1["condition_means"])
    expected = {row["condition"]: row["mean"] for row in a1["condition_means"]}
    for cond, expected_mean in expected.items():
        plotted = _as_float(df.loc[cond, "mean"], f"{cond} plotted mean")
        expected_value = _as_float(expected_mean, f"{cond} expected mean")
        if abs(plotted - expected_value) > _TOL:
            raise ValueError(
                f"A1 spot check failed for {cond}: plotted={plotted} expected={expected_value}"
            )


def render_a1_mean_ci(a1: dict, out_path: Path) -> None:
    _numerical_guard(a1)
    df = condition_means_df(a1["condition_means"])

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(CONDITION_ORDER))
    means = df["mean"].to_numpy()
    lo = means - df["ci_low"].to_numpy()
    hi = df["ci_high"].to_numpy() - means
    colors = [CONDITION_PALETTE[c] for c in CONDITION_ORDER]

    ax.errorbar(
        x, means, yerr=[lo, hi], fmt="none",
        ecolor="#444444", elinewidth=1.5, capsize=4, zorder=1,
    )
    ax.scatter(x, means, c=colors, s=90, zorder=2, edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Weber contrast (site-level mean)")
    ax.set_title("A1  ·  Weber contrast by condition (mean ± 95% CI)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.08)

    fig.text(0.5, 0.035, _omnibus_line(a1), ha="center", fontsize=9)
    fig.text(0.5, 0.005, _provenance_line(a1), ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def render_a1_posthoc_matrix(a1: dict, out_path: Path) -> None:
    by_pair = {}
    for row in a1["posthoc"]:
        by_pair[(row["a"], row["b"])] = row
        by_pair[(row["b"], row["a"])] = row

    n = len(CONDITION_ORDER)
    reject_grid = np.full((n, n), np.nan)
    for i, a in enumerate(CONDITION_ORDER):
        for j, b in enumerate(CONDITION_ORDER):
            if i == j:
                continue
            row = by_pair.get((a, b))
            if row is None:
                continue
            reject_grid[i, j] = 1.0 if row["reject_holm_at_0.05"] else 0.0

    fig, ax = plt.subplots(figsize=FIGSIZE)
    cmap = mcolors.ListedColormap(["#f0f0f0", "#4a6fa5"])
    ax.imshow(np.nan_to_num(reject_grid, nan=-1), cmap=cmap, vmin=0, vmax=1)

    for i, a in enumerate(CONDITION_ORDER):
        for j, b in enumerate(CONDITION_ORDER):
            if i == j:
                ax.text(j, i, "—", ha="center", va="center",
                        color="#888888", fontsize=11)
                continue
            row = by_pair.get((a, b))
            if row is None:
                continue
            p_holm = row["p_holm"]
            reject = row["reject_holm_at_0.05"]
            color = "white" if reject else "black"
            weight = "bold" if reject else "normal"
            label = "p<0.001" if p_holm < 0.001 else f"p={p_holm:.3f}"
            ax.text(j, i, label, ha="center", va="center",
                    color=color, fontsize=8.5, fontweight=weight)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(CONDITION_ORDER, rotation=20)
    ax.set_yticklabels(CONDITION_ORDER)
    ax.set_title("A1  ·  Holm-corrected pairwise comparisons")

    path_label = "RM-ANOVA paired t-tests" if a1["path"] == "rm_anova" else "Wilcoxon signed-rank"
    fig.text(
        0.5, 0.03,
        f"{path_label}  ·  Holm correction  ·  α={a1['omnibus']['alpha']}  ·  "
        "filled = reject H₀",
        ha="center", fontsize=9,
    )

    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
