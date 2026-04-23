"""Michelson — Sensitivity check alongside the A1 Weber result."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from config import (
    CONDITION_ORDER,
    CONDITION_PALETTE,
    DPI,
    FIGSIZE_WIDE,
    fmt_num,
    fmt_p,
)
from loader import condition_means_df


def _metric_line(payload: dict, label: str) -> str:
    omni = payload["omnibus"]
    path = payload["path"]
    if path == "rm_anova":
        gg = omni.get("sphericity_correction_applied", False)
        corr = " (GG)" if gg else ""
        return (
            f"{label}: RM-ANOVA F({omni['df1']:g},{omni['df2']:g})="
            f"{fmt_num(omni['F'], 2)}, {fmt_p(omni['p_reported'])}{corr}"
        )
    return (
        f"{label}: Friedman Q({omni['dof']:g})={fmt_num(omni['Q'], 2)}, {fmt_p(omni['p'])}"
    )


def _draw_panel(ax, df, title: str, ylabel: str) -> None:
    x = np.arange(len(CONDITION_ORDER))
    means = df["mean"].to_numpy()
    lo = means - df["ci_low"].to_numpy()
    hi = df["ci_high"].to_numpy() - means
    colors = [CONDITION_PALETTE[c] for c in CONDITION_ORDER]

    ax.errorbar(x, means, yerr=[lo, hi], fmt="none",
                ecolor="#444444", elinewidth=1.3, capsize=3, zorder=1)
    ax.scatter(x, means, c=colors, s=70, zorder=2, edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER, rotation=15)
    ax.set_xlabel("Condition")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.08)


def render_michelson_vs_weber(a1: dict, michelson: dict, out_path: Path) -> None:
    weber_df = condition_means_df(a1["condition_means"])
    mich_df = condition_means_df(michelson["condition_means"])

    fig, (ax_w, ax_m) = plt.subplots(1, 2, figsize=FIGSIZE_WIDE)
    _draw_panel(ax_w, weber_df,
                "Weber (primary)", "Weber contrast (site-level mean)")
    _draw_panel(ax_m, mich_df,
                "Michelson (sensitivity check)", "Michelson contrast (site-level mean)")

    fig.suptitle("Sensitivity check  ·  Weber vs. Michelson condition means (± 95% CI)",
                 fontsize=12)

    rc = michelson["ranking_comparison_vs_weber"]
    match_str = "rankings match" if rc["rankings_match"] else "rankings differ"
    caption = (
        f"{_metric_line(a1, 'Weber')}  ·  {_metric_line(michelson, 'Michelson')}  ·  "
        f"{match_str}"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=8.5, color="#333333")

    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
