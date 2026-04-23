"""A3 — Condition-level repeatability (CV with bootstrap CIs)."""

from pathlib import Path

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
from loader import bootstrap_cis_df


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


def render_a3_condition_cv(a3: dict, out_path: Path) -> None:
    df = bootstrap_cis_df(a3["bootstrap_cis"])

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(CONDITION_ORDER))
    cv = df["cv_point_estimate"].to_numpy()
    lo = cv - df["ci_low"].to_numpy()
    hi = df["ci_high"].to_numpy() - cv
    colors = [CONDITION_PALETTE[c] for c in CONDITION_ORDER]

    ax.errorbar(
        x, cv, yerr=[lo, hi], fmt="none",
        ecolor="#444444", elinewidth=1.5, capsize=4, zorder=1,
    )
    ax.scatter(x, cv, c=colors, s=90, zorder=2, edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Coefficient of variation (CV)")
    ax.set_title("A3  ·  Condition-level repeatability (point estimate ± bootstrap CI)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.08)
    ax.set_ylim(bottom=0)

    fig.text(0.5, 0.035, _stat_line(a3), ha="center", fontsize=9)
    fig.text(0.5, 0.005, _provenance_line(a3), ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
