"""A2 — Condition-by-site interaction profile (conditional v1 artifact).

A2 ships only when the go/no-go criteria in the implementation plan hold:

1. ``a2_results.json`` is present with an assessable ``interaction_primary``.
2. Per-(site, condition) Weber cell means can be derived from the raw dataset
   via a single groupby over ``(site, condition)``.
3. The figure renders cleanly at the module's default wide figure size.

If any criterion fails, ``render_a2_interaction`` logs a skip reason and
returns ``False`` so the runner can continue without failing the stage.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import (
    CONDITION_ORDER,
    DPI,
    FIGSIZE_WIDE,
    SITE_ORDER,
    SITE_PALETTE,
    fmt_num,
    fmt_p,
)


def _cell_means_from_raw(df: pd.DataFrame) -> pd.DataFrame:
    means = (
        df.groupby(["site", "condition"], as_index=False)["weber"]
        .mean()
        .pivot(index="site", columns="condition", values="weber")
    )
    return means.reindex(index=SITE_ORDER, columns=CONDITION_ORDER)


def render_a2_interaction(
    a2: dict | None,
    raw_df: pd.DataFrame,
    out_path: Path,
) -> tuple[bool, str]:
    if a2 is None:
        return False, "a2_results.json not present"

    interaction = a2.get("interaction_primary")
    if interaction is None or interaction.get("p") is None or interaction.get("F") is None:
        return False, "interaction_primary missing F/p fields"

    if not a2.get("interaction_primary_assessable"):
        return False, "interaction_primary_assessable is False"

    cell_means = _cell_means_from_raw(raw_df)
    if cell_means.isna().any().any():
        missing = cell_means.isna().sum().sum()
        return False, f"raw-data cell means have {missing} missing cells"

    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)
    x = np.arange(len(CONDITION_ORDER))
    for site in SITE_ORDER:
        y = cell_means.loc[site].to_numpy()
        ax.plot(
            x, y, marker="o", linewidth=1.6, markersize=6,
            color=SITE_PALETTE[site], label=site,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Weber contrast (cell mean, display-only)")
    ax.set_title("A2  ·  Condition × site interaction profile")
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.06)
    ax.legend(title="Site", loc="best", fontsize=8, title_fontsize=9)

    denom_df = a2["interaction_denominator_df"]
    stat = (
        f"Two-way fixed-effects ANOVA (trial-level, n={a2.get('n_obs')}):  "
        f"F({interaction['df']:g}, {denom_df:g}) = {fmt_num(interaction['F'], 2)},  "
        f"{fmt_p(interaction['p'])};  η²ₚ = {fmt_num(interaction['partial_eta_sq'], 3)}  ·  "
        f"α = {a2['alpha']}"
    )
    note = (
        "Lines show per-cell Weber means (display-only). "
        f"Inferential test: {a2['interaction_interpretation']}"
    )
    fig.text(0.5, 0.04, stat, ha="center", fontsize=9)
    fig.text(0.5, 0.010, note, ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return True, "ok"
