"""A2 — Condition-by-site interaction profile (conditional v1 artifact).

A2 ships only when the go/no-go criteria in the implementation plan hold:

1. ``a2_results.json`` is present with an assessable ``interaction_primary``.
2. ``a2_results.json::cell_means`` contains a complete descriptive
   site-condition matrix.
3. The figure renders cleanly at the module's default wide figure size.

If any criterion fails, ``render_a2_interaction`` logs a skip reason and
returns ``False`` so the runner can continue without failing the stage.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from config import (
    CONDITION_DISPLAY_LABELS,
    CONDITION_ORDER,
    DPI,
    FIGSIZE,
    SITE_ORDER,
    SITE_PALETTE,
    fmt_num,
    fmt_p,
)
from loader import ArtifactSchemaError, a2_cell_means_df, condition_means_df


_TOL = 1e-9
_A2_PRESENTATION_YMIN = 0.05
_A2_PRESENTATION_YMAX = 0.20

_SITE_LABELS = {
    "dorsal_hand_L": "DH-L",
    "dorsal_hand_R": "DH-R",
    "antecubital_L": "AC-L",
    "antecubital_R": "AC-R",
}

_SITE_MARKERS = {
    "dorsal_hand_L": "o",
    "dorsal_hand_R": "s",
    "antecubital_L": "^",
    "antecubital_R": "D",
}


def _condition_labels() -> list[str]:
    missing = [c for c in CONDITION_ORDER if c not in CONDITION_DISPLAY_LABELS]
    if missing:
        raise ValueError(f"Missing A2 condition display labels: {missing}")
    return [CONDITION_DISPLAY_LABELS[c] for c in CONDITION_ORDER]


def _pooled_a1_means(a1: dict) -> np.ndarray:
    if "condition_means" not in a1:
        raise ValueError("a1_results.json missing condition_means")
    df = condition_means_df(a1["condition_means"], condition_order=CONDITION_ORDER)
    if "mean" not in df.columns:
        raise ValueError("a1_results.json::condition_means missing mean column")
    if df["mean"].isna().any():
        missing = df.index[df["mean"].isna()].tolist()
        raise ValueError(
            f"a1_results.json::condition_means missing condition(s): {missing}"
        )
    try:
        values = df["mean"].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "a1_results.json::condition_means contains non-numeric mean values"
        ) from exc
    if not np.isfinite(values).all():
        raise ValueError("a1_results.json::condition_means contains non-finite means")
    return values


def _a2_cell_means(a2: dict):
    if "cell_means" not in a2:
        raise ValueError("a2_results.json missing cell_means")
    matrix = a2_cell_means_df(a2["cell_means"], condition_order=CONDITION_ORDER)
    expected_shape = (len(SITE_ORDER), len(CONDITION_ORDER))
    if matrix.shape != expected_shape:
        raise ValueError(
            f"a2_results.json::cell_means has shape {matrix.shape}; "
            f"expected {expected_shape}"
        )
    return matrix


def _guard_a1_a2_consistency(cell_means, pooled_means: np.ndarray) -> None:
    site_weighted_means = cell_means.mean(axis=0).to_numpy(dtype=float)
    for cond, a2_mean, a1_mean in zip(CONDITION_ORDER, site_weighted_means, pooled_means):
        if abs(a2_mean - a1_mean) > _TOL:
            raise ValueError(
                "A2/A1 mean consistency check failed for "
                f"{cond}: mean across A2 site cell means={a2_mean:.12g}, "
                f"official A1 pooled mean={a1_mean:.12g}. This strict guard "
                "assumes equal site weighting in the current balanced synthetic "
                "dataset; future unbalanced trial counts may require revisiting it."
            )


def _a2_ylim(cell_means, pooled_means: np.ndarray) -> tuple[float, float]:
    observed_min = min(
        float(cell_means.min().min()),
        float(np.min(pooled_means)),
    )
    observed_max = max(
        float(cell_means.max().max()),
        float(np.max(pooled_means)),
    )

    if observed_min >= _A2_PRESENTATION_YMIN and observed_max <= _A2_PRESENTATION_YMAX:
        return _A2_PRESENTATION_YMIN, _A2_PRESENTATION_YMAX

    span = max(observed_max - observed_min, 1e-9)
    lower = max(0.0, observed_min - 0.12 * span)
    upper = observed_max + 0.12 * span
    return lower, upper


def _fmt_df(value) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):g}"


def render_a2_interaction(
    a2: dict | None,
    a1: dict,
    out_path: Path,
) -> tuple[bool, str]:
    if a2 is None:
        return False, "a2_results.json not present"

    interaction = a2.get("interaction_primary")
    if interaction is None or interaction.get("p") is None or interaction.get("F") is None:
        return False, "interaction_primary missing F/p fields"

    if not a2.get("interaction_primary_assessable"):
        return False, "interaction_primary_assessable is False"

    try:
        xlabels = _condition_labels()
        cell_means = _a2_cell_means(a2)
        pooled_means = _pooled_a1_means(a1)
        _guard_a1_a2_consistency(cell_means, pooled_means)
    except (ArtifactSchemaError, KeyError, TypeError, ValueError) as exc:
        return False, str(exc)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(CONDITION_ORDER))
    reference_line, = ax.plot(
        x, pooled_means,
        color="#555555",
        linestyle=(0, (4, 3)),
        linewidth=1.35,
        alpha=0.68,
        label="Pooled A1 mean",
        zorder=1,
    )
    site_handles = []
    for site in SITE_ORDER:
        y = cell_means.loc[site].to_numpy()
        handle, = ax.plot(
            x, y,
            marker=_SITE_MARKERS[site],
            linewidth=1.8,
            markersize=7.5,
            markeredgecolor="white",
            markeredgewidth=0.7,
            color=SITE_PALETTE[site],
            label=_SITE_LABELS[site],
            zorder=3,
        )
        site_handles.append(handle)

    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Weber contrast (site-condition mean)")
    title = "A2 · Condition × site interaction"
    if a2.get("exploratory", True):
        title += " (exploratory)"
    ax.set_title(title)
    ax.set_ylim(*_a2_ylim(cell_means, pooled_means))
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.06)
    ax.legend(
        handles=site_handles + [reference_line],
        title="Site profiles / reference",
        loc="upper right",
        fontsize=8,
        title_fontsize=9,
    )

    denom_df = a2["interaction_denominator_df"]
    stat = (
        f"Two-way fixed-effects ANOVA (trial-level, n={a2.get('n_obs')}):  "
        f"F({_fmt_df(interaction['df'])}, {_fmt_df(denom_df)}) = {fmt_num(interaction['F'], 2)},  "
        f"{fmt_p(interaction['p'])};  η²ₚ = {fmt_num(interaction['partial_eta_sq'], 3)}  ·  "
        f"α = {a2['alpha']}"
    )
    note = (
        "Site profiles show descriptive A2 site-condition means; "
        "dashed gray line shows pooled A1 condition means for reference."
    )
    interpretation = f"Inferential test: {a2['interaction_interpretation']}"
    fig.text(0.5, 0.085, stat, ha="center", fontsize=9)
    fig.text(0.5, 0.050, note, ha="center", fontsize=8, color="#555555")
    fig.text(0.5, 0.020, interpretation, ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.16, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return True, "ok"
