"""A1 — Condition main effect on Weber contrast.

Produces two polished figures from ``a1_results.json``:

* ``a1_condition_means_ci.png`` — site-level mean ± 95% CI per condition.
* ``a1_posthoc_matrix.png`` — compact pairwise Holm-corrected matrix.
"""

from numbers import Real
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from config import (
    # A1_MAIN_FIGURE_CONDITION_ORDER,
    CONDITION_DISPLAY_LABELS,
    CONDITION_ORDER,
    CONDITION_PALETTE,
    DPI,
    FIGSIZE,
    SITE_ORDER,
    fmt_num,
    fmt_p,
)
from loader import condition_means_df, site_condition_means_df


_TOL = 1e-9
_POSTHOC_DIAGONAL_FILL = "#ededed"
_POSTHOC_EMPTY_FILL = "#ffffff"
_POSTHOC_REJECT_EDGE = "#265f85"
_POSTHOC_NONREJECT_EDGE = "#b8b8b8"
_POSTHOC_MAGNITUDE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "a1_posthoc_magnitude",
    ["#f7fbff", "#d8e8f4", "#9ecae1"],
)


def _as_float(value: object, label: str) -> float:
    if not isinstance(value, Real):
        raise TypeError(f"{label} must be real-valued, got {type(value).__name__}")
    return float(value)


def _fmt_p_named(name: str, p: float) -> str:
    if p is None:
        return f"{name}=n/a"
    if p < 0.001:
        return f"{name}<.001"
    return f"{name}={p:.3f}"


def _fmt_trimmed(x: float, digits: int = 3) -> str:
    text = f"{float(x):.{digits}f}"
    if text.startswith("0."):
        return text[1:]
    if text.startswith("-0."):
        return "-" + text[2:]
    return text


def _fmt_p_holm(p: object) -> str:
    if p is None:
        return "pH=n/a"
    p_float = _as_float(p, "Holm-adjusted p-value")
    if not np.isfinite(p_float):
        return "pH=n/a"
    if p_float < 0.001:
        return "pH<.001"
    return f"pH={_fmt_trimmed(p_float)}"


def _fmt_delta(value: float, label: str) -> str:
    return f"{label}={value:+.3f}"


def _posthoc_difference_field(a1: dict) -> tuple[str, str]:
    if a1["path"] == "rm_anova":
        return "mean_diff", "Δ"
    if a1["path"] == "friedman":
        return "median_diff", "Δmed"
    raise ValueError(f"Unsupported A1 post-hoc path: {a1['path']!r}")


def _unordered_pair(a: str, b: str) -> frozenset[str]:
    if a == b:
        raise ValueError(f"A1 post-hoc row compares {a!r} to itself")
    return frozenset((a, b))


def _build_posthoc_pair_lookup(posthoc_rows: list, condition_order: list[str]) -> dict:
    known_conditions = set(condition_order)
    expected_pairs = {
        _unordered_pair(a, b)
        for i, a in enumerate(condition_order)
        for b in condition_order[i + 1:]
    }
    by_pair = {}
    for row in posthoc_rows:
        a = row.get("a")
        b = row.get("b")
        unknown = sorted({a, b} - known_conditions)
        if unknown:
            raise ValueError(f"A1 post-hoc row uses unknown condition key(s): {unknown}")
        pair = _unordered_pair(a, b)
        if pair in by_pair:
            labels = sorted(pair)
            raise ValueError(f"A1 post-hoc has duplicate row for pair: {labels}")
        by_pair[pair] = row

    missing = expected_pairs - set(by_pair)
    extra = set(by_pair) - expected_pairs
    if missing:
        labels = [sorted(pair) for pair in sorted(missing, key=lambda p: sorted(p))]
        raise ValueError(f"A1 post-hoc missing expected pair(s): {labels}")
    if extra:
        labels = [sorted(pair) for pair in sorted(extra, key=lambda p: sorted(p))]
        raise ValueError(f"A1 post-hoc has unexpected pair(s): {labels}")
    return by_pair


def _oriented_posthoc_delta(
    row: dict,
    row_condition: str,
    column_condition: str,
    field: str,
) -> float:
    delta = _as_float(row[field], f"A1 post-hoc {field}")
    if row["a"] == row_condition and row["b"] == column_condition:
        return delta
    if row["a"] == column_condition and row["b"] == row_condition:
        return -delta
    raise ValueError(
        "A1 post-hoc row does not match requested comparison: "
        f"{row.get('a')} vs {row.get('b')} for {row_condition} vs {column_condition}"
    )


def _posthoc_caption(a1: dict, delta_label: str) -> str:
    n_sites = a1["n_sites"]
    if a1["path"] == "rm_anova":
        path_label = "Paired t-tests"
        metric = "in Weber contrast"
    else:
        path_label = "Wilcoxon signed-rank tests"
        metric = "from official post-hoc rows"
    return (
        f"{path_label} on site-level means, n={n_sites} sites; "
        f"Holm-adjusted p-values. Cells show {delta_label}(row-column) "
        f"{metric} and pH; shading tracks |{delta_label}|, "
        "border marks Holm rejection."
    )


def _posthoc_text_color(fill: str) -> str:
    r, g, b = mcolors.to_rgb(fill)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if luminance < 0.45 else "#1f1f1f"


def _posthoc_plot_rows(a1: dict) -> list[dict]:
    for condition in CONDITION_ORDER:
        if condition not in CONDITION_DISPLAY_LABELS:
            raise ValueError(f"Missing display label for A1 condition: {condition}")

    field, delta_label = _posthoc_difference_field(a1)
    by_pair = _build_posthoc_pair_lookup(a1["posthoc"], CONDITION_ORDER)
    rows = []
    n = len(CONDITION_ORDER)
    for i, row_condition in enumerate(CONDITION_ORDER):
        for j, column_condition in enumerate(CONDITION_ORDER):
            if i >= j:
                continue
            pair = _unordered_pair(row_condition, column_condition)
            row = by_pair[pair]
            for required in (field, "p_holm", "reject_holm_at_0.05"):
                if required not in row:
                    raise ValueError(
                        f"A1 post-hoc row {row['a']} vs {row['b']} "
                        f"missing {required!r}"
                    )
            delta = _oriented_posthoc_delta(row, row_condition, column_condition, field)
            p_holm = row["p_holm"]
            reject = row["reject_holm_at_0.05"]
            if not isinstance(reject, bool):
                raise TypeError(
                    "A1 post-hoc reject_holm_at_0.05 must be boolean for "
                    f"{row['a']} vs {row['b']}"
                )
            rows.append(
                {
                    "i": i,
                    "j": j,
                    "delta": delta,
                    "delta_label": delta_label,
                    "p_holm": p_holm,
                    "reject": reject,
                }
            )

    expected_count = n * (n - 1) // 2
    if len(rows) != expected_count:
        raise ValueError(
            f"A1 post-hoc matrix would draw {len(rows)} comparisons; "
            f"expected {expected_count}"
        )
    return rows


def _omnibus_line(a1: dict) -> str:
    omni = a1["omnibus"]
    es = a1["effect_size"]
    path = a1["path"]
    if path == "rm_anova":
        gg = omni.get("sphericity_correction_applied", False)
        if gg:
            eps_gg = _as_float(omni["eps_gg"], "A1 eps_gg")
            df1 = _as_float(omni["df1"], "A1 df1") * eps_gg
            df2 = _as_float(omni["df2"], "A1 df2") * eps_gg
            stat = (
                f"RM-ANOVA, GG-corrected: F({df1:.2f}, {df2:.2f})="
                f"{fmt_num(omni['F'], 2)}, {_fmt_p_named('p_GG', omni['p_gg'])}, "
                f"ε_GG={_fmt_trimmed(omni['eps_gg'])}"
            )
        else:
            stat = (
                f"RM-ANOVA: F({omni['df1']:g}, {omni['df2']:g})="
                f"{fmt_num(omni['F'], 2)}, {fmt_p(omni['p_reported'])}"
            )
    else:
        stat = (
            f"Friedman: Q({omni['dof']:g})={fmt_num(omni['Q'], 2)}, "
            f"{fmt_p(omni['p'])}"
        )
    if es.get("value") is not None:
        label = "η²p" if es["type"] == "partial_eta_sq" else "W"
        stat += f", {label}={_fmt_trimmed(es['value'])}"
    return stat


def _provenance_line(a1: dict) -> str:
    return (
        "Colored points: condition means ±95% CI across sites; "
        "gray lines: anonymous site means; "
        f"repeated-measures unit=site, n={a1['n_sites']}."
    )


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


def _site_profile_guard(a1: dict, profiles) -> None:
    expected_rows = int(a1["n_sites"]) * int(a1["n_conditions"])
    if len(a1["site_condition_means"]) != expected_rows:
        raise ValueError(
            "A1 site_condition_means must contain one site-level mean per "
            f"site-condition cell; got {len(a1['site_condition_means'])}, "
            f"expected {expected_rows}"
        )

    official_means = {
        row["condition"]: _as_float(row["mean"], f"{row['condition']} official mean")
        for row in a1["condition_means"]
    }
    profile_means = profiles.mean(axis=0)
    for cond, official_mean in official_means.items():
        profile_mean = _as_float(profile_means.loc[cond], f"{cond} site-profile mean")
        if abs(profile_mean - official_mean) > _TOL:
            raise ValueError(
                f"A1 site-profile check failed for {cond}: "
                f"profile_mean={profile_mean} official_mean={official_mean}"
            )


def _a1_ylim(summary_df, profiles) -> tuple[float, float]:
    upper = max(
        _as_float(summary_df["ci_high"].max(), "A1 max CI high"),
        _as_float(profiles.max().max(), "A1 max site profile mean"),
    )
    return 0.0, upper * 1.12


def render_a1_mean_ci(a1: dict, out_path: Path) -> None:
    _numerical_guard(a1)
    df = condition_means_df(
        a1["condition_means"],
        condition_order=CONDITION_ORDER,
    )
    profiles = site_condition_means_df(
        a1["site_condition_means"],
        condition_order=CONDITION_ORDER,
    )
    _site_profile_guard(a1, profiles)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.arange(len(CONDITION_ORDER))
    means = df["mean"].to_numpy()
    lo = means - df["ci_low"].to_numpy()
    hi = df["ci_high"].to_numpy() - means
    colors = [CONDITION_PALETTE[c] for c in CONDITION_ORDER]

    for site in SITE_ORDER:
        y = profiles.loc[site].to_numpy()
        ax.plot(
            x, y, color="#8a8a8a", linewidth=1.0, alpha=0.45, zorder=1,
        )
        ax.scatter(
            x, y, color="#8a8a8a", s=24, alpha=0.6, linewidth=0, zorder=2,
        )

    ax.errorbar(
        x, means, yerr=[lo, hi], fmt="none",
        ecolor="#333333", elinewidth=1.6, capsize=4, zorder=3,
    )
    ax.scatter(x, means, c=colors, s=95, zorder=4, edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [CONDITION_DISPLAY_LABELS[c] for c in CONDITION_ORDER]
    )
    ax.set_xlabel("Condition")
    ax.set_ylabel("Weber contrast (site-level mean)")
    ax.set_title("A1  ·  Weber contrast by condition")
    ax.set_ylim(*_a1_ylim(df, profiles))
    ax.grid(True, alpha=0.3, axis="y")
    ax.margins(x=0.08)

    fig.text(0.5, 0.045, _omnibus_line(a1), ha="center", fontsize=9)
    fig.text(0.5, 0.015, _provenance_line(a1), ha="center", fontsize=8, color="#555555")

    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def render_a1_posthoc_matrix(a1: dict, out_path: Path) -> None:
    plot_rows = _posthoc_plot_rows(a1)
    _, delta_label = _posthoc_difference_field(a1)
    n = len(CONDITION_ORDER)
    max_abs_delta = max(abs(row["delta"]) for row in plot_rows) if plot_rows else 0.0

    fig, ax = plt.subplots(figsize=FIGSIZE)

    for i in range(n):
        rect = Rectangle(
            (i - 0.5, i - 0.5), 1, 1,
            facecolor=_POSTHOC_DIAGONAL_FILL,
            edgecolor="#d0d0d0",
            linewidth=0.9,
        )
        ax.add_patch(rect)
        ax.text(
            i, i, "—",
            ha="center", va="center",
            color="#777777", fontsize=12,
        )

    for row in plot_rows:
        i = row["i"]
        j = row["j"]
        relative_delta = abs(row["delta"]) / max_abs_delta if max_abs_delta else 0.0
        fill = mcolors.to_hex(_POSTHOC_MAGNITUDE_CMAP(0.16 + 0.68 * relative_delta))
        edgecolor = _POSTHOC_REJECT_EDGE if row["reject"] else _POSTHOC_NONREJECT_EDGE
        linewidth = 1.6 if row["reject"] else 1.0
        rect = Rectangle(
            (j - 0.5, i - 0.5), 1, 1,
            facecolor=fill,
            edgecolor=edgecolor,
            linewidth=linewidth,
        )
        ax.add_patch(rect)
        ax.text(
            j, i,
            f"{_fmt_delta(row['delta'], row['delta_label'])}\n"
            f"{_fmt_p_holm(row['p_holm'])}",
            ha="center", va="center",
            color=_posthoc_text_color(fill),
            fontsize=8.3,
            linespacing=1.25,
            fontweight="bold" if row["reject"] else "normal",
        )

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(
        [CONDITION_DISPLAY_LABELS[c] for c in CONDITION_ORDER],
        rotation=20,
        ha="left",
    )
    ax.set_yticklabels([CONDITION_DISPLAY_LABELS[c] for c in CONDITION_ORDER])
    ax.xaxis.tick_top()
    ax.tick_params(
        axis="x", top=True, bottom=False, labeltop=True, labelbottom=False,
        length=0, pad=5,
    )
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect("equal")
    ax.set_facecolor(_POSTHOC_EMPTY_FILL)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("A1  ·  Holm-corrected post-hoc comparisons", pad=34)

    fig.text(
        0.5, 0.04,
        _posthoc_caption(a1, delta_label),
        ha="center", fontsize=8.3, color="#444444", wrap=True,
    )

    fig.tight_layout(rect=(0, 0.12, 1, 0.96))
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
