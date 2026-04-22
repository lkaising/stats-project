"""Plot builders for the reporting layer."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from tables import CONDITION_ORDER, SITE_ORDER


def _style_axes(ax, title, xlabel, ylabel):
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def weber_boxplot_by_condition(df, out_path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    data = [df.loc[df["condition"] == c, "weber"].values for c in CONDITION_ORDER]
    ax.boxplot(data, showmeans=True)
    ax.set_xticks(range(1, len(CONDITION_ORDER) + 1))
    ax.set_xticklabels(CONDITION_ORDER)
    _style_axes(ax, "Weber contrast by condition", "Condition", "Weber")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def weber_interaction_by_site(df, out_path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    means = (
        df.groupby(["site", "condition"])["weber"]
        .mean()
        .unstack("condition")
        .reindex(index=SITE_ORDER, columns=CONDITION_ORDER)
    )
    for site in SITE_ORDER:
        ax.plot(CONDITION_ORDER, means.loc[site].values, marker="o", label=site)
    ax.legend(fontsize=8, loc="best")
    _style_axes(ax, "Mean Weber by condition, per site", "Condition", "Mean Weber")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def cv_by_condition(df, out_path):
    cell_cv = (
        df.groupby(["site", "condition"])["weber"]
        .agg(lambda x: x.std(ddof=1) / x.mean())
        .unstack("condition")
        .reindex(index=SITE_ORDER, columns=CONDITION_ORDER)
    )
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, c in enumerate(CONDITION_ORDER):
        vals = cell_cv[c].values
        ax.scatter(
            [i] * len(vals), vals, alpha=0.7, color="tab:blue",
            label="per-cell CV" if i == 0 else None,
        )
        ax.plot(
            [i - 0.25, i + 0.25], [vals.mean()] * 2, color="tab:red", linewidth=2,
            label="per-condition mean" if i == 0 else None,
        )
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.legend(fontsize=8, loc="best")
    _style_axes(ax, "Per-cell Weber CV by condition", "Condition", "CV")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
