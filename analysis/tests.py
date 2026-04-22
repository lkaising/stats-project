"""A1, A2, A3, and the Michelson sensitivity check."""

import math
from itertools import combinations
from typing import Any, cast

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

from assumptions import site_level_matrix
from bootstrap import bootstrap_cv_cis, cv_matrix
from loader import CONDITION_ORDER, SITE_ORDER

OMNIBUS_ALPHA = 0.05
CI_LEVEL = 0.95


def _matrix_to_long(matrix):
    return (
        matrix.reset_index()
        .melt(id_vars="site", var_name="condition", value_name="value")
    )


def _float_or_nan(value):
    return float(value) if pd.notna(value) else float("nan")


def _row_value(row, candidates):
    for col in candidates:
        if col in row.index:
            return _float_or_nan(row[col])
    return float("nan")


def _source_row(table, source_name, fallback_index):
    if "Source" in table.columns:
        matches = table["Source"].astype(str).str.casefold() == source_name.casefold()
        if matches.any():
            return table.loc[matches].iloc[0]
    return table.iloc[fallback_index]


def _condition_means_with_ci(matrix):
    """t-based 95% CI on the per-condition mean across the 4 site-level values."""
    n = matrix.shape[0]
    tcrit = stats.t.ppf(0.5 + CI_LEVEL / 2.0, df=n - 1)
    out = []
    for cond in matrix.columns:
        col = matrix[cond].to_numpy()
        mean = float(np.mean(col))
        sd = float(np.std(col, ddof=1))
        se = sd / math.sqrt(n)
        out.append(
            {
                "condition": str(cond),
                "mean": mean,
                "sd": sd,
                "n_sites": int(n),
                "ci_low": mean - tcrit * se,
                "ci_high": mean + tcrit * se,
            }
        )
    return out


def _condition_cv_point_estimates(matrix):
    out = []
    for cond in matrix.columns:
        col = matrix[cond].dropna().to_numpy(dtype=float)
        point_estimate = float(np.mean(col)) if col.size else float("nan")
        out.append(
            {
                "condition": str(cond),
                "cv_point_estimate": point_estimate,
            }
        )
    return out


def _rm_anova_path(matrix, sphericity_correction_applied):
    long_df = _matrix_to_long(matrix)
    aov = pg.rm_anova(
        data=long_df,
        dv="value",
        within="condition",
        subject="site",
        correction=cast(Any, True),
        detailed=True,
        effsize="np2",
    )
    row = _source_row(aov, "condition", 0)
    error_row = _source_row(aov, "error", 1) if len(aov) > 1 else None
    F = _row_value(row, ("F",))
    ddof1 = _row_value(row, ("ddof1", "DF"))
    ddof2 = _row_value(row, ("ddof2",))
    if math.isnan(ddof2) and error_row is not None:
        ddof2 = _row_value(error_row, ("ddof2", "DF"))
    p_unc = _row_value(row, ("p_unc", "p-unc"))
    p_gg = _row_value(row, ("p_GG_corr", "p-GG-corr"))
    eps_gg = _row_value(row, ("eps",))
    np2 = _row_value(row, ("np2",))
    p_reported = p_gg if sphericity_correction_applied else p_unc
    significant = bool(p_reported < OMNIBUS_ALPHA)

    posthoc = []
    if significant:
        for a, b in combinations(matrix.columns, 2):
            x = matrix[a].to_numpy()
            y = matrix[b].to_numpy()
            t_res = stats.ttest_rel(x, y)
            diff = x - y
            n = diff.size
            mean_diff = float(np.mean(diff))
            sd_diff = float(np.std(diff, ddof=1))
            se = sd_diff / math.sqrt(n)
            tcrit = stats.t.ppf(0.5 + CI_LEVEL / 2.0, df=n - 1)
            posthoc.append(
                {
                    "a": str(a),
                    "b": str(b),
                    "t": float(t_res.statistic),
                    "df": n - 1,
                    "p_raw": float(t_res.pvalue),
                    "mean_diff": mean_diff,
                    "ci_low": mean_diff - tcrit * se,
                    "ci_high": mean_diff + tcrit * se,
                }
            )
        _apply_holm(posthoc)

    return {
        "path": "rm_anova",
        "omnibus": {
            "F": F,
            "df1": ddof1,
            "df2": ddof2,
            "p_unc": p_unc,
            "p_gg": p_gg,
            "eps_gg": eps_gg,
            "p_reported": p_reported,
            "sphericity_correction_applied": bool(sphericity_correction_applied),
            "significant": significant,
            "alpha": OMNIBUS_ALPHA,
        },
        "effect_size": {"type": "partial_eta_sq", "value": np2},
        "posthoc": posthoc,
    }


def _friedman_path(matrix):
    long_df = _matrix_to_long(matrix)
    res = pg.friedman(
        data=long_df,
        dv="value",
        within="condition",
        subject="site",
        method="chisq",
    )
    row = res.iloc[0]
    Q = float(row["Q"])
    dof = float(row["ddof1"]) if "ddof1" in res.columns else float(row["dof"])
    p = float(row["p_unc"]) if "p_unc" in res.columns else float(row["p-unc"])
    W = float(row["W"]) if "W" in res.columns else float("nan")
    significant = bool(p < OMNIBUS_ALPHA)

    posthoc = []
    if significant:
        for a, b in combinations(matrix.columns, 2):
            x = matrix[a].to_numpy()
            y = matrix[b].to_numpy()
            try:
                w_res = cast(Any, stats.wilcoxon(x, y, zero_method="wilcox", alternative="two-sided"))
                w_stat = float(w_res.statistic)
                w_p = float(w_res.pvalue)
            except ValueError as exc:
                w_stat = float("nan")
                w_p = float("nan")
                note = f"{type(exc).__name__}: {exc}"
            else:
                note = None
            diff = x - y
            entry = {
                "a": str(a),
                "b": str(b),
                "wilcoxon_W": w_stat,
                "p_raw": w_p,
                "median_diff": float(np.median(diff)),
            }
            if note:
                entry["note"] = note
            posthoc.append(entry)
        _apply_holm(posthoc)

    return {
        "path": "friedman",
        "omnibus": {
            "Q": Q,
            "dof": dof,
            "p": p,
            "significant": significant,
            "alpha": OMNIBUS_ALPHA,
        },
        "effect_size": {"type": "kendalls_w", "value": W},
        "posthoc": posthoc,
    }


def _apply_holm(posthoc_rows):
    """Holm-Bonferroni on `p_raw`; annotates rows in place with
    `p_holm` and `reject_holm_at_0.05`. Rows whose p_raw is NaN
    are preserved and marked as not-rejected."""
    if not posthoc_rows:
        return
    valid = [(i, r["p_raw"]) for i, r in enumerate(posthoc_rows) if not math.isnan(r["p_raw"])]
    valid.sort(key=lambda t: t[1])
    m = len(valid)
    running_max = 0.0
    adjusted = {}
    for rank, (orig_i, p) in enumerate(valid):
        adj = p * (m - rank)
        adj = min(adj, 1.0)
        adj = max(adj, running_max)
        running_max = adj
        adjusted[orig_i] = adj
    for i, row in enumerate(posthoc_rows):
        if i in adjusted:
            row["p_holm"] = float(adjusted[i])
            row["reject_holm_at_0.05"] = bool(adjusted[i] < OMNIBUS_ALPHA)
        else:
            row["p_holm"] = float("nan")
            row["reject_holm_at_0.05"] = False


def run_rm_pipeline(df, metric, assumption_block):
    matrix = site_level_matrix(df, metric)
    means = _condition_means_with_ci(matrix)
    fallback = bool(assumption_block["friedman_fallback_triggered"])
    if fallback:
        omnibus = _friedman_path(matrix)
    else:
        omnibus = _rm_anova_path(
            matrix,
            sphericity_correction_applied=assumption_block["sphericity_correction_applied"],
        )
    ranking = sorted(means, key=lambda d: d["mean"], reverse=True)
    return {
        "metric": metric,
        "n_sites": int(matrix.shape[0]),
        "n_conditions": int(matrix.shape[1]),
        "assumption_path": {
            "fallback_used": fallback,
            "mauchly_status": assumption_block["mauchly"]["status"],
            "sphericity_correction_applied": bool(
                assumption_block["sphericity_correction_applied"]
            ),
            "sphericity_reason": assumption_block["sphericity_reason"],
        },
        "condition_means": means,
        "ranking_desc_by_mean": [d["condition"] for d in ranking],
        **omnibus,
    }


def run_a1(df, assumptions):
    return run_rm_pipeline(df, "weber", assumptions["weber"])


def run_michelson(df, assumptions):
    out = run_rm_pipeline(df, "michelson", assumptions["michelson"])
    return out


def compare_rankings(a1_res, michelson_res):
    return {
        "weber_ranking": a1_res["ranking_desc_by_mean"],
        "michelson_ranking": michelson_res["ranking_desc_by_mean"],
        "rankings_match": a1_res["ranking_desc_by_mean"] == michelson_res["ranking_desc_by_mean"],
    }


def run_a2(df):
    aov = pg.anova(
        data=df.assign(
            condition=df["condition"].astype(str),
            site=df["site"].astype(str),
        ),
        dv="weber",
        between=["condition", "site"],
        ss_type=2,
        detailed=True,
        effsize="np2",
    )
    table = []
    for _, row in aov.iterrows():
        entry = {
            "source": str(row["Source"]),
            "ss": float(row["SS"]),
            "df": float(row["DF"]),
        }
        for col_candidates, key in [
            (("MS",), "ms"),
            (("F",), "F"),
            (("p_unc", "p-unc"), "p"),
            (("np2",), "partial_eta_sq"),
        ]:
            val = None
            for col in col_candidates:
                if col in aov.columns:
                    val = row[col]
                    break
            if val is None:
                entry[key] = float("nan")
            else:
                entry[key] = float(val) if pd.notna(val) else float("nan")
        table.append(entry)
    interaction = next((e for e in table if e["source"] == "condition * site"), None)
    return {
        "metric": "weber",
        "unit_of_analysis": "trial",
        "n_obs": int(len(df)),
        "model": "two-way fixed-effects ANOVA: condition + site + condition:site",
        "interaction_primary": interaction,
        "anova_table": table,
        "exploratory": True,
        "note": "A2 is exploratory given the single-subject design and 4 sites.",
    }


def run_a3(df):
    matrix, unstable = cv_matrix(df, metric="weber")
    any_unstable = len(unstable) > 0
    point_estimates = _condition_cv_point_estimates(matrix)
    point_map = {
        row["condition"]: row["cv_point_estimate"] for row in point_estimates
    }

    bootstrap = bootstrap_cv_cis(df, metric="weber")
    for row in bootstrap:
        row["cv_point_estimate"] = point_map[row["condition"]]

    if any_unstable:
        return {
            "metric": "weber_cv",
            "unit_of_analysis": "per-cell CV",
            "a3_mode": "descriptive",
            "fallback_reason": "unstable_cells_present",
            "unstable_cells": unstable,
            "per_condition_cv": point_estimates,
            "cv_matrix": matrix.to_dict(),
            "bootstrap_cis": bootstrap,
        }

    long_df = (
        matrix.reset_index()
        .melt(id_vars="site", var_name="condition", value_name="cv")
    )
    fried = pg.friedman(
        data=long_df,
        dv="cv",
        within="condition",
        subject="site",
        method="chisq",
    )
    row = fried.iloc[0]
    Q = float(row["Q"])
    dof = float(row["ddof1"]) if "ddof1" in fried.columns else float(row["dof"])
    p = float(row["p_unc"]) if "p_unc" in fried.columns else float(row["p-unc"])
    W = float(row["W"]) if "W" in fried.columns else float("nan")
    return {
        "metric": "weber_cv",
        "unit_of_analysis": "per-cell CV",
        "a3_mode": "inferential",
        "unstable_cells": unstable,
        "per_condition_cv": point_estimates,
        "friedman": {
            "Q": Q,
            "dof": dof,
            "p": p,
            "kendalls_w": W,
            "significant": bool(p < OMNIBUS_ALPHA),
            "alpha": OMNIBUS_ALPHA,
        },
        "cv_matrix": matrix.to_dict(),
        "bootstrap_cis": bootstrap,
    }
