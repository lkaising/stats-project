"""Normality and sphericity checks with locked fallback decisions."""

import math

import pingouin as pg
from scipy import stats

from loader import CONDITION_ORDER, SITE_ORDER

NORMALITY_ALPHA = 0.05
SPHERICITY_ALPHA = 0.05


def site_level_matrix(df, metric):
    pivot = (
        df.groupby(["site", "condition"], observed=True)[metric]
        .mean()
        .unstack("condition")
    )
    pivot = pivot.reindex(index=SITE_ORDER, columns=CONDITION_ORDER)
    return pivot


def _shapiro_per_condition(matrix):
    results = []
    any_violated = False
    for cond in matrix.columns:
        col = matrix[cond].to_numpy()
        stat, p = stats.shapiro(col)
        violated = bool(p < NORMALITY_ALPHA)
        if violated:
            any_violated = True
        results.append(
            {
                "condition": cond,
                "W": float(stat),
                "p": float(p),
                "violated": violated,
            }
        )
    return results, any_violated


def _mauchly_plausible(W, chi2, dof, pval):
    """Mauchly outputs are only meaningful when 0 ≤ W ≤ 1, chi² ≥ 0,
    dof > 0, and 0 ≤ p ≤ 1. With n_subjects ≤ k_conditions − 1 the
    contrast covariance is singular and pingouin returns numerically
    meaningless values instead of raising."""
    for x in (W, chi2, dof, pval):
        if math.isnan(x):
            return False
    if W < 0.0 or W > 1.0:
        return False
    if chi2 < 0.0:
        return False
    if dof <= 0:
        return False
    if pval < 0.0 or pval > 1.0:
        return False
    return True


def _mauchly(matrix):
    n_subjects, n_conditions = matrix.shape
    if n_subjects <= n_conditions - 1:
        return {
            "status": "noncomputable",
            "W": float("nan"),
            "chi2": float("nan"),
            "dof": float("nan"),
            "p": float("nan"),
            "violated": True,
            "reason": "mauchly_noncomputable",
            "note": (
                f"n_subjects={n_subjects} ≤ k_conditions-1={n_conditions - 1}; "
                "contrast covariance is singular and Mauchly is undefined."
            ),
        }

    long_df = (
        matrix.reset_index()
        .melt(id_vars="site", var_name="condition", value_name="value")
    )
    try:
        res = pg.sphericity(
            data=long_df,
            dv="value",
            within="condition",
            subject="site",
            method="mauchly",
        )
        W = float(res.W) if res.W is not None else float("nan")
        chi2 = float(res.chi2) if res.chi2 is not None else float("nan")
        dof = float(res.dof) if res.dof is not None else float("nan")
        pval = float(res.pval) if res.pval is not None else float("nan")
        spher = bool(res.spher)
        if not _mauchly_plausible(W, chi2, dof, pval):
            return {
                "status": "noncomputable",
                "W": W,
                "chi2": chi2,
                "dof": dof,
                "p": pval,
                "violated": True,
                "reason": "mauchly_noncomputable",
                "note": "Mauchly returned values outside the valid range.",
            }
        violated = bool(pval < SPHERICITY_ALPHA) or not spher
        return {
            "status": "computed",
            "W": W,
            "chi2": chi2,
            "dof": dof,
            "p": pval,
            "violated": violated,
            "reason": "mauchly_rejected" if violated else "mauchly_accepted",
        }
    except Exception as exc:
        return {
            "status": "noncomputable",
            "W": float("nan"),
            "chi2": float("nan"),
            "dof": float("nan"),
            "p": float("nan"),
            "violated": True,
            "reason": "mauchly_noncomputable",
            "error": f"{type(exc).__name__}: {exc}",
        }


def check_metric(df, metric):
    matrix = site_level_matrix(df, metric)
    shapiro, any_violated = _shapiro_per_condition(matrix)
    mauchly = _mauchly(matrix)
    return {
        "metric": metric,
        "matrix_shape": list(matrix.shape),
        "shapiro_alpha": NORMALITY_ALPHA,
        "shapiro_per_condition": shapiro,
        "normality_violated": any_violated,
        "friedman_fallback_triggered": any_violated,
        "sphericity_alpha": SPHERICITY_ALPHA,
        "mauchly": mauchly,
        "sphericity_correction_applied": bool(mauchly["violated"]),
        "sphericity_reason": mauchly["reason"],
    }


def run_assumption_checks(df):
    return {
        "weber": check_metric(df, "weber"),
        "michelson": check_metric(df, "michelson"),
    }
