"""Per-condition CV bootstrap CIs with the A3 stability guard."""

import math

import numpy as np
import pandas as pd

from loader import CONDITION_ORDER, SITE_ORDER

CV_MEAN_EPS = 1e-6
BOOTSTRAP_ITERS = 1000
BOOTSTRAP_CI = 0.95
BOOTSTRAP_SEED = 42


def _cell_cv(values):
    arr = np.asarray(values, dtype=float)
    mean = arr.mean()
    if abs(mean) < CV_MEAN_EPS:
        return float("nan")
    return float(arr.std(ddof=1) / mean)


def cv_matrix(df, metric="weber"):
    rows = []
    unstable = []
    for site in SITE_ORDER:
        row = {}
        for cond in CONDITION_ORDER:
            sel = df[(df["site"] == site) & (df["condition"] == cond)][metric].to_numpy()
            cv = _cell_cv(sel)
            row[cond] = cv
            if math.isnan(cv):
                unstable.append({"site": site, "condition": cond})
        rows.append(pd.Series(row, name=site))
    mat = pd.DataFrame(rows)[CONDITION_ORDER]
    mat.index.name = "site"
    return mat, unstable


def bootstrap_cv_cis(
    df,
    metric="weber",
    iterations=BOOTSTRAP_ITERS,
    ci=BOOTSTRAP_CI,
    seed=BOOTSTRAP_SEED,
):
    rng = np.random.default_rng(seed)
    cell_arrays = {}
    for site in SITE_ORDER:
        for cond in CONDITION_ORDER:
            arr = df[(df["site"] == site) & (df["condition"] == cond)][metric].to_numpy()
            cell_arrays[(site, cond)] = arr

    per_condition_draws = {cond: [] for cond in CONDITION_ORDER}
    for _ in range(iterations):
        for cond in CONDITION_ORDER:
            stable_cvs = []
            for site in SITE_ORDER:
                arr = cell_arrays[(site, cond)]
                idx = rng.integers(0, arr.size, size=arr.size)
                draw = arr[idx]
                mean = draw.mean()
                if abs(mean) < CV_MEAN_EPS:
                    continue
                stable_cvs.append(draw.std(ddof=1) / mean)
            if len(stable_cvs) >= 2:
                per_condition_draws[cond].append(float(np.mean(stable_cvs)))

    alpha = (1.0 - ci) / 2.0
    results = []
    for cond in CONDITION_ORDER:
        draws = np.asarray(per_condition_draws[cond], dtype=float)
        if draws.size == 0:
            results.append(
                {
                    "condition": cond,
                    "ci_low": float("nan"),
                    "ci_high": float("nan"),
                    "iterations_kept": 0,
                }
            )
            continue
        results.append(
            {
                "condition": cond,
                "ci_low": float(np.quantile(draws, alpha)),
                "ci_high": float(np.quantile(draws, 1.0 - alpha)),
                "iterations_kept": int(draws.size),
            }
        )
    return results
