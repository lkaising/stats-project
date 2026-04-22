"""Per-trial vein and background ROI intensity generation.

Order of operations for every (condition, site, trial):
  baseline -> + site offset -> + interaction -> + per-trial noise
  -> (+ registration noise if ratiometric).
"""

import pandas as pd


def _baseline_for_condition(params, condition):
    if condition in params["single_band_conditions"]:
        return params["mu_V"][condition], params["mu_B"][condition]
    lam_1, lam_2 = params["ratio_parents"][condition]
    base_V = 0.5 * (params["mu_V"][lam_1] + params["mu_V"][lam_2])
    base_B = 0.5 * (params["mu_B"][lam_1] + params["mu_B"][lam_2])
    return base_V, base_B


def _draw_trial(params, condition, site, rng):
    base_V, base_B = _baseline_for_condition(params, condition)
    delta_s = params["site_offsets"][site]
    eta_cs = params["interactions"][(site, condition)]

    pre_V = base_V + delta_s + eta_cs
    pre_B = base_B + delta_s + eta_cs

    eps_V = rng.normal(0.0, params["sigma_V"])
    eps_B = rng.normal(0.0, params["sigma_B"])

    if condition in params["ratio_conditions"]:
        eps_R_V = rng.normal(0.0, params["sigma_R"])
        eps_R_B = rng.normal(0.0, params["sigma_R"])
    else:
        eps_R_V = 0.0
        eps_R_B = 0.0

    return pre_V + eps_V + eps_R_V, pre_B + eps_B + eps_R_B


def build_dataset(params, rng):
    rows = []
    for condition in params["conditions"]:
        for site in params["sites"]:
            for trial in range(1, params["n_trials"] + 1):
                I_V, I_B = _draw_trial(params, condition, site, rng)
                rows.append({
                    "subject_id": params["subject_id"],
                    "site": site,
                    "condition": condition,
                    "trial": trial,
                    "I_V": I_V,
                    "I_B": I_B,
                })
    return pd.DataFrame(
        rows,
        columns=["subject_id", "site", "condition", "trial", "I_V", "I_B"],
    )
