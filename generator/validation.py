"""Post-generation sanity checks.

Hard checks raise ValidationError. Soft checks return ok=False but do not raise.
"""

from parameters import Parameters


class ValidationError(Exception):
    pass


def _hard_bounds(df):
    if not ((df["I_V"] > 0) & (df["I_V"] < 1)).all():
        return False, "I_V outside (0, 1)"
    if not ((df["I_B"] > 0) & (df["I_B"] < 1)).all():
        return False, "I_B outside (0, 1)"
    if not (df["I_V"] < df["I_B"]).all():
        return False, "I_V >= I_B in at least one row"
    if not ((df["weber"] > 0) & (df["weber"] < 1)).all():
        return False, "weber outside (0, 1)"
    if not ((df["michelson"] > 0) & (df["michelson"] < 1)).all():
        return False, "michelson outside (0, 1)"
    return True, "ok"


def _row_count(df, params: Parameters):
    expected = len(params.sites) * len(params.conditions) * params.n_trials
    if len(df) != expected:
        return False, f"row count {len(df)} != expected {expected}"
    counts = df.groupby(["site", "condition"]).size()
    if not (counts == params.n_trials).all():
        return False, "not all (site, condition) cells have the expected trial count"
    return True, "ok"


def _condition_ordering(df):
    pooled = df.groupby("condition")["weber"].mean()
    order = ["850", "850_940", "940", "940_1050", "1050"]
    values = [float(pooled[c]) for c in order]
    ok = all(values[i] > values[i + 1] for i in range(len(values) - 1))
    return ok, {c: v for c, v in zip(order, values)}


def _ratio_variability(df, params: Parameters):
    cell_sd = df.groupby(["site", "condition"])[["I_V", "I_B"]].std(ddof=1)
    per_site = {}
    overall_ok = True
    for site in params.sites:
        sb_cells = cell_sd.loc[site].loc[list(params.single_band_conditions)]
        rt_cells = cell_sd.loc[site].loc[list(params.ratio_conditions)]
        sb_sd_V = float(sb_cells["I_V"].mean())
        rt_sd_V = float(rt_cells["I_V"].mean())
        sb_sd_B = float(sb_cells["I_B"].mean())
        rt_sd_B = float(rt_cells["I_B"].mean())
        site_ok = rt_sd_V > sb_sd_V and rt_sd_B > sb_sd_B
        overall_ok = overall_ok and site_ok
        per_site[site] = {
            "sb_sd_V": sb_sd_V,
            "rt_sd_V": rt_sd_V,
            "sb_sd_B": sb_sd_B,
            "rt_sd_B": rt_sd_B,
            "ok": site_ok,
        }
    return overall_ok, per_site


def _site_modesty(df, params: Parameters):
    offsets_ok = all(abs(v) <= 0.012 + 1e-12 for v in params.site_offsets.values())
    per_site_weber = df.groupby("site")["weber"].mean()
    pooled = float(df["weber"].mean())
    deviations = {s: float(per_site_weber[s] - pooled) for s in params.sites}
    weber_ok = all(abs(d) < 0.010 for d in deviations.values())
    return (offsets_ok and weber_ok), {
        "offsets_ok": offsets_ok,
        "pooled_mean_weber": pooled,
        "per_site_deviation": deviations,
    }


def _interaction_modesty(params: Parameters):
    magnitudes_ok = all(abs(v) <= 0.006 + 1e-12 for v in params.interactions.values())
    per_condition_sums = {}
    sums_ok = True
    for c in params.conditions:
        s = sum(params.interactions[(site, c)] for site in params.sites)
        per_condition_sums[c] = s
        if abs(s) > 1e-9:
            sums_ok = False
    return (magnitudes_ok and sums_ok), {
        "magnitudes_ok": magnitudes_ok,
        "per_condition_sums": per_condition_sums,
    }


def run_checks(df, params: Parameters):
    bounds_ok, bounds_msg = _hard_bounds(df)
    count_ok, count_msg = _row_count(df, params)
    order_ok, order_info = _condition_ordering(df)
    ratio_ok, ratio_info = _ratio_variability(df, params)
    site_ok, site_info = _site_modesty(df, params)
    inter_ok, inter_info = _interaction_modesty(params)

    results = {
        "hard": {
            "bounds": {"ok": bounds_ok, "message": bounds_msg},
            "row_count": {"ok": count_ok, "message": count_msg},
        },
        "soft": {
            "condition_ordering": {"ok": order_ok, "mean_weber": order_info},
            "ratio_variability": {"ok": ratio_ok, "per_site": ratio_info},
            "site_modesty": {"ok": site_ok, "detail": site_info},
            "interaction_modesty": {"ok": inter_ok, "detail": inter_info},
        },
    }

    if not bounds_ok:
        raise ValidationError(f"Bounds check failed: {bounds_msg}")
    if not count_ok:
        raise ValidationError(f"Row count check failed: {count_msg}")

    return results
