"""Summary table builders for the reporting layer."""

import pandas as pd

CONDITION_ORDER = ["850", "850_940", "940", "940_1050", "1050"]
SITE_ORDER = ["dorsal_hand_L", "dorsal_hand_R", "antecubital_L", "antecubital_R"]


def cell_counts(df):
    counts = df.groupby(["site", "condition"]).size().unstack("condition")
    return counts.reindex(index=SITE_ORDER, columns=CONDITION_ORDER)


def condition_summary(df):
    g = df.groupby("condition")
    out = pd.DataFrame({
        "mean_weber": g["weber"].mean(),
        "sd_weber": g["weber"].std(ddof=1),
        "mean_michelson": g["michelson"].mean(),
    })
    out["cv_weber"] = out["sd_weber"] / out["mean_weber"]
    out = out[["mean_weber", "sd_weber", "cv_weber", "mean_michelson"]]
    return out.reindex(CONDITION_ORDER)


def site_condition_summary(df):
    g = df.groupby(["site", "condition"])
    out = pd.DataFrame({
        "mean_weber": g["weber"].mean(),
        "sd_weber": g["weber"].std(ddof=1),
    })
    out["cv_weber"] = out["sd_weber"] / out["mean_weber"]
    out = out.reset_index()
    out["site"] = pd.Categorical(out["site"], SITE_ORDER, ordered=True)
    out["condition"] = pd.Categorical(out["condition"], CONDITION_ORDER, ordered=True)
    out = out.sort_values(["site", "condition"]).reset_index(drop=True)
    out["site"] = out["site"].astype(str)
    out["condition"] = out["condition"].astype(str)
    return out
