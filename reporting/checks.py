"""Post-generation metric coherence recheck."""

import numpy as np


def metric_coherence(df, tol=1e-9):
    recomputed_weber = (df["I_B"] - df["I_V"]) / df["I_B"]
    recomputed_michelson = (df["I_B"] - df["I_V"]) / (df["I_B"] + df["I_V"])
    weber_residual = float(np.max(np.abs(recomputed_weber - df["weber"])))
    michelson_residual = float(np.max(np.abs(recomputed_michelson - df["michelson"])))
    ok = weber_residual <= tol and michelson_residual <= tol
    return {
        "ok": ok,
        "tolerance": tol,
        "weber_max_abs_residual": weber_residual,
        "michelson_max_abs_residual": michelson_residual,
    }
