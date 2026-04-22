"""Weber and Michelson contrast formulas."""


def weber(I_V, I_B):
    return (I_B - I_V) / I_B


def michelson(I_V, I_B):
    return (I_B - I_V) / (I_B + I_V)


def attach_metrics(df):
    df = df.copy()
    df["weber"] = weber(df["I_V"], df["I_B"])
    df["michelson"] = michelson(df["I_V"], df["I_B"])
    return df
