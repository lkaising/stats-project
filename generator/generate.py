"""CLI entry point for the NIR vein-contrast synthetic-data generator."""

import argparse
import hashlib
import json
import sys
from dataclasses import fields
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from parameters import Parameters, load_parameters
from model import build_dataset
from metrics import attach_metrics
from validation import run_checks


OUTPUT_DIR = _HERE / "output"


def _params_to_jsonable(params: Parameters) -> dict:
    out = {}
    for field in fields(params):
        key = field.name
        value = getattr(params, key)

        if key == "interactions":
            out[key] = {f"{site}|{cond}": val for (site, cond), val in value.items()}
        elif key == "ratio_parents":
            out[key] = {k: list(v) for k, v in value.items()}
        elif hasattr(value, "items"):
            out[key] = dict(value)
        elif isinstance(value, tuple):
            out[key] = list(value)
        else:
            out[key] = value
    return out


def _run_once(seed: int | None) -> tuple[Parameters, pd.DataFrame]:
    params = load_parameters(seed=seed) if seed is not None else load_parameters()
    rng = np.random.default_rng(params.seed)
    df = attach_metrics(build_dataset(params, rng))
    return params, df

def _format_log(df: pd.DataFrame, params: Parameters, results: dict) -> str:
    lines = [
        f"Run timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"Seed: {params.seed}",
        f"Rows written: {len(df)}",
        "",
        "Hard checks:",
    ]
    for name, res in results["hard"].items():
        lines.append(f"  {name}: {'OK' if res['ok'] else 'FAIL'} — {res['message']}")
    lines.append("")
    lines.append("Soft checks:")
    for name, res in results["soft"].items():
        lines.append(f"  {name}: {'OK' if res['ok'] else 'WARN'}")
    lines.append("")
    lines.append("Condition-pooled mean Weber:")
    for cond, val in results["soft"]["condition_ordering"]["mean_weber"].items():
        lines.append(f"  {cond}: {val:.6f}")
    lines.append("")
    lines.append("Per-site Weber deviation from pooled mean:")
    for site, dev in results["soft"]["site_modesty"]["detail"]["per_site_deviation"].items():
        lines.append(f"  {site}: {dev:+.6f}")
    lines.append("")
    lines.append("Ratio vs single-band trial-level SD (per site):")
    for site, info in results["soft"]["ratio_variability"]["per_site"].items():
        lines.append(
            f"  {site}: "
            f"I_V sb={info['sb_sd_V']:.5f} rt={info['rt_sd_V']:.5f} | "
            f"I_B sb={info['sb_sd_B']:.5f} rt={info['rt_sd_B']:.5f}"
        )
    return "\n".join(lines) + "\n"

def _write_outputs(df: pd.DataFrame, params: Parameters, results: dict, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_path = output_dir / "synthetic_dataset.csv"
    df.to_csv(dataset_path, index=False)

    params_path = output_dir / "parameters_used.json"
    with params_path.open("w") as f:
        json.dump(_params_to_jsonable(params), f, indent=2)

    log_path = output_dir / "run_log.txt"
    with log_path.open("w") as f:
        f.write(_format_log(df, params, results))

    return dataset_path, params_path, log_path


def main(seed: int | None = None, selftest: bool = False) -> None:
    params, df = _run_once(seed)
    results = run_checks(df, params)

    dataset_path, params_path, log_path = _write_outputs(df, params, results, OUTPUT_DIR)

    for name, res in results["soft"].items():
        if not res["ok"]:
            print(f"WARNING: soft check '{name}' failed.")

    print(f"Wrote {dataset_path}")
    print(f"Wrote {params_path}")
    print(f"Wrote {log_path}")

    if selftest:
        _, df2 = _run_once(seed)
        h1 = hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()
        h2 = hashlib.sha256(df2.to_csv(index=False).encode()).hexdigest()
        if h1 != h2:
            raise RuntimeError(
                "Self-test failed: regenerating with the same seed produced a different dataset."
            )
        print(f"Self-test OK (sha256={h1[:16]}...)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate the NIR vein-contrast synthetic dataset.")
    parser.add_argument("--seed", type=int, default=None, help="Override the locked seed (for testing only).")
    parser.add_argument("--selftest", action="store_true", help="Regenerate and verify byte-identical output.")
    args = parser.parse_args()
    main(seed=args.seed, selftest=args.selftest)
