"""Entry point for the Step 6 formal analysis pipeline."""

import json
import math
import sys

import numpy as np
import pandas as pd

from assumptions import run_assumption_checks
from bootstrap import BOOTSTRAP_CI, BOOTSTRAP_ITERS
from loader import ensure_output_dir, load_dataset, load_parameters
from report import render_report
from tests import compare_rankings, run_a1, run_a2, run_a3, run_michelson


def _sanitize_for_json(obj):
    if obj is None or isinstance(obj, (str, bool, int)):
        return obj
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize_for_json(v) for v in obj]
    if isinstance(obj, tuple):
        return [_sanitize_for_json(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        obj = float(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    return obj


def _write_json(path, payload):
    sanitized = _sanitize_for_json(payload)
    path.write_text(json.dumps(sanitized, indent=2, allow_nan=False))


def _posthoc_dataframe(res):
    rows = []
    for r in res["posthoc"]:
        row = {"a": r["a"], "b": r["b"], "path": res["path"]}
        row.update({k: v for k, v in r.items() if k not in {"a", "b"}})
        rows.append(row)
    return pd.DataFrame(rows)


def _resolve_bootstrap_seed(params):
    seed = params.get("seed")
    if isinstance(seed, bool) or seed is None:
        raise RuntimeError(
            "Analysis requires generator/output/parameters_used.json to contain a valid integer seed."
        )
    if isinstance(seed, int):
        return seed
    if isinstance(seed, float) and seed.is_integer():
        return int(seed)
    if isinstance(seed, str):
        try:
            return int(seed)
        except ValueError:
            pass
    raise RuntimeError(
        "Analysis requires generator/output/parameters_used.json to contain a valid integer seed."
    )


def main():
    out_dir = ensure_output_dir()
    df = load_dataset()
    params = load_parameters()
    bootstrap_seed = _resolve_bootstrap_seed(params)

    assumptions = run_assumption_checks(df)
    _write_json(out_dir / "assumption_checks.json", assumptions)

    a1 = run_a1(df, assumptions)
    _write_json(out_dir / "a1_results.json", a1)
    _posthoc_dataframe(a1).to_csv(out_dir / "a1_posthoc.csv", index=False)

    a2 = run_a2(df)
    _write_json(out_dir / "a2_results.json", a2)

    a3 = run_a3(df, bootstrap_seed=bootstrap_seed)
    _write_json(out_dir / "a3_results.json", a3)
    pd.DataFrame(a3["bootstrap_cis"]).to_csv(out_dir / "a3_bootstrap_cis.csv", index=False)

    michelson = run_michelson(df, assumptions)
    comparison = compare_rankings(a1, michelson)
    michelson_payload = {**michelson, "ranking_comparison_vs_weber": comparison}
    _write_json(out_dir / "michelson_results.json", michelson_payload)
    _posthoc_dataframe(michelson).to_csv(
        out_dir / "michelson_posthoc.csv", index=False
    )

    report_md = render_report(
        params=params,
        n_rows=len(df),
        assumptions=assumptions,
        a1=a1,
        a2=a2,
        a3=a3,
        michelson=michelson,
        comparison=comparison,
    )
    (out_dir / "analysis_report.md").write_text(report_md)

    print(
        f"analysis complete — {len(df)} trial rows; "
        f"A1 path={a1['path']}, significant={a1['omnibus'].get('significant')}; "
        f"A2 interaction p={a2['interaction_primary']['p']:.4g}; "
        f"A3 mode={a3['a3_mode']}; "
        f"bootstrap iters={BOOTSTRAP_ITERS}, CI={BOOTSTRAP_CI}, seed={a3['bootstrap']['seed']}."
    )


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1) from exc
