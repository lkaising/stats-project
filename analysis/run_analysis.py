"""Entry point for the Step 6 formal analysis pipeline."""

import json
import math
import sys

import pandas as pd

from assumptions import run_assumption_checks
from bootstrap import BOOTSTRAP_CI, BOOTSTRAP_ITERS, BOOTSTRAP_SEED
from loader import ensure_output_dir, load_dataset, load_parameters
from report import render_report
from tests import compare_rankings, run_a1, run_a2, run_a3, run_michelson


def _json_default(obj):
    if isinstance(obj, float):
        if math.isnan(obj):
            return None
        return obj
    raise TypeError(f"Not JSON-serializable: {type(obj).__name__}")


def _write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, default=_json_default))


def _posthoc_dataframe(res):
    rows = []
    for r in res["posthoc"]:
        row = {"a": r["a"], "b": r["b"], "path": res["path"]}
        row.update({k: v for k, v in r.items() if k not in {"a", "b"}})
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    out_dir = ensure_output_dir()
    df = load_dataset()
    params = load_parameters()

    assumptions = run_assumption_checks(df)
    _write_json(out_dir / "assumption_checks.json", assumptions)

    a1 = run_a1(df, assumptions)
    _write_json(out_dir / "a1_results.json", a1)
    _posthoc_dataframe(a1).to_csv(out_dir / "a1_posthoc.csv", index=False)

    a2 = run_a2(df)
    _write_json(out_dir / "a2_results.json", a2)

    a3 = run_a3(df)
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
        f"bootstrap iters={BOOTSTRAP_ITERS}, CI={BOOTSTRAP_CI}, seed={BOOTSTRAP_SEED}."
    )


if __name__ == "__main__":
    main()
