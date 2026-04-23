"""Entry point for the polished presentation layer.

Required inputs (fail loud if missing):
  - analysis/output/a1_results.json
  - analysis/output/a3_results.json
  - analysis/output/michelson_results.json
  - generator/output/synthetic_dataset.csv

Optional:
  - analysis/output/a2_results.json (enables A2 interaction figure)

Required outputs written to plots/output/:
  - a1_condition_means_ci.png
  - a1_posthoc_matrix.png
  - a3_condition_cv_ci.png
  - michelson_vs_weber.png
  - a2_interaction_profile.png (only if A2 go/no-go criteria hold)
"""

import sys

from config import ensure_output_dir
from loader import (
    MissingArtifactError,
    ArtifactSchemaError,
    load_a1,
    load_a2,
    load_a3,
    load_michelson,
    load_raw_dataset,
)
from a1 import render_a1_mean_ci, render_a1_posthoc_matrix
from a2 import render_a2_interaction
from a3 import render_a3_condition_cv
from michelson import render_michelson_vs_weber


def main() -> int:
    out_dir = ensure_output_dir()

    try:
        a1 = load_a1()
        a3 = load_a3()
        michelson = load_michelson()
        raw_df = load_raw_dataset()
    except (MissingArtifactError, ArtifactSchemaError) as exc:
        print(f"plots: required input error — {exc}", file=sys.stderr)
        return 1

    a2 = load_a2(required=False)

    p_a1 = out_dir / "a1_condition_means_ci.png"
    p_a1_post = out_dir / "a1_posthoc_matrix.png"
    p_a3 = out_dir / "a3_condition_cv_ci.png"
    p_mich = out_dir / "michelson_vs_weber.png"
    p_a2 = out_dir / "a2_interaction_profile.png"

    render_a1_mean_ci(a1, p_a1)
    render_a1_posthoc_matrix(a1, p_a1_post)
    render_a3_condition_cv(a3, p_a3)
    render_michelson_vs_weber(a1, michelson, p_mich)

    ok, reason = render_a2_interaction(a2, raw_df, p_a2)
    if ok:
        a2_summary = f"A2 interaction figure written: {p_a2.name}"
    else:
        a2_summary = f"A2 interaction figure skipped — {reason}"

    print("plots stage complete")
    print(f"  wrote: {p_a1.name}")
    print(f"  wrote: {p_a1_post.name}")
    print(f"  wrote: {p_a3.name}")
    print(f"  wrote: {p_mich.name}")
    print(f"  {a2_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
