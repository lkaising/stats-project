"""Entry point for post-generation reporting."""

import json
import sys

from loader import load_dataset, load_parameters, ensure_output_dir
from checks import metric_coherence
from tables import cell_counts, condition_summary, site_condition_summary
from plots import (
    weber_boxplot_by_condition,
    weber_interaction_by_site,
    cv_by_condition,
)


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6f}"
    return str(v)


def _df_to_markdown(df, index=True):
    if index:
        headers = [df.index.name or ""] + list(df.columns)
        rows = [[str(idx)] + [_fmt(v) for v in row] for idx, row in zip(df.index, df.values)]
    else:
        headers = list(df.columns)
        rows = [[_fmt(v) for v in row] for row in df.values]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def main():
    out_dir = ensure_output_dir()

    df = load_dataset()
    params = load_parameters()

    coherence = metric_coherence(df)
    (out_dir / "metric_coherence.json").write_text(json.dumps(coherence, indent=2))

    t1 = cell_counts(df)
    t2 = condition_summary(df)
    t3 = site_condition_summary(df)

    t1.to_csv(out_dir / "cell_counts.csv")
    t2.to_csv(out_dir / "condition_summary.csv")
    t3.to_csv(out_dir / "site_condition_summary.csv", index=False)

    p1 = out_dir / "weber_boxplot_by_condition.png"
    p2 = out_dir / "weber_interaction_by_site.png"
    p3 = out_dir / "cv_by_condition.png"
    weber_boxplot_by_condition(df, p1)
    weber_interaction_by_site(df, p2)
    cv_by_condition(df, p3)

    status = "PASS" if coherence["ok"] else "FAIL"
    report = [
        "# Post-Generation Validation Report",
        "",
        "- Input: `generator/output/synthetic_dataset.csv`",
        f"- Seed: `{params.get('seed')}`",
        f"- Rows: `{len(df)}`",
        "",
        "## Metric coherence",
        "",
        f"- Status: **{status}**",
        f"- Tolerance: `{coherence['tolerance']}`",
        f"- Max |weber residual|: `{coherence['weber_max_abs_residual']:.2e}`",
        f"- Max |michelson residual|: `{coherence['michelson_max_abs_residual']:.2e}`",
        "",
        "## T1 — Cell counts (site × condition)",
        "",
        _df_to_markdown(t1),
        "",
        "## T2 — Condition summary",
        "",
        _df_to_markdown(t2),
        "",
        "## T3 — Site × condition summary",
        "",
        _df_to_markdown(t3, index=False),
        "",
        "## Plots",
        "",
        f"![Weber boxplot by condition]({p1.name})",
        "",
        f"![Weber interaction by site]({p2.name})",
        "",
        f"![CV by condition]({p3.name})",
        "",
    ]
    (out_dir / "report.md").write_text("\n".join(report))

    if not coherence["ok"]:
        print(
            f"Metric coherence FAILED "
            f"(weber={coherence['weber_max_abs_residual']:.2e}, "
            f"michelson={coherence['michelson_max_abs_residual']:.2e})",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
