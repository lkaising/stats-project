# analysis/

Runs the locked A1, A2, A3, and Michelson
sensitivity-check analyses on the single synthetic dataset produced by
`generator/`, and writes a markdown report plus per-analysis artifacts.

## What it does

- Reads `generator/output/synthetic_dataset.csv` (240 trial-level rows) and
  `generator/output/parameters_used.json`.
- Runs assumption checks (Shapiro-Wilk per condition on site-aggregated means;
  Mauchly's sphericity on the 4 × 5 site-aggregated matrix) and records the
  fallback decisions.
- **A1** — one-way repeated-measures ANOVA on site-level mean Weber contrast
  (site = repeated-measures unit, condition = within-unit factor).
  - Greenhouse-Geisser correction if sphericity is violated or noncomputable
    (Mauchly is singular when n_subjects ≤ k_conditions − 1; the pipeline
    defaults to GG conservatively in that case).
  - Friedman test nonparametric fallback if any condition fails Shapiro-Wilk.
  - Post-hoc: paired t-tests with Holm (parametric path) or Wilcoxon
    signed-rank with Holm (Friedman path), run only when the omnibus is
    significant.
  - Effect size: partial η² (parametric) or Kendall's W (Friedman).
- **A2** — two-way fixed-effects ANOVA on trial-level Weber
  (`condition + site + condition:site`). The interaction F-test is the
  primary quantity; treated as exploratory.
- **A3** — per-cell Weber CV (from the 12 trials per site × condition cell),
  Friedman across conditions with site as blocking factor, plus per-condition
  bootstrap 95% CIs (1000 iterations, fixed seed). A numerical stability guard
  (`CV_MEAN_EPS = 1e-6`) excludes near-zero-mean cells and routes the analysis
  to descriptive mode if any cell is unstable.
- **Michelson sensitivity check** — same pipeline as A1 on Michelson contrast,
  with a Weber-vs-Michelson ranking comparison.

## What it does not do

- Does not regenerate or modify the dataset.
- Does not re-verify metric coherence, re-run generator structural checks,
  re-compute descriptive tables, or reproduce the descriptive plots already
  emitted by `reporting/`.

## Running

From the project root, after the generator has been run at least once:

```
python analysis/run_analysis.py
```

Outputs are written to `analysis/output/`.

## Outputs

Written to `analysis/output/`:

- `analysis_report.md` — human-readable summary of all four analyses, with
  synthetic-data framing.
- `assumption_checks.json` — Shapiro-Wilk + Mauchly outputs and fallback
  decisions (per metric).
- `a1_results.json` — A1 condition means with 95% CIs, omnibus (RM ANOVA or
  Friedman), effect size, and Holm-corrected post-hoc rows.
- `a1_posthoc.csv` — the A1 post-hoc rows in tabular form.
- `a2_results.json` — full two-way ANOVA table with partial η².
- `a3_results.json` — A3 Friedman result (or descriptive fallback), CV matrix,
  and bootstrap CIs.
- `a3_bootstrap_cis.csv` — per-condition CV bootstrap 95% CIs.
- `michelson_results.json` — Michelson sensitivity check, mirroring A1's
  structure, plus the Weber-vs-Michelson ranking comparison.
- `michelson_posthoc.csv` — Michelson post-hoc rows in tabular form.

## Framing

Every p-value, F-statistic, effect size, and CI in this folder is a pipeline
output computed on a single **synthetic** dataset. They are pipeline
demonstration outputs, not substantive claims about NIR imaging. The
`analysis_report.md` header and closing section restate this explicitly.

## Dependencies

`numpy`, `pandas`, `scipy`, `pingouin` (all listed in the project
`requirements.txt`; `pingouin` transitively pulls `statsmodels`).
