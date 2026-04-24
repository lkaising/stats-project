# `plots/`

Polished presentation layer for the project's completed analyses.

`plots/` is downstream of `analysis/`. It consumes analysis artifacts and
renders presentation-ready figures; it does **not** run new statistical tests,
recompute p-values, or replace any inferential decision made upstream.

## Role in the pipeline

```
generator/  →  reporting/  →  analysis/  →  plots/
```

`reporting/` owns fast dataset-facing QC visuals. `analysis/` owns the locked
inferential pipeline. `plots/` owns final presentation figures derived from
the analysis outputs.

## Inputs

Required (the stage fails loud if any are missing):

- `analysis/output/a1_results.json`
- `analysis/output/a3_results.json`
- `analysis/output/michelson_results.json`
- `generator/output/synthetic_dataset.csv`

Optional (enables additional figures if present):

- `analysis/output/a2_results.json` — gates the A2 interaction figure.

## Outputs

Written to `plots/output/`:

| File | Source |
| --- | --- |
| `a1_condition_means_ci.png` | `a1_results.json::condition_means` + exported site means |
| `a1_posthoc_matrix.png` | `a1_results.json::posthoc` |
| `a3_condition_cv_ci.png` | `a3_results.json::bootstrap_cis` + `a3_results.json::cv_matrix` |
| `michelson_vs_weber.png` | `a1_results.json` + `michelson_results.json` |
| `a2_interaction_profile.png` | `a2_results.json` + raw cell means (conditional) |

The A2 figure is only produced when `a2_results.json` is present, contains an
assessable interaction term, and the raw dataset yields a complete
`(site × condition)` cell-means matrix. If any of those fail, the runner logs
a skip reason and continues; the stage still exits 0.

## How to run

From the repository root:

```
python plots/run_plots.py
```

The plots stage is also wired into `run_pipeline.py` as the final stage.

## Conventions

- Canonical condition order (locked in v1): `["850", "850_940", "940", "940_1050", "1050"]`.
- A1 main-figure display order: `["850", "940", "1050", "850_940", "940_1050"]`.
- Canonical site order: `["dorsal_hand_L", "dorsal_hand_R", "antecubital_L", "antecubital_R"]`.
- Library: matplotlib only.
- Captions surface the analysis pathway (RM-ANOVA vs. Friedman fallback,
  Greenhouse-Geisser correction, A3 inferential vs. descriptive mode) so the
  presentation layer stays consistent with the upstream analysis.
- The A3 figure layers site-level CV context points from `cv_matrix` under the
  condition-level point estimate and bootstrap CI summary.
