# reporting/

External post-generation validation and reporting layer.

Reads `generator/output/synthetic_dataset.csv` and `generator/output/parameters_used.json`.
Writes summary tables, plots, and a single Markdown report to `reporting/output/`.

## What it does

- Re-derives Weber and Michelson from `I_V` and `I_B` and verifies they match the
  stored values within tolerance (`1e-9`).
- Builds three tables: cell counts (site × condition), condition summary,
  site × condition summary.
- Builds three plots: Weber boxplot by condition, Weber interaction by site,
  per-cell Weber CV by condition.
- Emits `report.md` indexing all of the above.

## What it does not do

- Does not regenerate or modify the dataset.
- Does not run any inferential tests (no A1/A2/A3, no p-values, no bootstrap CIs,
  no Friedman, no Holm).
- Does not re-run the structural checks already performed in
  `generator/validation.py`.

## Running

From the project root, after the generator has been run at least once:

```
python reporting/run_report.py
```

Outputs are written to `reporting/output/`. Exits non-zero only if the
metric-coherence recheck fails.

## Dependencies

`numpy`, `pandas`, `matplotlib` (all listed in the project `requirements.txt`).
