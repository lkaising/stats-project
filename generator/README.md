# NIR Vein-Contrast Synthetic-Data Generator

Produces the single synthetic dataset used by the locked A1, A2, A3, and Michelson
sensitivity analyses. Values are simulated at the ROI-intensity level for vein and
background, and Weber and Michelson contrasts are derived deterministically from
the same intensities.

## Locked Dataset Structure

- 1 simulated subject (`S01`)
- 4 anatomical sites: `dorsal_hand_L`, `dorsal_hand_R`, `antecubital_L`, `antecubital_R`
- 5 imaging conditions: `850`, `940`, `1050`, `850_940`, `940_1050`
- 12 trials per (site × condition) cell
- 240 total rows

## Primary Generated Fields

Each row carries the design keys (`subject_id`, `site`, `condition`, `trial`) plus:

- `I_V` — vein ROI intensity (primary simulated, 0–1 scale)
- `I_B` — background ROI intensity (primary simulated, 0–1 scale)
- `weber` — `(I_B − I_V) / I_B` (derived)
- `michelson` — `(I_B − I_V) / (I_B + I_V)` (derived)

## Dependencies

- Python 3.14+
- `numpy`
- `pandas`

Installed via `pip install -r requirements.txt` at the project root.

## Running the Generator

From the project root:

```
python generator/generate.py
```

Optional self-test (regenerate and verify the output is byte-identical under the
locked seed):

```
python generator/generate.py --selftest
```

## Outputs

Written to `generator/output/`:

- `synthetic_dataset.csv` — the 240-row dataset with all fields above
- `parameters_used.json` — exact parameter dump used for this run (including seed)
- `run_log.txt` — timestamp, seed, row counts, validation summary, pooled mean
  Weber per condition, per-site Weber deviations, within-cell trial-level SDs

## Validation Checks

Every run verifies:

- **Bounds (hard)** — every `I_V`, `I_B`, `weber`, and `michelson` is strictly in
  `(0, 1)` and `I_V < I_B` in every row.
- **Row count (hard)** — exactly 240 rows and exactly 12 trials per (site,
  condition) cell.
- **Condition ordering (soft)** — pooled mean Weber satisfies
  `850 > 850_940 > 940 > 940_1050 > 1050`.
- **Ratio-condition variability (soft)** — mean within-cell trial-level SD for
  both ROIs is higher in ratiometric conditions than in single-band conditions,
  reflecting the extra registration-noise term σ_R.
- **Site-effect modesty (soft)** — |δ_s| ≤ 0.012 and per-site mean Weber
  deviates from the pooled mean by less than ±0.010.
- **Interaction modesty (soft)** — |η_{c,s}| ≤ 0.006 and per-condition site
  sums of η are ≈ 0.
- **Reproducibility (hard, on `--selftest`)** — regenerating with the same seed
  produces a byte-identical dataset.

Hard checks raise; soft checks warn and are logged.
