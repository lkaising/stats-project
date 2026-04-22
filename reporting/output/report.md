# Post-Generation Validation Report

- Input: `generator/output/synthetic_dataset.csv`
- Seed: `123`
- Rows: `240`

## Metric coherence

- Status: **PASS**
- Tolerance: `1e-09`
- Max |weber residual|: `1.94e-16`
- Max |michelson residual|: `1.67e-16`

## T1 — Cell counts (site × condition)

| site | 850 | 850_940 | 940 | 940_1050 | 1050 |
| --- | --- | --- | --- | --- | --- |
| dorsal_hand_L | 12 | 12 | 12 | 12 | 12 |
| dorsal_hand_R | 12 | 12 | 12 | 12 | 12 |
| antecubital_L | 12 | 12 | 12 | 12 | 12 |
| antecubital_R | 12 | 12 | 12 | 12 | 12 |

## T2 — Condition summary

| condition | mean_weber | sd_weber | cv_weber | mean_michelson |
| --- | --- | --- | --- | --- |
| 850 | 0.175787 | 0.015257 | 0.086792 | 0.096439 |
| 850_940 | 0.150855 | 0.022688 | 0.150398 | 0.081740 |
| 940 | 0.123196 | 0.015041 | 0.122090 | 0.065708 |
| 940_1050 | 0.097441 | 0.024755 | 0.254052 | 0.051390 |
| 1050 | 0.074800 | 0.018657 | 0.249421 | 0.038949 |

## T3 — Site × condition summary

| site | condition | mean_weber | sd_weber | cv_weber |
| --- | --- | --- | --- | --- |
| dorsal_hand_L | 850 | 0.182684 | 0.017945 | 0.098227 |
| dorsal_hand_L | 850_940 | 0.161750 | 0.020929 | 0.129393 |
| dorsal_hand_L | 940 | 0.123693 | 0.019249 | 0.155619 |
| dorsal_hand_L | 940_1050 | 0.099341 | 0.022607 | 0.227573 |
| dorsal_hand_L | 1050 | 0.077371 | 0.021304 | 0.275345 |
| dorsal_hand_R | 850 | 0.175587 | 0.012817 | 0.072995 |
| dorsal_hand_R | 850_940 | 0.152673 | 0.028910 | 0.189361 |
| dorsal_hand_R | 940 | 0.123144 | 0.015756 | 0.127944 |
| dorsal_hand_R | 940_1050 | 0.103816 | 0.025226 | 0.242986 |
| dorsal_hand_R | 1050 | 0.071018 | 0.019442 | 0.273758 |
| antecubital_L | 850 | 0.174410 | 0.014853 | 0.085164 |
| antecubital_L | 850_940 | 0.152030 | 0.013263 | 0.087238 |
| antecubital_L | 940 | 0.122937 | 0.010161 | 0.082655 |
| antecubital_L | 940_1050 | 0.085104 | 0.025484 | 0.299440 |
| antecubital_L | 1050 | 0.071517 | 0.018911 | 0.264428 |
| antecubital_R | 850 | 0.170468 | 0.014240 | 0.083537 |
| antecubital_R | 850_940 | 0.136967 | 0.020117 | 0.146871 |
| antecubital_R | 940 | 0.123010 | 0.015628 | 0.127046 |
| antecubital_R | 940_1050 | 0.101501 | 0.024274 | 0.239150 |
| antecubital_R | 1050 | 0.079296 | 0.015529 | 0.195838 |

## Plots

![Weber boxplot by condition](weber_boxplot_by_condition.png)

![Weber interaction by site](weber_interaction_by_site.png)

![CV by condition](cv_by_condition.png)
