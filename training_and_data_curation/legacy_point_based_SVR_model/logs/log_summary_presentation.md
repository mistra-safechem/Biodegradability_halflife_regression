# SVR Half-Life Models — Summary (Plain Regression)

Generated: 2026-04-21 16:56 | Source: `logs/log_review_20260421_141215.md`

> **Note:** This pipeline treats half-life targets as point estimates on a log scale — it is a plain regression, not an interval-aware model. Metrics are absolute RMSE and MAE in days and R², not interval coverage or rank-correlation. RMSE values appear large because no forgiveness is given for predictions near-but-outside a persistence class boundary. Compare to the interval-aware pipeline (`WP3bioddeg_V2`) with care.

---

## Model Performance Overview (Reduced-Feature Final Model)

> CV–Test gap = CV R² minus full-feature held-out Test R². Values > 0.15 indicate overfitting or out-of-distribution test split (flagged in warnings). **Bold** gap values exceed the threshold.

| Compartment | Dataset | CV R² | Test R² (full) | CV–Test gap | RMSE (days) | MAE (days) | Features (→) | Test in AD |
|---|---|---|---|---|---|---|---|---|
| Air | HSBD | 0.763 ±0.0850 | 0.355 | **0.408** | 185.39 | 36.02 | 212→150 | 52.1% |
| Water | HSBD | 0.418 ±0.0768 | 0.415 | 0.003 | 210.74 | 53.18 | 246→144 | 72.2% |
| Water | VEGA | 0.672 ±0.0522 | 0.891 | **-0.219** | 340.03 | 118.94 | 198→59 | 67.9% |
| Soil | HSBD | 0.424 ±0.0412 | 0.132 | **0.292** | 3,085.17 | 345.23 | 281→147 | 72.6% |
| Soil | VEGA | 0.792 ±0.1204 | 0.567 | **0.225** | 394.49 | 154.3 | 199→129 | 40.7% |
| Sediment | HSBD | 0.472 ±0.0354 | 0.466 | 0.006 | 531.26 | 251.86 | 246→158 | 34.6% |
| Sediment | VEGA | 0.790 ±0.0556 | 0.564 | **0.226** | 550.37 | 301.22 | 188→131 | 43.4% |
| Water | COMBINED | 0.528 ±0.0567 | 0.505 | 0.023 | 245.67 | 56.84 | 250→114 | 78.3% |
| Soil | COMBINED | 0.490 ±0.1002 | 0.408 | 0.082 | 430.04 | 138.55 | 286→122 | 79.3% |
| Sediment | COMBINED | 0.602 ±0.0656 | 0.522 | 0.080 | 589.85 | 193.23 | 243→79 | 58.8% |

---

## Key Findings

- **Scope:** 10 SVR plain regression models across 4 compartments (air, sediment, soil, water) and 3 data sources (COMBINED, HSBD, VEGA). Feature reduction strategies used: A  top-SVR (importance>0); C  top-N pca99. Targets are treated as point estimates on a log scale — this is not an interval-aware pipeline.

- **Soil HSBD is the clear failure case.** Full-feature test R² = 0.13 and RMSE = 3,085.17 days — the model learns something in cross-validation (CV R² = 0.42) but generalises almost not at all to the held-out test set (CV–test gap = 0.29). More representative training data is the only viable fix.

- **4 models show material overfitting (CV–test gap > 0.15):** Air HSBD (gap = 0.41), Soil HSBD (gap = 0.29), Sediment VEGA (gap = 0.23), Soil VEGA (gap = 0.23). In all cases the test set is likely out of distribution relative to the training fold structure — stronger regularisation or rebalanced splits are warranted.

- **Water VEGA is an anomaly that should not be trusted.** The full-feature test R² (0.89) is implausibly higher than CV R² (0.67), and reduced-model RMSE is 340.03 days — the high test R² almost certainly reflects a favourable random split, not genuine generalisation. Learning curve validation R² signal: Water VEGA — learning curve validation R² = -4.70; severe overfitting on small training subsets..

- **RMSE on the raw day scale is dominated by high-persistence outliers and is not comparable across compartments.** Air HSBD RMSE = 185 days vs Soil HSBD RMSE = 3,085 days — this reflects the target distribution scale, not relative model quality. Log-scale RMSE is the more meaningful cross-compartment metric.

- **AD coverage is insufficient in 5 models** (< 60% of test compounds inside the leverage-based AD): Sediment HSBD (34.6%), Soil VEGA (40.7%), Sediment VEGA (43.4%), Air HSBD (52.1%), Sediment COMBINED (58.8%). The majority or near-majority of their predictions are structural extrapolations. Best AD coverage: Soil COMBINED (79.3%), Water COMBINED (78.3%), Soil HSBD (72.6%), Water HSBD (72.2%), Water VEGA (67.9%).

- **VEGA-source models have higher CV R² but raw RMSE is not straightforwardly better.** On log-scale RMSE (more honest cross-compartment comparison), VEGA models range 0.26–0.43 log₁₀ d, vs HSBD 0.41–0.65 log₁₀ d. The raw-day RMSE contrast is driven by extreme outliers in HSBD soil and sediment.

- **Combined-source models are the most practically reliable.** Water, Soil, and Sediment COMBINED all show small CV–test gaps (max 0.08) and the best AD coverage of any source group. They do not reach the CV R² of the best VEGA models, but they generalise more consistently.

- **Overall assessment:** Models with small CV–test gap and adequate AD coverage: Water HSBD, Water VEGA, Water COMBINED, Soil COMBINED. These are the only candidates for cautious screening use. All remaining models require data augmentation, regularisation tuning, or train/test rebalancing before applied use.

---

## Warnings

| Model | Warning |
|---|---|
| Air HSBD | reduced-model CV R² (0.763) vs full-feature held-out R² (0.355) gap > 0.15; test set may be out of distribution. |
| Water VEGA | learning curve validation R² = -4.70; severe overfitting on small training subsets. |
| Soil HSBD | very low full-feature test R² (0.132); model generalises poorly. |
| Soil HSBD | reduced-model CV R² (0.424) vs full-feature held-out R² (0.132) gap > 0.15; test set may be out of distribution. |
| Soil VEGA | only 40.7% of test set inside AD; many predictions are extrapolations. |
| Soil VEGA | reduced-model CV R² (0.792) vs full-feature held-out R² (0.567) gap > 0.15; test set may be out of distribution. |
| Sediment HSBD | only 34.6% of test set inside AD; many predictions are extrapolations. |
| Sediment VEGA | only 43.4% of test set inside AD; many predictions are extrapolations. |
| Sediment VEGA | learning curve validation R² = -6.09; severe overfitting on small training subsets. |
| Sediment VEGA | reduced-model CV R² (0.790) vs full-feature held-out R² (0.564) gap > 0.15; test set may be out of distribution. |
