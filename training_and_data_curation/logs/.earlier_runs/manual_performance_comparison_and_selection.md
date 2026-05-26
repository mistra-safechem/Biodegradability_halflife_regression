# Model Selection & Cross-Pipeline Comparison

Date: 2026-04-21 (based on earlier run, but due to (thanks to) reproducibility, the outcome is the same)
Pipelines compared:
- **Interval-aware SVR** — `logs/log_summary_presentation.md`
- **Plain regression SVR** — `logs/log_summary_presentation.md`

---

# Overview
This file was created by AI prompted summary of the results for decision support and should be taken with appropriate caution. The criteria applied here are based on standard QSAR publication norms and the specific context of half-life prediction, but they are ultimately subjective and domain expert review might come to different conclusions. 

The purpose of this document is to transparently outline the rationale for selecting certain models for publication and further use, while acknowledging the limitations of the data and methods.

---

## Publishability Criteria

The bar applied here requires: good rank correlation, CV R² ≥ 0.50, no severe overfitting, and AD coverage that is quantified and acknowledged. Large AD caveats are acceptable if transparent — QSAR publications routinely include models with partial AD coverage provided limitations are explicit.

---

## Interval-Aware Pipeline (WP3bioddeg_V2)

### Selected for publication

| Model | ρ | CV R² | Test coverage | AD outside | Rationale |
|---|---|---|---|---|---|
| **SEDIMENT / VEGA** | 0.921 | 0.790 | 73.6% | 26% | Best overall. All four key metrics are defensible. The only model in the full set that would pass reviewer scrutiny without major caveats. |
| **SOIL / VEGA** | 0.921 | 0.806 | 27.8% | 72% | Strong internally (ρ and CV R² both > 0.79). Coverage and AD are real weaknesses but publishable with an explicit AD section and framing as a screening tool. |
| **AIR / HSBD** | 0.846 | 0.769 | 23.3% | 77% | ρ and CV R² are strong. The AD problem is serious but the HSBD air dataset is the only available source for this endpoint — scarcity of alternatives is itself a publishable argument. Must be presented with a hard AD filter applied to predictions. |

### Not recommended for publication as-is

| Model | Primary disqualifier |
|---|---|
| WATER / hsbd | ρ = 0.70, R² = 0.41, 79% outside AD |
| SOIL / hsbd | ρ = 0.53 — near random for rank ordering |
| WATER / vega | 49% outside AD, moderate ρ |
| WATER / combined | No advantage over single-source; 80% outside AD |
| SOIL / combined | No advantage over single-source; 80% outside AD |
| SEDIMENT / combined | No advantage over single-source |

---

## Plain Regression Pipeline (WP3biodegradability)

The bar here is harder: there is no interval-awareness to justify wide errors. Test R² and CV–test gap are the primary filters.

### Selected for publication — conditionally

| Model | CV R² | Test R² | CV–Test gap | AD inside | Condition |
|---|---|---|---|---|---|
| **Water COMBINED** | 0.528 | 0.505 | 0.02 | 78.3% | Must explicitly acknowledge RMSE scale (246 days) and present alongside a naive baseline (e.g. mean predictor). Most internally consistent model in the set. |
| **Soil COMBINED** | 0.490 | 0.408 | 0.08 | 79.3% | Borderline R² (< 0.50) with RMSE of 430 days — not a useful screening tool on its own. Only publishable framed as a negative result or benchmark baseline. |

### Not recommended for publication

| Model | Primary disqualifier |
|---|---|
| Water VEGA | Anomalous split (test R² = 0.89 vs CV R² = 0.67); learning curve R² = -4.70 |
| Air HSBD | CV–test gap = 0.41 — cross-validation is not representative of generalisation |
| Soil HSBD | Test R² = 0.13 — model fails to generalise |
| Soil VEGA | CV–test gap = 0.23; only 40.7% of test set inside AD |
| Sediment VEGA | CV–test gap = 0.23; only 43.4% inside AD; learning curve R² = -6.09 |
| Sediment HSBD | Only 34.6% of test set inside AD |

---

## Cross-Pipeline Verdict

The interval-aware pipeline produces more publishable candidates because the methodology is epistemically honest about the data — treating quantised targets as intervals rather than pretending they are precise measurements. **SEDIMENT/VEGA (interval-aware) is the single strongest candidate across both pipelines.**

The plain regression pipeline's best outputs (Water COMBINED, Soil COMBINED) would only be publishable as a methodological comparison or negative result, not as standalone predictive tools.

### Recommended publication strategy

1. **Primary results:** interval-aware SEDIMENT/VEGA and SOIL/VEGA — present as the main models.
2. **Supporting result:** interval-aware AIR/HSBD — include with a hard AD filter and explicit extrapolation warning.
3. **Baseline comparison:** plain regression Water COMBINED and Soil COMBINED — include explicitly to justify why interval-aware modelling was necessary and what is gained by the approach.
4. All other models should be reported in supplementary material at most, with clear statements that they do not meet the quality threshold for applied use.
