# SVR Pipeline — Log Summary

Generated: 2026-04-30 12:04  
Parsed 10 run(s) from `logs/`

---

## 1. Main Performance (Final Model)

> Metrics from the final refitted model after feature reduction (`[final]` tag).
> Coverage = fraction of test compounds whose predicted log₁₀(T½) falls inside the target interval.

| Compartment | Source | N (raw) | Test coverage | MIL (log₁₀ d) | MIL (×days) | Spearman ρ | Kendall τ | Class acc. | CV R² | CV R² std | RMSE (×days) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AIR | hsbd | 308 | 17/73 (23.3%) | 0.5384 | 3.45× | 0.846 | 0.712 | 0.233 | 0.7688 ±0.0765 | 0.0765 | 2.56× |
| WATER | hsbd | 673 | 33/158 (20.9%) | 0.5143 | 3.27× | 0.699 | 0.525 | 0.209 | 0.4104 ±0.0686 | 0.0686 | 4.05× |
| WATER | vega | 223 | 27/53 (50.9%) | 0.1722 | 1.49× | 0.826 | 0.718 | 0.509 | 0.6718 ±0.0530 | 0.0530 | 2.98× |
| SOIL | hsbd | 672 | 34/157 (21.7%) | 0.6025 | 4.00× | 0.531 | 0.394 | 0.217 | 0.4882 ±0.0492 | 0.0492 | 4.16× |
| SOIL | vega | 226 | 15/54 (27.8%) | 0.5086 | 3.23× | 0.921 | 0.802 | 0.278 | 0.8063 ±0.1184 | 0.1184 | 2.19× |
| SEDIMENT | hsbd | 347 | 19/81 (23.5%) | 0.5154 | 3.28× | 0.808 | 0.648 | 0.235 | 0.5142 ±0.0674 | 0.0674 | 3.59× |
| SEDIMENT | vega | 221 | 39/53 (73.6%) | 0.0697 | 1.17× | 0.921 | 0.822 | 0.736 | 0.7903 ±0.0556 | 0.0556 | 2.12× |
| WATER | combined | 743 | 35/175 (20.0%) | 0.5297 | 3.39× | 0.664 | 0.495 | 0.200 | 0.5277 ±0.0567 | 0.0567 | 3.91× |
| SOIL | combined | 743 | 35/174 (20.1%) | 0.5720 | 3.73× | 0.709 | 0.527 | 0.201 | 0.4776 ±0.0758 | 0.0758 | 4.38× |
| SEDIMENT | combined | 414 | 25/97 (25.8%) | 0.5010 | 3.17× | 0.848 | 0.685 | 0.258 | 0.6564 ±0.0827 | 0.0827 | 3.11× |

---

## 2. Feature Reduction

> Winner strategy selected by highest (coverage, then ρ) on the test set.
> Initial = full descriptor set; Final = reduced feature set (refitted model).

| Compartment | Source | Winner strategy | N features | Init. CV R² | Final CV R² | Init. test ρ | Final test ρ | LC val R² (±std) |
|---|---|---|---|---|---|---|---|---|
| AIR | hsbd | C  top-N pca99 (N=66) | 66 | 0.7571 ±0.0885 | 0.7688 ±0.0765 | 0.851 | 0.846 | 0.6965 ±0.1379 |
| WATER | hsbd | C  top-N pca99 (N=87) | 87 | 0.4147 ±0.0847 | 0.4104 ±0.0686 | 0.604 | 0.699 | 0.2572 ±0.2134 |
| WATER | vega | C  top-N pca99 (N=56) | 56 | 0.7414 ±0.0513 | 0.6718 ±0.0530 | 0.847 | 0.826 | — |
| SOIL | hsbd | D  RFE (LinearSVR→RBF) | 147 | 0.4837 ±0.0593 | 0.4882 ±0.0492 | 0.533 | 0.531 | 0.4157 ±0.1981 |
| SOIL | vega | C  top-N pca99 (N=55) | 55 | 0.8024 ±0.1176 | 0.8063 ±0.1184 | 0.868 | 0.921 | 0.6666 ±0.1444 |
| SEDIMENT | hsbd | D  RFE (LinearSVR→RBF) | 158 | 0.5278 ±0.0646 | 0.5142 ±0.0674 | 0.821 | 0.808 | 0.2157 ±0.4628 |
| SEDIMENT | vega | A  top-SVR (importance>0) | 131 | 0.7948 ±0.0550 | 0.7903 ±0.0556 | 0.903 | 0.921 | — |
| WATER | combined | A  top-SVR (importance>0) | 114 | 0.5500 ±0.0710 | 0.5277 ±0.0567 | 0.657 | 0.664 | 0.3667 ±0.2267 |
| SOIL | combined | C  top-N pca99 (N=118) | 118 | 0.5028 ±0.0808 | 0.4776 ±0.0758 | 0.639 | 0.709 | 0.4226 ±0.1989 |
| SEDIMENT | combined | D  RFE (LinearSVR→RBF) | 132 | 0.6610 ±0.0733 | 0.6564 ±0.0827 | 0.800 | 0.848 | 0.5098 ±0.4157 |

---

## 3. Applicability Domain

> Leverage-based AD (Williams plot).  
> AD-outside = test compounds with leverage h > h★ (structurally outside the training space).  
> PCA coverage = % of test compounds within the 95% PCA ellipse of the training set.

| Compartment | Source | N test | AD-outside (n) | AD-outside (%) | PCA coverage (%) |
|---|---|---|---|---|---|
| AIR | hsbd | 73 | 56 | 76.7% | 94.5% |
| WATER | hsbd | 158 | 125 | 79.1% | 99.4% |
| WATER | vega | 53 | 26 | 49.1% | 96.2% |
| SOIL | hsbd | 157 | 123 | 78.3% | 100.0% |
| SOIL | vega | 54 | 39 | 72.2% | 96.3% |
| SEDIMENT | hsbd | 81 | 62 | 76.5% | 97.5% |
| SEDIMENT | vega | 53 | 14 | 26.4% | 98.1% |
| WATER | combined | 175 | 140 | 80.0% | 99.4% |
| SOIL | combined | 174 | 139 | 79.9% | 98.3% |
| SEDIMENT | combined | 97 | 72 | 74.2% | 100.0% |

---

## 4. Train / Test Gap

> High train ρ with low test ρ signals overfitting.

| Compartment | Source | Train coverage | Test coverage | Train ρ | Test ρ | RMSE train (×days) | RMSE test (×days) |
|---|---|---|---|---|---|---|---|
| AIR | hsbd | 50/219 (22.8%) | 14/73 (19.2%) | 0.983 | 0.851 | 2.63× | 2.63× |
| WATER | hsbd | 120/472 (25.4%) | 30/158 (19.0%) | 0.811 | 0.604 | 4.02× | 4.02× |
| WATER | vega | 123/158 (77.8%) | 33/53 (62.3%) | 0.980 | 0.847 | 2.63× | 2.63× |
| SOIL | hsbd | 104/471 (22.1%) | 32/157 (20.4%) | 0.868 | 0.533 | 4.19× | 4.19× |
| SOIL | vega | 34/160 (21.2%) | 12/54 (22.2%) | 0.982 | 0.868 | 2.20× | 2.20× |
| SEDIMENT | hsbd | 57/243 (23.5%) | 19/81 (23.5%) | 0.978 | 0.821 | 3.54× | 3.54× |
| SEDIMENT | vega | 127/157 (80.9%) | 37/53 (69.8%) | 0.968 | 0.903 | 2.11× | 2.11× |
| WATER | combined | 108/522 (20.7%) | 36/175 (20.6%) | 0.878 | 0.657 | 3.80× | 3.80× |
| SOIL | combined | 104/522 (19.9%) | 32/174 (18.4%) | 0.862 | 0.639 | 4.23× | 4.23× |
| SEDIMENT | combined | 61/290 (21.0%) | 23/97 (23.7%) | 0.981 | 0.800 | 3.09× | 3.09× |

---

## 5. Run Notes

### AIR / hsbd  `20260430_112714`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/hsbd_t_half_all.db`
- **Raw rows**: 308
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_air_hsbd_20260430_112800.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_air_hsbd_20260430_112800_ad.npz`

### WATER / hsbd  `20260430_112803`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/hsbd_t_half_all.db`
- **Raw rows**: 673
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_water_hsbd_20260430_113116.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_water_hsbd_20260430_113116_ad.npz`

### WATER / vega  `20260430_113123`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/vega_t_half_soil_water_sediment.db`
- **Raw rows**: 223
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_water_vega_20260430_113200.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_water_vega_20260430_113200_ad.npz`

### SOIL / hsbd  `20260430_113203`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/hsbd_t_half_all.db`
- **Raw rows**: 672
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_soil_hsbd_20260430_113859.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_soil_hsbd_20260430_113859_ad.npz`

### SOIL / vega  `20260430_113908`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/vega_t_half_soil_water_sediment.db`
- **Raw rows**: 226
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_soil_vega_20260430_113935.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_soil_vega_20260430_113935_ad.npz`

### SEDIMENT / hsbd  `20260430_113938`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/hsbd_t_half_all.db`
- **Raw rows**: 347
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_sediment_hsbd_20260430_114041.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_sediment_hsbd_20260430_114041_ad.npz`

### SEDIMENT / vega  `20260430_114046`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/vega_t_half_soil_water_sediment.db`
- **Raw rows**: 221
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_sediment_vega_20260430_114113.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_sediment_vega_20260430_114113_ad.npz`

### WATER / combined  `20260430_114117`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/combined_t_half_vega_hsbd_soil_water_sediment.db`
- **Raw rows**: 743
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_water_combined_20260430_114531.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_water_combined_20260430_114531_ad.npz`

### SOIL / combined  `20260430_114540`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/combined_t_half_vega_hsbd_soil_water_sediment.db`
- **Raw rows**: 743
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_soil_combined_20260430_115543.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_soil_combined_20260430_115543_ad.npz`

### SEDIMENT / combined  `20260430_115555`

- **Database**: `/Users/a/dev/WP3bioddeg_V2/processed_data/combined_t_half_vega_hsbd_soil_water_sediment.db`
- **Raw rows**: 414
- **Model**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_sediment_combined_20260430_115732.joblib`
- **AD artefact**: `/Users/a/dev/WP3bioddeg_V2/models/SVR_sediment_combined_20260430_115732_ad.npz`
