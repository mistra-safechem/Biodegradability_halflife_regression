# SVR Half-Life Models — Summary

Generated: 2026-04-30 12:05 | Source: `logs/log_summary.md`

---

## Model Performance Overview

| Compartment | Source | N (raw) | Test coverage | MIL (×days) | Spearman ρ | Kendall τ | Class acc. | CV R² | RMSE (×days) |
|---|---|---|---|---|---|---|---|---|---|
| AIR | hsbd | 308 | 17/73 (23.3%) | 3.45× | 0.846 | 0.712 | 0.233 | 0.7688 | 2.56× |
| WATER | hsbd | 673 | 33/158 (20.9%) | 3.27× | 0.699 | 0.525 | 0.209 | 0.4104 | 4.05× |
| WATER | vega | 223 | 27/53 (50.9%) | 1.49× | 0.826 | 0.718 | 0.509 | 0.6718 | 2.98× |
| SOIL | hsbd | 672 | 34/157 (21.7%) | 4.00× | 0.531 | 0.394 | 0.217 | 0.4882 | 4.16× |
| SOIL | vega | 226 | 15/54 (27.8%) | 3.23× | 0.921 | 0.802 | 0.278 | 0.8063 | 2.19× |
| SEDIMENT | hsbd | 347 | 19/81 (23.5%) | 3.28× | 0.808 | 0.648 | 0.235 | 0.5142 | 3.59× |
| SEDIMENT | vega | 221 | 39/53 (73.6%) | 1.17× | 0.921 | 0.822 | 0.736 | 0.7903 | 2.12× |
| WATER | combined | 743 | 35/175 (20.0%) | 3.39× | 0.664 | 0.495 | 0.200 | 0.5277 | 3.91× |
| SOIL | combined | 743 | 35/174 (20.1%) | 3.73× | 0.709 | 0.527 | 0.201 | 0.4776 | 4.38× |
| SEDIMENT | combined | 414 | 25/97 (25.8%) | 3.17× | 0.848 | 0.685 | 0.258 | 0.6564 | 3.11× |

---

## Key Findings

- **Scope:** 10 SVR models trained across 4 compartments (air, sediment, soil, water) and 3 data sources (combined, hsbd, vega); training sets range from 221–743 compounds.

- **Rank-order performance is strongly source-dependent.** VEGA models achieve Spearman ρ = 0.83–0.92; HSBD is notably weaker for water (ρ = 0.70), soil (ρ = 0.53) — indicating HSBD class structure or data quality limits predictability in those compartments.

- **Interval coverage is inadequate for most models.** WATER/vega (51%), SEDIMENT/vega (74%) reach useful coverage (≥50%); the remainder (AIR/hsbd (23%), WATER/hsbd (21%), SOIL/hsbd (22%), SOIL/vega (28%), SEDIMENT/hsbd (24%), WATER/combined (20%), SOIL/combined (20%), SEDIMENT/combined (26%)) have ~22% average coverage — roughly 78 in 100 test predictions fall outside the assigned persistence interval.

- **CV R² spans 0.41–0.81 across all models; 3 models are below R² = 0.5: WATER/hsbd (R² = 0.41), SOIL/hsbd (R² = 0.49), SOIL/combined (R² = 0.48). Models below this threshold should be treated as coarse classifiers rather than regressors.

- **Combining data sources does not improve performance.** All 3 combined models (water, soil, sediment) sit between their HSBD and VEGA counterparts rather than above both — dataset heterogeneity introduces conflicting label structure that the SVR cannot reconcile.

- **Overfitting is a material concern.** Train-to-test Spearman ρ drops of ≥ 0.15 observed in: WATER/hsbd (Δρ = 0.21); SOIL/hsbd (Δρ = 0.33); SEDIMENT/hsbd (Δρ = 0.16); WATER/combined (Δρ = 0.22); SOIL/combined (Δρ = 0.22); SEDIMENT/combined (Δρ = 0.18). Worst case is SOIL/hsbd (train ρ = 0.868 → test ρ = 0.533). VEGA models show smaller gaps, likely reflecting more homogeneous label distributions.

- **Applicability domain is the most critical limitation.** 8 models have ≥60% of test compounds outside the leverage threshold: AIR/hsbd (77%), WATER/hsbd (79%), SOIL/hsbd (78%), SOIL/vega (72%), SEDIMENT/hsbd (76%), WATER/combined (80%), SOIL/combined (80%), SEDIMENT/combined (74%). The majority of predictions from these models are technically extrapolations. Adequate AD coverage only in: WATER/vega (49%), SEDIMENT/vega (26%).

- **RMSE of 2.12–4.38× days is structurally expected but practically limiting.** Errors of this magnitude are consistent with quantised, interval-valued targets rather than kinetic noise; they preclude use of any of these models for point-accurate half-life estimation in a regulatory context.

- **Overall assessment:** The strongest models are SOIL/vega (ρ = 0.92, R² = 0.81); SEDIMENT/vega (ρ = 0.92, R² = 0.79); AIR/hsbd (ρ = 0.85, R² = 0.77). These are suitable for screening-level persistence classification. HSBD water and soil models, and all combined-source models, require further data curation before any regulatory or comparative use.
