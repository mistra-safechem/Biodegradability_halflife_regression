# METHODS.md — Interval-Aware SVR Pipeline for Biodegradability Half-Life Prediction

## Overview

This pipeline predicts environmental persistence categories (half-life) for
chemical compounds across four compartments: air, soil, water, sediment.
Targets are treated as **epistemically uncertain intervals**, not precise
kinetic measurements. Reported half-life values in the source databases
are quantised class representatives, not individual compound measurements.
The pipeline reflects this by mapping each target value to a log10 interval
and evaluating predictions accordingly.

---

## Data Sources

Three SQLite databases in `processed_data/`:

| File | Compartments |
|---|---|
| `hsbd_t_half_all.db` | air, soil, water, sediment |
| `vega_t_half_soil_water_sediment.db` | soil, water, sediment |
| `combined_t_half_vega_hsbd_soil_water_sediment.db` | soil, water, sediment |

Each table (`AirData`, `SoilData`, `WaterData`, `SedimentData`) contains:
- `T_half_days` — raw half-life class value (days)
- `T_half_class_days`, `T_half_lower_bound`, `T_half_upper_bound` — class metadata
- `T_half_log10_lower`, `T_half_log10_upper` — pre-computed log10 interval bounds
- `Canonical_smiles` — RDKit-canonical SMILES string
- RDKit physicochemical descriptors (~200 columns)
- MACCS fingerprint bits (MACCS_001 … MACCS_166)

Data is retrievable via 
- `src/db_utils.py:get_basic_data()` = only table reference, T_half_days, Canonical_smiles
- `src/db_utils.py:get_all_data()` = as basic_data() plus all columns, including all T_half & metadata, descriptors and MACCS bits
- `src/db_utils.py:get_selected_data()` = as basic_data() plus only T_half_days, descriptors and MACSS bits (no T_half metadata columns)

`get_selected_data()` is mainly used for backwards compatibility with the legacy point-based pipeline, which does not use the T_half metadata columns. 

The interval-aware pipeline uses `get_all_data()` to access the T_half metadata columns for interval construction and evaluation.

---
## Endpoint Interpretation

Half-life values in the databases are **class representatives** of persistence
categories (e.g. 1.4 d, 5.6 d, 14 d, …). Each unique class value C is
treated as the centroid of a half-order-of-magnitude interval in log10 space:

```
lower = log10(C) - half_step
upper = log10(C) + half_step
half_step = median(gaps between adjacent log10 class values) / 2
```

This mapping is computed in `src/ml_tools.py:build_log10_intervals()` and
exported to `logs/<compartment>_<ds>/<ts>/interval_table_<compartment>.csv`
for reproducibility.

**DB-stored bounds are preferred.** When `T_half_log10_lower` / `T_half_log10_upper`
are present and non-null in the database, those values are used directly.
`build_log10_intervals()` is only invoked as a fallback.

---

## Feature Engineering

Features for each compound are:

- **RDKit descriptors** (~200 physicochemical properties, minus selected descriptors) — computed by RDKit,
  stored in the DB, retrieved as columns.
- **MACCS fingerprints** (166-bit binary keys) — computed by RDKit, stored as
  `MACCS_001` … `MACCS_166`.

For details on the specific descriptors used and the feature selection process, see [`src/feature_selection_in_rdkit_tools.py.md`](src/feature_selection_in_rdkit_tools.py.md).

### Columns excluded from the feature matrix

The following column groups are **never** used as predictors:

| Group | Columns |
|---|---|
| Endpoint-derived | `T_half_class_days`, `T_half_lower_bound`, `T_half_upper_bound`, `T_half_log10_lower`, `T_half_log10_upper` |
| Row identifiers | `id`, `Canonical_smiles`, `reference` |

These are stashed in `Preprocessor.t_half_meta` before any dropping and
exported to `t_half_meta_<compartment>.csv` for traceability.

---

## Preprocessing Pipeline (`src/ml_tools.py:Preprocessor`)

Executed in order:

1. **Stash T_half metadata and SMILES** before any column removal.
2. **Drop irrelevant columns** (`drop_irrelevant_columns()`).
3. **Replace ±inf with NaN** — RDKit emits `inf` for degenerate molecular geometries.
4. **Drop rows** with any NaN in features or target; apply the same mask to
   metadata and SMILES series.
5. **Build interval table** from `T_half_days`; resolve per-sample bounds from
   DB or fallback.
6. **Outlier removal** via `IsolationForest` (contamination=0.05) on the
   feature matrix.
7. **Scale features** (`scale_features()`): replace ±inf → NaN, drop all-NaN
   columns, impute remaining NaNs with column median, apply `StandardScaler`
   to RDKit descriptor columns only. MACCS bits are left unscaled.
8. **Remove zero-variance and highly correlated features**
   (`remove_variance_and_correlation()`): Pearson |r| threshold = 0.95.

---

## Train / Test Split

- 75 % train / 25 % test, `random_state=42`.
- Any features that become constant within one split are dropped after splitting.
- SMILES strings are tracked per split and exported to
  `training_split_<compartment>.csv` (answers "was this compound in training?").


Regarding splitting: 75/25 is a standard ratio. Earlier iterations considered different splits for different datasets,
and while some different ratios gave somewhat better results, it ultimately did not impact the overall conclusions and 
was not worth the added complexity. 
The random state is fixed for reproducibility, and the same split is used across all feature reduction strategies to ensure a fair comparison.

---

## Model

**Support Vector Regression (SVR)** with RBF or linear kernel.

### Hyperparameter Search (`svr_grid_search()`)

Grid search over:

| Parameter | Values |
|---|---|
| `C` | 0.1, 1, 10, 100 |
| `epsilon` | 0.01, 0.1, 0.2, 0.5 |
| `kernel` | rbf, linear |
| `gamma` | scale, auto |

5-fold CV, scored on `neg_mean_squared_error` against the log10 class centroid.
This is a practical proxy; final evaluation uses interval metrics (below).

### Training target

`y = log10(T_half_days)` — the log10 class centroid, used as a point proxy
during grid search and CV. Interval bounds (`y_lower`, `y_upper`) are used
only in evaluation and feature reduction selection.

---

## Feature Reduction

After an initial full-feature fit, four reduction strategies are compared:

| Strategy | Description |
|---|---|
| A | Features with permutation importance > 0 |
| B | Top-N by importance (N = count of strategy A features) |
| C | Top-N by PCA explained variance ≥ 99 % |
| D | RFE via LinearSVR → retrain with RBF SVR |

Each strategy is evaluated on the test set with interval metrics (below).
**Winner** = highest interval coverage probability; Spearman ρ as tie-break.
The final model is retrained on the winning feature set.

---

## Evaluation Metrics

All metrics operate in **log10 space**; human-readable equivalents in days
are reported alongside every log10 figure.

### Primary metrics (interval-aware)

| Metric | Definition |
|---|---|
| **Interval coverage probability** | Fraction of test predictions falling inside [y_lower, y_upper] |
| **Mean interval loss (log10)** | Mean distance outside bounds; 0 when inside |
| **Mean interval loss (days, geom.)** | 10^(mean_interval_loss_log10) — geometric-mean day equivalent |
| **Class accuracy** | Same as coverage probability; named per OECD terminology |

### Rank-order metrics

| Metric | Definition |
|---|---|
| **Spearman ρ** | Rank correlation between log10 class centroids and predictions |
| **Kendall τ** | Concordance-based rank correlation |

### CV reference metric

5-fold CV logs `R²` and `RMSE(log10)` alongside `RMSE(days-geom) = 10^RMSE(log10)`.
These are **reference diagnostics only** — R² and RMSE on log10 values do not
carry kinetic interpretation and are not used to select models.

### Explicitly absent

RMSE on raw half-life values and R² interpreted as kinetic fidelity are
**not reported**. The data are persistence categories, not individual kinetic
measurements; such metrics would falsely imply compound-specific accuracy.

---

## Applicability Domain

### Leverage-based AD (Williams plot)

Computed in `applicability_domain_leverage()`:

1. Standardise train and test matrices with a `StandardScaler` fitted on
   X_train (identical scaler parameters stored in the AD artefact).
2. Compute the hat-matrix inverse: `(XᵀX)⁻¹` via pseudoinverse.
3. Per-sample leverage: `h_i = xᵢᵀ (XᵀX)⁻¹ xᵢ`.
4. AD boundary: `h* = 3p / n` (p = features, n = training samples).
5. `h ≤ h*` → compound is structurally within the training AD.

The Williams plot shows leverage (x-axis) vs. signed distance from the
nearest interval bound (y-axis; 0 when the prediction is inside the interval),
with a secondary y-axis in days.

### Chemical-space coverage (PCA)

2-component PCA on standardised features, fit on X_train, projected for
X_test. Coverage percentage = fraction of test points within the
axis-aligned bounding box of the training PCA scores.

### Morgan fingerprint PCA + Butina clustering

Morgan fingerprints (radius=2, 1024 bits) provide a structure-based view
independent of the descriptor feature set. AD-outside test compounds are
clustered with Butina (Tanimoto distance cutoff = 0.5) to identify whether
failures are structurally diverse or concentrated.

---

## Saved Artefacts

All paths to dependent files are recorded in the model JSON card so inference requires only
the JSON file.

### `models/` directory

| File | Contents |
|---|---|
| `SVR_<c>_<ds>_<ts>.joblib` | Fitted SVR model |
| `SVR_<c>_<ds>_<ts>.json` | Model card (params, metrics, all artefact paths) |
| `SVR_<c>_<ds>_<ts>_ad.npz` | AD artefact: `X_train`, `X_train_scaled`, `XtX_inv`, `X_train_mean`, `X_train_std`, `h_star`, `feature_cols` |
| `t_half_meta_<c>_<ds>_<ts>.csv` | Stashed T_half auxiliary columns, row-aligned with X/y |
| `training_split_<c>_<ds>_<ts>.csv` | SMILES + split label + y_log10 + y_lower + y_upper |


### `logs/<c>_<ds>/<ts>/` directory

| File | Contents |
|---|---|
| `SVR_interval_model.log` | Full run log |
| `interval_table_<c>.csv` | Class value → log10 bound mapping |
| `t_half_meta_<c>.csv` | Stashed T_half auxiliary columns, row-aligned with X/y |
| `training_split_<c>.csv` | SMILES + split label + y_log10 + y_lower + y_upper |
| `ad_membership_<c>.csv` | Per-compound leverage, h_star, AD flag, interval membership |
| `learning_curve_<c>.png` | Training vs CV R² as a function of training set size |
| `true_vs_pred_<c>.png` | True vs predicted — 2 panels: log10 and days |
| `residuals/qqplot_<c>.png` | Q-Q (log10) + day-scale histogram |
| `residuals/residual_hist_<c>.png` | 2-panel residual histogram |
| `ad_analysis/chemical_space_pca_<c>.png` | PCA coverage scatter |
| `ad_analysis/williams_plot_<c>.png` | Williams plot |
| `ad_analysis/morgan_pca_ad_<c>.png` | Morgan FP PCA coloured by AD membership |
| `ad_analysis/butina_cluster_sizes_<c>.png` | Butina cluster size bar chart |
| `ad_analysis/butina_cluster_pca_<c>.png` | Butina clusters overlaid on Morgan FP PCA |
|- `.earlier_runs/` | Previous iterations' logs, artefacts, and plots (not overwritten by new runs) |

---

## Inference (`SVR_interval_inference.py`)

Accepts a SMILES string (or file of SMILES) and a model JSON card.

Steps:

1. Load SVR (`.joblib`) and AD artefact (`.npz`) from paths in the JSON card.
2. Compute RDKit descriptors and MACCS bits for the input SMILES.
3. Select and order features to match `feature_columns` from the JSON card.
4. Replace any inf/NaN values with training-set column medians (from the
   stored `X_train` matrix in the AD artefact).
5. Scale using a `StandardScaler` refit on the stored `X_train` — guarantees
   exact reproducibility without a separate scaler file.
6. Predict log10 half-life centroid; convert to days.
7. Compute leverage `h = xᵀ (XᵀX)⁻¹ x`; compare to stored `h*`.
8. Check exact SMILES membership in `training_split_<c>.csv`.
9. Print human-readable report and/or write JSON output.

```
python SVR_interval_inference.py --smiles "CCCC" --model models/SVR_air_hsbd_<ts>.json
python SVR_interval_inference.py --smiles-file compounds.txt --model models/SVR_air_hsbd_<ts>.json
```

An example script for testing (or expanding into batch inference) is provided in `SVR_interval_inference_test.sh`.

---

## Entry Point

```
python SVR_interval_model_and_analysis.py --compartment <c> --data-source <ds>
```

| Argument | Choices |
|---|---|
| `--compartment` | air, soil, water, sediment |
| `--data-source` | hsbd, vega, combined |

Air is only available in `hsbd`. `vega` and `combined` cover soil, water,
sediment only.

---

## Code Layout

```
./
├── SVR_interval_model_and_analysis.py Model entry point and orchestration
├── SVR_interval_run_all_combos.sh     Script automates running all model compartment/data-source combos
├── SVR_interval_inference.py          Inference script
├── SVR_interval_inference_test.sh     Example script for testing or batch inference
├── METHODS.md                         This file
├── README.md                          Project overview and instructions
├── src/
│   ├── ml_tools.py             Interval construction, loss, evaluation,
│   │                           preprocessing, Preprocessor class, AD functions
│   ├── log_utils.py            Python logging wrapper (get_logger, _SectionLogger)
│   ├── db_schema.py            SQLAlchemy ORM models (unchanged)
│   ├── db_utils.py             DB query helpers (unchanged)
│   └── rdkit_tools.py          DESCRIPTOR_NAMES, MACCS_NAMES constants (unchanged)
├── models/                     Trained models and AD artefacts
├── logs/                       Per-run logs, CSVs, plots
└── processed_data/             SQLite databases
```

---
## Legacy model and modules
Model and inference from the point-based SVR (as well as early optimization attempts) scripts are found in `legacy_point_based_SVR_model/` directory, see also `legacy_point_based_SVR_model/readme.md` there for details.
The `src/legacy/` directory contains modules from the previous point-based API version. These are retained for reference and potential reuse but are not part of the current interval-aware pipeline. 

Layout in brief:
- `chemistry_analysis.py` — chemical structure analysis functions (e.g. Murcko scaffolds, BRICS decomposition) as separate utilities, not integrated into the main pipeline as in newer versions.
- `log_utils.py` — different logging utilities
- `ml_tools.py` — the original evaluation functions based on point predictions

