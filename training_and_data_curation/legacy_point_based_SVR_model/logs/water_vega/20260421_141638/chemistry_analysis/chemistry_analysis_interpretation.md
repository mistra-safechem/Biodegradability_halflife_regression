# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/water_vega/20260421_141638/chemistry_analysis/
```
---

## Section 1: Scaffold Analysis

### `scaffold_distribution.csv`

**Purpose:** Shows the distribution of Murcko scaffolds (core ring structures) across your dataset.

**Columns:**
| Column | Description | Type |
|--------|-------------|------|
| `scaffold_smiles` | SMILES representation of the scaffold | string |
| `count` | Number of molecules containing this scaffold | int |
| `coverage_pct` | Percentage of total dataset | float (0-100) |
| `cumulative_pct` | Cumulative percentage (sorted by count) | float (0-100) |

**Example data (top rows):**
```
scaffold_smiles,count,coverage_pct,cumulative_pct
c1ccccc1,16,30.19,30.19
,13,24.53,54.72
c1ccc(-c2ccccc2)cc1,8,15.09,69.81
c1ccc2c(c1)Oc1ccccc1O2,4,7.55,77.36
c1ncncn1,2,3.77,81.13
c1ccc(Cc2ccccc2)cc1,1,1.89,83.02
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (water_vega):**
- 53 total molecules, 15 unique scaffolds.
- Top scaffold: benzene (`c1ccccc1`) at 30.2%
- Top 3 scaffolds cover 69.8% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 10/15 scaffolds appear only once (67%) → long tail, limited generalisation per scaffold

---

### `scaffold_performance.csv` ⭐ **KEY FILE**

**Purpose:** Evaluates model prediction quality for each scaffold group.

**Columns:**
| Column | Description | Type | Unit |
|--------|-------------|------|------|
| `scaffold_smiles` | SMILES of the scaffold | string | - |
| `n_samples` | Number of molecules with this scaffold | int | - |
| `coverage_pct` | Percentage of dataset | float | % |
| `rmse` | Root Mean Square Error | float | days |
| `r2` | R² coefficient of determination | float | - |
| `mean_error` | Average prediction bias | float | days |

**Data (top 10 by coverage):**
```
scaffold_smiles,n_samples,coverage_pct,rmse,r2,mean_error
c1ccccc1,16,30.19,301.521,0.698,-76.804
,13,24.53,14.696,0.66,-5.105
c1ccc(-c2ccccc2)cc1,8,15.09,252.614,0.934,45.732
c1ccc2c(c1)Oc1ccccc1O2,4,7.55,4.019,0.742,-1.61
c1ncncn1,2,3.77,489.029,-1.036,-348.947
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (water_vega):**

```
⚠️ GOOD R²: scaffold=benzene (16 samples, 30.2%)
   R²=0.698, RMSE=302 days, mean_error=-76.8 days

✅ GOOD R²: scaffold=non-ring / acyclic (13 samples, 24.5%)
   R²=0.660, RMSE=15 days, mean_error=-5.1 days

⚠️ GOOD R²: scaffold=biphenyl (8 samples, 15.1%)
   R²=0.934, RMSE=253 days, mean_error=+45.7 days

✅ GOOD R²: scaffold=diphenyl ether (fused) (4 samples, 7.5%)
   R²=0.742, RMSE=4 days, mean_error=-1.6 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ncncn1 (2 samples, 3.8%)
   R²=-1.036, RMSE=489 days, mean_error=-348.9 days | systematic under-prediction (-349 days)

```

**Actionable insights:**
- **Trust:** non-ring / acyclic, diphenyl ether (fused) predictions
- **Question:** benzene, biphenyl, c1ncncn1 predictions

---

## Section 2: Substituent Analysis

### `rgroup_analysis.csv` (if generated)

**Purpose:** Identifies R-groups (substituents) attached to dominant scaffolds and their effect on predicted activity.

**Columns:**
| Column | Description | Type |
|--------|-------------|------|
| `scaffold_smiles` | Core scaffold | string |
| `rgroup_smiles` | R-group substituent | string |
| `position` | Attachment position on scaffold | int |
| `count` | Occurrences of this R-group | int |
| `mean_activity` | Average half-life | float |
| `activity_std` | Standard deviation | float |

**Interpretation:**

| Pattern | Chemical Meaning |
|---------|------------------|
| NO₂, CN → high activity | Electron-withdrawing groups increase persistence |
| OH, NH₂ → low activity | Electron-donating groups decrease persistence |
| Cl, Br → variable | Halogen effect depends on position/context |
| High activity_std | R-group effect is context-dependent |

**Note:** `rgroup_analysis.csv` was not generated for this run (requires multiple R-groups on the same scaffold with sufficient data).

---

## Section 3: Fragment Analysis

### `fragment_frequency.csv`

**Purpose:** Lists BRICS fragments and their associated biodegradation activity.

**Columns:**
| Column | Description | Type | Unit |
|--------|-------------|------|------|
| `fragment_smiles` | SMILES of fragment (with `[*]` attachment points) | string | - |
| `frequency` | Number of occurrences | int | - |
| `mean_activity` | Mean half-life for molecules containing fragment | float | days |
| `activity_std` | Standard deviation of activity | float | days |

**Example data (fragments with frequency ≥ 3, up to 10 shown):**
```
fragment_smiles,frequency,mean_activity,activity_std
[3*]O[3*],6,13.403,10.568
[4*]CCC,5,12.458,9.745
[5*]N[5*],4,365.625,395.725
[16*]c1ccccc1Cl,4,864.583,977.836
[4*]CC,4,190.312,345.428
[4*]CCCl,3,22.917,0.0
[4*]C(C)C,3,283.333,368.061
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Frequency** | ≥3 | 1-2 | Cannot generalise from singletons |
| **activity_std** | <100 days | >500 days | Fragment too context-dependent |
| **Trend clarity** | Chemically sensible | Contradicts literature | Check for data artifacts |

**Chemical patterns to expect:**
| Fragment type | Expected activity | Reason |
|---------------|-------------------|--------|
| Halogenated aromatics | High (persistent) | C-Cl bond stability, toxicity |
| Nitro groups | High (persistent) | Electron-withdrawing, recalcitrant |
| Esters, ethers | Low (labile) | Hydrolyzable bonds |
| Carboxylic acids | Low-moderate | Readily metabolized |
| Unsubstituted aromatics | Moderate | Ring cleavage required |

**Findings (water_vega):**

```
✅ GOOD: [3*]O[3*]
   freq=6, mean=13.4 days, std=10.6

✅ GOOD: [4*]CCC
   freq=5, mean=12.5 days, std=9.7

→: [5*]N[5*]
   freq=4, mean=365.6 days, std=395.7

⚠️ WARNING: [16*]c1ccccc1Cl
   freq=4, mean=864.6 days, std=977.8

→: [4*]CC
   freq=4, mean=190.3 days, std=345.4

✅ GOOD: [4*]CCCl
   freq=3, mean=22.9 days, std=0.0

→: [4*]C(C)C
   freq=3, mean=283.3 days, std=368.1

```

---

### `outlier_fragments.csv`

**Purpose:** Identifies fragments associated with outlier predictions (large prediction errors).

**Columns:**
| Column | Description | Type |
|--------|-------------|------|
| `fragment_smiles` | SMILES of fragment | string |
| `in_high_activity` | Present in high-activity molecules? | bool |
| `in_outliers` | Present in outlier predictions? | bool |

**Data:**
```
fragment_smiles,in_high_activity,in_outliers
ClC1=C(Cl)C2(Cl)C3C(Cl)C=CC3C1(Cl)C2(Cl)Cl,False,True
```

**Interpretation:**

| Pattern | Meaning |
|---------|---------|
| `in_outliers=True, in_high_activity=False` | Fragment confuses model (not just extreme activity) |
| `in_outliers=True, in_high_activity=True` | Fragment associated with extreme but predictable activity |
| `in_outliers=False, in_high_activity=True` | Fragment well-handled by model |

**Findings (water_vega):**

```
Fragment: ClC1=C(Cl)C2(Cl)C3C(Cl)C=CC3C1(Cl)C2(Cl)Cl
  in_outliers=True, in_high_activity=False
  → Model doesn't know how to handle this chemistry → likely extrapolation

```

**Actionable insights:**
- Outlier fragments define the **boundary of model applicability**
- Highly specific fragments (polychlorinated, fused cages) → likely extrapolation
- Consider adding similar structures to training data if these are important

---

## Section 4: Consistency Checks

### `descriptor_chemistry_map.json`

**Purpose:** Maps molecular descriptors to their feature importance.

**Structure:**
```json
{
  "DESCRIPTOR_NAME": {
    "importance": 0.00928,
    "interpretation": "Structural fingerprint bit"
  }
}
```

**Descriptor types:**

| Prefix | Type | Interpretability |
|--------|------|------------------|
| `fr_*` | Functional group count | High |
| `Num*` | Count descriptor | High |
| `PEOE_*` / `SlogP_*` | Partial charge / surface area | Medium |
| `Min/MaxPartialCharge` | Partial charge extremes | Medium |
| `NumH*` | H-bond counts | Medium |
| `MACCS_*` | MACCS fingerprint bit | Low (structural pattern) |

**Example data (all features, sorted by importance):**
```json
{
  "MACCS_080": {
    "importance": 0.002721429366295039,
    "interpretation": "Structural fingerprint bit"
  },
  "PEOE_VSA1": {
    "importance": 0.004041880280607868,
    "interpretation": "Molecular descriptor"
  },
  "EState_VSA10": {
    "importance": 0.004649009216339839,
    "interpretation": "Molecular descriptor"
  },
  "PEOE_VSA4": {
    "importance": 0.004843697309485734,
    "interpretation": "Molecular descriptor"
  },
  "FractionCSP3": {
    "importance": 0.005143310607856827,
    "interpretation": "3D character (sp3 fraction)"
  },
  "PEOE_VSA6": {
    "importance": 0.005188321553258204,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_128": {
    "importance": 0.005313027558287292,
    "interpretation": "Structural fingerprint bit"
  },
  "RingCount": {
    "importance": 0.008928840347649228,
    "interpretation": "Total ring count"
  }
}
```

**Findings (water_vega):**

```
✅ GOOD: RingCount is top feature
   → Total ring count

❌ BAD: Only MACCS fingerprint bits (2) — no interpretable functional group features
✅ GOOD: Physicochemical/charge descriptors present (6 features)
```

---

### `motif_stability.csv`

**Purpose:** Assesses prediction consistency within scaffold groups.

**Columns:**
| Column | Description | Type | Unit |
|--------|-------------|------|------|
| `scaffold_smiles` | Core scaffold | string | - |
| `n_samples` | Number of samples | int | - |
| `prediction_std` | Standard deviation of predictions | float | days |
| `prediction_variance` | Variance of predictions | float | days² |
| `mean_absolute_error` | MAE for this scaffold | float | days |
| `stable` | Stability flag (std < threshold) | bool | - |

**Data:**
```
scaffold_smiles,n_samples,prediction_std,prediction_variance,mean_absolute_error,stable
c1ccccc1,16,258.549,66847.628,86.009,False
,13,15.532,241.24,10.523,True
c1ccc2c(c1)Oc1ccccc1O2,4,5.597,31.33,3.854,True
c1ccc(-c2ccccc2)cc1,8,985.724,971651.566,184.949,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (water_vega):**

```
⚠️ MEDIUM: benzene (16 samples)
   std=258.5 days, MAE=86.0 days, stable=False

✅ GOOD: non-ring / acyclic (13 samples)
   std=15.5 days, MAE=10.5 days, stable=True

✅ GOOD: diphenyl ether (fused) (4 samples)
   std=5.6 days, MAE=3.9 days, stable=True

❌ BAD: biphenyl (8 samples)
   std=985.7 days, MAE=184.9 days, stable=False

```

**Key insight:** Compare with `scaffold_performance.csv`:
- **Both good:** Reliable scaffold
- **Good R², bad stability:** Model fits but predictions are sensitive
- **Both bad:** Avoid this chemistry

---

### `ad_scaffold_analysis.csv` ⭐ **CRITICAL FILE**

**Purpose:** Shows which scaffolds fall outside the model's Applicability Domain (AD).

**Columns:**
| Column | Description | Type | Unit |
|--------|-------------|------|------|
| `scaffold_smiles` | Core scaffold | string | - |
| `n_samples` | Total samples | int | - |
| `n_outside_ad` | Count outside AD | int | - |
| `pct_outside_ad` | Percentage outside AD | float | % |

**Data (sorted by n_samples):**
```
scaffold_smiles,n_samples,n_outside_ad,pct_outside_ad
c1ccccc1,16,6,37.5
,13,5,38.5
c1ccc(-c2ccccc2)cc1,8,0,0.0
c1ccc2c(c1)Oc1ccccc1O2,4,1,25.0
c1ncncn1,2,2,100.0
c1ccc(Cc2ccccc2)cc1,1,1,100.0
C1=CC2C3C=CC(C3)C2C1,1,1,100.0
c1ccc2ncccc2c1,1,1,100.0
c1ccc2c(c1)ccc1c3ccccc3ccc21,1,0,0.0
c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34,1,0,0.0
c1ccc2ccccc2c1,1,0,0.0
c1ccc2cc3c(ccc4ccccc43)cc2c1,1,0,0.0
c1ccc2c(c1)ccc1ccccc12,1,0,0.0
C1CCCCC1,1,0,0.0
c1ccc2c(c1)CCO2,1,0,0.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (water_vega):**

```
⚠️ CONCERNING: benzene (16 samples)
   6/16 outside AD (37.5%)

⚠️ CONCERNING: non-ring / acyclic (13 samples)
   5/13 outside AD (38.5%)

✅ GOOD: biphenyl (8 samples)
   0/8 outside AD (0.0%)

⚠️ CONCERNING: diphenyl ether (fused) (4 samples)
   1/4 outside AD (25.0%)

❌ BAD: c1ncncn1 (2 samples)
   2/2 outside AD (100.0%)

```

**Critical insight:**
```
High R² + High % outside AD = Model is LUCKY, not RELIABLE

Example pattern:
  R²>0.8 but 100% outside AD
  → Model may have memorized small sample, not learned generalizable patterns
  → Do NOT trust predictions for new molecules with this scaffold
```

---

## Visualizations

### `scaffold_bar_plot.png`

**What it shows:** Bar chart of top 10 scaffolds by count

**Use:** Quick visual of dataset composition

### `fragment_activity_plot.png`

**What it shows:** Scatter plot of fragment frequency vs. mean activity

**Use:** Identify fragments that correlate with high/low persistence

### `substituent_trends.png` (if generated)

**What it shows:** R-group effect on activity for dominant scaffolds

**Use:** Visualize substituent electronic/steric effects

---
