# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/soil_vega/20260421_142502/chemistry_analysis/
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
c1ccccc1,20,37.04,37.04
,13,24.07,61.11
c1ccc2c(c1)Oc1ccccc1O2,6,11.11,72.22
c1ccc(-c2ccccc2)cc1,5,9.26,81.48
c1ccc2c(c1)oc1ccccc12,2,3.7,85.18
C1=CCC=CC1,1,1.85,87.03
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (soil_vega):**
- 54 total molecules, 13 unique scaffolds.
- Top scaffold: benzene (`c1ccccc1`) at 37.0%
- Top 3 scaffolds cover 72.2% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 8/13 scaffolds appear only once (62%) → long tail, limited generalisation per scaffold

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
c1ccccc1,20,37.04,36.997,0.792,1.705
,13,24.07,146.035,-2.547,38.825
c1ccc2c(c1)Oc1ccccc1O2,6,11.11,656.106,0.191,-268.004
c1ccc(-c2ccccc2)cc1,5,9.26,933.325,-0.448,106.793
c1ccc2c(c1)oc1ccccc12,2,3.7,1177.129,-0.124,-836.119
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (soil_vega):**

```
✅ GOOD R²: scaffold=benzene (20 samples, 37.0%)
   R²=0.792, RMSE=37 days, mean_error=+1.7 days

❌ BAD R²: scaffold=non-ring / acyclic (13 samples, 24.1%)
   R²=-2.547, RMSE=146 days, mean_error=+38.8 days

❌ BAD R²: scaffold=diphenyl ether (fused) (6 samples, 11.1%)
   R²=0.191, RMSE=656 days, mean_error=-268.0 days | systematic under-prediction (-268 days)

❌ BAD R²: scaffold=biphenyl (5 samples, 9.3%)
   R²=-0.448, RMSE=933 days, mean_error=+106.8 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ccc2c(c1)oc1ccccc12 (2 samples, 3.7%)
   R²=-0.124, RMSE=1177 days, mean_error=-836.1 days | systematic under-prediction (-836 days)

```

**Actionable insights:**
- **Trust:** benzene predictions
- **Ignore:** non-ring / acyclic, diphenyl ether (fused), biphenyl, c1ccc2c(c1)oc1ccccc12 predictions

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
[16*]c1ccccc1,11,358.182,696.831
[3*]O[3*],6,22.917,0.0
[4*]CCCl,4,22.917,0.0
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

**Findings (soil_vega):**

```
⚠️ WARNING: [16*]c1ccccc1
   freq=11, mean=358.2 days, std=696.8

✅ GOOD: [3*]O[3*]
   freq=6, mean=22.9 days, std=0.0

✅ GOOD: [4*]CCCl
   freq=4, mean=22.9 days, std=0.0

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

_`outlier_fragments.csv` not found or empty._

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
  "MACCS_074": {
    "importance": 0.003720908833219902,
    "interpretation": "Structural fingerprint bit"
  },
  "Chi0n": {
    "importance": 0.0037547646017682622,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_122": {
    "importance": 0.003886096155735119,
    "interpretation": "Structural fingerprint bit"
  },
  "BalabanJ": {
    "importance": 0.004112089580239675,
    "interpretation": "Balaban connectivity index"
  },
  "MACCS_103": {
    "importance": 0.00421220810234395,
    "interpretation": "Structural fingerprint bit"
  },
  "Kappa3": {
    "importance": 0.004320288540674658,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_120": {
    "importance": 0.004844427028659293,
    "interpretation": "Structural fingerprint bit"
  },
  "SMR_VSA3": {
    "importance": 0.004982946697553889,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_098": {
    "importance": 0.005333266195463876,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_057": {
    "importance": 0.006240602091190603,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_129": {
    "importance": 0.006371344993237526,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_112": {
    "importance": 0.007047612740222356,
    "interpretation": "Structural fingerprint bit"
  },
  "VSA_EState2": {
    "importance": 0.007056688257565965,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_072": {
    "importance": 0.008816548859641029,
    "interpretation": "Structural fingerprint bit"
  },
  "VSA_EState7": {
    "importance": 0.009927391964474402,
    "interpretation": "Molecular descriptor"
  },
  "SlogP_VSA2": {
    "importance": 0.011274814775728633,
    "interpretation": "Molecular descriptor"
  },
  "fr_bicyclic": {
    "importance": 0.012960476514701502,
    "interpretation": "Molecular descriptor"
  }
}
```

**Findings (soil_vega):**

```
✅ GOOD: fr_bicyclic is top feature
   → Molecular descriptor

✅ GOOD: Functional group / count descriptors present (1 fr_*, 0 Num*)
⚠️ Note: 9 MACCS fingerprint bit(s) in feature set
✅ GOOD: Physicochemical/charge descriptors present (7 features)
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
c1ccccc1,20,87.788,7706.748,24.153,True
c1ccc(-c2ccccc2)cc1,5,1216.894,1480831.783,745.139,False
,13,182.055,33144.126,71.019,True
c1ccc2c(c1)Oc1ccccc1O2,6,165.966,27544.56,352.993,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (soil_vega):**

```
✅ GOOD: benzene (20 samples)
   std=87.8 days, MAE=24.2 days, stable=True

❌ BAD: biphenyl (5 samples)
   std=1216.9 days, MAE=745.1 days, stable=False

✅ GOOD: non-ring / acyclic (13 samples)
   std=182.1 days, MAE=71.0 days, stable=True

⚠️ MEDIUM: diphenyl ether (fused) (6 samples)
   std=166.0 days, MAE=353.0 days, stable=False

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
c1ccccc1,20,10,50.0
,13,12,92.3
c1ccc2c(c1)Oc1ccccc1O2,6,2,33.3
c1ccc(-c2ccccc2)cc1,5,0,0.0
c1ccc2c(c1)oc1ccccc12,2,2,100.0
C1=CCC=CC1,1,1,100.0
c1ccc(Oc2ccccc2)cc1,1,1,100.0
c1ccc2ncccc2c1,1,1,100.0
O=C(OCc1ccccc1)c1ccccc1,1,1,100.0
c1ccc(Nc2ccccc2)cc1,1,1,100.0
c1ccc(Cc2ccccc2)cc1,1,1,100.0
c1ccc2ccccc2c1,1,0,0.0
C1CCCCC1,1,0,0.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (soil_vega):**

```
⚠️ CONCERNING: benzene (20 samples)
   10/20 outside AD (50.0%)

❌ BAD: non-ring / acyclic (13 samples)
   12/13 outside AD (92.3%)

⚠️ CONCERNING: diphenyl ether (fused) (6 samples)
   2/6 outside AD (33.3%)

✅ GOOD: biphenyl (5 samples)
   0/5 outside AD (0.0%)

❌ BAD: c1ccc2c(c1)oc1ccccc12 (2 samples)
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
