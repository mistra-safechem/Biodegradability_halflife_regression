# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/water_hsbd/20260421_141308/chemistry_analysis/
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
,61,38.61,38.61
c1ccccc1,38,24.05,62.66
C1CCCCC1,6,3.8,66.46
c1ccc2c(c1)Oc1ccccc1O2,5,3.16,69.62
c1ccc(-c2ccccc2)cc1,4,2.53,72.15
c1ccc2c(c1)CCCC2,2,1.27,73.42
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (water_hsbd):**
- 158 total molecules, 45 unique scaffolds.
- Top scaffold: non-ring / acyclic (``) at 38.6%
- Top 3 scaffolds cover 66.5% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 36/45 scaffolds appear only once (80%) → long tail, limited generalisation per scaffold

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
,61,38.61,31.452,0.558,-4.329
c1ccccc1,38,24.05,186.44,-0.057,-41.528
C1CCCCC1,6,3.8,8.44,-0.011,-4.417
c1ccc2c(c1)Oc1ccccc1O2,5,3.16,75.79,0.148,-37.062
c1ccc(-c2ccccc2)cc1,4,2.53,580.261,0.611,-312.022
c1ccc2c(c1)CCCC2,2,1.27,1.412,0.578,1.412
c1ccc(Oc2ccccc2)cc1,2,1.27,29.11,-160.353,11.213
c1ccc2ccccc2c1,2,1.27,2.731,-1.717,-0.558
O=C(c1ccccc1)c1ccccc1,2,1.27,4.395,-4.908,-2.935
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (water_hsbd):**

```
✅ GOOD R²: scaffold=non-ring / acyclic (61 samples, 38.6%)
   R²=0.558, RMSE=31 days, mean_error=-4.3 days

❌ BAD R²: scaffold=benzene (38 samples, 24.1%)
   R²=-0.057, RMSE=186 days, mean_error=-41.5 days

❌ BAD R²: scaffold=cyclohexane (6 samples, 3.8%)
   R²=-0.011, RMSE=8 days, mean_error=-4.4 days

❌ BAD R²: scaffold=diphenyl ether (fused) (5 samples, 3.2%)
   R²=0.148, RMSE=76 days, mean_error=-37.1 days

✅ GOOD R²: scaffold=biphenyl (4 samples, 2.5%)
   R²=0.611, RMSE=580 days, mean_error=-312.0 days | systematic under-prediction (-312 days)

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ccc2c(c1)CCCC2 (2 samples, 1.3%)
   R²=0.578, RMSE=1 days, mean_error=+1.4 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=diphenyl ether (2 samples, 1.3%)
   R²=-160.353, RMSE=29 days, mean_error=+11.2 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=naphthalene (2 samples, 1.3%)
   R²=-1.717, RMSE=3 days, mean_error=-0.6 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=O=C(c1ccccc1)c1ccccc1 (2 samples, 1.3%)
   R²=-4.908, RMSE=4 days, mean_error=-2.9 days

```

**Actionable insights:**
- **Trust:** non-ring / acyclic predictions
- **Ignore:** benzene, cyclohexane, diphenyl ether (fused), biphenyl predictions
- **Question:** c1ccc2c(c1)CCCC2, diphenyl ether, naphthalene, O=C(c1ccccc1)c1ccccc1 predictions

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
[3*]O[3*],21,9.593,11.764
[16*]c1ccccc1,18,35.821,54.308
[4*]CC,9,35.014,60.684
[3*]OC,8,16.503,22.788
[4*]CCO,7,4.488,5.998
[4*]CC[4*],6,16.07,19.527
[1*]C([6*])=O,5,4.088,3.648
[8*]CC,5,222.647,271.03
[16*]c1ccc([16*])cc1,4,12.142,5.841
[5*]N[5*],3,13.328,16.464
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

**Findings (water_hsbd):**

```
✅ GOOD: [3*]O[3*]
   freq=21, mean=9.6 days, std=11.8

✅ GOOD: [16*]c1ccccc1
   freq=18, mean=35.8 days, std=54.3

✅ GOOD: [4*]CC
   freq=9, mean=35.0 days, std=60.7

✅ GOOD: [3*]OC
   freq=8, mean=16.5 days, std=22.8

✅ GOOD: [4*]CCO
   freq=7, mean=4.5 days, std=6.0

✅ GOOD: [4*]CC[4*]
   freq=6, mean=16.1 days, std=19.5

✅ GOOD: [1*]C([6*])=O
   freq=5, mean=4.1 days, std=3.6

→: [8*]CC
   freq=5, mean=222.6 days, std=271.0

✅ GOOD: [16*]c1ccc([16*])cc1
   freq=4, mean=12.1 days, std=5.8

✅ GOOD: [5*]N[5*]
   freq=3, mean=13.3 days, std=16.5

✅ GOOD: [8*]C(C)C
   freq=3, mean=15.3 days, std=10.6

✅ GOOD: [8*]C(F)(F)F
   freq=3, mean=41.0 days, std=43.9

✅ GOOD: [3*]OP(=S)(OC)OC
   freq=3, mean=27.1 days, std=17.1

✅ GOOD: [4*]CC(=O)O
   freq=3, mean=24.0 days, std=0.0

→: [16*]c1ccc(Cl)cc1
   freq=3, mean=93.0 days, std=118.1

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
  "MACCS_158": {
    "importance": 0.002547023922856212,
    "interpretation": "Structural fingerprint bit"
  },
  "PEOE_VSA9": {
    "importance": 0.0025651940703593835,
    "interpretation": "Molecular descriptor"
  },
  "EState_VSA6": {
    "importance": 0.0029424053523162685,
    "interpretation": "Molecular descriptor"
  },
  "VSA_EState9": {
    "importance": 0.002948508971834252,
    "interpretation": "Molecular descriptor"
  },
  "VSA_EState8": {
    "importance": 0.0029527515152261152,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_080": {
    "importance": 0.0029644887315885492,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_078": {
    "importance": 0.0030343496500283675,
    "interpretation": "Structural fingerprint bit"
  },
  "BCUT2D_LOGPHI": {
    "importance": 0.0031309383772040677,
    "interpretation": "Molecular descriptor"
  },
  "fr_thiophene": {
    "importance": 0.0032046142668832266,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_062": {
    "importance": 0.0033409214239930375,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_ketone": {
    "importance": 0.0034852467379372864,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_092": {
    "importance": 0.0035368190341355022,
    "interpretation": "Structural fingerprint bit"
  },
  "Chi3v": {
    "importance": 0.003609480122263936,
    "interpretation": "Molecular descriptor"
  },
  "NHOHCount": {
    "importance": 0.003842715372531913,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_110": {
    "importance": 0.004375072608756886,
    "interpretation": "Structural fingerprint bit"
  },
  "SMR_VSA1": {
    "importance": 0.004986609714443868,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_048": {
    "importance": 0.0050210971073886505,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_NH1": {
    "importance": 0.006653096476783704,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_029": {
    "importance": 0.008471026593894293,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_081": {
    "importance": 0.017056414117193996,
    "interpretation": "Structural fingerprint bit"
  }
}
```

**Findings (water_hsbd):**

```
⚠️ WARNING: MACCS_081 is top feature
   → Structural fingerprint bit — unclear chemical meaning

✅ GOOD: Functional group / count descriptors present (3 fr_*, 0 Num*)
⚠️ Note: 9 MACCS fingerprint bit(s) in feature set
✅ GOOD: Physicochemical/charge descriptors present (8 features)
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
c1ccccc1,38,48.933,2394.404,58.07,False
C1CCCCC1,6,2.466,6.082,5.815,True
,61,24.479,599.22,12.189,True
c1ccc2c(c1)Oc1ccccc1O2,5,19.374,375.338,39.722,True
c1ccc(-c2ccccc2)cc1,4,458.386,210117.286,349.218,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (water_hsbd):**

```
⚠️ MEDIUM: benzene (38 samples)
   std=48.9 days, MAE=58.1 days, stable=False

✅ GOOD: cyclohexane (6 samples)
   std=2.5 days, MAE=5.8 days, stable=True

✅ GOOD: non-ring / acyclic (61 samples)
   std=24.5 days, MAE=12.2 days, stable=True

✅ GOOD: diphenyl ether (fused) (5 samples)
   std=19.4 days, MAE=39.7 days, stable=True

⚠️ MEDIUM: biphenyl (4 samples)
   std=458.4 days, MAE=349.2 days, stable=False

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
,61,16,26.2
c1ccccc1,38,12,31.6
C1CCCCC1,6,1,16.7
c1ccc2c(c1)Oc1ccccc1O2,5,1,20.0
c1ccc(-c2ccccc2)cc1,4,0,0.0
O=C(c1ccccc1)c1ccccc1,2,1,50.0
c1ccc(Oc2ccccc2)cc1,2,1,50.0
c1ccc2c(c1)CCCC2,2,0,0.0
c1ccc2ccccc2c1,2,0,0.0
c1cncnc1,1,1,100.0
C(=Nc1ccccc1)NC1CCCCC1,1,1,100.0
O=c1[nH]c(=O)c2[nH]cnc2[nH]1,1,1,100.0
O=C(N=Cc1ccccc1)Nc1ccccc1,1,1,100.0
C1CCCCCCCCCCC1,1,1,100.0
C1C[C@H]2C[C@@H]3CCC[C@]3(C1)C2,1,1,100.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (water_hsbd):**

```
⚠️ CONCERNING: non-ring / acyclic (61 samples)
   16/61 outside AD (26.2%)

⚠️ CONCERNING: benzene (38 samples)
   12/38 outside AD (31.6%)

✅ GOOD: cyclohexane (6 samples)
   1/6 outside AD (16.7%)

⚠️ CONCERNING: diphenyl ether (fused) (5 samples)
   1/5 outside AD (20.0%)

✅ GOOD: biphenyl (4 samples)
   0/4 outside AD (0.0%)

⚠️ CONCERNING: O=C(c1ccccc1)c1ccccc1 (2 samples)
   1/2 outside AD (50.0%)

⚠️ CONCERNING: diphenyl ether (2 samples)
   1/2 outside AD (50.0%)

✅ GOOD: c1ccc2c(c1)CCCC2 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: naphthalene (2 samples)
   0/2 outside AD (0.0%)

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
