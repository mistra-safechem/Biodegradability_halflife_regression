# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/soil_hsbd/20260421_141719/chemistry_analysis/
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
c1ccccc1,41,26.11,26.11
,36,22.93,49.04
c1ccncc1,5,3.18,52.22
c1ccc(-c2ccccc2)cc1,5,3.18,55.4
c1ccc2ccccc2c1,3,1.91,57.31
c1ccc2c(c1)Oc1ccccc1O2,2,1.27,58.58
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (soil_hsbd):**
- 157 total molecules, Very high scaffold diversity — 66 unique scaffolds (many singletons expected).
- Top scaffold: benzene (`c1ccccc1`) at 26.1%
- Top 3 scaffolds cover 52.2% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 55/66 scaffolds appear only once (83%) → long tail, limited generalisation per scaffold

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
c1ccccc1,41,26.11,245.249,0.529,-57.288
,36,22.93,49.127,0.12,-7.478
c1ccncc1,5,3.18,567.428,-0.262,-239.998
c1ccc(-c2ccccc2)cc1,5,3.18,17270.072,-0.107,-7731.184
c1ccc2ccccc2c1,3,1.91,66.397,-296.454,65.917
c1ccc2c(c1)Oc1ccccc1O2,2,1.27,684.069,-7.152,492.954
c1ccc(Oc2ccccc2)cc1,2,1.27,101.418,-358.352,81.469
c1ccsc1,2,1.27,51.465,-2.28,13.454
O=C1CCCC(=O)C1C(=O)c1ccccc1,2,1.27,9.577,-747.756,9.521
O=C1C=CC(=O)c2ccccc21,2,1.27,10.044,-629.525,10.04
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (soil_hsbd):**

```
⚠️ GOOD R²: scaffold=benzene (41 samples, 26.1%)
   R²=0.529, RMSE=245 days, mean_error=-57.3 days

❌ BAD R²: scaffold=non-ring / acyclic (36 samples, 22.9%)
   R²=0.120, RMSE=49 days, mean_error=-7.5 days

❌ BAD R²: scaffold=pyridine (5 samples, 3.2%)
   R²=-0.262, RMSE=567 days, mean_error=-240.0 days | systematic under-prediction (-240 days)

❌ BAD R²: scaffold=biphenyl (5 samples, 3.2%)
   R²=-0.107, RMSE=17270 days, mean_error=-7731.2 days | systematic under-prediction (-7731 days)

❌ BAD R²: scaffold=naphthalene (3 samples, 1.9%)
   R²=-296.454, RMSE=66 days, mean_error=+65.9 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=diphenyl ether (fused) (2 samples, 1.3%)
   R²=-7.152, RMSE=684 days, mean_error=+493.0 days | systematic over-prediction (+493 days)

⚠️ UNRELIABLE (n≤2) R²: scaffold=diphenyl ether (2 samples, 1.3%)
   R²=-358.352, RMSE=101 days, mean_error=+81.5 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=thiophene (2 samples, 1.3%)
   R²=-2.280, RMSE=51 days, mean_error=+13.5 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=O=C1CCCC(=O)C1C(=O)c1ccccc1 (2 samples, 1.3%)
   R²=-747.756, RMSE=10 days, mean_error=+9.5 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=O=C1C=CC(=O)c2ccccc21 (2 samples, 1.3%)
   R²=-629.525, RMSE=10 days, mean_error=+10.0 days

```

**Actionable insights:**
- **Ignore:** non-ring / acyclic, pyridine, biphenyl, naphthalene, diphenyl ether (fused) predictions
- **Question:** benzene, diphenyl ether, thiophene, O=C1CCCC(=O)C1C(=O)c1ccccc1, O=C1C=CC(=O)c2ccccc21 predictions

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
[3*]O[3*],32,138.323,339.698
[4*]CC,23,77.846,86.886
[16*]c1ccccc1,19,130.339,325.328
[3*]OC,12,72.164,85.627
[1*]C([6*])=O,9,241.444,454.473
[4*]C(C)C,6,190.097,228.494
[4*]C[8*],6,248.283,569.229
[8*]C(F)(F)F,6,68.231,122.24
[16*]c1ccc(Cl)cc1,5,278.133,219.204
[4*]CC[4*],5,14.8,6.573
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

**Findings (soil_hsbd):**

```
→: [3*]O[3*]
   freq=32, mean=138.3 days, std=339.7

✅ GOOD: [4*]CC
   freq=23, mean=77.8 days, std=86.9

→: [16*]c1ccccc1
   freq=19, mean=130.3 days, std=325.3

✅ GOOD: [3*]OC
   freq=12, mean=72.2 days, std=85.6

→: [1*]C([6*])=O
   freq=9, mean=241.4 days, std=454.5

→: [4*]C(C)C
   freq=6, mean=190.1 days, std=228.5

⚠️ WARNING: [4*]C[8*]
   freq=6, mean=248.3 days, std=569.2

→: [8*]C(F)(F)F
   freq=6, mean=68.2 days, std=122.2

→: [16*]c1ccc(Cl)cc1
   freq=5, mean=278.1 days, std=219.2

✅ GOOD: [4*]CC[4*]
   freq=5, mean=14.8 days, std=6.6

✅ GOOD: [5*]N([5*])[5*]
   freq=5, mean=47.1 days, std=24.5

→: [8*]CC
   freq=5, mean=267.4 days, std=199.2

⚠️ WARNING: [4*]CC(=O)O
   freq=4, mean=369.0 days, std=622.7

✅ GOOD: [5*]N(C)C
   freq=4, mean=18.8 days, std=13.4

→: [8*]C(C)C
   freq=4, mean=106.0 days, std=101.4

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
  "MACCS_073": {
    "importance": 0.0020102723851308955,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_Imine": {
    "importance": 0.002026287378366065,
    "interpretation": "Molecular descriptor"
  },
  "fr_Ar_N": {
    "importance": 0.0021893981219321735,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_133": {
    "importance": 0.002236103799721512,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_125": {
    "importance": 0.0022932901962830154,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_nitro_arom_nonortho": {
    "importance": 0.002296529140041603,
    "interpretation": "Molecular descriptor"
  },
  "VSA_EState6": {
    "importance": 0.0023615111630362782,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_148": {
    "importance": 0.0023892833040604854,
    "interpretation": "Structural fingerprint bit"
  },
  "NumAromaticHeterocycles": {
    "importance": 0.0026901404641125116,
    "interpretation": "Molecular descriptor"
  },
  "Kappa2": {
    "importance": 0.0029303385202844034,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_112": {
    "importance": 0.003140888645758555,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_136": {
    "importance": 0.003359945856500923,
    "interpretation": "Structural fingerprint bit"
  },
  "MinEStateIndex": {
    "importance": 0.003552368229935854,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_075": {
    "importance": 0.0035801235291259004,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_121": {
    "importance": 0.0036148609999565153,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_091": {
    "importance": 0.004269447696388251,
    "interpretation": "Structural fingerprint bit"
  },
  "SMR_VSA4": {
    "importance": 0.005278730333610426,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_028": {
    "importance": 0.005483611765014096,
    "interpretation": "Structural fingerprint bit"
  },
  "SMR_VSA6": {
    "importance": 0.006011301508062261,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_036": {
    "importance": 0.007500340839739661,
    "interpretation": "Structural fingerprint bit"
  }
}
```

**Findings (soil_hsbd):**

```
⚠️ WARNING: MACCS_036 is top feature
   → Structural fingerprint bit — unclear chemical meaning

✅ GOOD: Functional group / count descriptors present (3 fr_*, 1 Num*)
⚠️ Note: 11 MACCS fingerprint bit(s) in feature set
✅ GOOD: Physicochemical/charge descriptors present (5 features)
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
c1ccc2ccccc2c1,3,11.492,132.055,65.917,True
,36,26.119,682.219,27.963,True
c1ccccc1,41,132.267,17494.616,85.927,True
c1ccncc1,5,33.978,1154.5,291.595,True
c1ccc(-c2ccccc2)cc1,5,1120.167,1254774.75,7960.357,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (soil_hsbd):**

```
✅ GOOD: naphthalene (3 samples)
   std=11.5 days, MAE=65.9 days, stable=True

✅ GOOD: non-ring / acyclic (36 samples)
   std=26.1 days, MAE=28.0 days, stable=True

✅ GOOD: benzene (41 samples)
   std=132.3 days, MAE=85.9 days, stable=True

✅ GOOD: pyridine (5 samples)
   std=34.0 days, MAE=291.6 days, stable=True

❌ BAD: biphenyl (5 samples)
   std=1120.2 days, MAE=7960.4 days, stable=False

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
c1ccccc1,41,5,12.2
,36,6,16.7
c1ccncc1,5,2,40.0
c1ccc(-c2ccccc2)cc1,5,0,0.0
c1ccc2ccccc2c1,3,0,0.0
O=C1CCCC(=O)C1C(=O)c1ccccc1,2,2,100.0
O=C1C=CC(=O)c2ccccc21,2,2,100.0
O=S(=O)(N=CN=c1ncnc[nH]1)c1ccccc1,2,0,0.0
c1ccc2c(c1)Oc1ccccc1O2,2,0,0.0
c1ccsc1,2,0,0.0
c1ccc(Oc2ccccc2)cc1,2,0,0.0
O=c1ccnc[nH]1,1,1,100.0
c1ccc2sncc2c1,1,1,100.0
c1ccc2[nH]c(-c3cscn3)nc2c1,1,1,100.0
C1=NN=CN(N=Cc2cccnc2)C1,1,1,100.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (soil_hsbd):**

```
✅ GOOD: benzene (41 samples)
   5/41 outside AD (12.2%)

✅ GOOD: non-ring / acyclic (36 samples)
   6/36 outside AD (16.7%)

⚠️ CONCERNING: pyridine (5 samples)
   2/5 outside AD (40.0%)

✅ GOOD: biphenyl (5 samples)
   0/5 outside AD (0.0%)

✅ GOOD: naphthalene (3 samples)
   0/3 outside AD (0.0%)

❌ BAD: O=C1CCCC(=O)C1C(=O)c1ccccc1 (2 samples)
   2/2 outside AD (100.0%)

❌ BAD: O=C1C=CC(=O)c2ccccc21 (2 samples)
   2/2 outside AD (100.0%)

✅ GOOD: O=S(=O)(N=CN=c1ncnc[nH]1)c1ccccc1 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: diphenyl ether (fused) (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: thiophene (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: diphenyl ether (2 samples)
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
