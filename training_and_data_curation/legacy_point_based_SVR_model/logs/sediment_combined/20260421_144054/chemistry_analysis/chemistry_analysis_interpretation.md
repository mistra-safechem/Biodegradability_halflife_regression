# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/sediment_combined/20260421_144054/chemistry_analysis/
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
,28,28.87,28.87
c1ccccc1,27,27.84,56.71
c1ccc(-c2ccccc2)cc1,8,8.25,64.96
c1ccc2ccccc2c1,4,4.12,69.08
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,3,3.09,72.17
c1ccc2c(c1)Oc1ccccc1O2,3,3.09,75.26
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (sediment_combined):**
- 97 total molecules, 28 unique scaffolds.
- Top scaffold: non-ring / acyclic (``) at 28.9%
- Top 3 scaffolds cover 65.0% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 20/28 scaffolds appear only once (71%) → long tail, limited generalisation per scaffold

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
,28,28.87,126.976,0.622,5.466
c1ccccc1,27,27.84,143.181,0.339,-14.158
c1ccc(-c2ccccc2)cc1,8,8.25,1795.856,0.211,-584.941
c1ccc2ccccc2c1,4,4.12,24.851,0.901,-12.141
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,3,3.09,74.054,-26.414,70.433
c1ccc2c(c1)Oc1ccccc1O2,3,3.09,861.339,-197742680332111.62,524.828
c1ccc(Cc2ccccc2)cc1,2,2.06,1411.396,-104952708805515.8,-1411.249
C=C(c1ccccc1)c1ccccc1,2,2.06,1566.805,0.0,-1561.528
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (sediment_combined):**

```
✅ GOOD R²: scaffold=non-ring / acyclic (28 samples, 28.9%)
   R²=0.622, RMSE=127 days, mean_error=+5.5 days

⚠️ MEDIUM R²: scaffold=benzene (27 samples, 27.8%)
   R²=0.339, RMSE=143 days, mean_error=-14.2 days

⚠️ MEDIUM R²: scaffold=biphenyl (8 samples, 8.2%)
   R²=0.211, RMSE=1796 days, mean_error=-584.9 days | systematic under-prediction (-585 days)

✅ GOOD R²: scaffold=naphthalene (4 samples, 4.1%)
   R²=0.901, RMSE=25 days, mean_error=-12.1 days

❌ BAD R²: scaffold=O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1 (3 samples, 3.1%)
   R²=-26.414, RMSE=74 days, mean_error=+70.4 days

❌ BAD R²: scaffold=diphenyl ether (fused) (3 samples, 3.1%)
   R²=-197742680332111.625, RMSE=861 days, mean_error=+524.8 days | systematic over-prediction (+525 days)

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ccc(Cc2ccccc2)cc1 (2 samples, 2.1%)
   R²=-104952708805515.797, RMSE=1411 days, mean_error=-1411.2 days | systematic under-prediction (-1411 days)

⚠️ UNRELIABLE (n≤2) R²: scaffold=C=C(c1ccccc1)c1ccccc1 (2 samples, 2.1%)
   R²=0.000, RMSE=1567 days, mean_error=-1561.5 days | systematic under-prediction (-1562 days)

```

**Actionable insights:**
- **Trust:** non-ring / acyclic, naphthalene predictions
- **Ignore:** biphenyl, O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1, diphenyl ether (fused), c1ccc(Cc2ccccc2)cc1, C=C(c1ccccc1)c1ccccc1 predictions
- **Question:** benzene predictions

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
[16*]c1ccccc1,7,221.059,333.559
[16*]c1ccc(Cl)cc1,7,1022.024,898.387
[3*]O[3*],6,205.889,218.761
[4*]CC,6,27.028,34.393
[3*]OC,5,51.667,26.245
[16*]c1ccccc1Cl,5,1245.833,974.546
[6*]C(=O)O,5,100.167,90.06
[16*]c1cc(Cl)ccc1Cl,3,2291.667,0.0
[4*]C(C)C,3,239.556,292.235
[4*]CCCC,3,180.0,0.0
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

**Findings (sediment_combined):**

```
→: [16*]c1ccccc1
   freq=7, mean=221.1 days, std=333.6

⚠️ WARNING: [16*]c1ccc(Cl)cc1
   freq=7, mean=1022.0 days, std=898.4

→: [3*]O[3*]
   freq=6, mean=205.9 days, std=218.8

✅ GOOD: [4*]CC
   freq=6, mean=27.0 days, std=34.4

✅ GOOD: [3*]OC
   freq=5, mean=51.7 days, std=26.2

⚠️ WARNING: [16*]c1ccccc1Cl
   freq=5, mean=1245.8 days, std=974.5

→: [6*]C(=O)O
   freq=5, mean=100.2 days, std=90.1

⚠️ HIGH PERSISTENCE: [16*]c1cc(Cl)ccc1Cl
   freq=3, mean=2291.7 days, std=0.0

→: [4*]C(C)C
   freq=3, mean=239.6 days, std=292.2

→: [4*]CCCC
   freq=3, mean=180.0 days, std=0.0

✅ GOOD: [4*]CCCl
   freq=3, mean=70.8 days, std=0.0

✅ GOOD: [16*]c1ccccc1[16*]
   freq=3, mean=54.9 days, std=27.7

→: [5*]N[5*]
   freq=3, mean=123.6 days, std=91.4

✅ GOOD: [16*]c1ccccc1C
   freq=3, mean=12.0 days, std=0.0

✅ GOOD: [3*]OP(=O)(O[3*])O[3*]
   freq=3, mean=19.2 days, std=17.3

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
  "EState_VSA9": {
    "importance": 0.0030480115697541334,
    "interpretation": "Molecular descriptor"
  },
  "BalabanJ": {
    "importance": 0.0031653384366877617,
    "interpretation": "Balaban connectivity index"
  },
  "fr_ArN": {
    "importance": 0.0034195499781586453,
    "interpretation": "Molecular descriptor"
  },
  "SMR_VSA4": {
    "importance": 0.0041188208959199776,
    "interpretation": "Molecular descriptor"
  },
  "fr_unbrch_alkane": {
    "importance": 0.004664099977568168,
    "interpretation": "Molecular descriptor"
  },
  "SMR_VSA1": {
    "importance": 0.005205256135885042,
    "interpretation": "Molecular descriptor"
  },
  "fr_Al_OH": {
    "importance": 0.005719879961563775,
    "interpretation": "Molecular descriptor"
  },
  "PEOE_VSA1": {
    "importance": 0.006566446659356717,
    "interpretation": "Molecular descriptor"
  },
  "SlogP_VSA8": {
    "importance": 0.008807559209917456,
    "interpretation": "Molecular descriptor"
  }
}
```

**Findings (sediment_combined):**

```
✅ GOOD: SlogP_VSA8 is top feature
   → Molecular descriptor

✅ GOOD: Functional group / count descriptors present (3 fr_*, 0 Num*)
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
c1ccc2c(c1)Oc1ccccc1O2,3,682.979,466460.327,524.828,False
c1ccccc1,27,127.086,16150.957,77.641,True
,28,212.85,45305.247,74.246,True
c1ccc(-c2ccccc2)cc1,8,591.637,350033.948,843.663,False
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,3,12.9,166.414,70.433,True
c1ccc2ccccc2c1,4,80.464,6474.474,20.314,True
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (sediment_combined):**

```
❌ BAD: diphenyl ether (fused) (3 samples)
   std=683.0 days, MAE=524.8 days, stable=False

✅ GOOD: benzene (27 samples)
   std=127.1 days, MAE=77.6 days, stable=True

⚠️ MEDIUM: non-ring / acyclic (28 samples)
   std=212.8 days, MAE=74.2 days, stable=True

❌ BAD: biphenyl (8 samples)
   std=591.6 days, MAE=843.7 days, stable=False

✅ GOOD: O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1 (3 samples)
   std=12.9 days, MAE=70.4 days, stable=True

✅ GOOD: naphthalene (4 samples)
   std=80.5 days, MAE=20.3 days, stable=True

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
,28,11,39.3
c1ccccc1,27,14,51.9
c1ccc(-c2ccccc2)cc1,8,1,12.5
c1ccc2ccccc2c1,4,0,0.0
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,3,3,100.0
c1ccc2c(c1)Oc1ccccc1O2,3,1,33.3
c1ccc(Cc2ccccc2)cc1,2,0,0.0
C=C(c1ccccc1)c1ccccc1,2,0,0.0
N=c1nc[nH]c(=N)[nH]1,1,1,100.0
c1ncncn1,1,1,100.0
c1ccc(Oc2ccccc2)cc1,1,1,100.0
O=c1ccn2c3c(cccc13)CCC2,1,1,100.0
O=C1c2ccccc2C(=O)c2ccccc21,1,1,100.0
O=c1ocnn1-c1ccccc1,1,1,100.0
c1ccc(OCn2cncn2)cc1,1,1,100.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (sediment_combined):**

```
⚠️ CONCERNING: non-ring / acyclic (28 samples)
   11/28 outside AD (39.3%)

❌ BAD: benzene (27 samples)
   14/27 outside AD (51.9%)

✅ GOOD: biphenyl (8 samples)
   1/8 outside AD (12.5%)

✅ GOOD: naphthalene (4 samples)
   0/4 outside AD (0.0%)

❌ BAD: O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1 (3 samples)
   3/3 outside AD (100.0%)

⚠️ CONCERNING: diphenyl ether (fused) (3 samples)
   1/3 outside AD (33.3%)

✅ GOOD: c1ccc(Cc2ccccc2)cc1 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: C=C(c1ccccc1)c1ccccc1 (2 samples)
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
