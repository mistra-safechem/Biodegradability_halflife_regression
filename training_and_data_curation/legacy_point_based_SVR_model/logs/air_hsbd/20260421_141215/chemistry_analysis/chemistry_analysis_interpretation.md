# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/air_hsbd/20260421_141215/chemistry_analysis/
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
c1ccccc1,29,39.73,39.73
,21,28.77,68.5
c1ccc(-c2ccccc2)cc1,6,8.22,76.72
c1ccc2c(c1)Oc1ccccc1O2,5,6.85,83.57
c1ccc2ccccc2c1,3,4.11,87.68
c1ccncc1,2,2.74,90.42
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (air_hsbd):**
- 73 total molecules, 13 unique scaffolds.
- Top scaffold: benzene (`c1ccccc1`) at 39.7%
- Top 3 scaffolds cover 76.7% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 7/13 scaffolds appear only once (54%) → long tail, limited generalisation per scaffold

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
c1ccccc1,29,39.73,4.486,0.903,-0.602
,21,28.77,175.221,0.272,-60.071
c1ccc(-c2ccccc2)cc1,6,8.22,737.435,0.198,-358.508
c1ccc2c(c1)Oc1ccccc1O2,5,6.85,2.385,0.0,0.984
c1ccc2ccccc2c1,3,4.11,1.265,-1.874,-0.177
c1ccncc1,2,2.74,14.659,0.0,-14.657
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (air_hsbd):**

```
✅ GOOD R²: scaffold=benzene (29 samples, 39.7%)
   R²=0.903, RMSE=4 days, mean_error=-0.6 days

⚠️ MEDIUM R²: scaffold=non-ring / acyclic (21 samples, 28.8%)
   R²=0.272, RMSE=175 days, mean_error=-60.1 days

❌ BAD R²: scaffold=biphenyl (6 samples, 8.2%)
   R²=0.198, RMSE=737 days, mean_error=-358.5 days | systematic under-prediction (-359 days)

❌ BAD R²: scaffold=diphenyl ether (fused) (5 samples, 6.8%)
   R²=0.000, RMSE=2 days, mean_error=+1.0 days

❌ BAD R²: scaffold=naphthalene (3 samples, 4.1%)
   R²=-1.874, RMSE=1 days, mean_error=-0.2 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=pyridine (2 samples, 2.7%)
   R²=0.000, RMSE=15 days, mean_error=-14.7 days

```

**Actionable insights:**
- **Trust:** benzene predictions
- **Ignore:** biphenyl, diphenyl ether (fused), naphthalene predictions
- **Question:** non-ring / acyclic, pyridine predictions

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
[16*]c1ccccc1,8,6.99,9.864
[6*]C(=O)O,4,2.292,0.0
[3*]O[3*],4,1.375,1.078
[16*]c1cc(Cl)c(Cl)cc1Cl,3,153.542,130.986
[4*]C(C)C,3,4.792,3.969
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

**Findings (air_hsbd):**

```
✅ GOOD: [16*]c1ccccc1
   freq=8, mean=7.0 days, std=9.9

✅ GOOD: [6*]C(=O)O
   freq=4, mean=2.3 days, std=0.0

✅ GOOD: [3*]O[3*]
   freq=4, mean=1.4 days, std=1.1

→: [16*]c1cc(Cl)c(Cl)cc1Cl
   freq=3, mean=153.5 days, std=131.0

✅ GOOD: [4*]C(C)C
   freq=3, mean=4.8 days, std=4.0

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
ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl,False,True
```

**Interpretation:**

| Pattern | Meaning |
|---------|---------|
| `in_outliers=True, in_high_activity=False` | Fragment confuses model (not just extreme activity) |
| `in_outliers=True, in_high_activity=True` | Fragment associated with extreme but predictable activity |
| `in_outliers=False, in_high_activity=True` | Fragment well-handled by model |

**Findings (air_hsbd):**

```
Fragment: ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl
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
  "MACCS_139": {
    "importance": 0.00486250735384588,
    "interpretation": "Structural fingerprint bit"
  },
  "SlogP_VSA10": {
    "importance": 0.005514082061394745,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_108": {
    "importance": 0.005677557702795458,
    "interpretation": "Structural fingerprint bit"
  },
  "SlogP_VSA11": {
    "importance": 0.005923494935154256,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_074": {
    "importance": 0.006240053934261788,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_aniline": {
    "importance": 0.006543268494457808,
    "interpretation": "Molecular descriptor"
  },
  "HallKierAlpha": {
    "importance": 0.006746149600909546,
    "interpretation": "Molecular descriptor"
  },
  "EState_VSA2": {
    "importance": 0.007289120961940459,
    "interpretation": "Molecular descriptor"
  },
  "fr_NH1": {
    "importance": 0.007386992454378363,
    "interpretation": "Molecular descriptor"
  },
  "fr_aryl_methyl": {
    "importance": 0.007886619298873862,
    "interpretation": "Molecular descriptor"
  },
  "BCUT2D_MWHI": {
    "importance": 0.00825859636888033,
    "interpretation": "Molecular descriptor"
  },
  "NOCount": {
    "importance": 0.008309445822377862,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_062": {
    "importance": 0.009041259021433528,
    "interpretation": "Structural fingerprint bit"
  },
  "BalabanJ": {
    "importance": 0.009613156995097408,
    "interpretation": "Balaban connectivity index"
  },
  "SlogP_VSA5": {
    "importance": 0.010236326292920418,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_107": {
    "importance": 0.012680866465279647,
    "interpretation": "Structural fingerprint bit"
  },
  "PEOE_VSA6": {
    "importance": 0.01796849502726258,
    "interpretation": "Molecular descriptor"
  },
  "fr_ketone": {
    "importance": 0.02027654410241008,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_119": {
    "importance": 0.041187803561149046,
    "interpretation": "Structural fingerprint bit"
  },
  "SMR_VSA6": {
    "importance": 0.05253089888936892,
    "interpretation": "Molecular descriptor"
  }
}
```

**Findings (air_hsbd):**

```
✅ GOOD: SMR_VSA6 is top feature
   → Molecular descriptor

✅ GOOD: Functional group / count descriptors present (4 fr_*, 0 Num*)
⚠️ Note: 6 MACCS fingerprint bit(s) in feature set
✅ GOOD: Physicochemical/charge descriptors present (10 features)
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
c1ccccc1,29,12.605,158.874,2.504,True
c1ccc(-c2ccccc2)cc1,6,180.87,32714.095,359.084,False
,21,43.302,1875.035,61.93,False
c1ccc2c(c1)Oc1ccccc1O2,5,2.173,4.72,1.949,True
c1ccc2ccccc2c1,3,0.631,0.398,1.033,True
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (air_hsbd):**

```
✅ GOOD: benzene (29 samples)
   std=12.6 days, MAE=2.5 days, stable=True

⚠️ MEDIUM: biphenyl (6 samples)
   std=180.9 days, MAE=359.1 days, stable=False

⚠️ MEDIUM: non-ring / acyclic (21 samples)
   std=43.3 days, MAE=61.9 days, stable=False

✅ GOOD: diphenyl ether (fused) (5 samples)
   std=2.2 days, MAE=1.9 days, stable=True

✅ GOOD: naphthalene (3 samples)
   std=0.6 days, MAE=1.0 days, stable=True

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
c1ccccc1,29,10,34.5
,21,18,85.7
c1ccc(-c2ccccc2)cc1,6,1,16.7
c1ccc2c(c1)Oc1ccccc1O2,5,0,0.0
c1ccc2ccccc2c1,3,1,33.3
c1ccncc1,2,2,100.0
C1CCCCC1,1,1,100.0
c1ccc(Cc2ccccc2)cc1,1,1,100.0
C1=CC2CC1[C@@H]1[C@H]2[C@H]2C=C[C@@H]1C2,1,1,100.0
c1ccc2c(c1)oc1ccccc12,1,0,0.0
C1CCCC1,1,0,0.0
c1cc2cccc3c4cccc5cccc(c(c1)c23)c54,1,0,0.0
c1cc2ccc3cccc4ccc(c1)c2c34,1,0,0.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (air_hsbd):**

```
⚠️ CONCERNING: benzene (29 samples)
   10/29 outside AD (34.5%)

❌ BAD: non-ring / acyclic (21 samples)
   18/21 outside AD (85.7%)

✅ GOOD: biphenyl (6 samples)
   1/6 outside AD (16.7%)

✅ GOOD: diphenyl ether (fused) (5 samples)
   0/5 outside AD (0.0%)

⚠️ CONCERNING: naphthalene (3 samples)
   1/3 outside AD (33.3%)

❌ BAD: pyridine (2 samples)
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
