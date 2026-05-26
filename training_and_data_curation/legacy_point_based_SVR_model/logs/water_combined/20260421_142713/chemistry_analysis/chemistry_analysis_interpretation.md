# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/water_combined/20260421_142713/chemistry_analysis/
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
,67,38.29,38.29
c1ccccc1,48,27.43,65.72
c1ccc(-c2ccccc2)cc1,8,4.57,70.29
c1ccc2c(c1)ccc1ccccc12,4,2.29,72.58
c1ccc2c(c1)Oc1ccccc1O2,4,2.29,74.87
C1CCCCC1,4,2.29,77.16
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (water_combined):**
- 175 total molecules, 40 unique scaffolds.
- Top scaffold: non-ring / acyclic (``) at 38.3%
- Top 3 scaffolds cover 70.3% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 29/40 scaffolds appear only once (72%) → long tail, limited generalisation per scaffold

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
,67,38.29,36.367,0.372,-1.216
c1ccccc1,48,27.43,325.768,-0.026,-57.383
c1ccc(-c2ccccc2)cc1,8,4.57,188.208,0.958,-89.575
c1ccc2c(c1)ccc1ccccc12,4,2.29,8.23,-0.031,7.211
c1ccc2c(c1)Oc1ccccc1O2,4,2.29,21.513,0.27,9.495
C1CCCCC1,4,2.29,4.753,0.662,-1.562
c1ccc2ccccc2c1,3,1.71,3.694,-2470107199949965.0,0.617
c1cc2ccc3cccc4ccc(c1)c2c34,2,1.14,27.448,-63.538,-26.872
c1ccc(Cc2ccccc2)cc1,2,1.14,27.381,-813.883,14.908
C1CCCC1,2,1.14,12.868,-1.265,-9.8
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (water_combined):**

```
⚠️ MEDIUM R²: scaffold=non-ring / acyclic (67 samples, 38.3%)
   R²=0.372, RMSE=36 days, mean_error=-1.2 days

❌ BAD R²: scaffold=benzene (48 samples, 27.4%)
   R²=-0.026, RMSE=326 days, mean_error=-57.4 days

✅ GOOD R²: scaffold=biphenyl (8 samples, 4.6%)
   R²=0.958, RMSE=188 days, mean_error=-89.6 days

❌ BAD R²: scaffold=c1ccc2c(c1)ccc1ccccc12 (4 samples, 2.3%)
   R²=-0.031, RMSE=8 days, mean_error=+7.2 days

⚠️ MEDIUM R²: scaffold=diphenyl ether (fused) (4 samples, 2.3%)
   R²=0.270, RMSE=22 days, mean_error=+9.5 days

✅ GOOD R²: scaffold=cyclohexane (4 samples, 2.3%)
   R²=0.662, RMSE=5 days, mean_error=-1.6 days

❌ BAD R²: scaffold=naphthalene (3 samples, 1.7%)
   R²=-2470107199949965.000, RMSE=4 days, mean_error=+0.6 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1cc2ccc3cccc4ccc(c1)c2c34 (2 samples, 1.1%)
   R²=-63.538, RMSE=27 days, mean_error=-26.9 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ccc(Cc2ccccc2)cc1 (2 samples, 1.1%)
   R²=-813.883, RMSE=27 days, mean_error=+14.9 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=C1CCCC1 (2 samples, 1.1%)
   R²=-1.265, RMSE=13 days, mean_error=-9.8 days

```

**Actionable insights:**
- **Trust:** biphenyl, cyclohexane predictions
- **Ignore:** benzene, c1ccc2c(c1)ccc1ccccc12, naphthalene predictions
- **Question:** non-ring / acyclic, diphenyl ether (fused), c1cc2ccc3cccc4ccc(c1)c2c34, c1ccc(Cc2ccccc2)cc1, C1CCCC1 predictions

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
[3*]O[3*],24,10.545,11.383
[16*]c1ccccc1,14,71.255,192.869
[4*]CC,10,79.05,221.377
[3*]OC,7,5.922,3.534
[1*]C([6*])=O,7,6.719,1.824
[16*]c1ccc([16*])cc1,6,8.455,7.315
[8*]CC,5,208.68,282.37
[4*]CC[4*],5,16.084,21.832
[6*]C(=O)O,4,2.375,3.31
[4*]CCC,4,14.758,9.944
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

**Findings (water_combined):**

```
✅ GOOD: [3*]O[3*]
   freq=24, mean=10.5 days, std=11.4

→: [16*]c1ccccc1
   freq=14, mean=71.3 days, std=192.9

→: [4*]CC
   freq=10, mean=79.0 days, std=221.4

✅ GOOD: [3*]OC
   freq=7, mean=5.9 days, std=3.5

✅ GOOD: [1*]C([6*])=O
   freq=7, mean=6.7 days, std=1.8

✅ GOOD: [16*]c1ccc([16*])cc1
   freq=6, mean=8.5 days, std=7.3

→: [8*]CC
   freq=5, mean=208.7 days, std=282.4

✅ GOOD: [4*]CC[4*]
   freq=5, mean=16.1 days, std=21.8

✅ GOOD: [6*]C(=O)O
   freq=4, mean=2.4 days, std=3.3

✅ GOOD: [4*]CCC
   freq=4, mean=14.8 days, std=9.9

✅ GOOD: [4*]CCCC
   freq=4, mean=5.3 days, std=2.1

✅ GOOD: [4*]CCO
   freq=4, mean=0.3 days, std=0.2

→: [5*]N[5*]
   freq=4, mean=362.4 days, std=399.7

✅ GOOD: [16*]c1ccc(N)cc1
   freq=4, mean=8.0 days, std=1.1

✅ GOOD: [16*]c1ccccc1[16*]
   freq=4, mean=11.3 days, std=7.8

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
ClC(Cl)Br,False,True
```

**Interpretation:**

| Pattern | Meaning |
|---------|---------|
| `in_outliers=True, in_high_activity=False` | Fragment confuses model (not just extreme activity) |
| `in_outliers=True, in_high_activity=True` | Fragment associated with extreme but predictable activity |
| `in_outliers=False, in_high_activity=True` | Fragment well-handled by model |

**Findings (water_combined):**

```
Fragment: ClC(Cl)Br
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
  "fr_NH1": {
    "importance": 0.003069265433320617,
    "interpretation": "Molecular descriptor"
  },
  "fr_ketone": {
    "importance": 0.0032911761015476217,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_024": {
    "importance": 0.0033503977017614595,
    "interpretation": "Structural fingerprint bit"
  },
  "Kappa3": {
    "importance": 0.0035303596294403326,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_021": {
    "importance": 0.003584682137388834,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_037": {
    "importance": 0.0038604434282516475,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_065": {
    "importance": 0.0042300393963789155,
    "interpretation": "Structural fingerprint bit"
  },
  "VSA_EState2": {
    "importance": 0.004308702704361772,
    "interpretation": "Molecular descriptor"
  },
  "BertzCT": {
    "importance": 0.004600580702807455,
    "interpretation": "Bertz complexity index"
  },
  "MaxEStateIndex": {
    "importance": 0.00462051577869344,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_101": {
    "importance": 0.0046562604535024255,
    "interpretation": "Structural fingerprint bit"
  },
  "VSA_EState5": {
    "importance": 0.005138683714010278,
    "interpretation": "Molecular descriptor"
  },
  "fr_oxime": {
    "importance": 0.00517827878360883,
    "interpretation": "Molecular descriptor"
  },
  "fr_COO": {
    "importance": 0.0060586352723148995,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_050": {
    "importance": 0.006545883825236387,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_079": {
    "importance": 0.010976383397335004,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_136": {
    "importance": 0.023624497289309283,
    "interpretation": "Structural fingerprint bit"
  }
}
```

**Findings (water_combined):**

```
⚠️ WARNING: MACCS_136 is top feature
   → Structural fingerprint bit — unclear chemical meaning

✅ GOOD: Functional group / count descriptors present (4 fr_*, 0 Num*)
⚠️ Note: 8 MACCS fingerprint bit(s) in feature set
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
,67,29.197,852.48,14.266,True
c1ccc2c(c1)Oc1ccccc1O2,4,42.273,1787.014,15.243,True
c1ccccc1,48,53.667,2880.114,73.355,False
c1ccc(-c2ccccc2)cc1,8,842.159,709232.302,110.182,False
c1ccc2ccccc2c1,3,3.642,13.264,3.111,True
c1ccc2c(c1)ccc1ccccc12,4,4.144,17.173,7.211,True
C1CCCCC1,4,3.809,14.507,3.056,True
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (water_combined):**

```
✅ GOOD: non-ring / acyclic (67 samples)
   std=29.2 days, MAE=14.3 days, stable=True

✅ GOOD: diphenyl ether (fused) (4 samples)
   std=42.3 days, MAE=15.2 days, stable=True

⚠️ MEDIUM: benzene (48 samples)
   std=53.7 days, MAE=73.4 days, stable=False

❌ BAD: biphenyl (8 samples)
   std=842.2 days, MAE=110.2 days, stable=False

✅ GOOD: naphthalene (3 samples)
   std=3.6 days, MAE=3.1 days, stable=True

✅ GOOD: c1ccc2c(c1)ccc1ccccc12 (4 samples)
   std=4.1 days, MAE=7.2 days, stable=True

✅ GOOD: cyclohexane (4 samples)
   std=3.8 days, MAE=3.1 days, stable=True

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
,67,11,16.4
c1ccccc1,48,14,29.2
c1ccc(-c2ccccc2)cc1,8,0,0.0
C1CCCCC1,4,0,0.0
c1ccc2c(c1)Oc1ccccc1O2,4,0,0.0
c1ccc2c(c1)ccc1ccccc12,4,0,0.0
c1ccc2ccccc2c1,3,0,0.0
c1ccncc1,2,0,0.0
C1CCCC1,2,0,0.0
c1ccc(Cc2ccccc2)cc1,2,0,0.0
c1cc2ccc3cccc4ccc(c1)c2c34,2,0,0.0
c1ncncn1,1,1,100.0
C1C[C@H]2C[C@@H]3CCC[C@]3(C1)C2,1,1,100.0
C=C1CC2CCC1C2,1,1,100.0
C1=C2CCCCC2C2CCc3ccccc3C2=C1,1,1,100.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (water_combined):**

```
✅ GOOD: non-ring / acyclic (67 samples)
   11/67 outside AD (16.4%)

⚠️ CONCERNING: benzene (48 samples)
   14/48 outside AD (29.2%)

✅ GOOD: biphenyl (8 samples)
   0/8 outside AD (0.0%)

✅ GOOD: cyclohexane (4 samples)
   0/4 outside AD (0.0%)

✅ GOOD: diphenyl ether (fused) (4 samples)
   0/4 outside AD (0.0%)

✅ GOOD: c1ccc2c(c1)ccc1ccccc12 (4 samples)
   0/4 outside AD (0.0%)

✅ GOOD: naphthalene (3 samples)
   0/3 outside AD (0.0%)

✅ GOOD: pyridine (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: C1CCCC1 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: c1ccc(Cc2ccccc2)cc1 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: c1cc2ccc3cccc4ccc(c1)c2c34 (2 samples)
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
