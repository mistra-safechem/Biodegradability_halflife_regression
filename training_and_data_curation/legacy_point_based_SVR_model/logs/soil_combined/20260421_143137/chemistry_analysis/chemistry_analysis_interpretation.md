# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/soil_combined/20260421_143137/chemistry_analysis/
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
,47,27.01,27.01
c1ccccc1,39,22.41,49.42
c1ccc(-c2ccccc2)cc1,6,3.45,52.87
c1ccc2ccccc2c1,6,3.45,56.32
c1ccncc1,5,2.87,59.19
c1ccc2c(c1)Oc1ccccc1O2,4,2.3,61.49
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (soil_combined):**
- 174 total molecules, Very high scaffold diversity — 68 unique scaffolds (many singletons expected).
- Top scaffold: non-ring / acyclic (``) at 27.0%
- Top 3 scaffolds cover 52.9% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 57/68 scaffolds appear only once (84%) → long tail, limited generalisation per scaffold

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
,47,27.01,102.14,0.211,-13.447
c1ccccc1,39,22.41,138.612,0.038,-46.994
c1ccc(-c2ccccc2)cc1,6,3.45,542.713,0.154,275.894
c1ccc2ccccc2c1,6,3.45,204.527,-0.474,-26.928
c1ccncc1,5,2.87,44.415,-1.408,26.868
c1ccc2c(c1)Oc1ccccc1O2,4,2.3,106.684,0.802,44.102
c1ccc(Cc2ccccc2)cc1,2,1.15,196.195,-0.186,120.581
N=c1nc[nH]c(=N)[nH]1,2,1.15,29.802,-0.553,28.733
c1cncnc1,2,1.15,147.317,-2.261,-109.285
c1ccc(Nc2ccccc2)cc1,2,1.15,46.264,0.255,-42.142
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (soil_combined):**

```
⚠️ MEDIUM R²: scaffold=non-ring / acyclic (47 samples, 27.0%)
   R²=0.211, RMSE=102 days, mean_error=-13.4 days

❌ BAD R²: scaffold=benzene (39 samples, 22.4%)
   R²=0.038, RMSE=139 days, mean_error=-47.0 days

❌ BAD R²: scaffold=biphenyl (6 samples, 3.5%)
   R²=0.154, RMSE=543 days, mean_error=+275.9 days | systematic over-prediction (+276 days)

❌ BAD R²: scaffold=naphthalene (6 samples, 3.5%)
   R²=-0.474, RMSE=205 days, mean_error=-26.9 days

❌ BAD R²: scaffold=pyridine (5 samples, 2.9%)
   R²=-1.408, RMSE=44 days, mean_error=+26.9 days

✅ GOOD R²: scaffold=diphenyl ether (fused) (4 samples, 2.3%)
   R²=0.802, RMSE=107 days, mean_error=+44.1 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ccc(Cc2ccccc2)cc1 (2 samples, 1.1%)
   R²=-0.186, RMSE=196 days, mean_error=+120.6 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=N=c1nc[nH]c(=N)[nH]1 (2 samples, 1.1%)
   R²=-0.553, RMSE=30 days, mean_error=+28.7 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1cncnc1 (2 samples, 1.1%)
   R²=-2.261, RMSE=147 days, mean_error=-109.3 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=c1ccc(Nc2ccccc2)cc1 (2 samples, 1.1%)
   R²=0.255, RMSE=46 days, mean_error=-42.1 days

```

**Actionable insights:**
- **Trust:** diphenyl ether (fused) predictions
- **Ignore:** benzene, biphenyl, naphthalene, pyridine predictions
- **Question:** non-ring / acyclic, c1ccc(Cc2ccccc2)cc1, N=c1nc[nH]c(=N)[nH]1, c1cncnc1, c1ccc(Nc2ccccc2)cc1 predictions

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
[3*]O[3*],28,152.87,288.361
[4*]CC,22,117.985,162.245
[16*]c1ccccc1,14,204.847,402.809
[4*]CCC,12,68.771,46.628
[8*]C(F)(F)F,10,106.942,126.45
[5*]N[5*],10,400.702,1017.402
[5*]N([5*])[5*],9,113.874,153.511
[8*]CC,9,263.912,283.402
[3*]OC,8,16.788,10.691
[1*]C([6*])=O,8,241.604,487.131
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

**Findings (soil_combined):**

```
→: [3*]O[3*]
   freq=28, mean=152.9 days, std=288.4

→: [4*]CC
   freq=22, mean=118.0 days, std=162.2

→: [16*]c1ccccc1
   freq=14, mean=204.8 days, std=402.8

✅ GOOD: [4*]CCC
   freq=12, mean=68.8 days, std=46.6

→: [8*]C(F)(F)F
   freq=10, mean=106.9 days, std=126.5

⚠️ WARNING: [5*]N[5*]
   freq=10, mean=400.7 days, std=1017.4

→: [5*]N([5*])[5*]
   freq=9, mean=113.9 days, std=153.5

→: [8*]CC
   freq=9, mean=263.9 days, std=283.4

✅ GOOD: [3*]OC
   freq=8, mean=16.8 days, std=10.7

→: [1*]C([6*])=O
   freq=8, mean=241.6 days, std=487.1

→: [4*]CCCC
   freq=7, mean=425.9 days, std=368.1

→: [16*]c1ccc(Cl)cc1
   freq=7, mean=301.6 days, std=273.1

⚠️ WARNING: [4*]C[8*]
   freq=7, mean=379.3 days, std=621.4

✅ GOOD: [8*]C(C)C
   freq=5, mean=49.8 days, std=24.3

→: [5*]N(C)C
   freq=5, mean=197.7 days, std=203.9

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
  "Kappa2": {
    "importance": 0.002931659894584887,
    "interpretation": "Molecular descriptor"
  },
  "fr_N_O": {
    "importance": 0.003479348218346495,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_139": {
    "importance": 0.0036083743272168956,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_hdrzone": {
    "importance": 0.003909782001282491,
    "interpretation": "Molecular descriptor"
  },
  "fr_lactone": {
    "importance": 0.00391767586896114,
    "interpretation": "Molecular descriptor"
  },
  "PEOE_VSA10": {
    "importance": 0.004008705498060922,
    "interpretation": "Molecular descriptor"
  },
  "EState_VSA9": {
    "importance": 0.0043325046055995575,
    "interpretation": "Molecular descriptor"
  },
  "PEOE_VSA8": {
    "importance": 0.0043724926895350015,
    "interpretation": "Molecular descriptor"
  },
  "Chi2v": {
    "importance": 0.004397439103344614,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_104": {
    "importance": 0.004406051141850059,
    "interpretation": "Structural fingerprint bit"
  },
  "PEOE_VSA3": {
    "importance": 0.004432636018769887,
    "interpretation": "Molecular descriptor"
  },
  "fr_ArN": {
    "importance": 0.005728579417577051,
    "interpretation": "Molecular descriptor"
  },
  "RingCount": {
    "importance": 0.006211632703400999,
    "interpretation": "Total ring count"
  },
  "MACCS_067": {
    "importance": 0.006464500610622881,
    "interpretation": "Structural fingerprint bit"
  },
  "SlogP_VSA8": {
    "importance": 0.006643073167140684,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_109": {
    "importance": 0.01174250281311861,
    "interpretation": "Structural fingerprint bit"
  }
}
```

**Findings (soil_combined):**

```
⚠️ WARNING: MACCS_109 is top feature
   → Structural fingerprint bit — unclear chemical meaning

✅ GOOD: Functional group / count descriptors present (4 fr_*, 0 Num*)
⚠️ Note: 4 MACCS fingerprint bit(s) in feature set
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
,47,34.72,1205.458,32.872,True
c1ccccc1,39,43.699,1909.581,60.87,True
c1ccncc1,5,23.329,544.221,37.346,True
c1ccc2c(c1)Oc1ccccc1O2,4,306.442,93906.936,102.064,False
c1ccc2ccccc2c1,6,74.665,5574.8,128.136,True
c1ccc(-c2ccccc2)cc1,6,944.863,892765.595,483.385,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (soil_combined):**

```
✅ GOOD: non-ring / acyclic (47 samples)
   std=34.7 days, MAE=32.9 days, stable=True

✅ GOOD: benzene (39 samples)
   std=43.7 days, MAE=60.9 days, stable=True

✅ GOOD: pyridine (5 samples)
   std=23.3 days, MAE=37.3 days, stable=True

⚠️ MEDIUM: diphenyl ether (fused) (4 samples)
   std=306.4 days, MAE=102.1 days, stable=False

✅ GOOD: naphthalene (6 samples)
   std=74.7 days, MAE=128.1 days, stable=True

❌ BAD: biphenyl (6 samples)
   std=944.9 days, MAE=483.4 days, stable=False

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
,47,10,21.3
c1ccccc1,39,2,5.1
c1ccc2ccccc2c1,6,1,16.7
c1ccc(-c2ccccc2)cc1,6,0,0.0
c1ccncc1,5,3,60.0
c1ccc2c(c1)Oc1ccccc1O2,4,0,0.0
N=c1nc[nH]c(=N)[nH]1,2,1,50.0
c1ccc(Nc2ccccc2)cc1,2,1,50.0
c1ccc2c(c1)oc1ccccc12,2,0,0.0
c1cncnc1,2,0,0.0
c1ccc(Cc2ccccc2)cc1,2,0,0.0
c1cc2c(c(-c3cc[nH]c3)c1)OCO2,1,1,100.0
O=C(N=Cc1ccccc1)Nc1ccccc1,1,1,100.0
C(=NNC=Nc1ccccc1)c1ccccn1,1,1,100.0
C1CCOCC1,1,1,100.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (soil_combined):**

```
⚠️ CONCERNING: non-ring / acyclic (47 samples)
   10/47 outside AD (21.3%)

✅ GOOD: benzene (39 samples)
   2/39 outside AD (5.1%)

✅ GOOD: naphthalene (6 samples)
   1/6 outside AD (16.7%)

✅ GOOD: biphenyl (6 samples)
   0/6 outside AD (0.0%)

❌ BAD: pyridine (5 samples)
   3/5 outside AD (60.0%)

✅ GOOD: diphenyl ether (fused) (4 samples)
   0/4 outside AD (0.0%)

⚠️ CONCERNING: N=c1nc[nH]c(=N)[nH]1 (2 samples)
   1/2 outside AD (50.0%)

⚠️ CONCERNING: c1ccc(Nc2ccccc2)cc1 (2 samples)
   1/2 outside AD (50.0%)

✅ GOOD: c1ccc2c(c1)oc1ccccc12 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: c1cncnc1 (2 samples)
   0/2 outside AD (0.0%)

✅ GOOD: c1ccc(Cc2ccccc2)cc1 (2 samples)
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
