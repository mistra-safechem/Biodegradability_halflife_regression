# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/sediment_hsbd/20260421_142534/chemistry_analysis/
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
,28,34.57,34.57
c1ccccc1,25,30.86,65.43
c1ccc(-c2ccccc2)cc1,4,4.94,70.37
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,2,2.47,72.84
c1ccc2c(c1)Cc1ccccc1-2,1,1.23,74.07
c1ccc2c(c1)ccc1c3ccccc3ccc21,1,1.23,75.3
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (sediment_hsbd):**
- 81 total molecules, 26 unique scaffolds.
- Top scaffold: non-ring / acyclic (``) at 34.6%
- Top 3 scaffolds cover 70.4% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 22/26 scaffolds appear only once (85%) → long tail, limited generalisation per scaffold

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
,28,34.57,142.211,0.596,-38.924
c1ccccc1,25,30.86,407.878,0.177,-75.504
c1ccc(-c2ccccc2)cc1,4,4.94,1091.932,-0.235,129.528
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,2,2.47,64.297,-587.696,63.782
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (sediment_hsbd):**

```
✅ GOOD R²: scaffold=non-ring / acyclic (28 samples, 34.6%)
   R²=0.596, RMSE=142 days, mean_error=-38.9 days

❌ BAD R²: scaffold=benzene (25 samples, 30.9%)
   R²=0.177, RMSE=408 days, mean_error=-75.5 days

❌ BAD R²: scaffold=biphenyl (4 samples, 4.9%)
   R²=-0.235, RMSE=1092 days, mean_error=+129.5 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1 (2 samples, 2.5%)
   R²=-587.696, RMSE=64 days, mean_error=+63.8 days

```

**Actionable insights:**
- **Trust:** non-ring / acyclic predictions
- **Ignore:** benzene, biphenyl predictions
- **Question:** O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1 predictions

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
[3*]OC,6,69.444,37.769
[6*]C(=O)O,5,100.167,90.06
[16*]c1ccccc1,5,548.933,977.752
[3*]O[3*],5,140.95,244.804
[4*]C(C)C,4,192.208,263.008
[5*]N([5*])[5*],4,103.096,103.357
[16*]c1ccc(C)cc1,3,6.7,0.0
[4*]CC,3,158.444,122.494
[4*]CC(C)O,3,14.3,0.0
[16*]c1ccccc1Cl,3,1604.167,1190.785
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

**Findings (sediment_hsbd):**

```
✅ GOOD: [3*]OC
   freq=6, mean=69.4 days, std=37.8

→: [6*]C(=O)O
   freq=5, mean=100.2 days, std=90.1

⚠️ WARNING: [16*]c1ccccc1
   freq=5, mean=548.9 days, std=977.8

→: [3*]O[3*]
   freq=5, mean=140.9 days, std=244.8

→: [4*]C(C)C
   freq=4, mean=192.2 days, std=263.0

→: [5*]N([5*])[5*]
   freq=4, mean=103.1 days, std=103.4

✅ GOOD: [16*]c1ccc(C)cc1
   freq=3, mean=6.7 days, std=0.0

→: [4*]CC
   freq=3, mean=158.4 days, std=122.5

✅ GOOD: [4*]CC(C)O
   freq=3, mean=14.3 days, std=0.0

⚠️ WARNING: [16*]c1ccccc1Cl
   freq=3, mean=1604.2 days, std=1190.8

✅ GOOD: [16*]c1ccccc1C
   freq=3, mean=12.0 days, std=0.0

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
  "MACCS_114": {
    "importance": 0.005917561932265744,
    "interpretation": "Structural fingerprint bit"
  },
  "MolLogP": {
    "importance": 0.006283034325223674,
    "interpretation": "Lipophilicity (partition coefficient)"
  },
  "BCUT2D_CHGLO": {
    "importance": 0.006375876689979199,
    "interpretation": "Molecular descriptor"
  },
  "fr_Ar_N": {
    "importance": 0.006466373057061489,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_141": {
    "importance": 0.006670935291555355,
    "interpretation": "Structural fingerprint bit"
  },
  "VSA_EState1": {
    "importance": 0.006757366690665584,
    "interpretation": "Molecular descriptor"
  },
  "fr_ketone_Topliss": {
    "importance": 0.007008346501965751,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_146": {
    "importance": 0.007279207115384381,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_Ar_OH": {
    "importance": 0.007310255177986005,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_163": {
    "importance": 0.007980522235819673,
    "interpretation": "Structural fingerprint bit"
  },
  "MACCS_096": {
    "importance": 0.008376065182105649,
    "interpretation": "Structural fingerprint bit"
  },
  "NumHeteroatoms": {
    "importance": 0.008692496903394778,
    "interpretation": "Heteroatom count"
  },
  "MACCS_099": {
    "importance": 0.008847711941639375,
    "interpretation": "Structural fingerprint bit"
  },
  "SlogP_VSA10": {
    "importance": 0.009038488520196512,
    "interpretation": "Molecular descriptor"
  },
  "NumRotatableBonds": {
    "importance": 0.009714248423605687,
    "interpretation": "Molecular flexibility"
  },
  "MACCS_165": {
    "importance": 0.00971461341111134,
    "interpretation": "Structural fingerprint bit"
  },
  "PEOE_VSA10": {
    "importance": 0.010856032032676228,
    "interpretation": "Molecular descriptor"
  },
  "EState_VSA9": {
    "importance": 0.01098224795665614,
    "interpretation": "Molecular descriptor"
  },
  "fr_Al_COO": {
    "importance": 0.011616295115302766,
    "interpretation": "Molecular descriptor"
  },
  "SMR_VSA1": {
    "importance": 0.013943041429516984,
    "interpretation": "Molecular descriptor"
  }
}
```

**Findings (sediment_hsbd):**

```
✅ GOOD: SMR_VSA1 is top feature
   → Molecular descriptor

✅ GOOD: Functional group / count descriptors present (4 fr_*, 2 Num*)
⚠️ Note: 7 MACCS fingerprint bit(s) in feature set
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
c1ccccc1,25,123.424,15233.379,172.465,True
,28,127.099,16154.23,73.002,True
c1ccc(-c2ccccc2)cc1,4,1519.552,2309037.145,927.484,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (sediment_hsbd):**

```
✅ GOOD: benzene (25 samples)
   std=123.4 days, MAE=172.5 days, stable=True

✅ GOOD: non-ring / acyclic (28 samples)
   std=127.1 days, MAE=73.0 days, stable=True

❌ BAD: biphenyl (4 samples)
   std=1519.6 days, MAE=927.5 days, stable=False

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
,28,17,60.7
c1ccccc1,25,20,80.0
c1ccc(-c2ccccc2)cc1,4,1,25.0
O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1,2,2,100.0
O=c1ocnn1-c1ccccc1,1,1,100.0
c1ccc(Oc2ccccc2)cc1,1,1,100.0
O=C(N=Cc1ccccc1)Nc1ccccc1,1,1,100.0
N=c1[nH]cc(Cc2ccccc2)c(=N)[nH]1,1,1,100.0
O=C1c2ccccc2C(=O)c2ccccc21,1,1,100.0
c1cncnc1,1,1,100.0
C1=CC2CC1[C@@H]1[C@H]2[C@H]2C=C[C@@H]1C2,1,1,100.0
C1CCOC1,1,1,100.0
c1ccc(Cc2ccccc2)cc1,1,1,100.0
c1ccc(C2CO2)cc1,1,1,100.0
O=c1ccn2c3c(cccc13)CCC2,1,1,100.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (sediment_hsbd):**

```
❌ BAD: non-ring / acyclic (28 samples)
   17/28 outside AD (60.7%)

❌ BAD: benzene (25 samples)
   20/25 outside AD (80.0%)

⚠️ CONCERNING: biphenyl (4 samples)
   1/4 outside AD (25.0%)

❌ BAD: O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1 (2 samples)
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
