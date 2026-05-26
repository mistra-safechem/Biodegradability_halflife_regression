# Chemistry Analysis Interpretation Guide
## File Structure
All chemistry analysis outputs are located in:
```
logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
```
For this run:
```
logs/sediment_vega/20260421_142642/chemistry_analysis/
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
c1ccccc1,19,35.85,35.85
,8,15.09,50.94
c1ccc(-c2ccccc2)cc1,7,13.21,64.15
c1ccc2c(c1)Oc1ccccc1O2,4,7.55,71.7
c1ccc2ccccc2c1,2,3.77,75.47
C1CCCCC1,2,3.77,79.24
```

**Interpretation:**

| Metric | Good | Bad | Action |
|--------|------|-----|--------|
| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |
| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |
| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |

**Findings (sediment_vega):**
- 53 total molecules, 17 unique scaffolds.
- Top scaffold: benzene (`c1ccccc1`) at 35.9%
- Top 3 scaffolds cover 64.2% of dataset
- Empty scaffold present (non-ring / acyclic molecules in dataset)
- ⚠️ 11/17 scaffolds appear only once (65%) → long tail, limited generalisation per scaffold

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
c1ccccc1,19,35.85,133.025,0.646,-22.558
,8,15.09,756.535,-0.092,-353.899
c1ccc(-c2ccccc2)cc1,7,13.21,575.536,0.353,-310.326
c1ccc2c(c1)Oc1ccccc1O2,4,7.55,315.252,0.837,143.74
c1ccc2ccccc2c1,2,3.77,129.544,0.0,129.52
C1CCCCC1,2,3.77,1222.445,-0.212,-857.134
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |
| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |
| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |
| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |

**Findings (sediment_vega):**

```
✅ GOOD R²: scaffold=benzene (19 samples, 35.9%)
   R²=0.646, RMSE=133 days, mean_error=-22.6 days

❌ BAD R²: scaffold=non-ring / acyclic (8 samples, 15.1%)
   R²=-0.092, RMSE=757 days, mean_error=-353.9 days | systematic under-prediction (-354 days)

⚠️ MEDIUM R²: scaffold=biphenyl (7 samples, 13.2%)
   R²=0.353, RMSE=576 days, mean_error=-310.3 days | systematic under-prediction (-310 days)

⚠️ GOOD R²: scaffold=diphenyl ether (fused) (4 samples, 7.5%)
   R²=0.837, RMSE=315 days, mean_error=+143.7 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=naphthalene (2 samples, 3.8%)
   R²=0.000, RMSE=130 days, mean_error=+129.5 days

⚠️ UNRELIABLE (n≤2) R²: scaffold=cyclohexane (2 samples, 3.8%)
   R²=-0.212, RMSE=1222 days, mean_error=-857.1 days | systematic under-prediction (-857 days)

```

**Actionable insights:**
- **Trust:** benzene predictions
- **Ignore:** non-ring / acyclic, biphenyl, cyclohexane predictions
- **Question:** diphenyl ether (fused), naphthalene predictions

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
[16*]c1ccccc1,6,1076.389,965.433
[3*]O[3*],5,70.833,0.0
[3*]OC,3,38.889,27.665
[16*]c1cc(Cl)ccc1Cl,3,1763.889,914.138
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

**Findings (sediment_vega):**

```
⚠️ WARNING: [16*]c1ccccc1
   freq=6, mean=1076.4 days, std=965.4

✅ GOOD: [3*]O[3*]
   freq=5, mean=70.8 days, std=0.0

✅ GOOD: [3*]OC
   freq=3, mean=38.9 days, std=27.7

⚠️ WARNING: [16*]c1cc(Cl)ccc1Cl
   freq=3, mean=1763.9 days, std=914.1

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

**Findings (sediment_vega):**

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
  "MACCS_129": {
    "importance": 0.002084393363456981,
    "interpretation": "Structural fingerprint bit"
  },
  "VSA_EState2": {
    "importance": 0.0021990654218941334,
    "interpretation": "Molecular descriptor"
  },
  "BCUT2D_LOGPHI": {
    "importance": 0.0022034594418437496,
    "interpretation": "Molecular descriptor"
  },
  "Ipc": {
    "importance": 0.002454937607403597,
    "interpretation": "Molecular descriptor"
  },
  "fr_nitro_arom_nonortho": {
    "importance": 0.002601555009652259,
    "interpretation": "Molecular descriptor"
  },
  "SlogP_VSA2": {
    "importance": 0.002790868837660369,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_125": {
    "importance": 0.0029579780574284995,
    "interpretation": "Structural fingerprint bit"
  },
  "EState_VSA10": {
    "importance": 0.003106295265058329,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_133": {
    "importance": 0.0031350937223280457,
    "interpretation": "Structural fingerprint bit"
  },
  "SMR_VSA5": {
    "importance": 0.0031418400572372452,
    "interpretation": "Molecular descriptor"
  },
  "BCUT2D_MWLOW": {
    "importance": 0.0032137122547828633,
    "interpretation": "Molecular descriptor"
  },
  "NHOHCount": {
    "importance": 0.003494253588057776,
    "interpretation": "Molecular descriptor"
  },
  "fr_nitro": {
    "importance": 0.003991558747576436,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_124": {
    "importance": 0.0041716022577692625,
    "interpretation": "Structural fingerprint bit"
  },
  "fr_aldehyde": {
    "importance": 0.0042450245596831445,
    "interpretation": "Molecular descriptor"
  },
  "fr_alkyl_halide": {
    "importance": 0.004466887673271363,
    "interpretation": "Molecular descriptor"
  },
  "PEOE_VSA4": {
    "importance": 0.004718991049184868,
    "interpretation": "Molecular descriptor"
  },
  "NumAliphaticHeterocycles": {
    "importance": 0.006300733991150738,
    "interpretation": "Molecular descriptor"
  },
  "PEOE_VSA1": {
    "importance": 0.008599295129175857,
    "interpretation": "Molecular descriptor"
  },
  "MACCS_127": {
    "importance": 0.009279711661289304,
    "interpretation": "Structural fingerprint bit"
  }
}
```

**Findings (sediment_vega):**

```
⚠️ WARNING: MACCS_127 is top feature
   → Structural fingerprint bit — unclear chemical meaning

✅ GOOD: Functional group / count descriptors present (4 fr_*, 1 Num*)
⚠️ Note: 5 MACCS fingerprint bit(s) in feature set
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
c1ccccc1,19,216.097,46697.909,74.341,True
c1ccc(-c2ccccc2)cc1,7,638.686,407919.173,456.444,False
c1ccc2c(c1)Oc1ccccc1O2,4,899.729,809512.317,251.693,False
,8,125.845,15837.03,370.125,False
```

**Interpretation:**

| Metric | Good | Medium | Bad |
|--------|------|--------|-----|
| **stable** | True | - | False |
| **prediction_std** | <200 days | 200-500 days | >500 days |
| **mean_absolute_error** | <100 days | 100-300 days | >300 days |

**Findings (sediment_vega):**

```
⚠️ MEDIUM: benzene (19 samples)
   std=216.1 days, MAE=74.3 days, stable=True

❌ BAD: biphenyl (7 samples)
   std=638.7 days, MAE=456.4 days, stable=False

❌ BAD: diphenyl ether (fused) (4 samples)
   std=899.7 days, MAE=251.7 days, stable=False

⚠️ MEDIUM: non-ring / acyclic (8 samples)
   std=125.8 days, MAE=370.1 days, stable=False

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
c1ccccc1,19,12,63.2
,8,8,100.0
c1ccc(-c2ccccc2)cc1,7,0,0.0
c1ccc2c(c1)Oc1ccccc1O2,4,1,25.0
C1CCCCC1,2,2,100.0
c1ccc2ccccc2c1,2,0,0.0
C1=CC2CC1C1C3C=CC(C3)C21,1,1,100.0
c1ccc2ncccc2c1,1,1,100.0
c1ccc(Nc2ccccc2)cc1,1,1,100.0
c1ccoc1,1,1,100.0
c1cc2c3c(cccc3c1)CC2,1,1,100.0
C1=CC2C3C=CC(C3)C2C1,1,1,100.0
O=C1CCCCC1,1,1,100.0
c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34,1,0,0.0
c1ccc2c(c1)ccc1ccccc12,1,0,0.0
```

**Interpretation:**

| Metric | Good | Medium | Bad | Action |
|--------|------|--------|-----|--------|
| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |
| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |

**Findings (sediment_vega):**

```
❌ BAD: benzene (19 samples)
   12/19 outside AD (63.2%)

❌ BAD: non-ring / acyclic (8 samples)
   8/8 outside AD (100.0%)

✅ GOOD: biphenyl (7 samples)
   0/7 outside AD (0.0%)

⚠️ CONCERNING: diphenyl ether (fused) (4 samples)
   1/4 outside AD (25.0%)

❌ BAD: cyclohexane (2 samples)
   2/2 outside AD (100.0%)

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
