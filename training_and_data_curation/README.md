# Biodegradability Regression Half-Life Modelling: Data Curation and Training Code
## MistraSafeChem WP3, Deliverable 3.2.2.2, Task 3.2.1.

This folder contains code for interval-based biodegradability half-life regression development. 
It includes data curation, descriptor calculation, database building, and model training and analysis components.

---
## References
For implementation details and data documentation, see:

- [METHODS.md](METHODS.md)
- [raw_data/readme.md](raw_data/readme.md)
- [processed_data/readme.md](processed_data/readme.md)
- [src/feature_selection_in_rdkit_tools.py.md](src/feature_selection_in_rdkit_tools.py.md)

---
## Reproducibility Overview
The full methodology is documented in [METHODS.md](METHODS.md).

Main workflow components:

- Data extraction, cleaning, descriptor calculation, and database building: [0_1_data_clean_descr_calcs_db_build/readme.md](0_1_data_clean_descr_calcs_db_build/readme.md)
- Interval-based model training and analysis: `SVR_interval_model_and_analysis.py`
- Interval-based inference: `SVR_interval_inference.py`
- Batch scripts for training and inference: `SVR_interval_model_run_all_combos.sh` and `SVR_interval_inference_test.sh` for running the training and inference pipelines.

Optional legacy point-based model:

- [legacy_point_based_SVR_model/readme.md](legacy_point_based_SVR_model/readme.md)

## Requirements
Python >= 3.11 (tested with 3.11, 3.12 and 3.13).

Using uv:
`uv sync`

Using pip:
`pip install .`

## Applicability
At its current stage, the model is not yet considered robust enough for practical deployment. 
This is of course absolutely debatable and up for experts to decide. The repository as is, is intended for research and development use.

The main limitation is data quality and structure: moderate dataset sizes, heterogeneous sources, and interval/discretized endpoint values constrain broad model generalization.

A potentially useful candidate is the interval-based SEDIMENT/Vega model, with:

- Test coverage: 39/53 (73.6%)
- MIL: 1.17x days
- Spearman rho: 0.921
- Kendall tau: 0.822
- Class accuracy: 0.736
- CV R2: 0.7903 +/- 0.0556
- Outside AD: 14/53 (26.4%)

Compared with other interval-based combinations, this model shows the strongest balance of rank performance, coverage, classification performance, 
and structural generalization. Other potentially viable combinations are SOIL/Vega, AIR/HSBD, and WATER/Vega.

Performance data for each compartment/dataset combination may be retrieved from respective combination in `logs/` folder or the corresponding JSON files in `models/` folder.

---
## Contributions
See the main [README](../README.md) for contributions and credits.

## License
This project is licensed under the MIT License - see the [LICENSE](../LICENSE.md) file for details.