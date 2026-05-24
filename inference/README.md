# Biodegradability Regression Half-Life Modelling
## MistraSafeChem WP3, Deliverable 3.2.2.2, Task 3.2.1.


# Standalone SVR Interval Inference
This folder contains code for running inference with the SVR Interval model in a standalone manner.

Code here is stripped down to inference only, with reduced overhead and 
dependencies, and excludes all training and database-building code.

---
## Model Date
The model files are from training date *20260516*. The timestamps were removed 
from the filenames for convenience (see optional model asset cleanup below), but
the original files with timestamps are still available in `training_and_data_curation/models` 
and can be copied back if needed.

---
## Installation - generic
For users running this inference package directly, install dependencies from
`pyproject.toml` with:
```bash
pip install .
```
(installs from `pyproject.toml`)

If you have UV package manager it will automatically create a venv and install dependencies there:
```bash
uv sync
# alt to select specific python version (should be 3.11 or later):
uv sync --python 3.11
```

## Installation - MistraSafeChem stack
For the MistraSafeChem admin toolbox setup, use the pinned environment file:
[requirements_msc_toolbox.yaml](requirements_msc_toolbox.yaml)

Create and activate the environment with conda/mamba:
```bash
conda env create -f requirements_msc_toolbox.yaml
conda activate msc-toolbox-py397
```

If the environment already exists and you need to update it:
```bash
conda env update -f requirements_msc_toolbox.yaml --prune
```

For this MSC toolbox setup, package installation from this folder is usually not required; 
it should be possible to run the inference script directly.


### Dependencies
Recommended: Python >=3.11 (3.11, 3.12, 3.13 tested).

The MistraSafeChem admin toolbox environment file can be used as-is.
If running with Python 3.10 or 3.9, inference may still work, but additional warnings can appear.

---
## Usage
- To run inference, use the `SVR_interval_inference.py` script, details on args are in the header section of that script.
- A bash script example for running inference (using UV) is provided in `SVR_interval_inference_test.sh` as well.

Example command:
```bash
python ./SVR_interval_inference.py \
  --smiles 'C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl' \
  --model models/SVR_air_hsbd.json \
  --json-out test.json
```

---
## Optional model asset cleanup
If you need or want to remove timestamp segments from the filenames under `inference/models`, use `optional_rename_inference_model_assets.py`. 
It only operates on that folder and also rewrites the path fields in the JSON files to match the renamed local assets.

For convenience, this script has been run and the model names are updated and are without timestamps.

---
## Contributions
See the main [README](../README.md) for contributions and credits.

---
## License
This project is licensed under the MIT License - see the [LICENSE](../LICENSE.md) file for details.