# Models directory
Contains model files .joblib with connected json file containing meta-data 
as well as .npz file with AD and other related data.

Contents:
------
A model consists of the following 5 files, named `*_<compartment>_<data_source>_<timestamp>.*`:
  - SVR_<c>_<ds>_<ts>.joblib           (model file)
  - SVR_<c>_<ds>_<ts>.json             (meta-data file)
  - SVR_<c>_<ds>_<ts>_ad.npz           (applicability domain artefact for inference)
  - t_half_meta_<c>_<ds>_<ts>.csv      (training-split T_half metadata for inference traceability)
  - training_split_<c>_<ds>_<ts>.csv   (SMILES + train/test split info for inference traceability)


The optional script `optional_rename_inference_model_assets.py` was run 
to remove timestamps from the inference model assets, the files under 
`inference/models` were thus renamed to `SVR_<c>_<ds>.joblib`, etc. for convenience.

If origin is important, the original files with timestamps are still available in
 `training_and_data_curation/models` and can be copied back if needed.

# Model performance overview and comparison
Models here are provided as is without interpretation or comparison.
A comparison/overview is provided in the repos `training_and_data_curation/logs` directory.

# Date of latest run
Final (re)run 2026-05-16 using updated relative paths instead of previous absolute paths;
won't impact the results since they should be reproducible.
