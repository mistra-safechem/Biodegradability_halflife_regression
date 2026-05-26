# Models directory
Contains model files .joblib with connected json file containing meta-data as well as .npz file with AD and other related data.

Contents:
------
A model consists of the following 5 files, named `*_<compartment>_<data_source>_<timestamp>.*`:
  - SVR_<c>_<ds>_<ts>.joblib           (model file)
  - SVR_<c>_<ds>_<ts>.json             (meta-data file)
  - SVR_<c>_<ds>_<ts>_ad.npz           (applicability domain artefact for inference)
  - t_half_meta_<c>_<ds>_<ts>.csv      (training-split T_half metadata for inference traceability)
  - training_split_<c>_<ds>_<ts>.csv   (SMILES + train/test split info for inference traceability)

Earlier runs are stored in `.earlier_runs/`.

# Date of last run
Final (re)run 2026-05-16 using updated relative paths instead of previous absolute paths;
won't impact the results since they should be reproducible.
An earlier run from 2026-04-30 is available in `.earlier_runs/` for comparison.