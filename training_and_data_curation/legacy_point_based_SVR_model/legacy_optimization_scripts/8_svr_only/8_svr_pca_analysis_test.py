"""
Further SVR Analysis – Chemical Space Coverage
2026-04-01 first version

Simplest approach: PCA on the training feature matrix.
- Fit PCA (2 components) on X_train.
- Project X_train, X_test, and (optionally) new query compounds into the same 2-D space.
- Plot the coverage to see how well the training set spans the chemical space and
  whether the test set / new compounds fall inside or outside the training domain.
"""

# paths are not correct for the subfolder - check the other scripts for how to set up paths and db access, and adapt as needed
import matplotlib

matplotlib.use("Agg")

import json
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sqlalchemy as sa
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sqlalchemy.orm import sessionmaker
from training_and_data_curation.src.db_schema import *
from training_and_data_curation.src.db_utils import get_all_data
from training_and_data_curation.src.log_utils import log_section, log_to_file
from training_and_data_curation.src.ml_tools import PreprocessedData, Preprocessor, chemical_space_pca, t_t_split

pd.set_option("display.max_columns", None)

SCRIPT_DIR = Path(__file__).parent.resolve()
DATA_DIR = SCRIPT_DIR / "processed_data"
DATABASE_FILE = DATA_DIR / "t_half_all.db"
ENGINE = sa.create_engine(f"sqlite:///{DATABASE_FILE}")
Session = sessionmaker(bind=ENGINE)


def main():
    compartment = "air"
    target_column = "T_half_days"
    use_outlier_removal = True
    test_size = 0.25
    random_state = 42
    to_drop = ["None"]

    working_dir = SCRIPT_DIR
    logs_dir = working_dir / "logs" / compartment / "chemical_space"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = str(logs_dir / "9_further_SVR_analysis.log")

    log_to_file("9 – Further SVR Analysis: Chemical Space Coverage", log_file)
    log_section("PCA Chemical Space", log_file)

    # ------------------------------------------------------------------ #
    # Load & preprocess (same pipeline as 8_svr_model_and_analysis.py)
    # ------------------------------------------------------------------ #
    data_to_use = get_all_data(compartment, Session)
    preprocessor = Preprocessor(target_column=target_column, to_drop=to_drop, remove_outliers=use_outlier_removal)
    data_clean = preprocessor.preprocess(data_to_use)
    preprocessed_data = PreprocessedData(
        name="AirData", df=data_clean, remove_outliers=use_outlier_removal, X=preprocessor.X, y_log=preprocessor.y_log
    )

    X = preprocessed_data.X
    y = preprocessed_data.y_log
    X_train, X_test, y_train, y_test = t_t_split(X, y, test_size=test_size, random_state=random_state)

    log_to_file(f"n_train={len(X_train)}, n_test={len(X_test)}, n_features={X_train.shape[1]}", log_file)

    # ------------------------------------------------------------------ #
    # Chemical space coverage via PCA
    # ------------------------------------------------------------------ #
    pca, scaler, train_pcs, test_pcs, coverage_pct = chemical_space_pca(X_train, X_test, logs_dir, compartment)
    log_to_file(f"Test coverage inside training bounding box: {coverage_pct:.1f}%", log_file)
    log_to_file(
        f"PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.3f}, PC2={pca.explained_variance_ratio_[1]:.3f}",
        log_file,
    )

    print("\nDone. Check logs/air/chemical_space/ for output.")


if __name__ == "__main__":
    main()

# -- ignore --# -- ignore --# -- ignore --# -- ignore --# -- ignore --
# -- ignore --# -- ignore --# -- ignore --# -- ignore --# -- ignore --
