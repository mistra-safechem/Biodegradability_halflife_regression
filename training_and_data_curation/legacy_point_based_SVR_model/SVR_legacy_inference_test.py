"""
# SVR Inference (legacy for point based model)
2026-02-19 SVR model loading and inference
responsibility cleanup and database check added
database check used to retrieve if compounds are in database since this model
did not create/store training data for such a check;
(vs how it is done to the new SVR model).

NOTE: this script is a simpler version of latter implemented one.
A list of three SMILES strings is hardcoded in the test data.

to run:
uv run SVR_legacy_inference_test.py --compartment water --data-source hsbd

no timestamp needs to be specified since the script will automatically load
the latest model file based on timestamp in the filename.
This is to simplify testing and avoid having to update the timestamp in the command every time a new model is trained.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

script_dir = Path(__file__).resolve().parent
src_dir = script_dir.parent
sys.path.append(str(src_dir))

import joblib
import numpy as np
import pandas as pd
import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker
from src.db_utils import get_selected_data
from src.legacy.log_utils import log_section, log_to_file
from src.legacy.ml_tools import scale_features
from src.rdkit_tools import calculate_descriptors

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 250)


# Use script directory for all relative paths
data_dir = src_dir / "processed_data"
model_dir = script_dir / "models"
# log_file - defined in main()

DATABASE_FILES = {
    "hsbd": "hsbd_t_half_all.db",
    "vega": "vega_t_half_soil_water_sediment.db",
    "combined": "combined_t_half_vega_hsbd_soil_water_sediment.db",
}


def load_model_and_params(model_name: Path) -> tuple[Any, dict]:
    # Function to load model and parameters
    """Load trained model and its parameters for inference"""
    model_path = model_name
    params_path = model_path.with_suffix(".json")

    # Load model
    model = joblib.load(model_path)

    # Load parameters
    with open(params_path, "r") as f:
        params_dict = json.load(f)

    print(f"Model loaded from: {model_path}")
    print(f"Parameters loaded from: {params_path}")

    return model, params_dict


def extract_scaled_features(input_df: pd.DataFrame, svr_params: dict):

    # Check which features are available
    print(f"Total features calculated: {len(input_df.columns) - 1}")  # -1 for SMILES column
    print(f"Model requires: {len(svr_params['feature_columns'])} features")

    # Scale features (entire descriptor set)
    scaled_features_all = scale_features(input_df.drop(columns=["Canonical_smiles"]))

    # Check for missing features
    required_features = svr_params["feature_columns"]
    missing_features = set(required_features) - set(scaled_features_all.columns)
    if missing_features:
        print(f"\nWARNING: Missing features: {missing_features}")
    else:
        print("\nAll required features are present.")

    # Select only the features used in training, in the correct order
    # This is critical - must match the exact order used during training
    scaled_features = scaled_features_all[required_features].copy()

    # Ensure index is reset to avoid any index-related issues
    scaled_features = scaled_features.reset_index(drop=True)

    print(f"Features prepared for inference: {scaled_features.shape}")
    print(f"Feature names match: {list(scaled_features.columns) == required_features}")

    return scaled_features


def predict_with_uncertainty(svr_model, scaled_features, svr_params):
    """
    Perform inference with the SVR model on scaled features and calculate uncertainty bounds based on CV RMSE.
    """
    # Perform inference - pass DataFrame with feature names
    predictions_log10 = svr_model.predict(scaled_features)

    # Convert from log10 to original scale (days)
    predictions_days = 10**predictions_log10

    # Get uncertainty estimate from CV RMSE (in log10 space)
    cv_rmse_log10 = svr_params["model_performance_scores"]["cv_rmse_log10_mean"]

    # Calculate uncertainty bounds in original scale (days)
    # Log10(T_half) +/- RMSE -> convert to days
    lower_bound_log10 = predictions_log10 - cv_rmse_log10
    upper_bound_log10 = predictions_log10 + cv_rmse_log10
    lower_bound_days = 10**lower_bound_log10
    upper_bound_days = 10**upper_bound_log10

    return (
        predictions_log10,
        predictions_days,
        lower_bound_days,
        upper_bound_days,
        cv_rmse_log10,
    )


def get_thalf_from_db(smiles: str, data_source, compartment, session) -> float | None:
    """Retrieve T_half value from the database for a given SMILES string."""
    table = f"{compartment}_data"
    result = session.execute(
        sa.select(sa.text("T_half_days")).select_from(sa.text(table)).where(sa.text("Canonical_smiles = :smiles")),
        {"smiles": smiles},
    ).fetchone()
    if result is not None:
        print(f"SMILES: {smiles} - T_half in database: {result[0]:.3f} days")
        return result[0]
    else:
        print(f"SMILES: {smiles} - Not found in database.")
        return None


def is_mol_in_db(smiles: str, data_source, compartment, session) -> float | None:
    """Check if a molecule with the given SMILES string is present in the database."""
    table = f"{compartment}_data"
    result = session.execute(
        sa.select(sa.func.count()).select_from(sa.text(table)).where(sa.text("Canonical_smiles = :smiles")),
        {"smiles": smiles},
    ).scalar()
    thalf: float | None = None
    if result is not None and result > 0:
        # retrieve where in db the molecule is found (for logging purposes)
        found_in = session.execute(
            sa.select(sa.text("id")).select_from(sa.text(table)).where(sa.text("Canonical_smiles = :smiles")),
            {"smiles": smiles},
        ).fetchall()
        print(f"SMILES: {smiles} - Found in {table} (IDs: {[row[0] for row in found_in]})")
        thalf = get_thalf_from_db(smiles, data_source, compartment, session)
    return None if result is None else thalf


def main(compartment: str, data_source: str):
    # ----------------------------------------------------------------
    # var definitions and validation
    if data_source == "vega" and compartment == "air":
        raise ValueError("VEGA dataset does not contain air compartment data.")

    # set directories and filenames, load database
    working_dir = src_dir
    # Database setup based on data source
    if data_source not in DATABASE_FILES:
        raise ValueError(f"Unknown data source: {data_source}")
    db_file = data_dir / DATABASE_FILES[data_source]

    ENGINE = sa.create_engine(f"sqlite:///{db_file}")
    Session = sessionmaker(bind=ENGINE)
    data_to_use = get_selected_data(compartment, Session)

    # logging setup (since we are testing, just write to same path for now)
    log_file = f"SVR_inference_{compartment}_{data_source}.log"
    if Path(log_file).exists():
        Path(log_file).unlink()  # remove existing log file to start fresh
    log_to_file("SVR Inference", log_file)

    # model loading
    model_name = f"SVR_{compartment}_{data_source}"
    # get latest model file (in case there are multiple with different timestamps)
    model_files = list(model_dir.glob(f"{model_name}_*.joblib"))
    if not model_files:
        raise FileNotFoundError(f"No model files found for {model_name} in {model_dir}")
    # Get the latest model file based on timestamp in the filename
    latest_model_file = max(model_files, key=lambda x: x.stat().st_mtime)
    svr_model, svr_params = load_model_and_params(model_name=latest_model_file)
    print(f"\nLoaded model: {model_name} from file: {latest_model_file}\n")

    log_to_file(
        f"Model and parameters for {compartment.capitalize()} {data_source.capitalize()} SVR loaded successfully.",
        log_file,
    )
    # ----------------------------------------------------------------

    # inference test - create dummy data and run through the same descriptor calculation and scaling steps as training
    # Create test data
    df = pd.DataFrame({"Canonical_smiles": ["CCO", "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl", "CCC"]})

    # check if test molecules are in db and retrieve T_half if available
    for smiles in df["Canonical_smiles"]:
        thalf = is_mol_in_db(smiles, data_source, compartment, Session())
    print("\n")

    # Calculate descriptors and prepare features for inference
    df_with_descriptors = calculate_descriptors(df)
    scaled_features = extract_scaled_features(input_df=df_with_descriptors, svr_params=svr_params)

    # Perform inference and calculate uncertainty bounds
    (
        predictions_log10,
        predictions_days,
        lower_bound_days,
        upper_bound_days,
        cv_rmse_log10,
    ) = predict_with_uncertainty(svr_model, scaled_features, svr_params)

    # Display results with uncertainty
    print("\nInference Results:")
    log_section("Inference Results", log_file)
    print("=" * 80)
    for i, (smiles, pred_log, pred_days, lower, upper) in enumerate(
        zip(
            df["Canonical_smiles"],
            predictions_log10,
            predictions_days,
            lower_bound_days,
            upper_bound_days,
        )
    ):
        uncertainty = (upper - lower) / 2  # Average uncertainty in days
        print(f"Molecule {i + 1}: {smiles}")
        print(f"  Predicted log10(T_half): {pred_log:.4f} ± {cv_rmse_log10:.4f}")
        print(f"  Predicted T_half: {pred_days:.2f} days (95% CI: {lower:.2f} - {upper:.2f} days)")
        print(f"  Uncertainty: ± {uncertainty:.2f} days")
        print()

        log_to_file(
            f"Molecule {i + 1}: {smiles}\n"
            f"  Predicted log10(T_half): {pred_log:.4f} ± {cv_rmse_log10:.4f}\n"
            f"  Predicted T_half: {pred_days:.2f} days (95% CI: {lower:.2f} - {upper:.2f} days)\n"
            f"  Uncertainty: ± {uncertainty:.2f} days\n",
            log_file,
        )

    print(f"Note: Uncertainty based on {svr_params['cross-validation_folds']}-fold CV RMSE = {cv_rmse_log10:.4f} (log10 space)")
    log_to_file(
        f"Note: Uncertainty based on {svr_params['cross-validation_folds']}-fold CV RMSE = {cv_rmse_log10:.4f} (log10 space)\n",
        log_file,
    )
    print(f"\nLogfile saved as: {log_file}")


if __name__ == "__main__":
    # maybe not so clean to have all this here, but it keeps the main() signature clean and focused on the core logic
    parser = argparse.ArgumentParser(description="SVR inference for biodegradability")
    parser.add_argument(
        "--compartment",
        required=True,
        choices=["air", "water", "sediment", "soil"],
        help="Compartment to analyze",
    )
    parser.add_argument(
        "--data-source",
        default="hsbd",
        choices=["hsbd", "vega"],
        help="Data source: hsbd or vega (default: hsbd)",
    )
    args = parser.parse_args()
    compartment = args.compartment
    data_source = args.data_source

    main(compartment, data_source)
