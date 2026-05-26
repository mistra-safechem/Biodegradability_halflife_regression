"""
SVR in depth analysis
2025-11-05 first version; Alexander Minidis

 -Support vector with optimization
 -IsolationForrest for the Outlier detection instead of LocalOutlier
 -drop constant features AFTER splitting;
 -proper x-validation, not split size, added, as well as a learning curve for general analysis quality
 -most preprocessing & model functions moved to ml_tools
 -model saving added
 -change to script format with logging to file instead of notebook cells; more structured and documented for reproducibility and future reference
"""

# paths are not correct for the subfolder - check the other scripts for how to set up paths and db access, and adapt as needed
import json
import sys

import joblib
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for script
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sqlalchemy as sa
from scipy import stats
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_validate, learning_curve
from sklearn.svm import SVR
from sqlalchemy.orm import sessionmaker

SCRIPT_DIR = Path(__file__).resolve().parent
TRAINING_DIR = SCRIPT_DIR.parents[2]
sys.path.append(str(TRAINING_DIR))

from src.db_schema import *
from src.db_utils import get_all_data
from src.legacy.log_utils import log_section, log_to_file
from src.legacy.ml_tools import (
    PreprocessedData,
    Preprocessor,
    output_metrics,
    output_metrics_w_return,
    svr_grid_search,
    t_t_split,
)
from src.rdkit_tools import MACCS_NAMES

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 250)

# Use script directory for all relative paths
DATA_DIR = TRAINING_DIR / "processed_data"
DATABASE_FILE = DATA_DIR / "hsbd_t_half_all.db"
ENGINE = sa.create_engine(f"sqlite:///{DATABASE_FILE}")
Session = sessionmaker(bind=ENGINE)


def main():
    # set directories and filenames, load database
    working_dir = SCRIPT_DIR
    compartment = "air"
    data_to_use = get_all_data(compartment, Session)

    # model settings
    target_column = "T_half_days"
    use_outlier_removal = True

    # logging setup
    log_time = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    logs_dir = working_dir / "logs" / compartment / log_time
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = str(logs_dir / f"8_svr_model_and_analysis.log")
    log_to_file("SVR Model and Analysis", log_file)
    log_section("Preprocessor class", log_file)

    # ## Preprocessing
    log_to_file(f"Dataset X and y shapes before preprocessing: X={data_to_use.shape[0]}, y={data_to_use.shape[0]}", log_file)
    to_drop = ["None"]
    log_to_file(f"Columns to drop: {to_drop}", log_file)

    preprocessor = Preprocessor(target_column=target_column, to_drop=to_drop, remove_outliers=use_outlier_removal)
    data_clean = preprocessor.preprocess(data_to_use)
    preprocessed_data = PreprocessedData(
        name="AirData", df=data_clean, remove_outliers=use_outlier_removal, X=preprocessor.X, y_log=preprocessor.y_log
    )
    for key, value in preprocessed_data.__dict__.items():
        if key == "name":
            log_to_file(f"{key}: {value}", log_file)
        if key == "remove_outliers":
            log_to_file(f"{key}: {value}", log_file)
    log_to_file(
        f"Dataset X and y shapes after preprocessing: X={preprocessed_data.X.shape[0]}, y_log={preprocessed_data.y_log.shape[0]}",
        log_file,
    )

    # ## Model
    X = preprocessed_data.X
    y = preprocessed_data.y_log
    log_section("Model section", log_file)

    # k-fold cross-validation for robust model evaluation
    # Configuration
    n_folds = 5
    test_size = 0.25
    random_state = 42

    print(f"Performing {n_folds}-fold cross-validation...")
    log_to_file(f"Performing {n_folds}-fold cross-validation", log_file, include_timestamp=True)

    # First do a train-test split to get final test set for model saving
    X_train, X_test, y_train, y_test = t_t_split(X, y, test_size=test_size, random_state=random_state)

    # Perform grid search on training data
    svr, svr_params = svr_grid_search(X_train, y_train)
    print(f"\nBest hyperparameters from grid search: {svr_params}")
    log_to_file(f"Best hyperparameters from grid search: {svr_params}", log_file, include_timestamp=True)

    # Cross-validation on training set (nested CV for robust evaluation)
    cv_scores = cross_validate(
        svr,
        X_train,
        y_train,
        cv=n_folds,
        scoring=["r2", "neg_mean_squared_error", "neg_mean_absolute_error"],
        return_train_score=True,
    )

    print(f"\n{n_folds}-Fold Cross-Validation Results (on training set):")
    print(f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})")
    print(
        f"  MSE (CV): {-cv_scores['test_neg_mean_squared_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_squared_error'].std():.4f})"
    )
    print(
        f"  MAE (CV): {-cv_scores['test_neg_mean_absolute_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_absolute_error'].std():.4f})"
    )

    log_to_file(f"\n{n_folds}-Fold Cross-Validation Results:", log_file, include_timestamp=True)
    log_to_file(f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})", log_file)
    log_to_file(
        f"  MSE (CV): {-cv_scores['test_neg_mean_squared_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_squared_error'].std():.4f})",
        log_file,
    )
    log_to_file(
        f"  MAE (CV): {-cv_scores['test_neg_mean_absolute_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_absolute_error'].std():.4f})",
        log_file,
    )

    # Train final model on full training set and evaluate on held-out test set
    svr.fit(X_train, y_train)
    y_pred_svr = svr.predict(X_test)

    # Inverse transform predictions and targets
    y_test_exp_svr = np.power(10, y_test)
    y_pred_exp_svr = np.power(10, y_pred_svr)

    print(f"\nFinal model performance on held-out test set (test_size={test_size}):")
    output_metrics(y_test_exp_svr, y_pred_exp_svr)

    # Store results in compatible format for next cell
    results = {
        test_size: {
            "y_test_exp_svr": y_test_exp_svr,
            "y_pred_exp_svr": y_pred_exp_svr,
            "params": svr_params,
            "cv_scores": cv_scores,
        }
    }

    # Extract results for logging and analysis
    # Since we now use proper CV, we have a single test_size with robust evaluation
    best_test_size = test_size
    params = results[best_test_size]["params"]
    cv_scores = results[best_test_size]["cv_scores"]

    print(f"\n{'=' * 60}")
    print(f"FINAL RESULTS SUMMARY")
    print(f"{'=' * 60}")
    print(f"Test size: {best_test_size}")
    print(f"SVR hyperparameters: {params}")
    print(f"\nCross-validation performance (training set, {n_folds} folds):")
    print(f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})")
    print(f"\nHeld-out test set performance:")

    # Extract predictions for analysis
    y_test_exp_svr = results[best_test_size]["y_test_exp_svr"]
    y_pred_exp_svr = results[best_test_size]["y_pred_exp_svr"]
    output_metrics(y_test_exp_svr, y_pred_exp_svr)
    print(f"{'=' * 60}")

    # Log to file
    log_to_file(f"\n{'=' * 60}", log_file, include_timestamp=True)
    log_to_file(f"FINAL RESULTS SUMMARY", log_file)
    log_to_file(f"{'=' * 60}", log_file)
    log_to_file(f"Test size: {best_test_size}", log_file)
    log_to_file(f"SVR hyperparameters: {params}", log_file)
    log_to_file(f"\nCross-validation performance ({n_folds} folds):", log_file)
    log_to_file(f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})", log_file)
    log_to_file(f"\nHeld-out test set performance:", log_file)
    log_to_file(f"{output_metrics_w_return(y_test_exp_svr, y_pred_exp_svr)}", log_file)
    log_to_file(f"{'=' * 60}", log_file)

    # ## Learning Curve Analysis
    # Evaluate how model performance changes with training set size to determine if more data would improve the model.

    # Perform learning curve analysis to assess training data requirements
    print("Performing learning curve analysis...")
    log_to_file("\nLearning Curve Analysis", log_file, include_timestamp=True)

    # Define training sizes to evaluate (10% to 100% of training data)
    train_sizes = np.linspace(0.1, 1.0, 10)

    # Compute learning curves
    # Note: This uses the full dataset X, y (not split) and internally does CV
    train_sizes_abs, train_scores, val_scores = learning_curve(
        svr,
        X,
        y,
        train_sizes=train_sizes,
        cv=n_folds,
        scoring="r2",
        random_state=random_state,
        n_jobs=-1,  # Use all CPU cores for faster computation
    )

    # Calculate mean and std for plotting
    train_scores_mean = train_scores.mean(axis=1)
    train_scores_std = train_scores.std(axis=1)
    val_scores_mean = val_scores.mean(axis=1)
    val_scores_std = val_scores.std(axis=1)

    # Log results
    print(f"\nLearning Curve Results:")
    print(f"Training sizes: {train_sizes_abs}")
    print(f"Validation R² scores: {val_scores_mean}")
    print(f"\nFinal performance with {train_sizes_abs[-1]} training samples:")
    print(f"  Training R²: {train_scores_mean[-1]:.4f} (+/- {train_scores_std[-1]:.4f})")
    print(f"  Validation R²: {val_scores_mean[-1]:.4f} (+/- {val_scores_std[-1]:.4f})")

    log_to_file(f"Training sizes evaluated: {train_sizes_abs.tolist()}", log_file)
    log_to_file(f"Final performance with {train_sizes_abs[-1]} training samples:", log_file)
    log_to_file(f"  Training R²: {train_scores_mean[-1]:.4f} (+/- {train_scores_std[-1]:.4f})", log_file)
    log_to_file(f"  Validation R²: {val_scores_mean[-1]:.4f} (+/- {val_scores_std[-1]:.4f})", log_file)

    # Check if curves are converging (would more data help?)
    if len(val_scores_mean) >= 3:
        recent_improvement = val_scores_mean[-1] - val_scores_mean[-3]
        if recent_improvement > 0.01:
            interpretation = "Curves are still improving - more data may help"
        elif abs(train_scores_mean[-1] - val_scores_mean[-1]) > 0.1:
            interpretation = "High variance - model may be overfitting"
        else:
            interpretation = "Curves have converged - more data unlikely to help significantly"
        print(f"\nInterpretation: {interpretation}")
        log_to_file(f"Interpretation: {interpretation}", log_file)

    # Plot learning curves
    plt.figure(figsize=(10, 6))

    # Plot training scores
    plt.plot(train_sizes_abs, train_scores_mean, "o-", color="royalblue", label="Training score", linewidth=2, markersize=6)
    plt.fill_between(
        train_sizes_abs,
        train_scores_mean - train_scores_std,
        train_scores_mean + train_scores_std,
        alpha=0.2,
        color="royalblue",
    )

    # Plot validation scores
    plt.plot(
        train_sizes_abs, val_scores_mean, "o-", color="darkorange", label="Cross-validation score", linewidth=2, markersize=6
    )
    plt.fill_between(
        train_sizes_abs, val_scores_mean - val_scores_std, val_scores_mean + val_scores_std, alpha=0.2, color="darkorange"
    )

    plt.xlabel("Training Set Size", fontsize=12)
    plt.ylabel("R² Score", fontsize=12)
    plt.title("Learning Curve: SVR Model Performance vs Training Data Size", fontsize=14)
    plt.legend(loc="lower right", fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    learning_curve_path = logs_dir / f"learning_curve_{compartment}.png"
    plt.savefig(learning_curve_path)
    print(f"Learning curve plot saved to {learning_curve_path}")
    log_to_file(f"Learning curve plot saved to {learning_curve_path}", log_file)
    plt.close()

    # Also plot the gap between training and validation scores to visualize overfitting
    plt.figure(figsize=(10, 6))
    gap = train_scores_mean - val_scores_mean
    plt.plot(train_sizes_abs, gap, "o-", color="crimson", linewidth=2, markersize=6)
    plt.fill_between(train_sizes_abs, 0, gap, alpha=0.3, color="crimson")
    plt.xlabel("Training Set Size", fontsize=12)
    plt.ylabel("Training-Validation Gap (R²)", fontsize=12)
    plt.title("Overfitting Analysis: Gap Between Training and Validation Scores", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0.05, color="gray", linestyle="--", label="5% gap threshold")
    plt.legend(loc="upper right", fontsize=10)
    plt.tight_layout()
    overfit_path = logs_dir / f"overfitting_gap_{compartment}.png"
    plt.savefig(overfit_path)
    print(f"Overfitting gap plot saved to {overfit_path}")
    log_to_file(f"Overfitting gap plot saved to {overfit_path}", log_file)
    plt.close()

    # ## Model Saving
    # Save the best SVR model and its parameters for later inference and analysis
    # possibly a bit of reduncancy with some calculations, but ensures we have all metrics documented for this final model
    models_dir = working_dir / "models"
    models_dir.mkdir(exist_ok=True)

    # Determine compartment name from data being used
    model_path = models_dir / f"SVR_{compartment}.joblib"
    joblib.dump(svr, model_path)
    print(f"Model saved to {model_path}")
    log_to_file(f"Model saved to {model_path}", log_file, include_timestamp=True)

    # Calculate additional metrics for model performance documentation
    # CV metrics are already in log10 space
    cv_rmse_log10 = np.sqrt(-cv_scores["test_neg_mean_squared_error"].mean())
    cv_rmse_log10_std = np.sqrt(cv_scores["test_neg_mean_squared_error"].std())

    # Test set metrics in log10 space
    test_mse_log10 = mean_squared_error(y_test, y_pred_svr)
    test_rmse_log10 = np.sqrt(test_mse_log10)
    test_mae_log10 = mean_absolute_error(y_test, y_pred_svr)

    # Test set metrics in original space (days)
    test_rmse_days = np.sqrt(mean_squared_error(y_test_exp_svr, y_pred_exp_svr))
    test_mae_days = mean_absolute_error(y_test_exp_svr, y_pred_exp_svr)

    print(f"\nAdditional metrics for model documentation:")
    print(f"  CV RMSE (log10): {cv_rmse_log10:.4f} (+/- {cv_rmse_log10_std:.4f})")
    print(f"  Test RMSE (log10): {test_rmse_log10:.4f}")
    print(f"  Test RMSE (days): {test_rmse_days:.2f}")
    print(f"  Test MAE (log10): {test_mae_log10:.4f}")
    print(f"  Test MAE (days): {test_mae_days:.2f}")

    # Verify feature counts
    print(f"\nFeature information:")
    print(f"  Original features (X): {X.shape[1]}")
    print(f"  Training features (X_train): {X_train.shape[1]}")
    print(f"  Features removed as constants: {X.shape[1] - X_train.shape[1]}")

    # Create model parameters dictionary and export to JSON
    model_parameters = {
        "model_name": "SVR",
        "compartment": compartment,
        "target_column": target_column,
        "model_params": params,
        "feature_columns": list(X_train.columns),  # CRITICAL: Use X_train columns (after constant removal)
        "y_value_transformation": "log10",
        "test_size": best_test_size,
        "cross-validation_folds": n_folds,
        "model_performance_scores": {
            "cv_r2_mean": cv_scores["test_r2"].mean(),
            "cv_r2_std": cv_scores["test_r2"].std(),
            "cv_mse_mean": -cv_scores["test_neg_mean_squared_error"].mean(),
            "cv_mse_std": cv_scores["test_neg_mean_squared_error"].std(),
            "cv_rmse_log10_mean": cv_rmse_log10,
            "cv_rmse_log10_std": cv_rmse_log10_std,
            "cv_mae_mean": -cv_scores["test_neg_mean_absolute_error"].mean(),
            "cv_mae_std": cv_scores["test_neg_mean_absolute_error"].std(),
            "test_rmse_log10": test_rmse_log10,
            "test_mae_log10": test_mae_log10,
            "test_rmse_days": test_rmse_days,
            "test_mae_days": test_mae_days,
        },
        "outlier_removal": use_outlier_removal,
        "preprocessing_drops": to_drop,
        "n_samples_train": len(X_train),
        "n_samples_test": len(X_test),
        "n_features": X_train.shape[1],
    }

    # Save to JSON
    params_path = models_dir / f"SVR_{compartment}_params.json"
    with open(params_path, "w") as f:
        json.dump(model_parameters, f, indent=2)
    print(f"Parameters saved to {params_path}")
    log_to_file(f"Parameters saved to {params_path}", log_file, include_timestamp=True)

    # # Model Analysis
    # SVR true vs predicted plot
    mae = mean_absolute_error(y_test_exp_svr, y_pred_exp_svr)
    plt.figure(figsize=(7, 7))
    plt.scatter(y_test_exp_svr, y_pred_exp_svr, alpha=0.6)
    plt.plot([min(y_test_exp_svr), max(y_test_exp_svr)], [min(y_test_exp_svr), max(y_test_exp_svr)], "r--", label="Ideal")
    plt.fill_between(
        [min(y_test_exp_svr), max(y_test_exp_svr)],
        [min(y_test_exp_svr) + mae, max(y_test_exp_svr) + mae],
        [min(y_test_exp_svr) - mae, max(y_test_exp_svr) - mae],
        color="gray",
        alpha=0.2,
        label=f"MAE ±{mae:.2f}",
    )
    plt.xlabel("True y (days)")
    plt.ylabel("Predicted y (days)")
    plt.title("SVR: True vs Predicted y with MAE region")
    plt.legend()
    plt.grid(True)
    true_vs_pred_path = logs_dir / f"true_vs_pred_{compartment}.png"
    plt.savefig(true_vs_pred_path)
    print(f"True vs Predicted plot saved to {true_vs_pred_path}")
    log_to_file(f"True vs Predicted plot saved to {true_vs_pred_path}", log_file)
    plt.close()

    # Add shaded region for MAE to SVR true vs predicted plot, limited to x/y <= 100
    mae = mean_absolute_error(y_test_exp_svr, y_pred_exp_svr)
    plt.figure(figsize=(7, 7))
    # Only plot points where both true and predicted are <= 100
    mask = (y_test_exp_svr <= 100) & (y_pred_exp_svr <= 100)
    plt.scatter(y_test_exp_svr[mask], y_pred_exp_svr[mask], alpha=0.6)
    plt.plot([0, 100], [0, 100], "r--", label="Ideal")
    plt.fill_between([0, 100], [mae, 100 + mae], [-mae, 100 - mae], color="gray", alpha=0.2, label=f"MAE ±{mae:.2f}")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.xlabel("True y (days)")
    plt.ylabel("Predicted y (days)")
    plt.title("SVR: True vs Predicted y with MAE region (y ≤ 100)")
    plt.legend()
    plt.grid(True)
    true_vs_pred_100_path = logs_dir / f"true_vs_pred_100_{compartment}.png"
    plt.savefig(true_vs_pred_100_path)
    print(f"True vs Predicted (≤100) plot saved to {true_vs_pred_100_path}")
    log_to_file(f"True vs Predicted (≤100) plot saved to {true_vs_pred_100_path}", log_file)
    plt.close()

    # retrieve best features by permutation importance for SVR
    result = permutation_importance(svr, X_test, y_test, n_repeats=10, random_state=42, scoring="neg_mean_squared_error")

    # plot top 20 features by permutation importance for SVR
    importances = result.importances_mean
    num_of_feats = 20
    indices = np.argsort(importances)[::-1][:num_of_feats]
    plt.figure(figsize=(10, 5))
    plt.bar(range(num_of_feats), importances[indices])
    plt.xticks(range(num_of_feats), X_test.columns[indices], rotation=90)
    plt.title(f"Top {num_of_feats} SVR Feature Importances (Permutation)")
    plt.ylabel("Mean Importance")
    plt.tight_layout()
    feat_imp_path = logs_dir / f"feature_importance_{compartment}.png"
    plt.savefig(feat_imp_path)
    print(f"Feature importance plot saved to {feat_imp_path}")
    log_to_file(f"Feature importance plot saved to {feat_imp_path}", log_file)
    plt.close()

    # analysis of SVR

    residuals = y_test_exp_svr - y_pred_exp_svr

    # Scatter plot: Residuals vs Predicted
    plt.figure(figsize=(8, 5))
    sns.scatterplot(x=y_pred_exp_svr, y=residuals)
    plt.axhline(0, color="red", linestyle="--")
    plt.xlabel("Predicted Values")
    plt.ylabel("Residuals")
    plt.title("Residuals vs Predicted")
    resid_vs_pred_path = logs_dir / f"residuals_vs_pred_{compartment}.png"
    plt.savefig(resid_vs_pred_path)
    print(f"Residuals vs Predicted plot saved to {resid_vs_pred_path}")
    log_to_file(f"Residuals vs Predicted plot saved to {resid_vs_pred_path}", log_file)
    plt.close()

    # Histogram of residuals
    plt.figure()
    sns.histplot(residuals, kde=True)
    plt.title("Residual Distribution")
    resid_hist_path = logs_dir / f"residual_hist_{compartment}.png"
    plt.savefig(resid_hist_path)
    print(f"Residual histogram plot saved to {resid_hist_path}")
    log_to_file(f"Residual histogram plot saved to {resid_hist_path}", log_file)
    plt.close()

    # Q-Q Plot using SciPy
    plt.figure()
    stats.probplot(residuals, dist="norm", plot=plt)
    plt.title("Q-Q Plot of Residuals")
    qqplot_path = logs_dir / f"qqplot_residuals_{compartment}.png"
    plt.savefig(qqplot_path)
    print(f"Q-Q plot of residuals saved to {qqplot_path}")
    log_to_file(f"Q-Q plot of residuals saved to {qqplot_path}", log_file)
    plt.close()

    # Residuals vs Top Features by Permutation Importance
    for col in X_test.columns[indices]:
        plt.figure(figsize=(6, 4))
        sns.scatterplot(x=X_test[col], y=residuals)
        plt.axhline(0, color="red", linestyle="--")
        plt.title(f"Residuals vs {col}")
        resid_vs_feat_path = logs_dir / f"residuals_vs_{col}_{compartment}.png"
        plt.savefig(resid_vs_feat_path)
        print(f"Residuals vs {col} plot saved to {resid_vs_feat_path}")
        log_to_file(f"Residuals vs {col} plot saved to {resid_vs_feat_path}", log_file)
        plt.close()

    # Or, print all residuals along with the feature value
    for val, res in zip(X_test["fr_term_acetylene"], residuals):
        print(f"fr_term_acetylene: {val}, residual: {res}")


if __name__ == "__main__":
    print("SVR model training and analysis complete. Check log file for details.")
    main()

# -- IGNORE --
# -- IGNORE --
# -- IGNORE --
# -- IGNORE --
