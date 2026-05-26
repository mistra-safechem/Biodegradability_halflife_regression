"""
SVR in depth analysis
2025-11-05 first version; Alexander Minidis
2026-04-15 latest edits; Alexander Minidis
"""

# paths are not correct for the subfolder - check the other scripts for how to set up paths and db access, and adapt as needed

import argparse
import json
import sys
import warnings

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
from sklearn.decomposition import PCA as _PCA
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import RFE
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_validate, learning_curve
from sklearn.svm import SVR, LinearSVR
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
    applicability_domain_leverage,
    chemical_space_morgan_pca,
    chemical_space_pca,
    output_metrics,
    output_metrics_w_return,
    svr_grid_search,
    t_t_split,
)

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 250)

# Paths and DB setup
data_dir = TRAINING_DIR / "processed_data"
DATABASE_FILES = {
    "hsbd": "hsbd_t_half_all.db",
    "vega": "vega_t_half_soil_water_sediment.db",
}


# ----------------------------------------------------------------
# ## Shared SVR fit + CV
# ----------------------------------------------------------------
def fit_svr(X_train, y_train, svr_params, n_folds, log_file):
    """
    Fit SVR and run cross-validation on X_train/y_train.

    Used by both train_model() (full features) and run_feature_reduction()
    (reduced features).  Grid search is NOT performed here.

    Returns
    -------
    svr        : fitted SVR instance
    cv_scores  : dict from cross_validate
    """
    svr = SVR(**svr_params)
    svr.fit(X_train, y_train)

    cv_scores = cross_validate(
        SVR(**svr_params),
        X_train,
        y_train,
        cv=n_folds,
        scoring=["r2", "neg_mean_squared_error", "neg_mean_absolute_error"],
        return_train_score=True,
    )

    log_to_file(f"\n{n_folds}-Fold Cross-Validation Results:", log_file, include_timestamp=True)
    log_to_file(
        f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})",
        log_file,
    )
    log_to_file(
        f"  MSE (CV): {-cv_scores['test_neg_mean_squared_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_squared_error'].std():.4f})",
        log_file,
    )
    log_to_file(
        f"  MAE (CV): {-cv_scores['test_neg_mean_absolute_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_absolute_error'].std():.4f})",
        log_file,
    )

    return svr, cv_scores


# ----------------------------------------------------------------
# ## Model training
# ----------------------------------------------------------------
def train_model(compartment, data_source, preprocessed_data, preprocessor, logs_dir, log_file):
    """
    Perform grid search, cross-validation, final fit, and model/JSON saving.

    Returns a model_artifacts dict containing all objects needed by run_analysis().
    """
    X = preprocessed_data.X
    y = preprocessed_data.y_log
    log_section("Model section", log_file)

    # Configuration
    n_folds = 5
    test_size = 0.25
    random_state = 42

    print(f"Performing {n_folds}-fold cross-validation...")
    log_to_file(f"Performing {n_folds}-fold cross-validation", log_file, include_timestamp=True)

    # Train/test split
    X_train, X_test, y_train, y_test = t_t_split(X, y, test_size=test_size, random_state=random_state)

    # Extract SMILES aligned with the split
    smiles_train = preprocessor.smiles.loc[X_train.index]
    smiles_test = preprocessor.smiles.loc[X_test.index]

    # Grid search on training data
    svr, svr_params = svr_grid_search(X_train, y_train)
    print(f"\nBest hyperparameters from grid search: {svr_params}")
    log_to_file(
        f"Best hyperparameters from grid search: {svr_params}",
        log_file,
        include_timestamp=True,
    )

    # Fit + cross-validate
    svr, cv_scores = fit_svr(X_train, y_train, svr_params, n_folds, log_file)

    print(f"\n{n_folds}-Fold Cross-Validation Results (on training set):")
    print(f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})")
    print(
        f"  MSE (CV): {-cv_scores['test_neg_mean_squared_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_squared_error'].std():.4f})"
    )
    print(
        f"  MAE (CV): {-cv_scores['test_neg_mean_absolute_error'].mean():.4f} (+/- {cv_scores['test_neg_mean_absolute_error'].std():.4f})"
    )

    y_pred_svr = svr.predict(X_test)
    y_pred_train = svr.predict(X_train)

    # Inverse transform predictions and targets
    y_test_exp_svr = np.power(10, y_test)
    y_pred_exp_svr = np.power(10, y_pred_svr)

    print(f"\nFinal model performance on held-out test set (test_size={test_size}):")
    output_metrics(y_test_exp_svr, y_pred_exp_svr)

    # Summary
    print(f"\n{'=' * 60}")
    print(f"FINAL RESULTS SUMMARY")
    print(f"{'=' * 60}")
    print(f"Test size: {test_size}")
    print(f"SVR hyperparameters: {svr_params}")
    print(f"\nCross-validation performance (training set, {n_folds} folds):")
    print(f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})")
    print(f"\nHeld-out test set performance:")
    output_metrics(y_test_exp_svr, y_pred_exp_svr)
    print(f"{'=' * 60}")

    log_to_file(f"\n{'=' * 60}", log_file, include_timestamp=True)
    log_to_file(f"FINAL RESULTS SUMMARY", log_file)
    log_to_file(f"{'=' * 60}", log_file)
    log_to_file(f"Test size: {test_size}", log_file)
    log_to_file(f"SVR hyperparameters: {svr_params}", log_file)
    log_to_file(f"\nCross-validation performance ({n_folds} folds):", log_file)
    log_to_file(
        f"  R² (CV): {cv_scores['test_r2'].mean():.4f} (+/- {cv_scores['test_r2'].std():.4f})",
        log_file,
    )
    log_to_file(f"\nHeld-out test set performance:", log_file)
    log_to_file(f"{output_metrics_w_return(y_test_exp_svr, y_pred_exp_svr)}", log_file)
    log_to_file(f"{'=' * 60}", log_file)

    # Additional metrics (printed for reference; saved with final reduced model)
    cv_rmse_log10 = np.sqrt(-cv_scores["test_neg_mean_squared_error"].mean())
    cv_rmse_log10_std = np.sqrt(cv_scores["test_neg_mean_squared_error"].std())

    test_mse_log10 = mean_squared_error(y_test, y_pred_svr)
    test_rmse_log10 = np.sqrt(test_mse_log10)
    test_mae_log10 = mean_absolute_error(y_test, y_pred_svr)

    test_rmse_days = np.sqrt(mean_squared_error(y_test_exp_svr, y_pred_exp_svr))
    test_mae_days = mean_absolute_error(y_test_exp_svr, y_pred_exp_svr)

    print(f"\nAdditional metrics for model documentation:")
    print(f"  CV RMSE (log10): {cv_rmse_log10:.4f} (+/- {cv_rmse_log10_std:.4f})")
    print(f"  Test RMSE (log10): {test_rmse_log10:.4f}")
    print(f"  Test RMSE (days): {test_rmse_days:.2f}")
    print(f"  Test MAE (log10): {test_mae_log10:.4f}")
    print(f"  Test MAE (days): {test_mae_days:.2f}")

    print(f"\nFeature information:")
    print(f"  Original features (X): {X.shape[1]}")
    print(f"  Training features (X_train): {X_train.shape[1]}")
    print(f"  Features removed as constants: {X.shape[1] - X_train.shape[1]}")

    return {
        "svr": svr,
        "X": X,
        "y": y,
        "X_train": X_train,
        "X_test": X_test,
        "y_train": y_train,
        "y_test": y_test,
        "y_pred_svr": y_pred_svr,
        "y_pred_train": y_pred_train,
        "y_test_exp_svr": y_test_exp_svr,
        "y_pred_exp_svr": y_pred_exp_svr,
        "smiles_train": smiles_train,
        "smiles_test": smiles_test,
        "cv_scores": cv_scores,
        "params": svr_params,
        "n_folds": n_folds,
        "test_size": test_size,
        "random_state": random_state,
        "data_source": data_source,
        "preprocessed_data": preprocessed_data,
    }


# ----------------------------------------------------------------
# ## Analysis — Section 1 & 2: Basic analysis
# ----------------------------------------------------------------
def analyse_basic(compartment, model_artifacts, logs_dir, log_file):
    """
    Sections 1-2: learning curve + overfitting gap, true vs predicted plots.

    Gate question: enough data? Model generalising? Headline prediction quality.
    """
    svr = model_artifacts["svr"]
    X = model_artifacts["X"]
    y = model_artifacts["y"]
    y_test = model_artifacts["y_test"]
    y_pred_svr = model_artifacts["y_pred_svr"]
    y_test_exp_svr = model_artifacts["y_test_exp_svr"]
    y_pred_exp_svr = model_artifacts["y_pred_exp_svr"]
    n_folds = model_artifacts["n_folds"]
    random_state = model_artifacts["random_state"]

    # ----------------------------------------------------------------
    # ## 1. Learning curve + overfitting gap
    # Gate question: is there enough data and is the model generalising?
    print("Performing learning curve analysis...")
    log_to_file("\nLearning Curve Analysis", log_file, include_timestamp=True)

    train_sizes = np.linspace(0.1, 1.0, 10)

    train_sizes_abs, train_scores, val_scores = learning_curve(
        svr,
        X,
        y,
        train_sizes=train_sizes,
        cv=n_folds,
        scoring="r2",
        random_state=random_state,
        n_jobs=-1,
    )

    train_scores_mean = train_scores.mean(axis=1)
    train_scores_std = train_scores.std(axis=1)
    val_scores_mean = val_scores.mean(axis=1)
    val_scores_std = val_scores.std(axis=1)

    print(f"\nLearning Curve Results:")
    print(f"Training sizes: {train_sizes_abs}")
    print(f"Validation R² scores: {val_scores_mean}")
    print(f"\nFinal performance with {train_sizes_abs[-1]} training samples:")
    print(f"  Training R²: {train_scores_mean[-1]:.4f} (+/- {train_scores_std[-1]:.4f})")
    print(f"  Validation R²: {val_scores_mean[-1]:.4f} (+/- {val_scores_std[-1]:.4f})")

    log_to_file(f"Training sizes evaluated: {train_sizes_abs.tolist()}", log_file)
    log_to_file(f"Final performance with {train_sizes_abs[-1]} training samples:", log_file)
    log_to_file(
        f"  Training R²: {train_scores_mean[-1]:.4f} (+/- {train_scores_std[-1]:.4f})",
        log_file,
    )
    log_to_file(
        f"  Validation R²: {val_scores_mean[-1]:.4f} (+/- {val_scores_std[-1]:.4f})",
        log_file,
    )

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

    plt.figure(figsize=(10, 6))
    plt.plot(
        train_sizes_abs,
        train_scores_mean,
        "o-",
        color="royalblue",
        label="Training score",
        linewidth=2,
        markersize=6,
    )
    plt.fill_between(
        train_sizes_abs,
        train_scores_mean - train_scores_std,
        train_scores_mean + train_scores_std,
        alpha=0.2,
        color="royalblue",
    )
    plt.plot(
        train_sizes_abs,
        val_scores_mean,
        "o-",
        color="darkorange",
        label="Cross-validation score",
        linewidth=2,
        markersize=6,
    )
    plt.fill_between(
        train_sizes_abs,
        val_scores_mean - val_scores_std,
        val_scores_mean + val_scores_std,
        alpha=0.2,
        color="darkorange",
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

    # ----------------------------------------------------------------
    # ## 2. True vs Predicted plots
    # Headline prediction quality on the held-out test set.
    mae = mean_absolute_error(y_test_exp_svr, y_pred_exp_svr)
    plt.figure(figsize=(7, 7))
    plt.scatter(y_test_exp_svr, y_pred_exp_svr, alpha=0.6)
    plt.plot(
        [min(y_test_exp_svr), max(y_test_exp_svr)],
        [min(y_test_exp_svr), max(y_test_exp_svr)],
        "r--",
        label="Ideal",
    )
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

    # True vs Predicted limited to y <= 100
    plt.figure(figsize=(7, 7))
    mask = (y_test_exp_svr <= 100) & (y_pred_exp_svr <= 100)
    plt.scatter(y_test_exp_svr[mask], y_pred_exp_svr[mask], alpha=0.6)
    plt.plot([0, 100], [0, 100], "r--", label="Ideal")
    plt.fill_between(
        [0, 100],
        [mae, 100 + mae],
        [-mae, 100 - mae],
        color="gray",
        alpha=0.2,
        label=f"MAE ±{mae:.2f}",
    )
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


# ----------------------------------------------------------------
# ## Analysis — Section 3 & 4: Coverage and domain
# ----------------------------------------------------------------
def analyse_coverage_domain(compartment, model_artifacts, logs_dir, log_file):
    """
    Sections 3-4: PCA chemical space coverage, applicability domain.

    Confirms test set sits inside training manifold; quantifies AD membership.
    """
    X_train = model_artifacts["X_train"]
    X_test = model_artifacts["X_test"]
    y_train = model_artifacts["y_train"]
    y_test = model_artifacts["y_test"]
    y_pred_svr = model_artifacts["y_pred_svr"]
    y_pred_train = model_artifacts["y_pred_train"]
    smiles_train = model_artifacts["smiles_train"]
    smiles_test = model_artifacts["smiles_test"]

    # ----------------------------------------------------------------
    # ## 3. PCA chemical space coverage
    # Confirm the test set sits inside the training descriptor manifold before
    # interpreting any errors — context for all residual plots that follow.
    pca, scaler, train_pcs, test_pcs, coverage_pct = chemical_space_pca(X_train, X_test, logs_dir, compartment)
    log_to_file(f"Test coverage inside training bounding box: {coverage_pct:.1f}%", log_file)
    log_to_file(
        f"PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.3f}, PC2={pca.explained_variance_ratio_[1]:.3f}",
        log_file,
    )

    # ----------------------------------------------------------------
    # ## 4. Applicability domain (Williams plot + Butina structural clustering)
    # Quantifies which test compounds are outside the AD; pairs with PCA above
    # (descriptor space → structural space, same question at two resolutions).
    ad_results = applicability_domain_leverage(
        X_train,
        X_test,
        y_train,
        y_test,
        y_pred_train,
        y_pred_svr,
        logs_dir / "ad_analysis",
        compartment,
    )
    log_section("Applicability Domain Analysis", log_file)
    log_to_file(f"h* threshold : {ad_results['h_star']:.4f}", log_file)
    log_to_file(f"Train inside AD : {ad_results['pct_train_inside']:.1f}%", log_file)
    log_to_file(f"Test  inside AD : {ad_results['pct_test_inside']:.1f}%", log_file)

    # Structural chemical space coloured by AD membership
    cluster_stats = chemical_space_morgan_pca(smiles_train, smiles_test, ad_results, logs_dir / "ad_analysis", compartment)
    log_to_file(f"AD-outside test compounds (n={cluster_stats['n_outside_ad']})", log_file)
    if "n_clusters" in cluster_stats:
        log_to_file(
            f"Butina clusters: {cluster_stats['n_clusters']}, "
            f"singletons: {cluster_stats['n_singletons']} ({cluster_stats['singleton_pct']:.1f}%)",
            log_file,
        )

    return ad_results


# ----------------------------------------------------------------
# ## Entry point
# ----------------------------------------------------------------
def analyse_features(compartment, model_artifacts, logs_dir, log_file, suffix=""):
    """
    Sections 5-10: permutation importance, residual analysis, feature reduction.

    Requires model_artifacts from train_model(). Uses indices from permutation
    importance (section 5) in per-feature residual plots (section 9).
    """
    svr = model_artifacts["svr"]
    X_train = model_artifacts["X_train"]
    X_test = model_artifacts["X_test"]
    y_test = model_artifacts["y_test"]
    y_pred_svr = model_artifacts["y_pred_svr"]

    # ----------------------------------------------------------------
    # ## 5. Permutation importance
    # Which features drive predictions globally — sets up the per-feature
    # residual plots in section 9.
    result = permutation_importance(
        svr,
        X_test,
        y_test,
        n_repeats=10,
        random_state=42,
        scoring="neg_mean_squared_error",
    )

    importances = result.importances_mean
    num_of_feats = 15
    indices = np.argsort(importances)[::-1]  # all features, descending importance

    plt.figure(figsize=(10, 5))
    plt.bar(range(num_of_feats), importances[indices[:num_of_feats]])
    plt.xticks(range(num_of_feats), X_test.columns[indices[:num_of_feats]], rotation=90)
    plt.title(f"Top {num_of_feats} SVR Feature Importances (Permutation)")
    plt.ylabel("Mean Importance")
    plt.tight_layout()
    feat_imp_path = logs_dir / f"feature_importance_{compartment}{suffix}.png"
    plt.savefig(feat_imp_path)
    print(f"Feature importance plot saved to {feat_imp_path}")
    log_to_file(f"Feature importance plot saved to {feat_imp_path}", log_file)
    plt.close()

    # ----------------------------------------------------------------
    # ## 6–8. Residual analysis (log10 space)
    # Both y_test and y_pred_svr are already in log10 space.
    residuals_log10 = y_test - y_pred_svr

    # Residual outputs go into a dedicated subfolder
    residuals_dir = logs_dir / "residuals"
    residuals_dir.mkdir(exist_ok=True)

    # --- 6. Scatter: Residuals vs Predicted ---
    # Model-level systematic bias check.
    res_vals_all = residuals_log10.values
    pred_vals_all = y_pred_svr.values if hasattr(y_pred_svr, "values") else y_pred_svr
    sp_r_pred, sp_p_pred = stats.spearmanr(pred_vals_all, res_vals_all)
    if abs(sp_r_pred) >= 0.3 and sp_p_pred < 0.05:
        direction_pred = "positive" if sp_r_pred > 0 else "negative"
        pred_comment = f"monotonic {direction_pred} trend (rho={sp_r_pred:.2f}, p={sp_p_pred:.3f})"
    else:
        pred_comment = f"no significant monotonic trend (rho={sp_r_pred:.2f}, p={sp_p_pred:.3f})"
    log_to_file(f"Residuals vs Predicted: {pred_comment}", log_file)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(pred_vals_all, res_vals_all, alpha=0.6, s=20, color="steelblue")
    ax.axhline(0, color="red", linestyle="--", linewidth=1)
    sort_idx_pred = np.argsort(pred_vals_all)
    window_pred = max(3, len(pred_vals_all) // 10)
    pv_sorted = pred_vals_all[sort_idx_pred]
    rv_sorted_pred = res_vals_all[sort_idx_pred]
    rv_smooth_pred = pd.Series(rv_sorted_pred).rolling(window_pred, center=True, min_periods=1).mean().values
    ax.plot(
        pv_sorted,
        rv_smooth_pred,
        color="darkorange",
        linewidth=1.5,
        label="rolling mean",
    )
    ax.set_xlabel("Predicted (log10 days)", fontsize=11)
    ax.set_ylabel("Residual (log10 days)", fontsize=11)
    ax.set_title(f"Residuals vs Predicted\n{pred_comment}", fontsize=9)
    ax.legend(fontsize=8)
    plt.tight_layout()
    resid_vs_pred_path = residuals_dir / f"residuals_vs_pred_{compartment}{suffix}.png"
    plt.savefig(resid_vs_pred_path, dpi=150)
    plt.close()
    print(f"Residuals vs Predicted plot saved to {resid_vs_pred_path}")
    log_to_file(f"Residuals vs Predicted plot saved to {resid_vs_pred_path}", log_file)

    # --- 7. Q-Q plot ---
    # Normality of residuals (distributional assumption check).
    plt.figure()
    stats.probplot(residuals_log10, dist="norm", plot=plt)
    plt.title("Q-Q Plot of Residuals (log10 days)")
    qqplot_path = residuals_dir / f"qqplot_residuals_{compartment}{suffix}.png"
    plt.savefig(qqplot_path)
    print(f"Q-Q plot of residuals saved to {qqplot_path}")
    log_to_file(f"Q-Q plot of residuals saved to {qqplot_path}", log_file)
    plt.close()

    # --- 8. Histogram of residuals ---
    # Complements Q-Q; shows shape and skew.
    plt.figure()
    sns.histplot(residuals_log10, kde=True)
    plt.xlabel("Residual (log10 days)")
    plt.title("Residual Distribution (log10 days)")
    resid_hist_path = residuals_dir / f"residual_hist_{compartment}{suffix}.png"
    plt.savefig(resid_hist_path)
    print(f"Residual histogram plot saved to {resid_hist_path}")
    log_to_file(f"Residual histogram plot saved to {resid_hist_path}", log_file)
    plt.close()

    # --- 9. Residuals vs Top Features (top 10 by permutation importance) ---
    # Feature-level error structure; meaningful now that we know which features
    # matter (section 5) and which compounds are in/out of AD (section 4).
    top_feat_indices = indices[:10]
    feat_stats = {}

    for col in X_test.columns[top_feat_indices]:
        feat_vals = X_test[col].values
        res_vals = residuals_log10.values

        spearman_r, spearman_p = stats.spearmanr(feat_vals, res_vals)

        median_val = np.median(feat_vals)
        group_low = res_vals[feat_vals <= median_val]
        group_high = res_vals[feat_vals > median_val]
        if len(group_low) >= 2 and len(group_high) >= 2:
            levene_stat, levene_p = stats.levene(group_low, group_high)
            hetero_flag = levene_p < 0.05
        else:
            levene_stat, levene_p = float("nan"), float("nan")
            hetero_flag = False

        trend_sig = abs(spearman_r) >= 0.3 and spearman_p < 0.05
        if trend_sig:
            direction = "positive" if spearman_r > 0 else "negative"
            trend_str = f"monotonic {direction} trend (rho={spearman_r:.2f}, p={spearman_p:.3f})"
        else:
            trend_str = f"no significant monotonic trend (rho={spearman_r:.2f}, p={spearman_p:.3f})"
        hetero_str = f"heteroscedastic (Levene p={levene_p:.3f})" if hetero_flag else f"homoscedastic (Levene p={levene_p:.3f})"
        comment_str = f"{trend_str}; {hetero_str}"

        feat_stats[col] = {
            "spearman_r": spearman_r,
            "spearman_p": spearman_p,
            "trend_sig": trend_sig,
            "hetero_flag": hetero_flag,
            "levene_p": levene_p,
        }

        log_to_file(f"  {col}: {comment_str}", log_file)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(feat_vals, res_vals, alpha=0.5, s=20, color="steelblue")
        ax.axhline(0, color="red", linestyle="--", linewidth=1)
        sort_idx = np.argsort(feat_vals)
        window = max(3, len(feat_vals) // 10)
        fv_sorted = feat_vals[sort_idx]
        rv_sorted = res_vals[sort_idx]
        rv_smooth = pd.Series(rv_sorted).rolling(window, center=True, min_periods=1).mean().values
        ax.plot(
            fv_sorted,
            rv_smooth,
            color="darkorange",
            linewidth=1.5,
            label="rolling mean",
        )
        ax.set_xlabel(col, fontsize=11)
        ax.set_ylabel("Residual (log10 days)", fontsize=11)
        ax.set_title(f"Residuals vs {col}\n{comment_str}", fontsize=9)
        ax.legend(fontsize=8)
        plt.tight_layout()
        resid_vs_feat_path = residuals_dir / f"residuals_vs_{col}_{compartment}{suffix}.png"
        plt.savefig(resid_vs_feat_path, dpi=150)
        plt.close()
        log_to_file(f"  Plot saved: {resid_vs_feat_path}", log_file)

    # --- 10. Feature reduction recommendations ---
    # Verdict and forward path — summarises all diagnostics above.
    log_to_file("", log_file)
    log_to_file("Feature Reduction Recommendations (top-10 permutation features):", log_file)
    log_to_file("-" * 70, log_file)
    for col, s in feat_stats.items():
        if s["trend_sig"] and s["hetero_flag"]:
            verdict = (
                "REVIEW — significant monotonic trend AND heteroscedasticity; "
                "model is systematically misusing this feature. "
                "Candidates: non-linear transform, binning, or removal."
            )
        elif s["trend_sig"]:
            verdict = (
                f"TREND (rho={s['spearman_r']:.2f}) — residuals shift with feature value; "
                "model may under-utilise its range. "
                "Consider interaction term or kernel/transform."
            )
        elif s["hetero_flag"]:
            verdict = (
                f"HETEROSCEDASTIC (Levene p={s['levene_p']:.3f}) — error variance changes across "
                "feature range; consider robust scaling or binning."
            )
        else:
            verdict = "OK — no systematic trend or heteroscedasticity detected; retain as-is."
        log_to_file(f"  {col}: {verdict}", log_file)

    log_to_file("", log_file)

    # --- 11. PCA explained-variance diagnostic ---
    # Fits PCA on training features to estimate the intrinsic dimensionality of
    # the descriptor space.  n_components_95 is stored in feat_stats so the
    # caller can use it to slice the permutation-importance-ranked feature list
    # (top-N named features) during feature reduction — keeping interpretability
    # instead of switching to anonymous PC axes.

    pca_diag = _PCA()
    pca_diag.fit(X_train)
    cumvar = pd.Series(pca_diag.explained_variance_ratio_.cumsum())
    n_components_95 = int((cumvar < 0.95).sum()) + 1  # first index that reaches ≥ 95%
    n_components_99 = int((cumvar < 0.99).sum()) + 1  # first index that reaches ≥ 99%

    log_to_file("", log_file)
    log_to_file(
        f"PCA dimensionality estimate: {n_components_95} components explain ≥95% variance "
        f"(out of {X_train.shape[1]} features).",
        log_file,
    )
    log_to_file(
        f"PCA dimensionality estimate: {n_components_99} components explain ≥99% variance "
        f"(out of {X_train.shape[1]} features).",
        log_file,
    )

    feat_stats["_pca95"] = {
        "n_components_95": n_components_95,
        "n_features_total": X_train.shape[1],
    }
    feat_stats["_pca99"] = {
        "n_components_99": n_components_99,
        "n_features_total": X_train.shape[1],
    }

    return feat_stats, indices, importances


# ----------------------------------------------------------------
# ## Model saving
# ----------------------------------------------------------------
def save_model(
    svr,
    feature_cols,
    cv_scores,
    X_train,
    X_test,
    y_test,
    y_pred,
    compartment,
    data_source,
    preprocessed_data,
    n_folds,
    test_size,
    feature_reduction_strategy,
    log_file,
):
    """
    Save fitted SVR to .joblib and model parameters to .json.

    Both files are timestamped: SVR_{compartment}_{data_source}_{YYYYMMDD_HHMMSS}.joblib/.json
    Overwrites any existing file at the same path.

    Parameters
    ----------
    feature_reduction_strategy : str or None
        Winner strategy label; None if no reduction was applied.
    """
    cv_rmse_log10 = np.sqrt(-cv_scores["test_neg_mean_squared_error"].mean())
    cv_rmse_log10_std = np.sqrt(cv_scores["test_neg_mean_squared_error"].std())
    test_rmse_log10 = np.sqrt(mean_squared_error(y_test, y_pred))
    test_mae_log10 = mean_absolute_error(y_test, y_pred)
    y_test_exp = np.power(10, y_test)
    y_pred_exp = np.power(10, y_pred)
    test_rmse_days = np.sqrt(mean_squared_error(y_test_exp, y_pred_exp))
    test_mae_days = mean_absolute_error(y_test_exp, y_pred_exp)

    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    models_dir = SCRIPT_DIR / "models"
    models_dir.mkdir(exist_ok=True)
    model_file = models_dir / f"SVR_{compartment}_{data_source}_{timestamp}.joblib"
    params_file = model_file.with_suffix(".json")

    joblib.dump(svr, model_file)
    log_to_file(f"Model saved to {model_file}", log_file, include_timestamp=True)

    model_parameters = {
        "model_name": "SVR",
        "compartment": compartment,
        "target_column": "T_half_days",
        "model_params": svr.get_params(),
        "feature_reduction_strategy": feature_reduction_strategy,
        "feature_columns": list(feature_cols),
        "y_value_transformation": "log10",
        "test_size": test_size,
        "cross-validation_folds": n_folds,
        "model_performance_scores": {
            "cv_r2_mean": float(cv_scores["test_r2"].mean()),
            "cv_r2_std": float(cv_scores["test_r2"].std()),
            "cv_mse_mean": float(-cv_scores["test_neg_mean_squared_error"].mean()),
            "cv_mse_std": float(cv_scores["test_neg_mean_squared_error"].std()),
            "cv_rmse_log10_mean": float(cv_rmse_log10),
            "cv_rmse_log10_std": float(cv_rmse_log10_std),
            "cv_mae_mean": float(-cv_scores["test_neg_mean_absolute_error"].mean()),
            "cv_mae_std": float(cv_scores["test_neg_mean_absolute_error"].std()),
            "test_rmse_log10": float(test_rmse_log10),
            "test_mae_log10": float(test_mae_log10),
            "test_rmse_days": float(test_rmse_days),
            "test_mae_days": float(test_mae_days),
        },
        "outlier_removal": preprocessed_data.remove_outliers,
        "preprocessing_drops": ["None"],
        "n_samples_train": len(X_train),
        "n_samples_test": len(X_test),
        "n_features": len(feature_cols),
    }

    with open(params_file, "w") as f:
        json.dump(model_parameters, f, indent=2)
    log_to_file(f"Parameters saved to {params_file}", log_file, include_timestamp=True)

    return model_file


# ----------------------------------------------------------------
# ## Feature reduction — strategy comparison
# ----------------------------------------------------------------
def run_feature_reduction(compartment, model_artifacts, feat_stats, indices, importances, logs_dir, log_file):
    """
    Compare four feature-reduction strategies using the same train/test split
    as the baseline model.  All evaluation is log-only — no plots, no model saves.

    Strategies
    ----------
    A  top-SVR      : all features with permutation importance > 0
    B  top-X        : top X = len(strategy A) by permutation importance
    C  top-N pca99  : top N = n_components_99 by permutation importance
    D  RFE          : LinearSVR inside RFE to select features, then re-eval with RBF SVR

    Uses the same svr_params, X_train, X_test, y_train, y_test as the baseline.
    """

    svr_params = model_artifacts["params"]
    X_train = model_artifacts["X_train"]
    X_test = model_artifacts["X_test"]
    y_train = model_artifacts["y_train"]
    y_test = model_artifacts["y_test"]
    feature_names = X_train.columns

    # Baseline metrics (log10 space)
    baseline_pred = model_artifacts["y_pred_svr"]
    baseline_r2 = 1 - np.sum((y_test - baseline_pred) ** 2) / np.sum((y_test - y_test.mean()) ** 2)
    baseline_rmse = np.sqrt(mean_squared_error(y_test, baseline_pred))
    baseline_mae = mean_absolute_error(y_test, baseline_pred)

    n_components_99 = feat_stats["_pca99"]["n_components_99"]

    # Strategy definitions: (label, feature_indices_array)
    # indices is the full descending-importance sort of all features
    pos_mask = importances[indices] > 0
    indices_pos = indices[pos_mask]  # strategy A
    top_svr_n = len(indices_pos)

    strategies = [
        ("A  top-SVR (importance>0)", indices_pos),
        ("B  top-X (N={})".format(top_svr_n), indices[:top_svr_n]),
        ("C  top-N pca99 (N={})".format(n_components_99), indices[:n_components_99]),
    ]

    log_section("Feature Reduction — Strategy Comparison", log_file)
    log_to_file(
        f"Baseline (full feature set, {X_train.shape[1]} features): "
        f"R²={baseline_r2:.4f}  RMSE={baseline_rmse:.4f}  MAE={baseline_mae:.4f}",
        log_file,
    )
    log_to_file("-" * 70, log_file)

    results = []
    cols_by_label = {}  # label → feature columns (Index)

    # --- Strategies A, B, C ---
    for label, feat_idx in strategies:
        cols = feature_names[feat_idx]
        n_feats = len(cols)
        Xtr = X_train[cols]
        Xte = X_test[cols]

        svr_red = SVR(**svr_params)
        svr_red.fit(Xtr, y_train)
        y_pred_red = svr_red.predict(Xte)

        r2 = 1 - np.sum((y_test - y_pred_red) ** 2) / np.sum((y_test - y_test.mean()) ** 2)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred_red))
        mae = mean_absolute_error(y_test, y_pred_red)

        delta_r2 = r2 - baseline_r2
        delta_rmse = rmse - baseline_rmse
        delta_mae = mae - baseline_mae

        log_to_file(
            f"Strategy {label}:\n"
            f"  Features selected ({n_feats}): {list(cols)}\n"
            f"  R²={r2:.4f} (Δ{delta_r2:+.4f})  "
            f"RMSE={rmse:.4f} (Δ{delta_rmse:+.4f})  "
            f"MAE={mae:.4f} (Δ{delta_mae:+.4f})",
            log_file,
        )
        results.append((label, n_feats, r2, rmse, mae))
        cols_by_label[label] = cols

    # --- Strategy D: RFE via LinearSVR, evaluate with RBF SVR ---
    log_to_file("", log_file)
    log_to_file("Strategy D  RFE (LinearSVR selector → RBF SVR evaluator):", log_file)

    rfe_n = top_svr_n  # same footprint as strategy A/B
    lin_svr = LinearSVR(max_iter=20000, random_state=42)
    rfe = RFE(estimator=lin_svr, n_features_to_select=rfe_n, step=1)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        rfe.fit(X_train, y_train)
    rfe_mask = rfe.support_
    rfe_cols = feature_names[rfe_mask]

    svr_rfe = SVR(**svr_params)
    svr_rfe.fit(X_train[rfe_cols], y_train)
    y_pred_rfe = svr_rfe.predict(X_test[rfe_cols])

    r2_rfe = 1 - np.sum((y_test - y_pred_rfe) ** 2) / np.sum((y_test - y_test.mean()) ** 2)
    rmse_rfe = np.sqrt(mean_squared_error(y_test, y_pred_rfe))
    mae_rfe = mean_absolute_error(y_test, y_pred_rfe)

    delta_r2_rfe = r2_rfe - baseline_r2
    delta_rmse_rfe = rmse_rfe - baseline_rmse
    delta_mae_rfe = mae_rfe - baseline_mae

    log_to_file(
        f"  RFE selector: LinearSVR, n_features_to_select={rfe_n}\n"
        f"  Features selected ({rfe_n}): {list(rfe_cols)}\n"
        f"  R²={r2_rfe:.4f} (Δ{delta_r2_rfe:+.4f})  "
        f"RMSE={rmse_rfe:.4f} (Δ{delta_rmse_rfe:+.4f})  "
        f"MAE={mae_rfe:.4f} (Δ{delta_mae_rfe:+.4f})",
        log_file,
    )
    rfe_label = "D  RFE (LinearSVR→RBF)"
    results.append((rfe_label, rfe_n, r2_rfe, rmse_rfe, mae_rfe))
    cols_by_label[rfe_label] = rfe_cols

    # --- Summary table ---
    log_to_file("", log_file)
    log_to_file("Summary (sorted by R², descending):", log_file)
    log_to_file(f"  {'Strategy':<35} {'N feats':>7} {'R²':>8} {'RMSE':>8} {'MAE':>8}", log_file)
    log_to_file(
        f"  {'Baseline (all features)':<35} {X_train.shape[1]:>7} {baseline_r2:>8.4f} {baseline_rmse:>8.4f} {baseline_mae:>8.4f}",
        log_file,
    )
    for label, n_feats, r2, rmse, mae in sorted(results, key=lambda x: x[2], reverse=True):
        log_to_file(f"  {label:<35} {n_feats:>7} {r2:>8.4f} {rmse:>8.4f} {mae:>8.4f}", log_file)

    # --- Best strategy selection ---
    # Highest R²; if A and B tie, prefer A
    best = max(results, key=lambda x: x[2])
    a_entry = next(r for r in results if r[0].startswith("A"))
    b_entry = next(r for r in results if r[0].startswith("B"))
    if best[0].startswith("B") and a_entry[2] == b_entry[2]:
        best = a_entry

    best_label, best_n, best_r2, best_rmse, best_mae = best
    best_cols = cols_by_label[best_label]

    log_to_file("", log_file)
    log_to_file(
        f"Winner: {best_label}  (N={best_n}, R²={best_r2:.4f}, RMSE={best_rmse:.4f}, MAE={best_mae:.4f})",
        log_file,
    )

    # --- Re-train final model on best feature set ---
    log_section("Feature Reduction — Final Model", log_file)

    n_folds = model_artifacts["n_folds"]
    test_size = model_artifacts["test_size"]
    preprocessed_data = model_artifacts["preprocessed_data"]
    data_source = model_artifacts["data_source"]

    svr_final, cv_scores_red = fit_svr(X_train[best_cols], y_train, svr_params, n_folds, log_file)

    y_pred_final = svr_final.predict(X_test[best_cols])

    cv_rmse_log10 = np.sqrt(-cv_scores_red["test_neg_mean_squared_error"].mean())
    cv_rmse_log10_std = np.sqrt(cv_scores_red["test_neg_mean_squared_error"].std())
    test_rmse_log10 = np.sqrt(mean_squared_error(y_test, y_pred_final))
    test_mae_log10 = mean_absolute_error(y_test, y_pred_final)
    y_test_exp = np.power(10, y_test)
    y_pred_final_exp = np.power(10, y_pred_final)
    test_rmse_days = np.sqrt(mean_squared_error(y_test_exp, y_pred_final_exp))
    test_mae_days = mean_absolute_error(y_test_exp, y_pred_final_exp)

    log_to_file(
        f"Final model CV ({n_folds}-fold on training set, reduced features):\n"
        f"  R² (CV): {cv_scores_red['test_r2'].mean():.4f} (+/- {cv_scores_red['test_r2'].std():.4f})\n"
        f"  RMSE (log10, CV): {cv_rmse_log10:.4f} (+/- {cv_rmse_log10_std:.4f})\n"
        f"  MAE (CV): {-cv_scores_red['test_neg_mean_absolute_error'].mean():.4f}",
        log_file,
    )
    log_to_file(
        f"Final model test set:\n"
        f"  RMSE (log10): {test_rmse_log10:.4f}  MAE (log10): {test_mae_log10:.4f}\n"
        f"  RMSE (days):  {test_rmse_days:.2f}   MAE (days):  {test_mae_days:.2f}",
        log_file,
    )

    save_model(
        svr_final,
        best_cols,
        cv_scores_red,
        X_train[best_cols],
        X_test[best_cols],
        y_test,
        y_pred_final,
        compartment,
        data_source,
        preprocessed_data,
        n_folds,
        test_size,
        best_label,
        log_file,
    )

    return {
        **model_artifacts,
        "svr": svr_final,
        "X_train": X_train[best_cols],
        "X_test": X_test[best_cols],
        "y_pred_svr": y_pred_final,
        "y_pred_train": svr_final.predict(X_train[best_cols]),
        "cv_scores": cv_scores_red,
    }


# ----------------------------------------------------------------
# ## Entry point
# ----------------------------------------------------------------
def main(compartment: str, data_source: str):
    # Validation
    if data_source == "vega" and compartment == "air":
        raise ValueError("VEGA dataset does not contain air compartment data.")

    working_dir = SCRIPT_DIR

    DATABASE_FILE = data_dir / DATABASE_FILES[data_source]
    ENGINE = sa.create_engine(f"sqlite:///{DATABASE_FILE}")
    Session = sessionmaker(bind=ENGINE)
    data_to_use = get_all_data(compartment, Session)

    target_column = "T_half_days"
    use_outlier_removal = True

    # Logging setup
    log_time = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    logs_dir = working_dir / "logs" / f"{compartment}_{data_source}" / log_time
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = str(logs_dir / f"9_SVR.log")
    log_to_file("SVR Model and Analysis", log_file)
    log_section("Preprocessor class", log_file)

    # Preprocessing
    log_to_file(
        f"Dataset X and y shapes before preprocessing: X={data_to_use.shape[0]}, y={data_to_use.shape[0]}",
        log_file,
    )
    to_drop = ["None"]
    log_to_file(f"Columns to drop: {to_drop}", log_file)

    preprocessor = Preprocessor(
        target_column=target_column,
        to_drop=to_drop,
        remove_outliers=use_outlier_removal,
    )
    data_clean = preprocessor.preprocess(data_to_use)
    preprocessed_data = PreprocessedData(
        name="AirData",
        df=data_clean,
        remove_outliers=use_outlier_removal,
        X=preprocessor.X,
        y_log=preprocessor.y_log,
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

    # Train
    model_artifacts = train_model(compartment, data_source, preprocessed_data, preprocessor, logs_dir, log_file)

    # Feature analysis on full-feature model (drives reduction strategy selection)
    log_section("Model Analysis", log_file)
    feat_stats, indices, importances = analyse_features(compartment, model_artifacts, logs_dir, log_file)

    # Feature reduction — returns updated artifacts for the final reduced model
    reduced_artifacts = run_feature_reduction(
        compartment,
        model_artifacts,
        feat_stats,
        indices,
        importances,
        logs_dir,
        log_file,
    )

    # Post-reduction diagnostics on final model
    analyse_basic(compartment, reduced_artifacts, logs_dir, log_file)
    analyse_features(compartment, reduced_artifacts, logs_dir, log_file, suffix="_reduced")
    analyse_coverage_domain(compartment, reduced_artifacts, logs_dir, log_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SVR analysis for biodegradability")
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

    main(args.compartment, args.data_source)
    print("SVR model training and analysis complete. Check log file for details.")
