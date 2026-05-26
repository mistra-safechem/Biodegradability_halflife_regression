"""Interval-aware SVR pipeline for quantised persistence endpoints.

Targets are treated as log10 intervals rather than precise kinetic measurements.
Evaluation uses interval coverage probability and rank-order correlation;
RMSE / R² on raw half-life values are intentionally absent.

Usage
-----
    uv run SVR_interval_model_and_analysis.py --compartment soil --data-source hsbd
    uv run SVR_interval_model_and_analysis.py --compartment water --data-source vega
    uv run SVR_interval_model_and_analysis.py --compartment sediment --data-source combined

See ``--help`` for all options.

Output
------
All artefacts are written to logs/<compartment>_<data_source>/<timestamp>/:
  - SVR_interval.log                        main log
  - interval_table_<c>.csv                  class → log10 bounds mapping
  - learning_curve_<c>.png
  - true_vs_pred_<c>.png                    interval-coloured true vs predicted (log10 + days panels)
  - feature_importance_<c>.png              permutation importances (full feature set)
  - feature_importance_<c>_reduced.png      permutation importances (reduced feature set)
  - ad_membership_<c>.csv                   per-compound leverage / interval membership
  - residuals/
      qqplot_<c>.png
      residual_hist_<c>.png
  - ad_analysis/                            (morgan PCA + Butina clustering)
      chemical_space_pca_<c>.png
      williams_plot_<c>.png
      …

All model and companion files are written to models/ with the same timestamp for traceability:
  - SVR_<c>_<ds>_<ts>.joblib           (model file)
  - SVR_<c>_<ds>_<ts>.json             (meta-data file)
  - SVR_<c>_<ds>_<ts>_ad.npz           (applicability domain artefact for inference)
  - t_half_meta_<c>_<ds>_<ts>.csv      (training-split T_half metadata for inference traceability)
  - training_split_<c>_<ds>_<ts>.csv   (SMILES + train/test split info for inference traceability)
"""

import argparse
import json
import os
import warnings
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")

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
from sklearn.model_selection import cross_validate, learning_curve
from sklearn.svm import SVR, LinearSVR
from sqlalchemy.orm import sessionmaker

from src.db_schema import *  # noqa: F401,F403  — registers ORM models
from src.db_utils import get_all_data
from src.log_utils import get_logger
from src.ml_tools import (
    PreprocessedData,
    Preprocessor,
    applicability_domain_leverage,
    chemical_space_morgan_pca,
    chemical_space_pca,
    log_evaluation_metrics,
    svr_grid_search,
    t_t_split,
)

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 250)

SCRIPT_DIR = Path(__file__).parent.resolve()
data_dir = SCRIPT_DIR / "processed_data"
log_base_dir = SCRIPT_DIR / "logs"
models_dir = SCRIPT_DIR / "models"
DATABASE_FILES = {
    "hsbd": "hsbd_t_half_all.db",
    "vega": "vega_t_half_soil_water_sediment.db",
    "combined": "combined_t_half_vega_hsbd_soil_water_sediment.db",
}


# ---------------------------------------------------------------------------
# Shared SVR fit + cross-validation
# ---------------------------------------------------------------------------


def fit_svr(X_train, y_train, svr_params, n_folds, log):
    """Fit SVR and run k-fold cross-validation.  Returns (svr, cv_scores)."""
    svr = SVR(**svr_params)
    svr.fit(X_train, y_train)

    cv = cross_validate(
        SVR(**svr_params),
        X_train,
        y_train,
        cv=n_folds,
        scoring=["r2", "neg_mean_squared_error", "neg_mean_absolute_error"],
        return_train_score=True,
    )
    rmse_log10 = np.sqrt(-cv["test_neg_mean_squared_error"].mean())
    log.info(
        "%d-fold CV on training set: R²=%.4f (±%.4f)  RMSE(log10)=%.4f  RMSE(days-geom)=%.2f×",
        n_folds,
        cv["test_r2"].mean(),
        cv["test_r2"].std(),
        rmse_log10,
        10**rmse_log10,
    )
    return svr, cv


# ---------------------------------------------------------------------------
# Model training
# ---------------------------------------------------------------------------


def train_model(compartment, data_source, preprocessed_data, preprocessor, logs_dir, log):
    """Grid search, CV, final fit, and interval-metric evaluation."""
    X = preprocessed_data.X
    y_log10 = preprocessed_data.y_log10
    y_lower = preprocessed_data.y_lower
    y_upper = preprocessed_data.y_upper

    log.section("Model Training")

    n_folds = 5
    test_size = 0.25
    random_state = 42

    X_train, X_test, y_log10_train, y_log10_test, y_lower_train, y_lower_test, y_upper_train, y_upper_test = t_t_split(
        X,
        y_log10,
        y_lower,
        y_upper,
        test_size=test_size,
        random_state=random_state,
    )

    smiles_train, smiles_test = preprocessor.get_smiles_for_split(X_train, X_test)

    # Grid search (uses log10 centroid as proxy target)
    svr, svr_params = svr_grid_search(X_train, y_log10_train)

    # CV
    svr, cv_scores = fit_svr(X_train, y_log10_train, svr_params, n_folds, log)

    y_pred_test = svr.predict(X_test)
    y_pred_train = svr.predict(X_train)

    # Interval evaluation — primary metrics
    log.section("Interval Evaluation — Test Set")
    test_metrics = log_evaluation_metrics(
        y_pred_test,
        y_lower_test.values,
        y_upper_test.values,
        y_log10_test.values,
        logger=log,
        prefix="  [test]  ",
    )

    log.section("Interval Evaluation — Train Set (reference)")
    train_metrics = log_evaluation_metrics(
        y_pred_train,
        y_lower_train.values,
        y_upper_train.values,
        y_log10_train.values,
        logger=log,
        prefix="  [train] ",
    )

    return {
        "svr": svr,
        "X": X,
        "X_train": X_train,
        "X_test": X_test,
        "y_log10": y_log10,
        "y_log10_train": y_log10_train,
        "y_log10_test": y_log10_test,
        "y_lower_train": y_lower_train,
        "y_upper_train": y_upper_train,
        "y_lower_test": y_lower_test,
        "y_upper_test": y_upper_test,
        "y_pred_test": y_pred_test,
        "y_pred_train": y_pred_train,
        "smiles_train": smiles_train,
        "smiles_test": smiles_test,
        "cv_scores": cv_scores,
        "params": svr_params,
        "n_folds": n_folds,
        "test_size": test_size,
        "random_state": random_state,
        "data_source": data_source,
        "preprocessed_data": preprocessed_data,
        "train_metrics": train_metrics,
        "test_metrics": test_metrics,
    }


# ---------------------------------------------------------------------------
# Analysis — learning curve and diagnostic plots
# ---------------------------------------------------------------------------


def analyse_basic(compartment, ma, logs_dir, log):
    """Learning curve, interval-coloured true-vs-predicted, residuals."""
    svr = ma["svr"]
    X = ma["X"]
    y_log10 = ma["y_log10"]
    y_log10_test = ma["y_log10_test"]
    y_pred_test = ma["y_pred_test"]
    y_lower_test = ma["y_lower_test"]
    y_upper_test = ma["y_upper_test"]
    n_folds = ma["n_folds"]
    random_state = ma["random_state"]

    log.section("Learning Curve")

    sizes, tr_sc, val_sc = learning_curve(
        svr,
        X,
        y_log10,
        train_sizes=np.linspace(0.1, 1.0, 10),
        cv=n_folds,
        scoring="r2",
        random_state=random_state,
        n_jobs=-1,
    )
    tr_mean, tr_std = tr_sc.mean(axis=1), tr_sc.std(axis=1)
    val_mean, val_std = val_sc.mean(axis=1), val_sc.std(axis=1)
    log.info("Learning curve val R² (final): %.4f ± %.4f", val_mean[-1], val_std[-1])

    plt.figure(figsize=(10, 6))
    plt.plot(sizes, tr_mean, "o-", color="royalblue", label="Training R²")
    plt.fill_between(sizes, tr_mean - tr_std, tr_mean + tr_std, alpha=0.2, color="royalblue")
    plt.plot(sizes, val_mean, "o-", color="darkorange", label="CV R²")
    plt.fill_between(sizes, val_mean - val_std, val_mean + val_std, alpha=0.2, color="darkorange")
    plt.xlabel("Training Set Size")
    plt.ylabel("R² (log10 centroid proxy)")
    plt.title(f"Learning Curve — {compartment}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(logs_dir / f"learning_curve_{compartment}.png")
    plt.close()

    # True-vs-predicted coloured by interval membership — 2-panel: log10 + days
    residuals_log10 = y_log10_test.values - y_pred_test
    inside = (y_pred_test >= y_lower_test.values) & (y_pred_test <= y_upper_test.values)
    colors = np.where(inside, "royalblue", "crimson")

    # Back-convert to days for the right panel
    true_days = 10**y_log10_test.values
    pred_days = 10**y_pred_test

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(color="royalblue", label=f"Inside interval (n={inside.sum()})"),
        Patch(color="crimson", label=f"Outside interval (n={(~inside).sum()})"),
    ]

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 7))

    # Left panel — log10 space
    ax_l.scatter(y_log10_test, y_pred_test, c=colors, alpha=0.6, s=25)
    lims_l = [min(y_log10_test.min(), y_pred_test.min()), max(y_log10_test.max(), y_pred_test.max())]
    ax_l.plot(lims_l, lims_l, "k--", linewidth=1, label="1:1")
    ax_l.legend(handles=legend_handles, fontsize=9)
    ax_l.set_xlabel("True (log10 days)", fontsize=11)
    ax_l.set_ylabel("Predicted (log10 days)", fontsize=11)
    ax_l.set_title(f"True vs Predicted — {compartment}\n(log10 scale)")

    # Right panel — days space
    ax_r.scatter(true_days, pred_days, c=colors, alpha=0.6, s=25)
    lims_r = [min(true_days.min(), pred_days.min()), max(true_days.max(), pred_days.max())]
    ax_r.plot(lims_r, lims_r, "k--", linewidth=1, label="1:1")
    ax_r.legend(handles=legend_handles, fontsize=9)
    ax_r.set_xlabel("True (days)", fontsize=11)
    ax_r.set_ylabel("Predicted (days)", fontsize=11)
    ax_r.set_title(f"True vs Predicted — {compartment}\n(days scale)")

    fig.suptitle("colour = interval membership", fontsize=10, y=0.02)
    fig.tight_layout()
    fig.savefig(logs_dir / f"true_vs_pred_{compartment}.png", dpi=150)
    plt.close(fig)

    # Q-Q and histogram of residuals
    residuals_dir = logs_dir / "residuals"
    residuals_dir.mkdir(exist_ok=True)

    # Signed days residuals: predicted − true (positive = over-prediction)
    residuals_days = pred_days - true_days

    # Q-Q plot (log10 space) + days histogram side-by-side
    fig_qq, (ax_qq, ax_dhist) = plt.subplots(1, 2, figsize=(14, 5))

    stats.probplot(residuals_log10, dist="norm", plot=ax_qq)
    ax_qq.set_title(f"Q-Q Plot of Residuals — {compartment}\n(log10 days)")

    sns.histplot(residuals_days, kde=True, ax=ax_dhist, color="steelblue")
    ax_dhist.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax_dhist.set_xlabel("Residual (days)  [predicted − true]", fontsize=11)
    ax_dhist.set_title(f"Residual Distribution — {compartment}\n(days scale)")

    fig_qq.tight_layout()
    fig_qq.savefig(residuals_dir / f"qqplot_{compartment}.png", dpi=150)
    plt.close(fig_qq)

    # Residual histogram — 2-panel: log10 + days
    fig_rh, (ax_rh_l, ax_rh_r) = plt.subplots(1, 2, figsize=(14, 5))

    sns.histplot(residuals_log10, kde=True, ax=ax_rh_l, color="royalblue")
    ax_rh_l.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax_rh_l.set_xlabel("Residual (log10 days)  [true − predicted]", fontsize=11)
    ax_rh_l.set_title(f"Residual Distribution — {compartment}\n(log10 scale)")

    sns.histplot(residuals_days, kde=True, ax=ax_rh_r, color="steelblue")
    ax_rh_r.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax_rh_r.set_xlabel("Residual (days)  [predicted − true]", fontsize=11)
    ax_rh_r.set_title(f"Residual Distribution — {compartment}\n(days scale)")

    fig_rh.tight_layout()
    fig_rh.savefig(residuals_dir / f"residual_hist_{compartment}.png", dpi=150)
    plt.close(fig_rh)

    log.info("Basic diagnostic plots written to %s", logs_dir)


# ---------------------------------------------------------------------------
# Analysis — applicability domain
# ---------------------------------------------------------------------------


def analyse_coverage_domain(compartment, ma, logs_dir, log):
    """PCA chemical-space coverage + Williams plot AD assessment."""
    log.section("Applicability Domain")

    pca, scaler, tr_pcs, te_pcs, cov_pct = chemical_space_pca(ma["X_train"], ma["X_test"], logs_dir, compartment)
    log.info("PCA coverage: %.1f%%", cov_pct)

    ad_results = applicability_domain_leverage(
        ma["X_train"],
        ma["X_test"],
        ma["y_lower_train"],
        ma["y_upper_train"],
        ma["y_lower_test"],
        ma["y_upper_test"],
        ma["y_pred_train"],
        ma["y_pred_test"],
        logs_dir / "ad_analysis",
        compartment,
    )

    cluster_stats = chemical_space_morgan_pca(
        ma["smiles_train"],
        ma["smiles_test"],
        ad_results,
        logs_dir / "ad_analysis",
        compartment,
    )
    log.info("AD-outside test compounds: %d", cluster_stats["n_outside_ad"])

    return ad_results, cluster_stats


# ---------------------------------------------------------------------------
# Feature importance and reduction
# ---------------------------------------------------------------------------


def analyse_features(compartment, ma, logs_dir, log, suffix=""):
    """Permutation importance, per-feature residuals, PCA dimensionality."""
    log.section(f"Feature Analysis{' (reduced)' if suffix else ''}")

    svr = ma["svr"]
    X_train = ma["X_train"]
    X_test = ma["X_test"]
    y_log10_test = ma["y_log10_test"]
    y_pred_test = ma["y_pred_test"]

    # Permutation importance
    result = permutation_importance(
        svr,
        X_test,
        y_log10_test,
        n_repeats=10,
        random_state=42,
        scoring="neg_mean_squared_error",
    )
    imp = result.importances_mean
    indices = np.argsort(imp)[::-1]
    top_n = 15

    plt.figure(figsize=(10, 5))
    plt.bar(range(top_n), imp[indices[:top_n]])
    plt.xticks(range(top_n), X_test.columns[indices[:top_n]], rotation=90)
    plt.title(f"Top {top_n} Permutation Importances — {compartment}{suffix}")
    plt.ylabel("Mean importance")
    plt.tight_layout()
    plt.savefig(logs_dir / f"feature_importance_{compartment}{suffix}.png")
    plt.close()

    # PCA dimensionality diagnostic
    pca_d = _PCA()
    pca_d.fit(X_train)
    cumvar = pd.Series(pca_d.explained_variance_ratio_.cumsum())
    n95 = int((cumvar < 0.95).sum()) + 1
    n99 = int((cumvar < 0.99).sum()) + 1
    log.info("PCA dimensionality: %d components ≥95%%, %d components ≥99%% (of %d total).", n95, n99, X_train.shape[1])

    feat_stats: dict = {
        "_pca95": {"n_components_95": n95, "n_features_total": X_train.shape[1]},
        "_pca99": {"n_components_99": n99, "n_features_total": X_train.shape[1]},
    }
    return feat_stats, indices, imp


# ---------------------------------------------------------------------------
# Feature reduction strategies
# ---------------------------------------------------------------------------


def run_feature_reduction(compartment, ma, feat_stats, indices, importances, logs_dir, log, timestamp=None):
    """Compare four feature-reduction strategies; retrain and save best."""
    log.section("Feature Reduction")

    svr_params = ma["params"]
    X_train, X_test = ma["X_train"], ma["X_test"]
    y_log10_train, y_log10_test = ma["y_log10_train"], ma["y_log10_test"]
    y_lower_train, y_upper_train = ma["y_lower_train"], ma["y_upper_train"]
    y_lower_test, y_upper_test = ma["y_lower_test"], ma["y_upper_test"]
    feat_names = X_train.columns

    n99 = feat_stats["_pca99"]["n_components_99"]
    pos_mask = importances[indices] > 0
    indices_pos = indices[pos_mask]
    top_svr_n = len(indices_pos)

    strategies = [
        ("A  top-SVR (importance>0)", indices_pos),
        ("B  top-X (N={})".format(top_svr_n), indices[:top_svr_n]),
        ("C  top-N pca99 (N={})".format(n99), indices[:n99]),
    ]

    results = []
    cols_by_label = {}

    def _eval(X_tr, X_te, y_lo_tr, y_up_tr, y_lo_te, y_up_te, y_log10_te, label):
        svr_r = SVR(**svr_params)
        svr_r.fit(X_tr, y_log10_train)
        pred = svr_r.predict(X_te)
        cov = log_evaluation_metrics(
            pred,
            y_lo_te.values,
            y_up_te.values,
            y_log10_te.values,
            logger=log,
            prefix=f"  [{label}] ",
        )
        return svr_r, pred, cov

    for label, feat_idx in strategies:
        cols = feat_names[feat_idx]
        svr_r, pred, cov = _eval(
            X_train[cols],
            X_test[cols],
            y_lower_train,
            y_upper_train,
            y_lower_test,
            y_upper_test,
            y_log10_test,
            label,
        )
        results.append((label, len(cols), cov["coverage_probability"], cov["mean_interval_loss"], cov["spearman_r"]))
        cols_by_label[label] = cols
        log.info(
            "Strategy %s: n=%d  coverage=%.3f  MIL=%.4f  rho=%.3f",
            label,
            len(cols),
            cov["coverage_probability"],
            cov["mean_interval_loss"],
            cov["spearman_r"],
        )

    # Strategy D — RFE via LinearSVR
    rfe_n = top_svr_n
    lin_svr = LinearSVR(max_iter=20_000, random_state=42)
    rfe = RFE(estimator=lin_svr, n_features_to_select=rfe_n, step=1)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        rfe.fit(X_train, y_log10_train)
    rfe_cols = feat_names[rfe.support_]
    d_label = "D  RFE (LinearSVR→RBF)"
    svr_d, pred_d, cov_d = _eval(
        X_train[rfe_cols],
        X_test[rfe_cols],
        y_lower_train,
        y_upper_train,
        y_lower_test,
        y_upper_test,
        y_log10_test,
        d_label,
    )
    results.append((d_label, rfe_n, cov_d["coverage_probability"], cov_d["mean_interval_loss"], cov_d["spearman_r"]))
    cols_by_label[d_label] = rfe_cols

    # Select best by coverage probability (primary), then Spearman rho (tie-break)
    best = max(results, key=lambda r: (r[2], r[4]))
    best_label, best_n, best_cov, best_mil, best_rho = best
    best_cols = cols_by_label[best_label]
    log.info("Winner: %s  (n=%d, coverage=%.3f, rho=%.3f)", best_label, best_n, best_cov, best_rho)

    # Retrain final model on best feature set
    log.section("Feature Reduction — Final Model")
    n_folds = ma["n_folds"]
    svr_final, cv_final = fit_svr(X_train[best_cols], y_log10_train, svr_params, n_folds, log)
    y_pred_final = svr_final.predict(X_test[best_cols])

    log.section("Final Model — Interval Evaluation")
    final_metrics = log_evaluation_metrics(
        y_pred_final,
        y_lower_test.values,
        y_upper_test.values,
        y_log10_test.values,
        logger=log,
        prefix="  [final] ",
    )

    model_file = save_model(
        svr_final,
        best_cols,
        cv_final,
        X_train[best_cols],
        X_test[best_cols],
        y_lower_test,
        y_upper_test,
        y_log10_test,
        y_pred_final,
        compartment,
        ma["data_source"],
        ma["preprocessed_data"],
        n_folds,
        ma["test_size"],
        best_label,
        final_metrics,
        log,
        timestamp=timestamp,
    )

    return {
        **ma,
        "svr": svr_final,
        "X_train": X_train[best_cols],
        "X_test": X_test[best_cols],
        "y_pred_test": y_pred_final,
        "y_pred_train": svr_final.predict(X_train[best_cols]),
        "cv_scores": cv_final,
        "feature_reduction_strategy": best_label,
        "model_file": model_file,
    }


# ---------------------------------------------------------------------------
# Model saving
# ---------------------------------------------------------------------------


def save_model(
    svr,
    feature_cols,
    cv_scores,
    X_train,
    X_test,
    y_lower_test,
    y_upper_test,
    y_log10_test,
    y_pred,
    compartment,
    data_source,
    preprocessed_data,
    n_folds,
    test_size,
    feature_reduction_strategy,
    final_metrics,
    log,
    extra_meta: dict = None,
    timestamp: str = None,
):
    """Save fitted SVR to .joblib and parameters to .json.

    Parameters
    ----------
    extra_meta:
        Optional dict of additional keys merged into the JSON model card.
        Used to inject artefact paths (ad_artefact_file, etc.) that are
        only known after this function returns.  Call
        :func:`patch_model_json` instead if the paths are computed later.
    timestamp:
        Timestamp string (``%Y%m%d_%H%M%S``).  When provided by the caller
        the same stamp is shared across the model file and all companion
        artefacts (CSV files, AD .npz) so filenames stay in sync.
        Defaults to the current time if omitted.
    """
    if timestamp is None:
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    models_dir.mkdir(exist_ok=True)
    model_file = models_dir / f"SVR_{compartment}_{data_source}_{timestamp}.joblib"
    params_file = model_file.with_suffix(".json")

    joblib.dump(svr, model_file)
    log.info("Model saved to %s", model_file.relative_to(SCRIPT_DIR))

    params = {
        "model_name": "SVR",
        "compartment": compartment,
        "data_source": data_source,
        "target_interpretation": "interval-valued (interval-aware regression)",
        "target_column": "T_half_days",
        "y_value_transformation": "log10",
        "interval_strategy": "half_order_of_magnitude",
        "model_params": svr.get_params(),
        "feature_reduction_strategy": feature_reduction_strategy,
        "feature_columns": list(feature_cols),
        "test_size": test_size,
        "cross_validation_folds": n_folds,
        "model_performance_scores": {
            "cv_r2_mean": float(cv_scores["test_r2"].mean()),
            "cv_r2_std": float(cv_scores["test_r2"].std()),
            "interval_coverage_probability": final_metrics["coverage_probability"],
            "mean_interval_loss_log10": final_metrics["mean_interval_loss"],
            "mean_interval_loss_days": final_metrics["mean_interval_loss_days"],
            "class_accuracy": final_metrics["class_accuracy"],
            "spearman_r": final_metrics["spearman_r"],
            "spearman_p": final_metrics["spearman_p"],
            "kendall_tau": final_metrics["kendall_tau"],
            "kendall_p": final_metrics["kendall_p"],
        },
        "outlier_removal": preprocessed_data.remove_outliers,
        "n_samples_train": len(X_train),
        "n_samples_test": len(X_test),
        "n_features": len(feature_cols),
        "note": (
            "Predictions represent interval-consistent categorisation, "
            "not compound-specific kinetic accuracy. "
            "Performance reflects persistence categorisation ability."
        ),
    }

    if extra_meta:
        params.update(extra_meta)

    with open(params_file, "w") as f:
        json.dump(params, f, indent=2)
    log.info("Parameters saved to %s", params_file.relative_to(SCRIPT_DIR))
    return model_file


def patch_model_json(model_file: Path, extra: dict) -> None:
    """Merge *extra* keys (with relative paths) into an already-written model JSON card.

    All Path values in *extra* are converted to paths relative to the directory
    that contains the JSON file before writing, so model cards are portable
    across machines and directory layouts.
    """
    params_file = model_file.with_suffix(".json")
    json_dir = params_file.parent
    with open(params_file) as f:
        params = json.load(f)
    relative_extra = {
        k: os.path.relpath(v, start=json_dir) if isinstance(v, (str, Path)) and Path(v).is_absolute() else v
        for k, v in extra.items()
    }
    params.update(relative_extra)
    with open(params_file, "w") as f:
        json.dump(params, f, indent=2)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(compartment: str, data_source: str) -> None:
    if (data_source in ("vega", "combined")) and compartment == "air":
        raise ValueError("VEGA and combined datasets do not contain air compartment data.")

    if data_source not in DATABASE_FILES:
        raise ValueError(f"Unknown data source: {data_source}")
    db_file = data_dir / DATABASE_FILES[data_source]

    engine = sa.create_engine(f"sqlite:///{db_file}")
    Session = sessionmaker(bind=engine)
    df_raw = get_all_data(compartment, Session)

    # Logging setup — timestamp shared with model artefacts for consistent naming
    log_time = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    timestamp = log_time  # reused for model filename, split CSV, and meta CSV
    logs_dir = log_base_dir / f"{compartment}_{data_source}" / log_time
    logs_dir.mkdir(parents=True, exist_ok=True)
    log = get_logger(
        f"{compartment}_{data_source}",
        log_dir=logs_dir,
        log_filename="SVR_interval.log",
    )

    log.section(f"Interval-Aware SVR — {compartment.upper()} / {data_source}")
    log.info("Database: %s", db_file)
    log.info("Raw dataset: %d rows", len(df_raw))

    # Preprocessing
    preprocessor = Preprocessor(
        target_column="T_half_days",
        to_drop={"None"},
        remove_outliers=True,
    )
    df_clean = preprocessor.preprocess(df_raw)

    # Export interval table for reproducibility (§11)
    interval_csv = logs_dir / f"interval_table_{compartment}.csv"
    preprocessor.interval_table.to_csv(interval_csv, index=False)
    log.info("Interval mapping table saved to %s", interval_csv)

    # Export stashed T_half auxiliary columns (row-aligned with X/y) for
    # traceability and sanity-checking — NOT used as features.
    # Stored in models/ alongside the model files to avoid filename clashes
    # across runs.  Actual writing is deferred until models_dir is created
    # inside save_model(); we store the DataFrame here and write it later.
    _t_half_meta_df = (
        preprocessor.t_half_meta if preprocessor.t_half_meta is not None and not preprocessor.t_half_meta.empty else None
    )

    preprocessed_data = PreprocessedData(
        name=compartment,
        df=df_clean,
        remove_outliers=True,
        X=preprocessor.X,
        y_log10=preprocessor.y_log10,
        y_lower=preprocessor.y_lower,
        y_upper=preprocessor.y_upper,
        interval_table=preprocessor.interval_table,
    )

    # Train
    model_artifacts = train_model(compartment, data_source, preprocessed_data, preprocessor, logs_dir, log)

    # Feature analysis
    feat_stats, indices, importances = analyse_features(compartment, model_artifacts, logs_dir, log)

    # Feature reduction → final model
    reduced = run_feature_reduction(
        compartment, model_artifacts, feat_stats, indices, importances, logs_dir, log, timestamp=timestamp
    )

    # Post-reduction diagnostics
    analyse_basic(compartment, reduced, logs_dir, log)
    analyse_features(compartment, reduced, logs_dir, log, suffix="_reduced")
    ad_results, cluster_stats = analyse_coverage_domain(compartment, reduced, logs_dir, log)

    # ------------------------------------------------------------------
    # Export C: training split CSV (models/ — answers "was this in training?")
    # Stored alongside the model files with timestamped name to avoid clashes.
    # ------------------------------------------------------------------
    split_rows = []
    for smi, y, lo, hi in zip(
        reduced["smiles_train"],
        reduced["y_log10_train"],
        reduced["y_lower_train"],
        reduced["y_upper_train"],
    ):
        split_rows.append({"smiles": smi, "split": "train", "y_log10": y, "y_lower": lo, "y_upper": hi})
    for smi, y, lo, hi in zip(
        reduced["smiles_test"],
        reduced["y_log10_test"],
        reduced["y_lower_test"],
        reduced["y_upper_test"],
    ):
        split_rows.append({"smiles": smi, "split": "test", "y_log10": y, "y_lower": lo, "y_upper": hi})
    models_dir.mkdir(exist_ok=True)
    split_csv = models_dir / f"training_split_{compartment}_{data_source}_{timestamp}.csv"
    pd.DataFrame(split_rows).to_csv(split_csv, index=False)
    log.info("Training split CSV saved to %s", split_csv.relative_to(SCRIPT_DIR))

    # ------------------------------------------------------------------
    # Export T_half meta CSV (models/ — alongside model files)
    # ------------------------------------------------------------------
    if _t_half_meta_df is not None:
        meta_csv = models_dir / f"t_half_meta_{compartment}_{data_source}_{timestamp}.csv"
        _t_half_meta_df.to_csv(meta_csv, index=False)
        log.info("T_half meta columns saved to %s", meta_csv.relative_to(SCRIPT_DIR))
    else:
        meta_csv = None

    # ------------------------------------------------------------------
    # Export A: AD artefact .npz (models/ — same stem as .joblib + _ad)
    # Stored alongside the model so inference only needs the models/ dir.
    # ------------------------------------------------------------------
    model_file = reduced["model_file"]
    ad_npz_path = model_file.with_name(model_file.stem + "_ad.npz")
    np.savez_compressed(
        ad_npz_path,
        X_train=ad_results["X_train"].values.astype(np.float64),
        X_train_scaled=ad_results["X_train_scaled"].astype(np.float64),
        XtX_inv=ad_results["XtX_inv"].astype(np.float64),
        X_train_mean=ad_results["X_train_mean"].astype(np.float64),
        X_train_std=ad_results["X_train_std"].astype(np.float64),
        h_star=np.array([ad_results["h_star"]], dtype=np.float64),
        feature_cols=np.array(list(reduced["X_train"].columns)),
    )
    log.info("AD artefact saved to %s", ad_npz_path.relative_to(SCRIPT_DIR))

    # ------------------------------------------------------------------
    # Export B: per-compound AD membership CSV (logs/)
    # ------------------------------------------------------------------
    ad_rows = []
    for smi, h, inside, il in zip(
        reduced["smiles_train"],
        ad_results["h_train"],
        ad_results["inside_train"],
        ad_results["interval_loss_train"],
    ):
        ad_rows.append(
            {
                "smiles": smi,
                "split": "train",
                "leverage_h": float(h),
                "h_star": float(ad_results["h_star"]),
                "in_ad_leverage": bool(h <= ad_results["h_star"]),
                "in_interval": bool(inside),
                "interval_loss": float(il),
            }
        )
    for smi, h, inside, il in zip(
        reduced["smiles_test"],
        ad_results["h_test"],
        ad_results["inside_test"],
        ad_results["interval_loss_test"],
    ):
        ad_rows.append(
            {
                "smiles": smi,
                "split": "test",
                "leverage_h": float(h),
                "h_star": float(ad_results["h_star"]),
                "in_ad_leverage": bool(h <= ad_results["h_star"]),
                "in_interval": bool(inside),
                "interval_loss": float(il),
            }
        )
    ad_csv = logs_dir / f"ad_membership_{compartment}.csv"
    pd.DataFrame(ad_rows).to_csv(ad_csv, index=False)
    log.info("AD membership CSV saved to %s", ad_csv)

    # ------------------------------------------------------------------
    # Patch JSON model card with artefact paths (stored as relative paths)
    # ------------------------------------------------------------------
    extra_patch = {
        "ad_artefact_file": str(ad_npz_path),
        "training_split_file": str(split_csv),
    }
    if meta_csv is not None:
        extra_patch["t_half_meta_file"] = str(meta_csv)
    patch_model_json(model_file, extra_patch)
    log.info("Model JSON patched with artefact paths.")

    log.section("Run complete")
    log.info("All outputs in %s", logs_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interval-aware SVR for biodegradability half-life prediction")
    parser.add_argument(
        "--compartment",
        required=True,
        choices=["air", "water", "sediment", "soil"],
    )
    parser.add_argument(
        "--data-source",
        default="hsbd",
        choices=["hsbd", "vega", "combined"],
    )
    args = parser.parse_args()
    main(args.compartment, args.data_source)
    print("Done. Check log file for details.")
