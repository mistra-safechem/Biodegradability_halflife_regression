"""Machine-learning tools for quantised persistence endpoints.

Primary framing is interval-aware modelling: targets are treated as
*epistemically uncertain intervals*, not precise kinetic measurements.
Methodological details are documented in METHODS.md.

Design principles
-----------------
- Each class value is mapped to a half-order-of-magnitude log10 interval
    [log10(C) - 0.5*step, log10(C) + 0.5*step] where *step* is the median
    gap between adjacent class values in log10 space.
- Training typically uses interval-aware objectives (or practical proxies)
    and interval-based evaluation.
- Evaluation emphasizes interval coverage probability and rank-order
    correlation (Spearman / Kendall); RMSE / R² on raw half-life values are
    not the primary metrics for the interval workflow.

For functions used by legacy scripts, see src/legacy/ and the readme.md in that folder.

Module layout
-------------
    Interval construction   build_log10_intervals(), attach_intervals()
    Loss functions          interval_loss(), quantile_interval_loss()
    Evaluation              interval_coverage(), rank_correlation_metrics(),
                                class_accuracy(), log_evaluation_metrics()
    Preprocessing           drop_irrelevant_columns(), scale_features(),
                                decorrelate(), remove_variance_and_correlation(),
                                detect_and_remove_outliers()
    Pipeline classes        Preprocessor, PreprocessedData
    Modelling helpers       t_t_split(), svr_grid_search()
    Applicability domain    chemical_space_pca(), applicability_domain_leverage(),
                                chemical_space_morgan_pca()
"""

import matplotlib

matplotlib.use("Agg")

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from src.rdkit_tools import DESCRIPTOR_NAMES, MACCS_NAMES

log = logging.getLogger(__name__)

# Columns that encode the endpoint or are derived from it.
# They must NEVER enter the feature matrix; they are stashed separately
# in Preprocessor.t_half_meta for later retrieval / sanity-checking.
T_HALF_META_COLS: list[str] = [
    "T_half_class_days",
    "T_half_lower_bound",
    "T_half_upper_bound",
    "T_half_log10_lower",
    "T_half_log10_upper",
]

# ---------------------------------------------------------------------------
# 1. Interval construction
# ---------------------------------------------------------------------------


def build_log10_intervals(
    y: pd.Series,
    strategy: Literal["half_order"] = "half_order",
) -> pd.DataFrame:
    """Map each class value to a log10 interval using half-order-of-magnitude bounds.

    The target column contains *quantised* class values (e.g. 1.4, 5.6, 14 days).
    Each unique class value C is mapped to

        lower = log10(C) - 0.5 * step
        upper = log10(C) + 0.5 * step

    where *step* is the median gap between adjacent class values in log10 space.
    This is conservative (symmetric, reproducible, independent of sparse edges).

    Parameters
    ----------
    y:
        Series of raw half-life class values (days, positive).
    strategy:
        Only ``"half_order"`` is currently supported.

    Returns
    -------
    DataFrame with columns ``["class_value", "log10_class",
    "log10_lower", "log10_upper"]``, one row per unique class value,
    sorted ascending.

    Notes
    -----
    - Caller is responsible for ensuring ``y > 0`` before calling.
    - The mapping is computed once on the full dataset and applied uniformly
      to train and test — there is no data leakage risk.
    """
    if strategy != "half_order":
        raise NotImplementedError(f"Strategy '{strategy}' is not implemented.")

    unique_vals = np.sort(np.unique(y.dropna().values))
    if len(unique_vals) < 2:
        raise ValueError("Need at least two distinct class values to compute step size.")

    log_vals = np.log10(unique_vals)
    gaps = np.diff(log_vals)
    half_step = np.median(gaps) / 2.0

    rows = []
    for cv, lv in zip(unique_vals, log_vals):
        rows.append(
            {
                "class_value": cv,
                "log10_class": lv,
                "log10_lower": lv - half_step,
                "log10_upper": lv + half_step,
            }
        )

    return pd.DataFrame(rows)


def attach_intervals(y: pd.Series, interval_table: pd.DataFrame) -> pd.DataFrame:
    """Attach log10 lower/upper bounds to each sample.

    Parameters
    ----------
    y:
        Target series (raw class values, same index as feature matrix).
    interval_table:
        Output of :func:`build_log10_intervals`.

    Returns
    -------
    DataFrame with columns ``["y_raw", "log10_lower", "log10_upper"]``,
    index aligned with ``y``.
    """
    mapping = interval_table.set_index("class_value")[["log10_lower", "log10_upper"]]
    result = pd.DataFrame({"y_raw": y})
    result = result.join(mapping, on="y_raw")
    return result


# ---------------------------------------------------------------------------
# 2. Loss functions
# ---------------------------------------------------------------------------


def interval_loss(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
) -> np.ndarray:
    """Per-sample interval loss.

    Zero penalty when the prediction falls inside [lower, upper].
    Linear distance penalty outside.

    Parameters
    ----------
    y_pred:
        Model predictions in log10 space, shape (n,).
    y_lower, y_upper:
        Interval bounds in log10 space, shape (n,).

    Returns
    -------
    Per-sample loss values, shape (n,).
    """
    below = np.maximum(0.0, y_lower - y_pred)
    above = np.maximum(0.0, y_pred - y_upper)
    return below + above


def mean_interval_loss(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
) -> float:
    """Mean interval loss across all samples."""
    return float(np.mean(interval_loss(y_pred, y_lower, y_upper)))


def quantile_interval_loss(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
    alpha_lower: float = 0.1,
    alpha_upper: float = 0.9,
) -> np.ndarray:
    """Asymmetric quantile (pinball) loss targeting the interval endpoints.

    Models lower bound with quantile ``alpha_lower`` and upper bound with
    ``alpha_upper`` — useful when training two separate quantile regressors.

    Parameters
    ----------
    y_pred:
        Scalar predictions in log10 space.
    y_lower, y_upper:
        Interval bounds in log10 space.
    alpha_lower:
        Quantile for the lower bound (default 0.1).
    alpha_upper:
        Quantile for the upper bound (default 0.9).

    Returns
    -------
    Per-sample combined pinball loss.
    """
    res_lo = y_lower - y_pred
    loss_lo = np.where(res_lo >= 0, alpha_lower * res_lo, (alpha_lower - 1) * res_lo)

    res_hi = y_upper - y_pred
    loss_hi = np.where(res_hi >= 0, alpha_upper * res_hi, (alpha_upper - 1) * res_hi)

    return loss_lo + loss_hi


# ---------------------------------------------------------------------------
# 3. Evaluation — recommended interval metrics (see METHODS.md)
# ---------------------------------------------------------------------------


def interval_coverage(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
) -> dict:
    """Compute interval coverage probability and related statistics.

    A prediction is *covered* when it falls inside [lower, upper].

    Parameters
    ----------
    y_pred:
        Model predictions in log10 space.
    y_lower, y_upper:
        Interval bounds in log10 space (same length as y_pred).

    Returns
    -------
    dict with keys:
        coverage_probability  — fraction of predictions inside the interval
        mean_interval_loss    — average distance outside bounds (0 if inside)
        n_covered             — count inside
        n_total               — total samples
    """
    inside = (y_pred >= y_lower) & (y_pred <= y_upper)
    return {
        "coverage_probability": float(inside.mean()),
        "mean_interval_loss": mean_interval_loss(y_pred, y_lower, y_upper),
        "n_covered": int(inside.sum()),
        "n_total": len(y_pred),
    }


def rank_correlation_metrics(
    y_true_log10: np.ndarray,
    y_pred: np.ndarray,
) -> dict:
    """Rank-order correlation between true log10 class centroids and predictions.

    Uses Spearman's rho and Kendall's tau — appropriate for ordinal / quantised
    targets where absolute distance is not meaningful.

    Parameters
    ----------
    y_true_log10:
        True values in log10 space (class centroids, i.e. log10(class_value)).
    y_pred:
        Model predictions in log10 space.

    Returns
    -------
    dict with keys: spearman_r, spearman_p, kendall_tau, kendall_p
    """
    from scipy.stats import kendalltau, spearmanr

    sp_r, sp_p = spearmanr(y_true_log10, y_pred)
    kt, kt_p = kendalltau(y_true_log10, y_pred)
    return {
        "spearman_r": float(sp_r),
        "spearman_p": float(sp_p),
        "kendall_tau": float(kt),
        "kendall_p": float(kt_p),
    }


def class_accuracy(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
) -> float:
    """Fraction of predictions whose rounded class assignment matches the true class.

    Equivalent to :func:`interval_coverage` ``coverage_probability`` but named
    explicitly as class-level accuracy.
    """
    return float(((y_pred >= y_lower) & (y_pred <= y_upper)).mean())


def log_evaluation_metrics(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
    y_true_log10: np.ndarray,
    logger: Optional[logging.Logger] = None,
    prefix: str = "",
) -> dict:
    """Compute and log all recommended evaluation metrics.

    Parameters
    ----------
    y_pred:
        Predictions in log10 space.
    y_lower, y_upper:
        Per-sample interval bounds in log10 space.
    y_true_log10:
        True log10 class centroids (for rank correlation).
    logger:
        Logger instance.  Prints to stdout if None.
    prefix:
        Optional string prepended to each log line (e.g. compartment name).

    Returns
    -------
    dict containing all metric values.
    """
    cov = interval_coverage(y_pred, y_lower, y_upper)
    rc = rank_correlation_metrics(y_true_log10, y_pred)
    acc = class_accuracy(y_pred, y_lower, y_upper)

    mil_log10 = cov["mean_interval_loss"]
    mil_days = float(10**mil_log10)  # geometric-mean distance in days
    metrics = {**cov, **rc, "class_accuracy": acc, "mean_interval_loss_days": mil_days}

    lines = [
        f"{prefix}Interval coverage probability : {cov['coverage_probability']:.3f} ({cov['n_covered']}/{cov['n_total']})",
        f"{prefix}Mean interval loss            : {mil_log10:.4f} log10 days  (~{mil_days:.2f} days geom. equiv.)",
        f"{prefix}Class accuracy                : {acc:.3f}",
        f"{prefix}Spearman rho                  : {rc['spearman_r']:.3f}  (p={rc['spearman_p']:.4f})",
        f"{prefix}Kendall tau                   : {rc['kendall_tau']:.3f}  (p={rc['kendall_p']:.4f})",
    ]
    for line in lines:
        if logger:
            logger.info(line)
        else:
            print(line)

    return metrics


# ---------------------------------------------------------------------------
# 4. Preprocessing
# ---------------------------------------------------------------------------


def drop_irrelevant_columns(df: pd.DataFrame, to_drop: set[str] = frozenset()) -> pd.DataFrame:
    """Drop non-feature columns and optionally entire feature sets.

    Parameters
    ----------
    df:
        Raw dataframe from the database.
    to_drop:
        Set of feature-set labels to drop.  Recognised values:
        ``"maccs"`` (drops MACCS_* columns) and ``"rdkit"`` (drops RDKit
        descriptor columns).  Pass an empty set to keep all features.

    Returns
    -------
    Dataframe with irrelevant columns removed.
    """
    # Always drop: row identifiers, SMILES string, source reference, and all
    # T_half auxiliary columns that are derived from / encode the endpoint.
    meta_cols = ["id", "Canonical_smiles", "reference"] + T_HALF_META_COLS
    df = df.drop(columns=meta_cols, errors="ignore")

    if "maccs" in to_drop:
        df = df.drop(columns=MACCS_NAMES, errors="ignore")
    if "rdkit" in to_drop:
        df = df.drop(columns=DESCRIPTOR_NAMES, errors="ignore")

    if df.shape[1] == 0:
        raise ValueError("No feature columns remain after dropping.  Check the 'to_drop' argument.")
    return df


def detect_and_remove_outliers(X: pd.DataFrame, y: pd.Series, contamination: float = 0.05) -> tuple[pd.DataFrame, pd.Series]:
    """Remove outliers from the feature matrix using IsolationForest.

    Parameters
    ----------
    X:
        Feature matrix.
    y:
        Target series (aligned index with X).
    contamination:
        Expected proportion of outliers (default 0.05).

    Returns
    -------
    (X_clean, y_clean) with outliers removed, index reset.
    """
    iso = IsolationForest(contamination=contamination, random_state=42)
    mask = iso.fit_predict(X) == 1
    n_removed = (~mask).sum()
    log.info("IsolationForest removed %d outliers (contamination=%.2f).", n_removed, contamination)
    return X[mask].reset_index(drop=True), y[mask].reset_index(drop=True)


def scale_features(X: pd.DataFrame) -> pd.DataFrame:
    """StandardScaler on RDKit descriptors; MACCS bits are left unscaled.

    When both descriptor and MACCS columns are present, only the descriptors
    are scaled.  If only MACCS bits are present, the dataframe is returned
    unchanged (binary features do not benefit from scaling).

    Infinite values (which RDKit occasionally produces for degenerate
    molecular geometries) are replaced with NaN before scaling so that
    sklearn does not raise a ValueError.  Columns that are *entirely* NaN
    after this substitution are dropped; remaining NaNs are imputed with the
    column median.

    Returns a copy with reset index.
    """
    # Replace ±inf with NaN so sklearn validators don't reject the matrix
    X = X.replace([np.inf, -np.inf], np.nan)

    scaler = StandardScaler()
    has_maccs = any(c in X.columns for c in MACCS_NAMES)
    has_rdkit = any(c in X.columns for c in DESCRIPTOR_NAMES)

    if has_rdkit and has_maccs:
        maccs_df = X[MACCS_NAMES].reset_index(drop=True)
        rdkit_df = X.drop(columns=MACCS_NAMES).reset_index(drop=True)
        # Drop columns that are entirely NaN; impute remaining NaNs with median
        rdkit_df = rdkit_df.dropna(axis=1, how="all")
        rdkit_df = rdkit_df.fillna(rdkit_df.median())
        rdkit_scaled = pd.DataFrame(scaler.fit_transform(rdkit_df), columns=rdkit_df.columns)
        return pd.concat([rdkit_scaled, maccs_df], axis=1)

    if has_rdkit:
        X_rdkit = X.dropna(axis=1, how="all")
        X_rdkit = X_rdkit.fillna(X_rdkit.median())
        return pd.DataFrame(scaler.fit_transform(X_rdkit), columns=X_rdkit.columns).reset_index(drop=True)

    # MACCS only — no scaling
    return X.copy().reset_index(drop=True)


def decorrelate(X: pd.DataFrame, target_column: str, threshold: float = 0.95) -> list[str]:
    """Return a list of columns to drop to bring pairwise correlation below threshold.

    The target column is excluded from the correlation calculation.
    When two features are correlated above *threshold*, the second one
    (by column order) is marked for removal.
    """
    X_num = X.select_dtypes(include=[np.number]).copy()
    if target_column in X_num.columns:
        X_num = X_num.drop(columns=[target_column])

    corr = X_num.corr().abs()
    columns = X_num.columns
    to_drop: set[str] = set()

    for i in range(len(columns) - 1):
        for j in range(i + 1, len(columns)):
            if corr.at[columns[i], columns[j]] > threshold and columns[i] not in to_drop:
                to_drop.add(columns[j])
    return list(to_drop)


def remove_variance_and_correlation(X: pd.DataFrame, target_column: str, threshold: float = 0.95) -> pd.DataFrame:
    """Drop zero-variance columns, then highly correlated columns.

    Parameters
    ----------
    X:
        Scaled feature matrix (target column may or may not be present).
    target_column:
        Name of the target column to exclude from correlation analysis.
    threshold:
        Pearson |r| above which one of a correlated pair is dropped.

    Returns
    -------
    Pruned dataframe.
    """
    zero_std = X.columns[X.std() == 0]
    X = X.drop(columns=zero_std)
    log.debug("Dropped %d zero-variance columns.", len(zero_std))

    to_drop = decorrelate(X, target_column, threshold=threshold)
    X = X.drop(columns=to_drop)
    log.info(
        "After variance + correlation pruning: %d features, %d samples.",
        X.shape[1],
        X.shape[0],
    )
    return X


# ---------------------------------------------------------------------------
# 5. Preprocessor pipeline
# ---------------------------------------------------------------------------


class Preprocessor:
    """End-to-end preprocessing pipeline for interval-aware modelling.

    Produces per-sample log10 interval bounds (``y_lower``, ``y_upper``) and
    the log10 class centroid (``y_log10``) alongside the feature matrix.

    Interval bounds are taken from the DB-stored columns
    ``T_half_log10_lower`` / ``T_half_log10_upper`` when they are present and
    non-null.  If those columns are absent or all-null, bounds are recomputed
    from scratch via :func:`build_log10_intervals` and a warning is logged.

    All T_half auxiliary columns (``T_half_class_days``, ``T_half_lower_bound``,
    ``T_half_upper_bound``, ``T_half_log10_lower``, ``T_half_log10_upper``) are
    stripped from the feature matrix before training.  They are preserved in
    ``self.t_half_meta`` (row-aligned with ``X`` and ``y``) for later
    retrieval, sanity-checking, or CSV export.

    Parameters
    ----------
    to_drop:
        Feature sets to remove.  ``"maccs"`` and/or ``"rdkit"`` are
        recognised; pass ``{"None"}`` or empty set to keep everything.
    target_column:
        Name of the raw half-life column (days).
    remove_outliers:
        Whether to apply IsolationForest outlier removal.
    corr_threshold:
        Pearson |r| threshold for the decorrelation step.

    Attributes set after :meth:`preprocess`
    ----------------------------------------
    X            : pd.DataFrame        — final feature matrix (index 0..N-1)
    y            : pd.Series           — raw target values (days), aligned
    y_log10      : pd.Series           — log10 class centroids
    y_lower      : pd.Series           — log10 interval lower bounds
    y_upper      : pd.Series           — log10 interval upper bounds
    interval_table : pd.DataFrame      — mapping from class_value to bounds
    smiles       : pd.Series | None    — Canonical_smiles, aligned with X
    t_half_meta  : pd.DataFrame        — stashed T_half auxiliary columns,
                                         row-aligned with X/y; may be saved
                                         as CSV for traceability
    """

    def __init__(
        self,
        to_drop: set[str] = frozenset({"None"}),
        target_column: str = "T_half_days",
        remove_outliers: bool = True,
        corr_threshold: float = 0.95,
    ) -> None:
        self.to_drop = to_drop
        self.target_column = target_column
        self.remove_outliers = remove_outliers
        self.corr_threshold = corr_threshold

        # Set after preprocess()
        self.X: Optional[pd.DataFrame] = None
        self.y: Optional[pd.Series] = None
        self.y_log10: Optional[pd.Series] = None
        self.y_lower: Optional[pd.Series] = None
        self.y_upper: Optional[pd.Series] = None
        self.interval_table: Optional[pd.DataFrame] = None
        self.smiles: Optional[pd.Series] = None
        self.t_half_meta: Optional[pd.DataFrame] = None

    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        """Run the full preprocessing pipeline.

        Parameters
        ----------
        df:
            Raw dataframe as returned by :func:`src.db_utils.get_all_data`.

        Returns
        -------
        DataFrame containing the final feature columns plus the target column.
        All attributes listed in the class docstring are populated.
        """
        # --- Stash auxiliary T_half columns BEFORE any dropping ---
        # These encode the endpoint; they must not enter the feature matrix.
        present_meta = [c for c in T_HALF_META_COLS if c in df.columns]
        t_half_meta = df[present_meta].copy().reset_index(drop=True)

        # Preserve SMILES before column dropping
        smiles = df["Canonical_smiles"].copy().reset_index(drop=True) if "Canonical_smiles" in df.columns else None

        # --- Feature selection (T_HALF_META_COLS dropped inside here) ---
        df_clean = drop_irrelevant_columns(df, self.to_drop)

        X = df_clean.drop(columns=[self.target_column])
        y = df_clean[self.target_column].reset_index(drop=True)
        log.info("Raw dataset: %d features, %d samples.", X.shape[1], X.shape[0])

        # Replace ±inf with NaN (RDKit can emit inf for degenerate molecules)
        X = X.reset_index(drop=True).replace([np.inf, -np.inf], np.nan)

        # Drop rows with NaN in features or target; apply mask to meta/smiles too
        nan_mask = ~(X.isna().any(axis=1) | y.isna())
        n_dropped = int((~nan_mask).sum())
        X = X[nan_mask].reset_index(drop=True)
        y = y[nan_mask].reset_index(drop=True)
        t_half_meta = t_half_meta[nan_mask].reset_index(drop=True)
        if smiles is not None:
            smiles = smiles[nan_mask].reset_index(drop=True)
        log.info(
            "After NaN/inf removal: %d features, %d samples (%d rows dropped).",
            X.shape[1],
            X.shape[0],
            n_dropped,
        )

        # --- Interval construction ---
        # Always build interval_table from T_half_days for the CSV export /
        # sanity-check, regardless of which bounds source is used.
        interval_table = build_log10_intervals(y)
        half_step = float((interval_table["log10_upper"] - interval_table["log10_lower"]).median() / 2)
        log.info(
            "Interval table built: %d unique class values, half-step=%.4f log10 days (~%.2f× in linear scale).",
            len(interval_table),
            half_step,
            10**half_step,
        )

        # Prefer DB-stored log10 bounds; fall back to recomputed bounds.
        db_lower = t_half_meta.get("T_half_log10_lower")
        db_upper = t_half_meta.get("T_half_log10_upper")
        use_db_bounds = db_lower is not None and db_upper is not None and db_lower.notna().any() and db_upper.notna().any()
        if use_db_bounds:
            log.info("Using DB-stored log10 interval bounds (T_half_log10_lower/upper).")
        else:
            log.warning("DB-stored log10 bounds absent or all-null; recomputing from T_half_days via build_log10_intervals().")

        # --- Outlier removal ---
        if self.remove_outliers:
            X, y = detect_and_remove_outliers(X, y)
            outlier_mask = X.index  # integer positions kept after IsolationForest
            t_half_meta = t_half_meta.loc[outlier_mask].reset_index(drop=True)
            if smiles is not None:
                smiles = smiles.loc[outlier_mask].reset_index(drop=True)
            log.info("After outlier removal: %d features, %d samples.", X.shape[1], X.shape[0])
        else:
            log.info("Outlier removal skipped.")

        # --- Scaling and decorrelation ---
        X_scaled = scale_features(X)
        X_final = remove_variance_and_correlation(X_scaled, self.target_column, self.corr_threshold)
        log.info(
            "Final preprocessed dataset: %d features, %d samples.",
            X_final.shape[1],
            X_final.shape[0],
        )

        # --- Resolve interval bounds ---
        y_reset = y.reset_index(drop=True)
        if use_db_bounds:
            # Re-index t_half_meta after all row filters
            y_lower = t_half_meta["T_half_log10_lower"].reset_index(drop=True).rename("log10_lower")
            y_upper = t_half_meta["T_half_log10_upper"].reset_index(drop=True).rename("log10_upper")
            # Fill any remaining NaN cells using the recomputed table as fallback
            if y_lower.isna().any() or y_upper.isna().any():
                fallback = attach_intervals(y_reset, interval_table)
                y_lower = y_lower.fillna(fallback["log10_lower"].reset_index(drop=True))
                y_upper = y_upper.fillna(fallback["log10_upper"].reset_index(drop=True))
                log.warning("Some rows had null DB bounds; filled with recomputed values.")
        else:
            intervals = attach_intervals(y_reset, interval_table)
            y_lower = intervals["log10_lower"].reset_index(drop=True)
            y_upper = intervals["log10_upper"].reset_index(drop=True)

        self.X = X_final.reset_index(drop=True)
        self.y = y_reset
        self.y_log10 = np.log10(y_reset).rename("y_log10")
        self.y_lower = y_lower
        self.y_upper = y_upper
        self.interval_table = interval_table
        self.smiles = smiles.reset_index(drop=True) if smiles is not None else None
        self.t_half_meta = t_half_meta.reset_index(drop=True)

        df_out = self.X.copy()
        df_out[self.target_column] = self.y
        return df_out

    def get_smiles_for_split(self, X_train: pd.DataFrame, X_test: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        """Return SMILES series aligned with train/test index splits."""
        return self.smiles.loc[X_train.index], self.smiles.loc[X_test.index]


# ---------------------------------------------------------------------------
# 6. Data container
# ---------------------------------------------------------------------------


@dataclass
class PreprocessedData:
    """Container for preprocessed data ready for modelling.

    Attributes
    ----------
    name:
        Compartment label (e.g. ``"soil"``).
    df:
        Full preprocessed dataframe (features + target).
    remove_outliers:
        Whether outlier removal was applied.
    X:
        Feature matrix.
    y_log10:
        Log10 class centroids (point target; used only for rank correlation).
    y_lower:
        Per-sample log10 interval lower bounds.
    y_upper:
        Per-sample log10 interval upper bounds.
    interval_table:
        DataFrame mapping class_value → log10 bounds.
    """

    name: str
    df: pd.DataFrame
    remove_outliers: bool = False
    X: Optional[pd.DataFrame] = None
    y_log10: Optional[pd.Series] = None
    y_lower: Optional[pd.Series] = None
    y_upper: Optional[pd.Series] = None
    interval_table: Optional[pd.DataFrame] = field(default=None)


# ---------------------------------------------------------------------------
# 7. Train / test split
# ---------------------------------------------------------------------------


def t_t_split(
    X: pd.DataFrame,
    y_log10: pd.Series,
    y_lower: pd.Series,
    y_upper: pd.Series,
    test_size: float = 0.2,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, pd.Series, pd.Series, pd.Series, pd.Series]:
    """Split data into train / test sets, dropping any constant features.

    Parameters
    ----------
    X:
        Feature matrix.
    y_log10, y_lower, y_upper:
        Target series (all aligned with X by integer index).
    test_size:
        Fraction of samples for the test set.
    random_state:
        Random seed.

    Returns
    -------
    X_train, X_test, y_log10_train, y_log10_test,
    y_lower_train, y_lower_test, y_upper_train, y_upper_test
    """
    idx = X.index
    train_idx, test_idx = train_test_split(idx, test_size=test_size, random_state=random_state)

    X_train, X_test = X.loc[train_idx], X.loc[test_idx]

    # Remove constant features (can appear after splitting)
    const = [c for c in X_train.columns if X_train[c].nunique() == 1 or X_test[c].nunique() == 1]
    if const:
        log.debug("Removing %d constant features after split: %s", len(const), const)
    X_train = X_train.drop(columns=const)
    X_test = X_test.drop(columns=const)

    return (
        X_train,
        X_test,
        y_log10.loc[train_idx],
        y_log10.loc[test_idx],
        y_lower.loc[train_idx],
        y_lower.loc[test_idx],
        y_upper.loc[train_idx],
        y_upper.loc[test_idx],
    )


# ---------------------------------------------------------------------------
# 8. SVR grid search
# ---------------------------------------------------------------------------


def svr_grid_search(
    X_train: pd.DataFrame,
    y_train: pd.Series,
) -> tuple[SVR, dict]:
    """Grid search over SVR hyperparameters.

    Scoring is ``neg_mean_squared_error`` on the log10 class centroid target.
    This is a practical proxy; the final model is evaluated with interval
    metrics rather than MSE.

    Returns
    -------
    (best_svr, best_params)
    """
    param_grid = {
        "C": [0.1, 1, 10, 100],
        "epsilon": [0.01, 0.1, 0.2, 0.5],
        "kernel": ["rbf", "linear"],
        "gamma": ["scale", "auto"],
    }
    grid = GridSearchCV(SVR(), param_grid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1)
    grid.fit(X_train, y_train)
    log.info("SVR grid search best params: %s", grid.best_params_)
    return grid.best_estimator_, grid.best_params_


# ---------------------------------------------------------------------------
# 9. Applicability domain
# ---------------------------------------------------------------------------


def chemical_space_pca(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    output_dir: Path,
    compartment: str,
) -> tuple[PCA, StandardScaler, np.ndarray, np.ndarray, float]:
    """2-component PCA chemical-space coverage plot.

    Fits PCA on X_train, projects X_test, plots scatter with bounding box,
    and returns coverage percentage.

    Returns
    -------
    (pca, scaler, train_pcs, test_pcs, coverage_pct)
    """
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    pca = PCA(n_components=2, random_state=42)
    train_pcs = pca.fit_transform(X_train_s)
    test_pcs = pca.transform(X_test_s)

    var = pca.explained_variance_ratio_
    log.info(
        "PCA explained variance: PC1=%.3f, PC2=%.3f (total=%.3f)",
        var[0],
        var[1],
        sum(var),
    )

    pc1_min, pc1_max = train_pcs[:, 0].min(), train_pcs[:, 0].max()
    pc2_min, pc2_max = train_pcs[:, 1].min(), train_pcs[:, 1].max()
    inside = (
        (test_pcs[:, 0] >= pc1_min) & (test_pcs[:, 0] <= pc1_max) & (test_pcs[:, 1] >= pc2_min) & (test_pcs[:, 1] <= pc2_max)
    )
    coverage_pct = float(inside.mean() * 100)
    log.info("Test coverage (PCA bounding box): %.1f%%", coverage_pct)

    output_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 6))
    plt.scatter(train_pcs[:, 0], train_pcs[:, 1], alpha=0.4, s=20, color="royalblue", label=f"Train (n={len(X_train)})")
    plt.scatter(
        test_pcs[:, 0], test_pcs[:, 1], alpha=0.7, s=30, color="darkorange", marker="^", label=f"Test (n={len(X_test)})"
    )
    rect = plt.Rectangle(
        (pc1_min, pc2_min),
        pc1_max - pc1_min,
        pc2_max - pc2_min,
        linewidth=1.5,
        edgecolor="royalblue",
        facecolor="none",
        linestyle="--",
        label="Train bounding box",
    )
    plt.gca().add_patch(rect)
    plt.xlabel(f"PC1 ({var[0] * 100:.1f}% var)", fontsize=11)
    plt.ylabel(f"PC2 ({var[1] * 100:.1f}% var)", fontsize=11)
    plt.title(f"Chemical Space Coverage — PCA ({compartment})", fontsize=13)
    plt.legend(fontsize=9)
    plt.tight_layout()
    out_path = output_dir / f"chemical_space_pca_{compartment}.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    log.info("PCA chemical space plot saved to %s", out_path)

    return pca, scaler, train_pcs, test_pcs, coverage_pct


def applicability_domain_leverage(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_lower_train: pd.Series,
    y_upper_train: pd.Series,
    y_lower_test: pd.Series,
    y_upper_test: pd.Series,
    y_pred_train: np.ndarray,
    y_pred_test: np.ndarray,
    output_dir: Path,
    compartment: str,
) -> dict:
    """Leverage-based Applicability Domain assessment (Williams plot).

    AD boundary:
        h* = 3p / n   (p = features, n = training samples)
        prediction inside interval  (replaces ±3 std-residual criterion)

    The Williams plot shows leverage (x) vs whether the prediction is inside
    or outside the target interval (y-axis: signed distance from nearest bound,
    zero when inside).

    Parameters
    ----------
    X_train, X_test:
        Unscaled feature matrices.
    y_lower_train, y_upper_train:
        Interval bounds for training samples in log10 space.
    y_lower_test, y_upper_test:
        Interval bounds for test samples in log10 space.
    y_pred_train, y_pred_test:
        SVR predictions in log10 space.
    output_dir:
        Directory for the Williams plot image.
    compartment:
        Label for plot title and filename.

    Returns
    -------
    dict with keys:

    h_star              float        AD leverage boundary (3p/n)
    h_train             np.ndarray   per-training-sample leverage
    h_test              np.ndarray   per-test-sample leverage
    inside_train        np.ndarray   bool — prediction inside interval (train)
    inside_test         np.ndarray   bool — prediction inside interval (test)
    pct_train_inside    float        % training predictions inside interval
    pct_test_inside     float        % test predictions inside interval
    interval_loss_train np.ndarray   per-training-sample interval loss
    interval_loss_test  np.ndarray   per-test-sample interval loss
    XtX_inv             np.ndarray   (p×p) hat-matrix inverse for new-compound leverage
    X_train_scaled      np.ndarray   standardised training feature matrix (n×p)
    X_train_mean        np.ndarray   column means used by the internal StandardScaler
    X_train_std         np.ndarray   column stds  used by the internal StandardScaler
    X_train             pd.DataFrame original (unscaled) training feature matrix
    """
    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_train)
    X_te_s = scaler.transform(X_test)

    n, p = X_tr_s.shape
    h_star = 3.0 * p / n

    XtX_inv = np.linalg.pinv(X_tr_s.T @ X_tr_s)
    h_train = (X_tr_s @ XtX_inv * X_tr_s).sum(axis=1)
    h_test = (X_te_s @ XtX_inv * X_te_s).sum(axis=1)

    # AD membership: inside interval
    il_train = interval_loss(y_pred_train, y_lower_train.values, y_upper_train.values)
    il_test = interval_loss(y_pred_test, y_lower_test.values, y_upper_test.values)
    inside_train = il_train == 0.0
    inside_test = il_test == 0.0

    pct_train = float(inside_train.mean() * 100)
    pct_test = float(inside_test.mean() * 100)

    log.info(
        "AD (Williams) [%s] — h*=%.4f  Train inside=%.1f%%  Test inside=%.1f%%",
        compartment,
        h_star,
        pct_train,
        pct_test,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    # Signed distance from nearest bound (0 when inside)
    def signed_dist(pred, lo, hi):
        below = lo - pred
        above = pred - hi
        return np.where(pred < lo, -below, np.where(pred > hi, above, 0.0))

    sd_train = signed_dist(y_pred_train, y_lower_train.values, y_upper_train.values)
    sd_test = signed_dist(y_pred_test, y_lower_test.values, y_upper_test.values)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(h_train, sd_train, alpha=0.5, s=20, color="royalblue", label=f"Train (n={n})")
    ax.scatter(h_test, sd_test, alpha=0.7, s=30, color="darkorange", marker="^", label=f"Test (n={len(X_test)})")
    ax.axvline(h_star, color="gray", linestyle="--", linewidth=1.5, label=f"h* = {h_star:.3f}")
    ax.axhline(0, color="green", linestyle="-", linewidth=0.8, label="Interval boundary (=0 inside)")
    ax.set_xlabel("Leverage (h)", fontsize=11)
    ax.set_ylabel("Signed distance from interval (log10 days)", fontsize=11)
    ax.set_title(f"Williams Plot — Applicability Domain ({compartment})", fontsize=13)
    ax.legend(fontsize=9)

    # Secondary Y-axis in days (non-linear mapping: days = 10^(log10_dist))
    # We annotate with representative day-equivalent ticks derived from the
    # primary axis range so the reader can interpret distances in familiar units.
    y_lo, y_hi = ax.get_ylim()
    # Build secondary axis via a twin that shares the same x but maps y → days
    ax2 = ax.twinx()
    ax2.set_ylim(y_lo, y_hi)
    # Choose representative log10 tick positions from the primary axis
    primary_ticks = ax.get_yticks()
    valid_ticks = primary_ticks[(primary_ticks >= y_lo) & (primary_ticks <= y_hi)]
    ax2.set_yticks(valid_ticks)
    ax2.set_yticklabels(
        [f"{10**t:.2g} d" if t != 0 else "0 d" for t in valid_ticks],
        fontsize=8,
    )
    ax2.set_ylabel("Signed distance (days equiv.)", fontsize=10)

    fig.tight_layout()
    out_path = output_dir / f"williams_plot_{compartment}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Williams plot saved to %s", out_path)

    return {
        "h_star": h_star,
        "h_train": h_train,
        "h_test": h_test,
        "inside_train": inside_train,
        "inside_test": inside_test,
        "pct_train_inside": pct_train,
        "pct_test_inside": pct_test,
        "interval_loss_train": il_train,
        "interval_loss_test": il_test,
        # --- inference artefacts ---
        "XtX_inv": XtX_inv,
        "X_train_scaled": X_tr_s,
        "X_train_mean": scaler.mean_,
        "X_train_std": scaler.scale_,
        "X_train": X_train,
    }


def chemical_space_morgan_pca(
    smiles_train: pd.Series,
    smiles_test: pd.Series,
    ad_results: dict,
    output_dir: Path,
    compartment: str,
) -> dict:
    """PCA on Morgan fingerprints coloured by AD membership + Butina clustering.

    Unchanged from the previous version in intent; uses the updated
    ``inside_train`` / ``inside_test`` boolean arrays from
    :func:`applicability_domain_leverage`.

    Parameters
    ----------
    smiles_train, smiles_test:
        SMILES series aligned with X_train / X_test by index.
    ad_results:
        Dict returned by :func:`applicability_domain_leverage`.
    output_dir:
        Directory for output images.
    compartment:
        Label for plot titles and filenames.

    Returns
    -------
    dict with Butina cluster statistics for AD-outside test compounds.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    from rdkit.ML.Cluster import Butina

    def _to_fp(smi: str):
        mol = Chem.MolFromSmiles(smi)
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024) if mol else None

    def _fp_arr(fp) -> np.ndarray:
        arr = np.zeros((1024,), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    fps_tr_raw = [_to_fp(s) for s in smiles_train]
    fps_te_raw = [_to_fp(s) for s in smiles_test]
    valid_tr = [i for i, fp in enumerate(fps_tr_raw) if fp is not None]
    valid_te = [i for i, fp in enumerate(fps_te_raw) if fp is not None]
    fps_tr = [fps_tr_raw[i] for i in valid_tr]
    fps_te = [fps_te_raw[i] for i in valid_te]
    inside_train = ad_results["inside_train"][valid_tr]
    inside_test = ad_results["inside_test"][valid_te]

    mat_tr = np.array([_fp_arr(fp) for fp in fps_tr])
    mat_te = np.array([_fp_arr(fp) for fp in fps_te])

    pca = PCA(n_components=2, random_state=42)
    tr_pcs = pca.fit_transform(mat_tr)
    te_pcs = pca.transform(mat_te)
    var = pca.explained_variance_ratio_

    output_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(9, 7))
    plt.scatter(
        tr_pcs[inside_train, 0],
        tr_pcs[inside_train, 1],
        alpha=0.3,
        s=18,
        color="royalblue",
        label=f"Train in AD (n={inside_train.sum()})",
    )
    plt.scatter(
        tr_pcs[~inside_train, 0],
        tr_pcs[~inside_train, 1],
        alpha=0.7,
        s=35,
        color="steelblue",
        marker="s",
        label=f"Train out AD (n={(~inside_train).sum()})",
    )
    plt.scatter(
        te_pcs[inside_test, 0],
        te_pcs[inside_test, 1],
        alpha=0.5,
        s=30,
        color="darkorange",
        marker="^",
        label=f"Test in AD (n={inside_test.sum()})",
    )
    plt.scatter(
        te_pcs[~inside_test, 0],
        te_pcs[~inside_test, 1],
        alpha=0.9,
        s=55,
        color="crimson",
        marker="^",
        label=f"Test out AD (n={(~inside_test).sum()})",
    )
    plt.xlabel(f"PC1 ({var[0] * 100:.1f}% var)", fontsize=11)
    plt.ylabel(f"PC2 ({var[1] * 100:.1f}% var)", fontsize=11)
    plt.title(f"Chemical Space — Morgan FP PCA — AD membership ({compartment})", fontsize=13)
    plt.legend(fontsize=9)
    plt.tight_layout()
    pca_path = output_dir / f"morgan_pca_ad_{compartment}.png"
    plt.savefig(pca_path, dpi=150)
    plt.close()
    log.info("Morgan FP PCA (AD coloured) saved to %s", pca_path)

    fps_outside = [fps_te[i] for i in range(len(fps_te)) if not inside_test[i]]
    n_outside = len(fps_outside)
    cluster_stats: dict[str, Any] = {"n_outside_ad": n_outside}

    if n_outside >= 2:
        dists = []
        for i in range(1, n_outside):
            sims = DataStructs.BulkTanimotoSimilarity(fps_outside[i], fps_outside[:i])
            dists.extend([1 - x for x in sims])
        clusters = Butina.ClusterData(dists, n_outside, 0.5, isDistData=True)
        n_singletons = sum(1 for c in clusters if len(c) == 1)
        cluster_sizes = sorted([len(c) for c in clusters], reverse=True)
        singleton_pct = n_singletons / n_outside * 100

        log.info(
            "Butina (cutoff=0.5) — %d AD-outside test compounds: %d clusters, %d singletons (%.1f%%)",
            n_outside,
            len(clusters),
            n_singletons,
            singleton_pct,
        )
        cluster_stats.update(
            {
                "n_clusters": len(clusters),
                "n_singletons": n_singletons,
                "singleton_pct": singleton_pct,
                "cluster_sizes": cluster_sizes,
            }
        )

        from matplotlib.patches import Patch

        colors_bar = ["#d62728" if s == 1 else "#1f77b4" for s in cluster_sizes]
        fig, ax = plt.subplots(figsize=(max(6, len(clusters) * 0.4 + 2), 4))
        ax.bar(range(len(cluster_sizes)), cluster_sizes, color=colors_bar, edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Cluster index (sorted by size)", fontsize=10)
        ax.set_ylabel("Cluster size (# compounds)", fontsize=10)
        ax.set_title(
            f"Butina cluster sizes — AD-outside test compounds ({compartment})\n"
            f"Tanimoto distance cutoff=0.5 | {len(clusters)} clusters | "
            f"{n_singletons} singletons ({singleton_pct:.1f}%)",
            fontsize=11,
        )
        ax.legend(
            handles=[
                Patch(color="#1f77b4", label="Multi-member"),
                Patch(color="#d62728", label="Singleton"),
            ],
            fontsize=9,
        )
        ax.set_xticks(range(len(cluster_sizes)))
        ax.set_xticklabels([str(i + 1) for i in range(len(cluster_sizes))], fontsize=7)
        plt.tight_layout()
        bar_path = output_dir / f"butina_cluster_sizes_{compartment}.png"
        plt.savefig(bar_path, dpi=150)
        plt.close()

        outside_idx = [i for i in range(len(fps_te)) if not inside_test[i]]
        cluster_label = np.full(n_outside, -1, dtype=int)
        for c_idx, c in enumerate(clusters):
            for m in c:
                cluster_label[m] = c_idx
        pcs_out = te_pcs[np.array(outside_idx)]
        n_multi = len(clusters) - n_singletons
        cmap = plt.get_cmap("tab20")

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(tr_pcs[:, 0], tr_pcs[:, 1], s=12, alpha=0.15, color="grey", label="Train")
        ax.scatter(te_pcs[inside_test, 0], te_pcs[inside_test, 1], s=18, alpha=0.25, color="darkorange", label="Test in AD")
        singleton_mask = np.array([cluster_label[i] >= n_multi for i in range(n_outside)])
        ax.scatter(
            pcs_out[singleton_mask, 0],
            pcs_out[singleton_mask, 1],
            s=70,
            alpha=0.9,
            color="crimson",
            marker="x",
            linewidths=1.5,
            label=f"Singleton ({n_singletons})",
        )
        for c_idx in range(n_multi):
            mask = cluster_label == c_idx
            if not mask.any():
                continue
            ax.scatter(
                pcs_out[mask, 0],
                pcs_out[mask, 1],
                s=70,
                alpha=0.85,
                color=cmap(c_idx % 20),
                marker="^",
                label=f"Cluster {c_idx + 1} (n={mask.sum()})",
            )
        ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}% var)", fontsize=10)
        ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}% var)", fontsize=10)
        ax.set_title(f"Butina Clusters — AD-outside Test Compounds ({compartment})", fontsize=12)
        ax.legend(fontsize=8, markerscale=1.2)
        plt.tight_layout()
        cluster_pca_path = output_dir / f"butina_cluster_pca_{compartment}.png"
        plt.savefig(cluster_pca_path, dpi=150)
        plt.close()
        log.info("Butina cluster PCA saved to %s", cluster_pca_path)
    else:
        log.info("Too few AD-outside test compounds for Butina clustering (n=%d).", n_outside)

    return cluster_stats
