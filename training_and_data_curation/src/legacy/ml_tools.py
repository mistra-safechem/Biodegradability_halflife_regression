"""_summary_

This module contains all functions related to machine learning,
including preprocessing, outlier detection, feature selection, and modeling.
The main class is `Preprocessor`, which handles all preprocessing steps in a single method.
There are also functions for outlier detection using IsolationForest and LocalOutlierFactor,
as well as a function for decorrelating features based on correlation threshold.
The module also includes functions for splitting data into train and test sets,
outputting regression metrics, and performing grid search for SVR hyperparameters.

Additionaly modules after March 2026 might have ended up in separate scripts due to developer inconsistency.

LEGACY: copied to src/legacy for backwardscompatibility, see also the readme.md in this folder.

"""

import matplotlib

matplotlib.use("Agg")

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from sklearn.metrics import explained_variance_score, mean_absolute_error, mean_squared_error, r2_score, root_mean_squared_error
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.neighbors import LocalOutlierFactor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from src.rdkit_tools import DESCRIPTOR_NAMES, MACCS_NAMES


def drop_irrelevant_columns(df, to_drop: set[str] = {"maccs", "rdkit"}) -> pd.DataFrame:
    # drop columns that are not relevant for ML
    cols_to_drop = ["id", "Canonical_smiles", "reference"]
    df = df.drop(columns=cols_to_drop, errors="ignore")
    # drop specified feature sets
    if "maccs" in to_drop:
        cols_to_drop = MACCS_NAMES
        df = df.drop(columns=cols_to_drop, errors="ignore")
    if "rdkit" in to_drop:
        cols_to_drop = DESCRIPTOR_NAMES
        df = df.drop(columns=cols_to_drop, errors="ignore")
    # raise error if no features left
    if df.shape[1] == 0:
        raise ValueError("No features left - perhaps 'to_drop' parameter forgotten or used incorrectly?")
    return df


# --------------------------------------------------------------
# Outlier detection methods (not used, using sklearn instead)
# --------------------------------------------------------------
def detect_outliers_iqr(df, factor=1.5) -> dict:
    """Outlier detection using IQR and Z-score methods
    added for poc, not used, using sklearn instead
    """
    outlier_indices = {}
    for col in df.select_dtypes(include=[np.number]).columns:
        Q1 = df[col].quantile(0.25)
        Q3 = df[col].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - factor * IQR
        upper_bound = Q3 + factor * IQR
        outliers = df[(df[col] < lower_bound) | (df[col] > upper_bound)].index.tolist()
        if outliers:
            outlier_indices[col] = outliers
    return outlier_indices


def detect_outliers_zscore(df, threshold=3) -> dict:
    """Outlier detection using IQR and Z-score methods
    added for poc, not used, using sklearn instead
    """
    outlier_indices = {}
    for col in df.select_dtypes(include=[np.number]).columns:
        z_scores = np.abs((df[col] - df[col].mean()) / df[col].std())
        outliers = df[z_scores > threshold].index.tolist()
        if outliers:
            outlier_indices[col] = outliers
    return outlier_indices


# # Example usage:
# iqr_outliers = detect_outliers_iqr(X)
# zscore_outliers = detect_outliers_zscore(X)
# iqr_df = pd.DataFrame.from_dict(iqr_outliers, orient='index').reset_index()
# zscore_df = pd.DataFrame.from_dict(zscore_outliers, orient='index').reset_index()
# print("IQR outliers DataFrame: Outliers are values outside [Q1 - 1.5 x IQR, Q3 + 1.5 x IQR]")
# print("Z-score outliers DataFrame: Outliers are values with Z-score > 3")

# --------------------------------------------------------------
# Preprocessing functions
# --------------------------------------------------------------


def detect_and_remove_outliers(X: pd.DataFrame, y: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    """Detect and remove outliers from dataset using IsolationForest and LocalOutlierFactor
    It seems IsolationsForest works better ?

    :param X: feature matrix
    :type X: pd.DataFrame
    :param y: target vector
    :type y: np.ndarray
    :return: cleaned feature matrix and target vector
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    # IsolationForest
    iso = IsolationForest(contamination=0.05, random_state=42)
    outlier_pred_iso = iso.fit_predict(X)
    iso_outliers = np.where(outlier_pred_iso == -1)[0]
    print(f"IsolationForest detected {len(iso_outliers)} outliers.")

    mask_inliers = outlier_pred_iso == 1
    X_clean = X[mask_inliers].copy()
    y_clean = y[mask_inliers].copy()
    return X_clean, y_clean


def detect_and_remove_outliers_alternate(X: pd.DataFrame, y: pd.Series) -> tuple[pd.DataFrame, pd.Series]:
    """Detect and remove outliers from dataset using IsolationForest and LocalOutlierFactor

    :param X: feature matrix
    :type X: pd.DataFrame
    :param y: target vector
    :type y: np.ndarray
    :return: cleaned feature matrix and target vector
    :rtype: tuple[np.ndarray, np.ndarray]
    """

    lof = LocalOutlierFactor(n_neighbors=20, contamination=0.05)
    outlier_pred_lof = lof.fit_predict(X)
    lof_outliers = np.where(outlier_pred_lof == -1)[0]
    print(f"LocalOutlierFactor detected {len(lof_outliers)} outliers.")
    mask_inliers = outlier_pred_lof == 1
    X_clean = X[mask_inliers].copy()
    y_clean = y[mask_inliers].copy()
    return X_clean, y_clean


# Example usage:
# X_clean, y_clean = detect_and_remove_outliers(X, y)
# print(f"Original dataset size: {X.shape[0]}, Cleaned dataset size: {X_clean.shape[0]}"


def scale_features(X: pd.DataFrame) -> pd.DataFrame:
    """will not scale if only maccs are present"""
    scaler = StandardScaler()
    if any(name in X.columns for name in MACCS_NAMES) and any(name in X.columns for name in DESCRIPTOR_NAMES):
        X_MACCS = pd.DataFrame(X[MACCS_NAMES])
        X_rdkit = pd.DataFrame(X.drop(columns=MACCS_NAMES))
        X_scaled = pd.DataFrame(scaler.fit_transform(X_rdkit), columns=X_rdkit.columns)
        X_scaled = X_scaled.reset_index(drop=True)
        X_MACCS = X_MACCS.reset_index(drop=True)
        X_scaled = pd.concat([X_scaled, X_MACCS], axis=1)
    elif any(name in X.columns for name in DESCRIPTOR_NAMES):
        X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)
        X_scaled = X_scaled.reset_index(drop=True)
    elif any(name in X.columns for name in MACCS_NAMES):
        # actually a simle else would do
        X_scaled = X.copy()
    return X_scaled


def decorrelate(X: pd.DataFrame, target_column: str, threshold=0.95) -> list:
    """Decorrelate dataframe by finding which columns shall \
    be removed to achieve correlation level below threshold

    :param X: input dataframe
    :type X: pandas.core.frame.DataFrame
    :param threshold: maximum correlation allowed in the frame X
    :type threshold: float
    :return: list of columns to prune in order to achieve desired decorrelation level
    """
    # Use only numeric columns for correlation
    X_num = X.select_dtypes(include=[np.number]).copy()
    if target_column in X_num.columns:
        X_num = X_num.drop(columns=[target_column])

    N = X_num.shape[1]
    corr = X_num.corr().abs()

    to_drop = set()
    columns = X_num.columns
    for i in range(N - 1):
        for j in range(i + 1, N):
            if corr.at[columns[i], columns[j]] > threshold and columns[i] not in to_drop:
                to_drop.add(columns[j])
    return list(to_drop)


def remove_variance_and_correlation(X: pd.DataFrame, target_column: str) -> pd.DataFrame:
    """two steps in one function, not clean, but since this is always used together...
    :param X: input dataframe, should preferably have scaled features
    :type X: pd.DataFrame
    :param target_column: name of the target column to exclude from decorrelation
    :type target_column: str
    :return: dataframe with removed zero variance columns and highly correlated columns
    :rtype: pd.DataFrame
    """
    # remove zero std columns (no variance)
    zero_std_cols = X.columns[X.std() == 0]
    X = X.drop(columns=zero_std_cols)
    print(f"Number of features: {X.shape[1]}, number of samples: {X.shape[0]}")

    # drop columns hihgly correlated to some others
    cols_to_drop = decorrelate(X, target_column, threshold=0.95)
    X_decorrelated = X.drop(columns=cols_to_drop)
    print(f"Number of features: {X_decorrelated.shape[1]}, number of samples: {X_decorrelated.shape[0]}")
    return X_decorrelated


class Preprocessor:
    """
    Class to handle all preprocessing steps. Automatically drops irrelevant columns (inkl. smiles, etc), scales features, removes outliers, and decorrelates features.


    :param to_drop: set of feature sets to drop, defaults to {"maccs", "rdkit"}
    :type to_drop: set[str], optional

    :param target_column: name of the target column, defaults to "t_half"
    :type target_column: str, optional

    :param remove_outliers: whether to remove outliers, defaults to True
    :type remove_outliers: bool, optional

    :return: preprocessed dataframe with features and target column
    :rtype: pd.DataFrame

    Attributes:
        X (pd.DataFrame): Preprocessed feature matrix
        y (pd.Series): Target vector
        y_log (pd.Series): Log-transformed target vector

    """

    def __init__(
        self,
        to_drop: set[str] = {"maccs", "rdkit", "None"},
        target_column: str = "t_half",
        remove_outliers: bool = True,
    ) -> None:
        self.to_drop = to_drop
        self.target_column = target_column
        self.remove_outliers = remove_outliers
        self.X = None
        self.y = None
        self.y_log = None
        self.smiles = None

    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        # Save SMILES before any column dropping so they can be aligned with the
        # manual selected feature removal added
        # final feature matrix for structural visualisation downstream.
        smiles = df["Canonical_smiles"].copy() if "Canonical_smiles" in df.columns else None

        # feature selection ---------------------------------
        df_clean = drop_irrelevant_columns(df, self.to_drop)

        # -------------------------------------------------
        X = df_clean.drop(columns=[self.target_column])
        y = df_clean[self.target_column]
        print(f"Dataset: {X.shape[1]} features, {X.shape[0]} samples.")

        nan_mask = ~(X.isna().any(axis=1) | y.isna())
        X = X[nan_mask].reset_index(drop=True)
        y = y[nan_mask].reset_index(drop=True)
        if smiles is not None:
            smiles = smiles[nan_mask].reset_index(drop=True)
        print(f"After dropping NaNs: {X.shape[1]} features, {X.shape[0]} samples.")

        if self.remove_outliers:
            X_clean, y_clean = detect_and_remove_outliers(X, y)
            if smiles is not None:
                # X_clean has non-contiguous index (subset of 0..N-1); mirror on smiles
                smiles = smiles.loc[X_clean.index].reset_index(drop=True)
            print(f"Outlier cleaned dataset: {X_clean.shape[1]} features, {X_clean.shape[0]} samples.")
        else:
            print("Outlier removal not applied.")
            X_clean, y_clean = X, y

        X_scaled = scale_features(X_clean)  # resets index to 0..M-1 internally
        X_final = remove_variance_and_correlation(X_scaled, self.target_column)
        print(f"Final preprocessed dataset: {X_final.shape[1]} features, {X_final.shape[0]} samples.")

        # Set attributes for direct access
        self.X = X_final.copy()
        self.y = y_clean.reset_index(drop=True)
        self.y_log = np.log10(self.y) if (self.y is not None) else None  # log10! (bug before was log-e)
        self.smiles = smiles  # index 0..M-1, aligned with self.X

        df_final = self.X.copy()
        df_final[self.target_column] = self.y
        return df_final

    def get_smiles_for_split(self, X_train: pd.DataFrame, X_test: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        """Get SMILES aligned with train/test splits.

        Parameters
        ----------
        X_train : pd.DataFrame
            Training feature matrix
        X_test : pd.DataFrame
            Test feature matrix

        Returns
        -------
        tuple[pd.Series, pd.Series]
            SMILES for training and test sets
        """
        return self.smiles.loc[X_train.index], self.smiles.loc[X_test.index]


@dataclass
class PreprocessedData:
    """class containing all relevant preprocessed data for a compartament ready for modeling"""

    name: str
    df: pd.DataFrame  # not really required; mainly for reference/debugging; remove if too much memory use
    remove_outliers: bool = False
    X: pd.DataFrame = None
    y_log: pd.Series = None


# --------------------------------------------------------------
#  Modelling functions
# --------------------------------------------------------------


def t_t_split(X: pd.DataFrame, y: pd.Series, test_size: float = 0.2, random_state: int = 42) -> Any:
    """Split data into train and test sets, ensuring no constant features in either set."""

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    # Find constant features in train and test
    const_train = [col for col in X_train.columns if X_train[col].nunique() == 1]
    const_test = [col for col in X_test.columns if X_test[col].nunique() == 1]
    const_any = set(const_train) | set(const_test)
    # Remove from both sets
    X_train = X_train.drop(columns=const_any)
    X_test = X_test.drop(columns=const_any)
    return X_train, X_test, y_train, y_test


def output_metrics(y_true: Any, y_pred: Any) -> None:
    """print the output metrics for regression tasks"""
    print(f"R2: {r2_score(y_true, y_pred):.3f}")
    print(f"MAE: {mean_absolute_error(y_true, y_pred):.3f}")
    print(f"MSE: {mean_squared_error(y_true, y_pred):.3f}")
    print(f"RMSE: {root_mean_squared_error(y_true, y_pred):.3f}")
    print(f"Explained Variance: {explained_variance_score(y_true, y_pred):.3f}")


def output_metrics_w_return(y_true: Any, y_pred: Any) -> dict:
    """second function so that I don't have to refactor existing code"""
    # print(f"R2: {r2_score(y_true, y_pred):.3f}")
    # print(f"MAE: {mean_absolute_error(y_true, y_pred):.3f}")
    # print(f"MSE: {mean_squared_error(y_true, y_pred):.3f}")
    # print(f"RMSE: {root_mean_squared_error(y_true, y_pred):.3f}")
    # print(f"Explained Variance: {explained_variance_score(y_true, y_pred):.3f}")

    return {
        "R2": r2_score(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "MSE": mean_squared_error(y_true, y_pred),
        "RMSE": root_mean_squared_error(y_true, y_pred),
        "Explained Variance": explained_variance_score(y_true, y_pred),
    }


def svr_grid_search(X_train: pd.DataFrame, y_train: pd.Series) -> tuple[SVR, dict]:
    """Perform grid search to find the best SVR model.

    Parameters
    ----------
    X_train : pd.DataFrame
        _description_
    y_train : pd.Series
        _description_

    Returns
    -------
    SVR
        _description_
    """

    param_grid = {
        "C": [0.1, 1, 10, 100],
        "epsilon": [0.01, 0.1, 0.2, 0.5],
        "kernel": ["rbf", "linear"],
        "gamma": ["scale", "auto"],
    }

    svr = SVR()
    grid_search = GridSearchCV(svr, param_grid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1)
    grid_search.fit(X_train, y_train)

    print("Best parameters:", grid_search.best_params_)
    best_svr_estimate = grid_search.best_estimator_
    return best_svr_estimate, grid_search.best_params_


def chemical_space_pca(
    X_train: pd.DataFrame, X_test: pd.DataFrame, output_dir: Path, compartment: str
) -> tuple[PCA, StandardScaler, np.ndarray, np.ndarray, float]:
    """
    Fit a 2-component PCA on X_train, project X_train and X_test, then save a scatter plot.
    The explained variance and simple coverage statistics are printed and logged.
    """
    # Scale first – PCA is sensitive to feature magnitudes
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    pca = PCA(n_components=2, random_state=42)
    train_pcs = pca.fit_transform(X_train_scaled)
    test_pcs = pca.transform(X_test_scaled)

    var_explained = pca.explained_variance_ratio_
    print(f"\nPCA explained variance: PC1={var_explained[0]:.3f}, PC2={var_explained[1]:.3f} (total={sum(var_explained):.3f})")

    # ------------------------------------------------------------------ #
    # Simple coverage check: what fraction of test points falls inside
    # the bounding box of the training set in PCA space?
    # ------------------------------------------------------------------ #
    pc1_min, pc1_max = train_pcs[:, 0].min(), train_pcs[:, 0].max()
    pc2_min, pc2_max = train_pcs[:, 1].min(), train_pcs[:, 1].max()

    inside = (
        (test_pcs[:, 0] >= pc1_min) & (test_pcs[:, 0] <= pc1_max) & (test_pcs[:, 1] >= pc2_min) & (test_pcs[:, 1] <= pc2_max)
    )
    coverage_pct = inside.mean() * 100
    print(f"Test set coverage (bounding box): {coverage_pct:.1f}% of test compounds fall inside training domain")

    # ------------------------------------------------------------------ #
    # Plot
    # ------------------------------------------------------------------ #
    plt.figure(figsize=(8, 6))
    plt.scatter(
        train_pcs[:, 0],
        train_pcs[:, 1],
        alpha=0.4,
        s=20,
        color="royalblue",
        label=f"Train (n={len(X_train)})",
    )
    plt.scatter(
        test_pcs[:, 0],
        test_pcs[:, 1],
        alpha=0.7,
        s=30,
        color="darkorange",
        marker="^",
        label=f"Test (n={len(X_test)})",
    )

    # Draw training bounding box
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

    plt.xlabel(f"PC1 ({var_explained[0] * 100:.1f}% var)", fontsize=11)
    plt.ylabel(f"PC2 ({var_explained[1] * 100:.1f}% var)", fontsize=11)
    plt.title(f"Chemical Space Coverage via PCA – {compartment}", fontsize=13)
    plt.legend(fontsize=9)
    plt.tight_layout()
    out_path = output_dir / f"chemical_space_pca_{compartment}.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"PCA chemical space plot saved to {out_path}")

    return pca, scaler, train_pcs, test_pcs, coverage_pct


def applicability_domain_leverage(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train: pd.Series,
    y_test: pd.Series,
    y_pred_train: np.ndarray,
    y_pred_test: np.ndarray,
    output_dir: Path,
    compartment: str,
) -> dict:
    """
    Leverage-based Applicability Domain assessment (Williams plot).

    Computes hat-matrix leverages h_i = x_i^T (X_train^T X_train)^{-1} x_i on
    StandardScaler-scaled features together with standardised residuals.

    AD boundary:
        h* = 3p / n   (p = features, n = training samples)
        |std. residual| <= 3

    Compounds outside either boundary are flagged as outside the AD.

    Parameters
    ----------
    X_train, X_test : feature matrices (unscaled DataFrames)
    y_train, y_test : true targets in log10 space
    y_pred_train, y_pred_test : model predictions in log10 space
    output_dir : directory for the Williams plot image
    compartment : label used in plot title and filename

    Returns
    -------
    dict with keys: h_star, h_train, h_test, std_res_train, std_res_test,
                    inside_train, inside_test, pct_train_inside, pct_test_inside
    """
    # Scale features (leverage is sensitive to scale)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    n, p = X_train_scaled.shape
    h_star = 3.0 * p / n

    # Hat-matrix diagonal: h_i = x_i^T (X^T X)^{-1} x_i
    # Vectorised: h = rowsum( (X @ pinv(X^T X)) * X )
    XtX_inv = np.linalg.pinv(X_train_scaled.T @ X_train_scaled)
    h_train = (X_train_scaled @ XtX_inv * X_train_scaled).sum(axis=1)
    h_test = (X_test_scaled @ XtX_inv * X_test_scaled).sum(axis=1)

    # Standardised residuals in log10 space
    res_train = np.asarray(y_train) - y_pred_train
    res_test = np.asarray(y_test) - y_pred_test
    std_res = res_train.std(ddof=1)
    std_res_train = res_train / std_res
    std_res_test = res_test / std_res

    # AD membership
    inside_train = (h_train <= h_star) & (np.abs(std_res_train) <= 3)
    inside_test = (h_test <= h_star) & (np.abs(std_res_test) <= 3)
    pct_train_inside = inside_train.mean() * 100
    pct_test_inside = inside_test.mean() * 100

    print(f"\nApplicability Domain (Williams plot) — {compartment}")
    print(f"  h* threshold : {h_star:.4f}  (3 × {p} features / {n} training samples)")
    print(f"  Train inside AD : {pct_train_inside:.1f}%  ({inside_train.sum()}/{n})")
    print(f"  Test  inside AD : {pct_test_inside:.1f}%  ({inside_test.sum()}/{len(X_test)})")

    # Williams plot
    plt.figure(figsize=(8, 6))
    plt.scatter(
        h_train,
        std_res_train,
        alpha=0.5,
        s=20,
        color="royalblue",
        label=f"Train (n={n})",
    )
    plt.scatter(
        h_test,
        std_res_test,
        alpha=0.7,
        s=30,
        color="darkorange",
        marker="^",
        label=f"Test (n={len(X_test)})",
    )
    plt.axvline(h_star, color="gray", linestyle="--", linewidth=1.5, label=f"h* = {h_star:.3f}")
    plt.axhline(3, color="red", linestyle="--", linewidth=1.0, label="±3 std residual")
    plt.axhline(-3, color="red", linestyle="--", linewidth=1.0)
    plt.xlabel("Leverage (h)", fontsize=11)
    plt.ylabel("Standardised Residual", fontsize=11)
    plt.title(f"Williams Plot — Applicability Domain ({compartment})", fontsize=13)
    plt.legend(fontsize=9)
    plt.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"williams_plot_{compartment}.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Williams plot saved to {out_path}")

    return {
        "h_star": h_star,
        "h_train": h_train,
        "h_test": h_test,
        "std_res_train": std_res_train,
        "std_res_test": std_res_test,
        "inside_train": inside_train,
        "inside_test": inside_test,
        "pct_train_inside": pct_train_inside,
        "pct_test_inside": pct_test_inside,
    }


def chemical_space_morgan_pca(
    smiles_train: pd.Series,
    smiles_test: pd.Series,
    ad_results: dict,
    output_dir: Path,
    compartment: str,
) -> dict:
    """
    PCA on Morgan fingerprints coloured by AD membership + Butina clustering of
    AD-outside test compounds.

    Intent: reveal whether AD-outside test compounds are structurally novel (form
    their own singletons/clusters in SMILES space) or are descriptor-space outliers
    for other reasons.

    Parameters
    ----------
    smiles_train, smiles_test : SMILES series aligned with X_train / X_test by index
    ad_results : dict returned by applicability_domain_leverage()
    output_dir : directory for output images
    compartment : label for plot titles and filenames

    Returns
    -------
    dict with Butina cluster statistics for AD-outside test compounds
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    from rdkit.ML.Cluster import Butina

    def _smiles_to_fp(smi: str):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)

    def _fp_to_array(fp) -> np.ndarray:
        arr = np.zeros((1024,), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    # Build fingerprints, track valid indices (drop invalid SMILES)
    fps_train_raw = [_smiles_to_fp(s) for s in smiles_train]
    fps_test_raw = [_smiles_to_fp(s) for s in smiles_test]

    valid_tr = [i for i, fp in enumerate(fps_train_raw) if fp is not None]
    valid_te = [i for i, fp in enumerate(fps_test_raw) if fp is not None]
    fps_train = [fps_train_raw[i] for i in valid_tr]
    fps_test = [fps_test_raw[i] for i in valid_te]

    inside_train = ad_results["inside_train"][valid_tr]
    inside_test = ad_results["inside_test"][valid_te]

    mat_train = np.array([_fp_to_array(fp) for fp in fps_train])
    mat_test = np.array([_fp_to_array(fp) for fp in fps_test])

    # PCA fitted on training fingerprints only
    pca = PCA(n_components=2, random_state=42)
    train_pcs = pca.fit_transform(mat_train)
    test_pcs = pca.transform(mat_test)
    var_explained = pca.explained_variance_ratio_

    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 4-colour plot: train/test × inside/outside AD
    # ------------------------------------------------------------------ #
    plt.figure(figsize=(9, 7))
    plt.scatter(
        train_pcs[inside_train, 0],
        train_pcs[inside_train, 1],
        alpha=0.3,
        s=18,
        color="royalblue",
        label=f"Train in AD (n={inside_train.sum()})",
    )
    plt.scatter(
        train_pcs[~inside_train, 0],
        train_pcs[~inside_train, 1],
        alpha=0.7,
        s=35,
        color="steelblue",
        marker="s",
        label=f"Train out AD (n={(~inside_train).sum()})",
    )
    plt.scatter(
        test_pcs[inside_test, 0],
        test_pcs[inside_test, 1],
        alpha=0.5,
        s=30,
        color="darkorange",
        marker="^",
        label=f"Test in AD (n={inside_test.sum()})",
    )
    plt.scatter(
        test_pcs[~inside_test, 0],
        test_pcs[~inside_test, 1],
        alpha=0.9,
        s=55,
        color="crimson",
        marker="^",
        label=f"Test out AD (n={(~inside_test).sum()})",
    )
    plt.xlabel(f"PC1 ({var_explained[0] * 100:.1f}% var)", fontsize=11)
    plt.ylabel(f"PC2 ({var_explained[1] * 100:.1f}% var)", fontsize=11)
    plt.title(f"Chemical Space — Morgan FP PCA — AD membership ({compartment})", fontsize=13)
    plt.legend(fontsize=9)
    plt.tight_layout()
    pca_path = output_dir / f"morgan_pca_ad_{compartment}.png"
    plt.savefig(pca_path, dpi=150)
    plt.close()
    print(f"Morgan FP PCA (AD coloured) saved to {pca_path}")

    # ------------------------------------------------------------------ #
    # Butina clustering of AD-outside test compounds
    # ------------------------------------------------------------------ #
    fps_outside = [fps_test[i] for i in range(len(fps_test)) if not inside_test[i]]
    n_outside = len(fps_outside)
    cluster_stats = {"n_outside_ad": n_outside}

    if n_outside >= 2:
        dists = []
        for i in range(1, n_outside):
            sims = DataStructs.BulkTanimotoSimilarity(fps_outside[i], fps_outside[:i])
            dists.extend([1 - x for x in sims])
        # cutoff=0.5 → Tanimoto similarity ≥ 0.5 — captures "same chemical series" band
        # for Morgan2 (bits); anything below 0.5 is weakly related or noise on this FP scale.
        clusters = Butina.ClusterData(dists, n_outside, 0.5, isDistData=True)

        n_singletons = sum(1 for c in clusters if len(c) == 1)
        cluster_sizes = sorted([len(c) for c in clusters], reverse=True)
        singleton_pct = n_singletons / n_outside * 100

        print(f"\nButina clustering (Tanimoto distance cutoff=0.5, similarity ≥ 0.5) of {n_outside} AD-outside test compounds:")
        print(f"  Clusters  : {len(clusters)}")
        print(f"  Singletons: {n_singletons} ({singleton_pct:.1f}%) — structurally unique outliers")
        print(f"  Cluster sizes (largest first): {cluster_sizes[:10]}")

        cluster_stats.update(
            {
                "n_clusters": len(clusters),
                "n_singletons": n_singletons,
                "singleton_pct": singleton_pct,
                "cluster_sizes": cluster_sizes,
            }
        )

        # ------------------------------------------------------------------ #
        # Plot 1: cluster size bar chart
        # ------------------------------------------------------------------ #
        from matplotlib.patches import Patch

        fig, ax = plt.subplots(figsize=(max(6, len(clusters) * 0.4 + 2), 4))
        colors_bar = ["#d62728" if s == 1 else "#1f77b4" for s in cluster_sizes]
        ax.bar(
            range(len(cluster_sizes)),
            cluster_sizes,
            color=colors_bar,
            edgecolor="white",
            linewidth=0.5,
        )
        ax.set_xlabel("Cluster index (sorted by size)", fontsize=10)
        ax.set_ylabel("Cluster size (# compounds)", fontsize=10)
        ax.set_title(
            f"Butina cluster sizes — AD-outside test compounds ({compartment})\n"
            f"Tanimoto distance cutoff=0.5 (sim ≥ 0.5) | {len(clusters)} clusters | {n_singletons} singletons ({singleton_pct:.1f}%)",
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
        print(f"Butina cluster size bar chart saved to {bar_path}")

        # ------------------------------------------------------------------ #
        # Plot 2: Morgan PCA coloured by Butina cluster assignment
        # ------------------------------------------------------------------ #
        outside_indices = [i for i in range(len(fps_test)) if not inside_test[i]]
        cluster_label = np.full(n_outside, -1, dtype=int)
        for c_idx, c in enumerate(clusters):
            for member_idx in c:
                cluster_label[member_idx] = c_idx

        pcs_outside = test_pcs[np.array(outside_indices)]
        cmap = plt.get_cmap("tab20")
        n_multi = len(clusters) - n_singletons

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(
            train_pcs[:, 0],
            train_pcs[:, 1],
            s=12,
            alpha=0.15,
            color="grey",
            label="Train",
        )
        ax.scatter(
            test_pcs[inside_test, 0],
            test_pcs[inside_test, 1],
            s=18,
            alpha=0.25,
            color="darkorange",
            label="Test in AD",
        )
        # singletons — Butina returns clusters sorted largest first, so singletons are at the tail
        singleton_mask = np.array([cluster_label[i] >= n_multi for i in range(n_outside)])
        ax.scatter(
            pcs_outside[singleton_mask, 0],
            pcs_outside[singleton_mask, 1],
            s=70,
            alpha=0.9,
            color="crimson",
            marker="x",
            linewidths=1.5,
            label=f"Singleton ({n_singletons})",
        )
        for c_idx in range(n_multi):
            mask = cluster_label == c_idx
            if mask.sum() == 0:
                continue
            ax.scatter(
                pcs_outside[mask, 0],
                pcs_outside[mask, 1],
                s=70,
                alpha=0.85,
                color=cmap(c_idx % 20),
                marker="^",
                label=f"Cluster {c_idx + 1} (n={mask.sum()})",
            )
        ax.set_xlabel(f"PC1 ({var_explained[0] * 100:.1f}% var)", fontsize=10)
        ax.set_ylabel(f"PC2 ({var_explained[1] * 100:.1f}% var)", fontsize=10)
        ax.set_title(f"Butina Clusters — AD-outside Test Compounds ({compartment})", fontsize=12)
        ax.legend(fontsize=8, markerscale=1.2)
        plt.tight_layout()
        cluster_pca_path = output_dir / f"butina_cluster_pca_{compartment}.png"
        plt.savefig(cluster_pca_path, dpi=150)
        plt.close()
        print(f"Butina cluster PCA plot saved to {cluster_pca_path}")

    else:
        print(f"Too few AD-outside test compounds for Butina clustering (n={n_outside})")

    return cluster_stats
