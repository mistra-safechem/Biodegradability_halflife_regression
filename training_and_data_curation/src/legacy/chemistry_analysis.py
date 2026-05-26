"""
Chemical structure analysis for QSAR model interpretation.

LEGACY: module copied to src/legacy for backwardscompatibility, see also the readme.md in this folder.

"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from src.legacy.log_utils import log_section, log_to_file


def get_murcko_scaffold(smiles: str) -> str:
    """Extract Murcko scaffold from SMILES.

    Parameters
    ----------
    smiles : str
        Input molecule SMILES

    Returns
    -------
    str
        Murcko scaffold SMILES, or empty string if invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)


def get_brics_fragments(smiles: str) -> List[str]:
    """Fragment molecule using BRICS algorithm.

    Parameters
    ----------
    smiles : str
        Input molecule SMILES

    Returns
    -------
    List[str]
        List of BRICS fragment SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    fragments = BRICS.BreakBRICSBonds(mol)
    frag_smiles = [Chem.MolToSmiles(frag) for frag in Chem.GetMolFrags(fragments, asMols=True)]
    return frag_smiles


def get_rgroups(scaffold_smi: str, molecule_smi: str) -> Dict[str, str]:
    """Identify R-group substituents on scaffold.

    Parameters
    ----------
    scaffold_smi : str
        Murcko scaffold SMILES
    molecule_smi : str
        Full molecule SMILES

    Returns
    -------
    Dict[str, str]
        Dictionary mapping attachment points to substituent SMILES
    """
    scaffold_mol = Chem.MolFromSmiles(scaffold_smi)
    molecule_mol = Chem.MolFromSmiles(molecule_smi)

    if scaffold_mol is None or molecule_mol is None:
        return {}

    # Find attachment points on scaffold (atoms not in scaffold but in molecule)
    scaffold_atoms = set(scaffold_mol.GetAtoms())
    rgroups = {}

    # Simple approach: identify substituents by difference
    # More sophisticated: use RDKit's R-group decomposition
    try:
        from rdkit.Chem import rdRGroupDecomposition

        rgroups_result = rdRGroupDecomposition.RGroupDecompose([scaffold_mol], [molecule_mol], asRows=True)
        if rgroups_result and len(rgroups_result) > 0:
            for i, (label, rgroup_mol) in enumerate(rgroups_result[0].items()):
                if rgroup_mol and not rgroup_mol.IsNone():
                    rgroups[label] = Chem.MolToSmiles(rgroup_mol)
    except Exception:
        pass

    return rgroups


def analyse_scaffolds(
    smiles: pd.Series,
    y_true: pd.Series,
    y_pred: np.ndarray,
    output_dir: Path,
    top_n: int = 10,
) -> Dict[str, Any]:
    """Section 1: Core (scaffold) analysis.

    Murcko scaffold decomposition with performance metrics per scaffold.

    Parameters
    ----------
    smiles : pd.Series
        SMILES strings aligned with y_true and y_pred
    y_true : pd.Series
        True target values (days)
    y_pred : np.ndarray
        Predicted target values (days)
    output_dir : Path
        Directory for output files
    top_n : int
        Number of top scaffolds to report

    Returns
    -------
    Dict[str, Any]
        Scaffold analysis results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Compute scaffolds
    scaffold_map = smiles.apply(get_murcko_scaffold)

    # Scaffold distribution
    scaffold_counts = scaffold_map.value_counts()
    total = scaffold_counts.sum()
    scaffold_dist = pd.DataFrame(
        {
            "scaffold_smiles": scaffold_counts.index,
            "count": scaffold_counts.values,
            "coverage_pct": (scaffold_counts.values / total * 100).round(2),
        }
    )

    # Cumulative coverage
    scaffold_dist["cumulative_pct"] = scaffold_dist["coverage_pct"].cumsum().round(2)

    # Save scaffold distribution
    scaffold_dist.to_csv(output_dir / "scaffold_distribution.csv", index=False)

    # Performance per scaffold
    df = pd.DataFrame(
        {
            "smiles": smiles.values,
            "scaffold": scaffold_map.values,
            "y_true": y_true.values,
            "y_pred": y_pred,
        }
    )

    scaffold_perf = []
    for scaffold in scaffold_counts.index[:top_n]:
        subset = df[df["scaffold"] == scaffold]
        if len(subset) < 2:
            continue

        y_t = subset["y_true"].values
        y_p = subset["y_pred"].values

        rmse = np.sqrt(np.mean((y_t - y_p) ** 2))
        ss_res = np.sum((y_t - y_p) ** 2)
        ss_tot = np.sum((y_t - np.mean(y_t)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        mean_err = np.mean(y_p - y_t)

        scaffold_perf.append(
            {
                "scaffold_smiles": scaffold,
                "n_samples": len(subset),
                "coverage_pct": scaffold_dist[scaffold_dist["scaffold_smiles"] == scaffold]["coverage_pct"].values[0],
                "rmse": round(rmse, 3),
                "r2": round(r2, 3),
                "mean_error": round(mean_err, 3),
            }
        )

    scaffold_perf_df = pd.DataFrame(scaffold_perf)
    scaffold_perf_df.to_csv(output_dir / "scaffold_performance.csv", index=False)

    # Flag single-scaffold dominance
    top_scaffold_cov = scaffold_dist["coverage_pct"].iloc[0] if len(scaffold_dist) > 0 else 0
    dominance_flag = top_scaffold_cov > 50

    # Plot top scaffolds
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        top_scaffolds = scaffold_dist.head(top_n)
        ax.barh(range(len(top_scaffolds)), top_scaffolds["count"], color="steelblue")
        ax.set_yticks(range(len(top_scaffolds)))
        ax.set_yticklabels([f"#{i + 1}" for i in range(len(top_scaffolds))], fontsize=9)
        ax.set_xlabel("Count", fontsize=11)
        ax.set_title(f"Top {top_n} Murcko Scaffolds", fontsize=13)
        ax.invert_yaxis()

        # Add coverage labels
        for i, (count, cov) in enumerate(zip(top_scaffolds["count"], top_scaffolds["coverage_pct"])):
            ax.text(count + 0.5, i, f"{cov:.1f}%", va="center", fontsize=8)

        plt.tight_layout()
        plt.savefig(output_dir / "scaffold_bar_plot.png", dpi=150)
        plt.close()
    except Exception:
        pass

    # Summary statistics
    n_scaffolds = len(scaffold_counts)
    n_singletons = (scaffold_counts == 1).sum()
    cum_90_idx = (scaffold_dist["cumulative_pct"] >= 90).idxmax()
    scaffolds_for_90 = cum_90_idx + 1 if cum_90_idx >= 0 else top_n

    return {
        "n_scaffolds_total": n_scaffolds,
        "n_singletons": n_singletons,
        "top_scaffold_coverage": top_scaffold_cov,
        "dominance_flag": dominance_flag,
        "scaffolds_for_90_pct": scaffolds_for_90,
        "scaffold_distribution": scaffold_dist,
        "scaffold_performance": scaffold_perf_df,
    }


def analyse_substituents(
    smiles: pd.Series,
    y_true: pd.Series,
    y_pred: np.ndarray,
    scaffold_results: Dict[str, Any],
    output_dir: Path,
) -> Dict[str, Any]:
    """Section 2: Substituent analysis.

    R-group analysis on dominant cores.

    Parameters
    ----------
    smiles : pd.Series
        SMILES strings
    y_true : pd.Series
        True target values
    y_pred : np.ndarray
        Predicted target values
    scaffold_results : Dict[str, Any]
        Results from analyse_scaffolds()
    output_dir : Path
        Directory for output files

    Returns
    -------
    Dict[str, Any]
        Substituent analysis results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get dominant scaffold
    if scaffold_results["scaffold_distribution"].empty:
        return {"error": "No scaffolds found"}

    dominant_scaffold = scaffold_results["scaffold_distribution"].iloc[0]["scaffold_smiles"]

    # Filter molecules with dominant scaffold
    scaffold_map = smiles.apply(get_murcko_scaffold)
    dominant_mask = scaffold_map == dominant_scaffold
    dominant_smiles = smiles[dominant_mask]
    dominant_y = y_true[dominant_mask]
    dominant_y_pred = y_pred[dominant_mask]

    if len(dominant_smiles) < 5:
        return {"error": "Too few molecules with dominant scaffold"}

    # Extract R-groups
    rgroup_data = []
    for smi, y_t, y_p in zip(dominant_smiles.values, dominant_y.values, dominant_y_pred):
        rgroups = get_rgroups(dominant_scaffold, smi)
        for label, rgroup_smi in rgroups.items():
            rgroup_data.append(
                {
                    "r_group_label": label,
                    "r_group_smiles": rgroup_smi,
                    "y_true": y_t,
                    "y_pred": y_p,
                    "error": y_p - y_t,
                }
            )

    if not rgroup_data:
        return {"error": "No R-groups identified"}

    rgroup_df = pd.DataFrame(rgroup_data)

    # Aggregate by R-group
    rgroup_agg = (
        rgroup_df.groupby(["r_group_label", "r_group_smiles"])
        .agg({"y_true": ["count", "mean", "std"], "error": ["mean", "std"]})
        .round(3)
    )
    rgroup_agg.columns = ["_".join(col).strip("_") for col in rgroup_agg.columns]
    rgroup_agg = rgroup_agg.reset_index()

    rgroup_agg.to_csv(output_dir / "rgroup_analysis.csv", index=False)

    # Activity vs substituent class trends
    # Group by R-group label and compute mean activity
    rgroup_trends = rgroup_df.groupby("r_group_label").agg({"y_true": "mean", "error": "mean"}).round(3)
    rgroup_trends.to_csv(output_dir / "substituent_trends.csv")

    # Plot if matplotlib available
    try:
        import matplotlib.pyplot as plt

        if len(rgroup_agg) > 0:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Left: R-group frequency
            rgroup_counts = rgroup_agg.groupby("r_group_label").size()
            axes[0].bar(range(len(rgroup_counts)), rgroup_counts.values, color="steelblue")
            axes[0].set_xlabel("R-group label")
            axes[0].set_ylabel("Count")
            axes[0].set_title("R-group Frequency")
            axes[0].set_xticks(range(len(rgroup_counts)))
            axes[0].set_xticklabels(rgroup_counts.index, rotation=45, fontsize=8)

            # Right: Mean error by R-group
            rgroup_err = rgroup_agg.groupby("r_group_label")["error_mean"].mean()
            colors = ["red" if e > 0 else "green" if e < 0 else "gray" for e in rgroup_err.values]
            axes[1].bar(range(len(rgroup_err)), rgroup_err.values, color=colors)
            axes[1].axhline(0, color="black", linestyle="--", linewidth=0.8)
            axes[1].set_xlabel("R-group label")
            axes[1].set_ylabel("Mean Error")
            axes[1].set_title("Prediction Error by R-group")
            axes[1].set_xticks(range(len(rgroup_err)))
            axes[1].set_xticklabels(rgroup_err.index, rotation=45, fontsize=8)

            plt.tight_layout()
            plt.savefig(output_dir / "substituent_trends.png", dpi=150)
            plt.close()
    except Exception:
        pass

    return {
        "dominant_scaffold": dominant_scaffold,
        "n_molecules_analyzed": len(dominant_smiles),
        "n_rgroups_identified": len(rgroup_df),
        "unique_rgroups": rgroup_df["r_group_smiles"].nunique(),
        "rgroup_data": rgroup_df,
        "rgroup_aggregated": rgroup_agg,
    }


def analyse_fragments(
    smiles: pd.Series,
    y_true: pd.Series,
    y_pred: np.ndarray,
    feature_importances: np.ndarray,
    feature_names: List[str],
    output_dir: Path,
) -> Dict[str, Any]:
    """Section 3: Fragment / frequency analysis.

    BRICS fragmentation with frequency vs activity analysis.

    Parameters
    ----------
    smiles : pd.Series
        SMILES strings
    y_true : pd.Series
        True target values
    y_pred : np.ndarray
        Predicted target values
    feature_importances : np.ndarray
        Feature importance scores (e.g., from permutation)
    feature_names : List[str]
        Feature names
    output_dir : Path
        Directory for output files

    Returns
    -------
    Dict[str, Any]
        Fragment analysis results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Fragment all molecules
    all_fragments = []
    for smi, y_t in zip(smiles.values, y_true.values):
        frags = get_brics_fragments(smi)
        for frag in frags:
            all_fragments.append({"fragment_smiles": frag, "y_true": y_t})

    if not all_fragments:
        return {"error": "No fragments generated"}

    frag_df = pd.DataFrame(all_fragments)

    # Fragment frequency
    frag_counts = frag_df.groupby("fragment_smiles").agg({"y_true": ["count", "mean", "std"]})
    frag_counts.columns = ["frequency", "mean_activity", "activity_std"]
    frag_counts = frag_counts.round(3).reset_index()
    frag_counts = frag_counts.sort_values("frequency", ascending=False)

    frag_counts.to_csv(output_dir / "fragment_frequency.csv", index=False)

    # Frequency vs mean activity plot
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 6))

        # Plot common fragments (frequency > 1)
        common = frag_counts[frag_counts["frequency"] > 1]
        if len(common) > 0:
            ax.scatter(
                common["frequency"],
                common["mean_activity"],
                alpha=0.6,
                s=50,
                color="steelblue",
            )

            # Label top 10 by frequency
            for _, row in common.head(10).iterrows():
                ax.annotate(
                    row["fragment_smiles"][:20] + "...",
                    (row["frequency"], row["mean_activity"]),
                    fontsize=6,
                    ha="left",
                )

        ax.set_xlabel("Fragment Frequency")
        ax.set_ylabel("Mean Activity (days)")
        ax.set_title("Fragment Frequency vs Mean Activity")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_dir / "fragment_activity_plot.png", dpi=150)
        plt.close()
    except Exception:
        pass

    # Fragment enrichment in high-activity compounds
    high_activity_threshold = y_true.quantile(0.75)
    high_activity_mask = y_true >= high_activity_threshold

    # Fragments from high-activity compounds
    high_act_fragments = set()
    for smi in smiles[high_activity_mask].values:
        frags = get_brics_fragments(smi)
        high_act_fragments.update(frags)

    # Fragments from outliers (high error)
    errors = np.abs(y_pred - y_true.values)
    outlier_threshold = np.percentile(errors, 90)
    outlier_mask = errors >= outlier_threshold

    outlier_fragments = set()
    for smi in smiles[outlier_mask].values:
        frags = get_brics_fragments(smi)
        outlier_fragments.update(frags)

    # Save outlier-enriched fragments
    outlier_enriched = list(outlier_fragments - high_act_fragments)
    outlier_df = pd.DataFrame(
        {
            "fragment_smiles": outlier_enriched,
            "in_high_activity": [f in high_act_fragments for f in outlier_enriched],
            "in_outliers": [True] * len(outlier_enriched),
        }
    )
    outlier_df.to_csv(output_dir / "outlier_fragments.csv", index=False)

    return {
        "n_fragments_total": len(frag_df),
        "n_unique_fragments": frag_df["fragment_smiles"].nunique(),
        "fragment_frequency": frag_counts,
        "n_high_activity_fragments": len(high_act_fragments),
        "n_outlier_fragments": len(outlier_fragments),
        "outlier_fragments": outlier_df,
    }


def analyse_chemistry_consistency(
    smiles: pd.Series,
    X_test: pd.DataFrame,
    y_true: pd.Series,
    y_pred: np.ndarray,
    residuals: np.ndarray,
    feature_importances: np.ndarray,
    feature_names: List[str],
    ad_inside: np.ndarray,
    output_dir: Path,
) -> Dict[str, Any]:
    """Section 4: Model–chemistry consistency check.

    Descriptor sanity check, prediction stability, AD violations analysis.

    Parameters
    ----------
    smiles : pd.Series
        SMILES strings
    X_test : pd.DataFrame
        Test features
    y_true : pd.Series
        True target values
    y_pred : np.ndarray
        Predicted target values
    residuals : np.ndarray
        Prediction residuals (y_true - y_pred)
    feature_importances : np.ndarray
        Feature importance scores
    feature_names : List[str]
        Feature names
    ad_inside : np.ndarray
        Boolean array indicating compounds inside AD
    output_dir : Path
        Directory for output files

    Returns
    -------
    Dict[str, Any]
        Chemistry consistency analysis results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Map important descriptors to chemistry
    descriptor_chemistry = {}
    top_20_idx = np.argsort(feature_importances)[-20:]

    for idx in top_20_idx:
        if idx >= len(feature_names):
            continue
        fname = feature_names[idx]
        importance = feature_importances[idx]

        # Interpret descriptor based on name
        interpretation = _interpret_descriptor(fname)

        descriptor_chemistry[fname] = {
            "importance": float(importance),
            "interpretation": interpretation,
        }

    # Save descriptor map
    with open(output_dir / "descriptor_chemistry_map.json", "w") as f:
        json.dump(descriptor_chemistry, f, indent=2)

    # Prediction stability per motif (scaffold)
    scaffold_map = smiles.apply(get_murcko_scaffold)

    motif_stability = []
    for scaffold in scaffold_map.unique():
        mask = scaffold_map == scaffold
        if mask.sum() < 3:
            continue

        scaffold_pred = y_pred[mask]
        scaffold_true = y_true[mask].values
        scaffold_res = residuals[mask]

        pred_var = np.var(scaffold_pred)
        pred_std = np.std(scaffold_pred)
        mean_abs_err = np.mean(np.abs(scaffold_res))

        # Stability flag: low variance + low error = stable
        stable = pred_std < np.std(y_pred) * 0.5 and mean_abs_err < np.mean(np.abs(residuals))

        motif_stability.append(
            {
                "scaffold_smiles": scaffold,
                "n_samples": mask.sum(),
                "prediction_std": round(pred_std, 3),
                "prediction_variance": round(pred_var, 3),
                "mean_absolute_error": round(mean_abs_err, 3),
                "stable": stable,
            }
        )

    motif_stability_df = pd.DataFrame(motif_stability)
    motif_stability_df.to_csv(output_dir / "motif_stability.csv", index=False)

    # AD violations scaffold-driven analysis
    ad_outside_mask = ~ad_inside

    ad_scaffold_analysis = []
    for scaffold in scaffold_map.unique():
        mask = scaffold_map == scaffold
        n_total = mask.sum()
        n_outside = (mask & ad_outside_mask).sum()
        pct_outside = (n_outside / n_total * 100) if n_total > 0 else 0

        ad_scaffold_analysis.append(
            {
                "scaffold_smiles": scaffold,
                "n_samples": n_total,
                "n_outside_ad": n_outside,
                "pct_outside_ad": round(pct_outside, 1),
            }
        )

    ad_scaffold_df = pd.DataFrame(ad_scaffold_analysis)
    ad_scaffold_df = ad_scaffold_df.sort_values("pct_outside_ad", ascending=False)
    ad_scaffold_df.to_csv(output_dir / "ad_scaffold_analysis.csv", index=False)

    # Summary: are AD violations scaffold-driven?
    # Check if certain scaffolds have disproportionately high AD-outside rates
    high_ad_violation = ad_scaffold_df[ad_scaffold_df["pct_outside_ad"] > 50]
    scaffold_driven = len(high_ad_violation) > 0

    return {
        "descriptor_map": descriptor_chemistry,
        "motif_stability": motif_stability_df,
        "ad_scaffold_analysis": ad_scaffold_df,
        "scaffold_driven_ad_violations": scaffold_driven,
        "n_scaffolds_with_high_ad_violations": len(high_ad_violation),
    }


def _interpret_descriptor(descriptor_name: str) -> str:
    """Interpret RDKit descriptor name into chemistry meaning.

    Parameters
    ----------
    descriptor_name : str
        RDKit descriptor name

    Returns
    -------
    str
        Human-readable interpretation
    """
    interpretations = {
        "MolLogP": "Lipophilicity (partition coefficient)",
        "MolWt": "Molecular weight",
        "NumRotatableBonds": "Molecular flexibility",
        "NumHDonors": "Hydrogen bond donor count",
        "NumHAcceptors": "Hydrogen bond acceptor count",
        "TPSA": "Polar surface area",
        "NumAromaticRings": "Aromaticity",
        "NumAliphaticRings": "Aliphatic ring count",
        "RingCount": "Total ring count",
        "FractionCSP3": "3D character (sp3 fraction)",
        "NumHeteroatoms": "Heteroatom count",
        "NumSaturatedRings": "Saturated ring count",
        "HeavyAtomCount": "Heavy atom count",
        "NumRadicalElectrons": "Radical character",
        "FormalCharge": "Net molecular charge",
        "MaxPartialCharge": "Maximum partial charge",
        "MinPartialCharge": "Minimum partial charge",
        "MaxAbsPartialCharge": "Maximum absolute partial charge",
        "AvgIpc": "Information content index",
        "BalabanJ": "Balaban connectivity index",
        "BertzCT": "Bertz complexity index",
    }

    # Check exact match
    if descriptor_name in interpretations:
        return interpretations[descriptor_name]

    # Check partial match
    for key, value in interpretations.items():
        if key.lower() in descriptor_name.lower():
            return value
        # could log if descriptor name isn't included in feature set

    # MACCS keys
    if descriptor_name.startswith("MACCS"):
        return "Structural fingerprint bit"

    # Default
    return "Molecular descriptor"


def analyse_chemistry_all(compartment: str, model_artifacts: Dict[str, Any], logs_dir: Path, log_file: str) -> Dict[str, Any]:
    """Run all chemistry analyses (sections 1-4).

    Parameters
    ----------
    compartment : str
        Compartment name (air, water, soil, sediment)
    model_artifacts : Dict[str, Any]
        Model artifacts from training
    logs_dir : Path
        Base logs directory
    log_file : str
        Log file path

    Returns
    -------
    Dict[str, Any]
        Combined chemistry analysis results
    """

    # Setup output directory
    chemistry_dir = logs_dir / "chemistry_analysis"
    chemistry_dir.mkdir(parents=True, exist_ok=True)

    log_section("Chemical Structure Analysis (Sections 1-4)", log_file)

    # Extract data from model_artifacts
    smiles_train = model_artifacts.get("smiles_train", pd.Series())
    smiles_test = model_artifacts.get("smiles_test", pd.Series())
    y_train = model_artifacts.get("y_train", pd.Series())
    y_test = model_artifacts.get("y_test", pd.Series())
    y_pred = model_artifacts.get("y_pred_svr", np.array([]))
    y_pred_train = model_artifacts.get("y_pred_train", np.array([]))
    X_test = model_artifacts.get("X_test", pd.DataFrame())

    # Inverse transform predictions if log-transformed
    y_test_exp = model_artifacts.get("y_test_exp_svr", np.power(10, y_test))
    y_pred_exp = model_artifacts.get("y_pred_exp_svr", np.power(10, y_pred))

    # Feature importances
    importances = model_artifacts.get("importances", np.array([]))
    feature_names = list(X_test.columns) if X_test is not None else []

    # AD results
    ad_results = model_artifacts.get("ad_results", {})
    ad_inside_test = ad_results.get("inside_test", np.ones(len(y_test), dtype=bool))

    results = {}

    # Section 1: Scaffold analysis
    log_to_file("Section 1: Scaffold Analysis", log_file)
    scaffold_results = analyse_scaffolds(
        smiles=smiles_test,
        y_true=y_test_exp,
        y_pred=y_pred_exp,
        output_dir=chemistry_dir,
        top_n=10,
    )
    results["scaffold"] = scaffold_results

    if "error" not in scaffold_results:
        log_to_file(f"  Total scaffolds: {scaffold_results['n_scaffolds_total']}", log_file)
        log_to_file(f"  Singletons: {scaffold_results['n_singletons']}", log_file)
        log_to_file(
            f"  Top scaffold coverage: {scaffold_results['top_scaffold_coverage']:.1f}%",
            log_file,
        )
        log_to_file(f"  Dominance flag: {scaffold_results['dominance_flag']}", log_file)
        log_to_file(
            f"  Scaffolds for 90% coverage: {scaffold_results['scaffolds_for_90_pct']}",
            log_file,
        )
    else:
        log_to_file(f"  Error: {scaffold_results['error']}", log_file)

    # Section 2: Substituent analysis
    log_to_file("Section 2: Substituent Analysis", log_file)
    substituent_results = analyse_substituents(
        smiles=smiles_test,
        y_true=y_test_exp,
        y_pred=y_pred_exp,
        scaffold_results=scaffold_results,
        output_dir=chemistry_dir,
    )
    results["substituent"] = substituent_results

    if "error" not in substituent_results:
        log_to_file(
            f"  Dominant scaffold: {substituent_results['dominant_scaffold'][:50]}...",
            log_file,
        )
        log_to_file(
            f"  Molecules analyzed: {substituent_results['n_molecules_analyzed']}",
            log_file,
        )
        log_to_file(f"  Unique R-groups: {substituent_results['unique_rgroups']}", log_file)
    else:
        log_to_file(f"  Error: {substituent_results['error']}", log_file)

    # Section 3: Fragment analysis
    log_to_file("Section 3: Fragment Analysis", log_file)
    fragment_results = analyse_fragments(
        smiles=smiles_test,
        y_true=y_test_exp,
        y_pred=y_pred_exp,
        feature_importances=importances if len(importances) > 0 else np.zeros(len(feature_names)),
        feature_names=feature_names,
        output_dir=chemistry_dir,
    )
    results["fragment"] = fragment_results

    if "error" not in fragment_results:
        log_to_file(f"  Total fragments: {fragment_results['n_fragments_total']}", log_file)
        log_to_file(f"  Unique fragments: {fragment_results['n_unique_fragments']}", log_file)
        log_to_file(
            f"  Outlier fragments: {len(fragment_results['outlier_fragments'])}",
            log_file,
        )
    else:
        log_to_file(f"  Error: {fragment_results['error']}", log_file)

    # Section 4: Chemistry consistency
    log_to_file("Section 4: Chemistry Consistency Check", log_file)
    residuals = y_test_exp - y_pred_exp
    consistency_results = analyse_chemistry_consistency(
        smiles=smiles_test,
        X_test=X_test,
        y_true=y_test_exp,
        y_pred=y_pred_exp,
        residuals=residuals,
        feature_importances=importances if len(importances) > 0 else np.zeros(len(feature_names)),
        feature_names=feature_names,
        ad_inside=ad_inside_test,
        output_dir=chemistry_dir,
    )
    results["consistency"] = consistency_results

    log_to_file(
        f"  Descriptor map entries: {len(consistency_results['descriptor_map'])}",
        log_file,
    )
    log_to_file(f"  Motifs analyzed: {len(consistency_results['motif_stability'])}", log_file)
    log_to_file(
        f"  Scaffold-driven AD violations: {consistency_results['scaffold_driven_ad_violations']}",
        log_file,
    )

    log_to_file(f"Chemistry analysis complete. Output: {chemistry_dir}", log_file)

    return results
