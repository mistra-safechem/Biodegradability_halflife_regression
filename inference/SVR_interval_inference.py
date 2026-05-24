"""Inference Script — predict biodegradability half-lifes

Given a SMILES string (or file list) and a trained model JSON card, this script:

1. Loads the SVR model (.joblib) and AD artefact (.npz) from the paths
   recorded in the model JSON card (as well as the training-split CSV).
2. Computes RDKit descriptors and MACCS fingerprints for the input SMILES.
3. Selects and orders features to match the training feature set exactly.
4. Scales features using the training-set statistics stored in the AD
   artefact (refits a StandardScaler on the saved X_train — identical to
   the scaler used during training).
5. Predicts the log10 half-life centroid.
6. Computes leverage h = xᵀ (XᵀX)⁻¹ x and compares to h* to decide AD.
7. Looks up the training-split CSV to check if the exact SMILES was in
   the training set.
8. Prints a human-readable prediction report.

Usage
-----
    uv run SVR_interval_inference.py --smiles "CCCC" --model models/SVR_air_hsbd_<ts>.json
    (or: Python .... with same args)

    # Multiple SMILES from a file (one per line):
    uv run SVR_interval_inference.py --smiles-file compounds.txt --model models/SVR_air_hsbd_<ts>.json

Notes
-----
- The script re-derives the StandardScaler from X_train stored in the AD
  artefact. This guarantees exact reproducibility without storing the
  scaler as a separate file (stems from initial point based SVR model,
  which is now legacy but the AD artefact structure was retained for consistency).
- Leverage-based AD reflects *structural* similarity to the training set.
  It (obviously) does NOT guarantee that the model's prediction is correct for a
  given compound.
- Training-set membership check is an exact SMILES string match. If the
  same compound was stored under a different SMILES representation it
  will not be detected. Consider canonical SMILES normalisation upstream.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys
from src.rdkit_tools import DESCRIPTOR_NAMES, MACCS_NAMES  # noqa: E402

SCRIPT_DIR = Path(__file__).parent


# ---------------------------------------------------------------------------
# Feature computation
# ---------------------------------------------------------------------------


def compute_features(smiles: str) -> pd.DataFrame | None:
    """Return a single-row DataFrame of RDKit descriptors + MACCS bits.

    Returns None if the SMILES cannot be parsed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    desc_vals = [func(mol) for _, func in Descriptors.descList if _ in DESCRIPTOR_NAMES]
    # Ensure order matches DESCRIPTOR_NAMES exactly
    desc_dict = {}
    for name, func in Descriptors.descList:
        if name in DESCRIPTOR_NAMES:
            try:
                desc_dict[name] = func(mol)
            except Exception:
                desc_dict[name] = float("nan")

    maccs = MACCSkeys.GenMACCSKeys(mol)
    maccs_dict = {f"MACCS_{i:03d}": int(maccs[i]) for i in range(1, 167)}

    row = {**{n: desc_dict.get(n, float("nan")) for n in DESCRIPTOR_NAMES}, **maccs_dict}
    return pd.DataFrame([row])


# ---------------------------------------------------------------------------
# Scaling (refit from stored X_train)
# ---------------------------------------------------------------------------


def make_scaler(X_train: np.ndarray, feature_cols: list[str]):
    """Return a StandardScaler fitted on the stored training matrix.

    Only RDKit descriptor columns are scaled; MACCS bits are left unchanged —
    matching the behaviour of ml_tools.scale_features().
    """
    from sklearn.preprocessing import StandardScaler

    rdkit_mask = [c in DESCRIPTOR_NAMES for c in feature_cols]
    rdkit_idx = [i for i, m in enumerate(rdkit_mask) if m]

    scaler = StandardScaler()
    if rdkit_idx:
        scaler.fit(X_train[:, rdkit_idx])
    return scaler, rdkit_idx


def scale_new_compound(x_raw: np.ndarray, scaler, rdkit_idx: list[int]) -> np.ndarray:
    """Apply the training scaler to a new compound's feature vector."""
    x = x_raw.copy()
    if rdkit_idx:
        x[:, rdkit_idx] = scaler.transform(x[:, rdkit_idx])
    return x


# ---------------------------------------------------------------------------
# Leverage computation
# ---------------------------------------------------------------------------


def compute_leverage(x_scaled: np.ndarray, XtX_inv: np.ndarray) -> float:
    """Compute hat-matrix diagonal element h = xᵀ (XᵀX)⁻¹ x."""
    return float(x_scaled @ XtX_inv @ x_scaled.T)


# ---------------------------------------------------------------------------
# Training-set membership check
# ---------------------------------------------------------------------------


def check_training_membership(smiles: str, split_csv_path: Path) -> dict:
    """Return membership info from the training-split CSV.

    Returns a dict with keys:
        in_training_set : bool
        split           : "train" | "test" | None
        y_log10         : float | None
        y_lower         : float | None
        y_upper         : float | None
    """
    result = {"in_training_set": False, "split": None, "y_log10": None, "y_lower": None, "y_upper": None}

    if not split_csv_path.exists():
        return result

    df = pd.read_csv(split_csv_path)
    match = df[df["smiles"] == smiles]
    if match.empty:
        return result

    row = match.iloc[0]
    result["in_training_set"] = True
    result["split"] = row.get("split")
    result["y_log10"] = row.get("y_log10")
    result["y_lower"] = row.get("y_lower")
    result["y_upper"] = row.get("y_upper")
    return result


# ---------------------------------------------------------------------------
# Core prediction
# ---------------------------------------------------------------------------


def predict(smiles: str, model_json_path: Path) -> dict:
    """Run the full inference pipeline for a single SMILES.

    Returns a dict with all prediction outputs and metadata.
    """
    # --- Load model card ---
    with open(model_json_path) as f:
        card = json.load(f)

    model_path = Path(card.get("model_file", model_json_path.with_suffix(".joblib")))
    if not model_path.is_absolute():
        model_path = model_json_path.parent / model_path
    if not model_path.exists():
        # Fall back: assume .joblib lives next to the .json
        model_path = model_json_path.with_suffix(".joblib")

    # abit of monkeypatching since earlier modelbuilding used absolute paths.
    ad_npz_raw = card["ad_artefact_file"]
    ad_npz_path = Path(ad_npz_raw) if Path(ad_npz_raw).is_absolute() else model_json_path.parent / ad_npz_raw
    # abit of monkeypatching since earlier modelbuilding used absolute paths.
    split_raw = card.get("training_split_file", "")
    split_csv_path = (
        (Path(split_raw) if Path(split_raw).is_absolute() else model_json_path.parent / split_raw) if split_raw else Path("")
    )
    feature_cols = card["feature_columns"]
    compartment = card["compartment"]

    # --- Load SVR ---
    svr = joblib.load(model_path)

    # --- Load AD artefact ---
    npz = np.load(ad_npz_path, allow_pickle=True)
    X_train_raw = npz["X_train"]  # (n_train × p) unscaled
    XtX_inv = npz["XtX_inv"]  # (p × p)
    h_star = float(npz["h_star"][0])

    # --- Compute features ---
    feat_df = compute_features(smiles)
    if feat_df is None:
        return {"error": f"RDKit could not parse SMILES: {smiles!r}"}

    # Select and order to match training feature set
    missing = [c for c in feature_cols if c not in feat_df.columns]
    if missing:
        return {"error": f"Missing features for SMILES (RDKit mismatch): {missing[:5]}..."}

    x_raw = feat_df[feature_cols].values  # shape (1, p)

    # Replace inf/NaN with training column medians
    col_medians = np.nanmedian(X_train_raw, axis=0)
    x_raw = np.where(~np.isfinite(x_raw), col_medians, x_raw)

    # --- Scale using training statistics ---
    scaler, rdkit_idx = make_scaler(X_train_raw, feature_cols)
    x_scaled = scale_new_compound(x_raw, scaler, rdkit_idx)

    # --- Predict ---
    y_log10_pred = float(svr.predict(x_scaled)[0])
    y_days_pred = 10**y_log10_pred

    # --- AD check ---
    h = compute_leverage(x_scaled[0], XtX_inv)
    in_ad = h <= h_star

    # --- Training-set membership ---
    membership = check_training_membership(smiles, split_csv_path)

    return {
        "smiles": smiles,
        "compartment": compartment,
        "predicted_log10_half_life": round(y_log10_pred, 4),
        "predicted_half_life_days": round(y_days_pred, 2),
        "applicability_domain": {
            "in_ad": in_ad,
            "leverage_h": round(h, 6),
            "h_star": round(h_star, 6),
            "note": (
                "Leverage-based structural AD. "
                "h > h* indicates the compound is structurally dissimilar "
                "to the training set — prediction reliability is reduced."
            ),
        },
        "training_set_membership": membership,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def print_report(result: dict) -> None:
    """Print a human-readable prediction report to stdout."""
    if "error" in result:
        print(f"ERROR: {result['error']}")
        return

    ad = result["applicability_domain"]
    mem = result["training_set_membership"]

    ad_flag = "IN AD" if ad["in_ad"] else "OUTSIDE AD"
    ad_note = "" if ad["in_ad"] else "  [!] Reduced reliability — structurally dissimilar to training set"

    if mem["in_training_set"]:
        split_label = mem["split"].upper()
        mem_note = (
            f"  This exact SMILES was in the {split_label} set  "
            f"(y_log10={mem['y_log10']:.3f}, "
            f"interval=[{mem['y_lower']:.3f}, {mem['y_upper']:.3f}] log10 days)"
        )
    else:
        mem_note = "  Not found in training or test set (novel compound)"

    print()
    print("=" * 60)
    print(f"  Prediction — {result['compartment'].upper()} compartment")
    print("=" * 60)
    print(f"  SMILES                 : {result['smiles']}")
    print(
        f"  Predicted half-life    : {result['predicted_half_life_days']:.2f} days"
        f"  (log10 = {result['predicted_log10_half_life']:.4f})"
    )
    print(f"  Applicability domain   : {ad_flag}  (h={ad['leverage_h']:.4f}, h*={ad['h_star']:.4f}){ad_note}")
    print(f"  Training set           :{mem_note}")
    print("=" * 60)
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict half-life interval for one or more SMILES strings.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--smiles",
        "-s",
        help="Single SMILES string",
    )
    group.add_argument(
        "--smiles-file",
        "-f",
        type=Path,
        help="Path to a plain-text file with one SMILES per line, no header.",
    )
    parser.add_argument(
        "--model",
        "-m",
        required=True,
        type=Path,
        help="Path to the model JSON card (e.g. models/SVR_air_hsbd_<ts>.json).",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional path to write results as JSON (one object per line).",
    )
    args = parser.parse_args()

    smiles_list: list[str] = []
    if args.smiles:
        print(f"Predicting for single SMILES: {args.smiles}")
        smiles_list = [args.smiles.strip().strip("'")]  # necessary for windows, not linux but doesn't hurt.
        print(smiles_list)
    else:
        raw = args.smiles_file.read_text().splitlines()
        smiles_list = [s.strip() for s in raw if s.strip() and not s.startswith("#")]

    # Validate model JSON path
    if not args.model.suffix:
        args.model = args.model.with_suffix(".json")
    elif args.model.suffix != ".json":
        args.model = args.model.with_suffix(".json")

    # ---------------------------------------------------------------------------
    # The actual prediction loop
    # ---------------------------------------------------------------------------
    results = []
    for smi in smiles_list:
        result = predict(smi, args.model)
        print_report(result)
        results.append(result)

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(results, fh, indent=2)
        print(f"Results written to {args.json_out}")


if __name__ == "__main__":
    main()
