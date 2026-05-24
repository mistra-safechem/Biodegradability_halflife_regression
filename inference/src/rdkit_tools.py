"""
Module for calculating RDKit molecular descriptors and MACCS fingerprints.

SMILES standardizer as used in MistraSafeChem (MSC) project
Two versions, old and new (better handling of pandas SettingWithCopyWarning)

In case rdkit2025 features were calculated, this will be reduced to the
2022 descriptors as used by the standard pipeline in Mistrasafechem project
with one manual addition, solubility logS.
"""

import os
from typing import Any, List

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.DataStructs import ConvertToNumpyArray
from tqdm import tqdm

RDLogger.logger().setLevel(RDLogger.CRITICAL)  # disable rdkit warnings

CURRENT_SRC_DIR = os.path.dirname(os.path.abspath(__file__))

_descriptor_names = []
for line in open(os.path.join(CURRENT_SRC_DIR, "rdkitdescriptors2022.txt"), "r"):
    line = line.strip()
    if line:
        _descriptor_names.append(line)
_descriptor_names.append("ESOL_logS")

# remove some manually curated descriptors since deemed redundant or not useful for the current task
# for original database curation, these are already excluded, but for inference on new data this is best kept here for consistency.
# alt: edit the rdkitdescriptors2022.txt file (call it e.g. "descriptors,txt"), but this way it's more explicit which descriptors are dropped and why
manual_drop_features = ["ExactMolWt", "HeavyAtomMolWt", "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3"]
DESCRIPTOR_NAMES = [name for name in _descriptor_names if name not in manual_drop_features]
MACCS_NAMES = [f"MACCS_{i:03d}" for i in range(1, 167)]  # MACCS keys names

descriptor_names_short = [
    name for name in DESCRIPTOR_NAMES if name != "ESOL_logS"
]  # remove ESOL_logS for descriptor calculation


def _compute_MACCS(smiles: str) -> dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {f"MACCS_{i:03d}": float("nan") for i in range(1, 167)}

    fp = MACCSkeys.GenMACCSKeys(mol)  # 167 bits; bit 0 is unused

    arr = np.zeros((fp.GetNumBits(),), dtype=int)
    ConvertToNumpyArray(fp, arr)
    return {f"MACCS_{i:03d}": arr[i] for i in range(1, 167)}


def _compute_descriptors(smiles: str) -> dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {name: float("nan") for name, _ in descriptor_names_short + [("ESOL_logS", None)]}

    descr = {}
    for name in descriptor_names_short:
        func = getattr(Descriptors, name)
        descr[name] = func(mol)

    _esol_logS = (
        0.16
        - 0.64 * descr["MolLogP"]
        - 0.0062 * descr["MolWt"]
        + 0.066 * descr["NumRotatableBonds"]
        - 0.74 * (sum(1 for a in mol.GetAtoms() if a.GetIsAromatic()) / (mol.GetNumHeavyAtoms() or 1))
    )
    descr["ESOL_logS"] = _esol_logS
    return descr


def calculate_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate RDKit molecular descriptors for a dataframe with a 'Canonical_smiles' column.
    and concatenate the results to the original dataframe.
    """

    try:
        df["Canonical_smiles"] = df["Canonical_smiles"].astype(str)
    except KeyError:
        raise KeyError("DataFrame must contain a 'Canonical_smiles' column.")

    descriptor_data = []
    for smi in tqdm(df["Canonical_smiles"], desc="Calculating descriptors"):
        descriptor_data.append(_compute_descriptors(smi))
    descriptor_df = pd.DataFrame(descriptor_data)

    maccs_data = []
    for smi in tqdm(df["Canonical_smiles"], desc="Calculating MACCS fingerprints"):
        maccs_data.append(_compute_MACCS(smi))
    maccs_df = pd.DataFrame(maccs_data)
    return pd.concat([df.reset_index(drop=True), descriptor_df, maccs_df], axis=1)


def _smiles_standardizer_msc_old_pd(smiles: str) -> List[str]:
    """Smiles standardizer as used in MistraSafeChem (MSC) project
    Adapted for single SMILES input
    a bit overkill using df but don't want to change the original function too much

    this was written for older pandas versions where SettingWithCopyWarning was not an issue
    keep it for reference
    """

    smiles_list = [smiles]
    _df = pd.DataFrame(smiles_list, columns=["smiles"])
    _df_st = pd.DataFrame(columns=["smiles"], index=range(0, len(_df)))

    for i in range(0, len(_df)):
        if "\\" in _df.smiles[i]:
            _df.smiles[i] = _df.smiles[i].replace("\\", "")
        if "/" in _df.smiles[i]:
            _df.smiles[i] = _df.smiles[i].replace("/", "")
        m = Chem.MolFromSmiles(_df.smiles[i], sanitize=False)
        m_sanitize = Chem.MolFromSmiles(_df.smiles[i], sanitize=True)
        if m_sanitize is not None:
            _df_st.loc[i, ["smiles"]] = Chem.MolToSmiles(m_sanitize)
        else:
            m.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(
                m,
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True,
            )
            _df_st.loc[i, ["smiles"]] = Chem.MolToSmiles(m)
    return _df_st["smiles"].tolist()


def smiles_standardizer_msc(smiles: str) -> List[str]:
    """Smiles standardizer as used in MistraSafeChem (MSC) project
    Adapted for single SMILES input
    a bit overkill using df but don't want to change the original function too much

    """

    smiles_list = [smiles]
    _df = pd.DataFrame(smiles_list, columns=["smiles"])
    _df_st = pd.DataFrame(columns=["smiles"], index=range(len(_df)))

    for i in range(len(_df)):
        s = str(_df.loc[i, "smiles"])
        if "\\" in s:
            s = s.replace("\\", "")
        if "/" in s:
            s = s.replace("/", "")
        _df.loc[i, "smiles"] = s
        m = Chem.MolFromSmiles(s, sanitize=False)
        m_sanitize = Chem.MolFromSmiles(s, sanitize=True)
        if m_sanitize is not None:
            _df_st.loc[i, "smiles"] = Chem.MolToSmiles(m_sanitize)
        else:
            m.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(
                m,
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True,
            )
            _df_st.loc[i, "smiles"] = Chem.MolToSmiles(m)
    return _df_st["smiles"].tolist()
