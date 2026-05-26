from pathlib import Path

import pandas as pd
from pandas import read_csv

print("Current working directory:")
work = Path.cwd()
print("Script directory:")
work = Path(__file__).resolve().parent

f1 = work / "missing_smiles_air.txt"
f2 = work / "missing_smiles_soil.txt"
f3 = work / "missing_smiles_water.txt"
f4 = work / "missing_smiles_sediment.txt"

f1_df = pd.read_csv(f1, header=None, sep="\t", names=["abc"])
f2_df = pd.read_csv(f2, header=None, sep="\t", names=["abc"])
f3_df = pd.read_csv(f3, header=None, sep="\t", names=["abc"])
f4_df = pd.read_csv(f4, header=None, sep="\t", names=["abc"])

corrected_names_df = pd.concat([f1_df, f2_df, f3_df, f4_df], axis=0)
# sort df by abc column
corrected_names_df = corrected_names_df.sort_values(by="abc").reset_index(drop=True)
# remove duplicate rows
corrected_names_df = corrected_names_df.drop_duplicates().reset_index(drop=True)

output_path = work / "missing_names_combo_all_compartments.csv"
print(f"Output will be written to: {output_path}")
corrected_names_df.to_csv(output_path, index=False, sep="\t", header=False)
print(corrected_names_df.head())
