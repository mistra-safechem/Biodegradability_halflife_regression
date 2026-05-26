"""For manual analysis / debug, outputfile not actually used"""

from pathlib import Path

import pandas as pd
from pandas import read_csv

print("Current working directory:")
work = Path.cwd()
print("Script directory:")
work = Path(__file__).resolve().parent

f1 = work / "missing_names_combo_manual_corrected_ready2merge.csv"
f2 = work / "missing_names_combo_all_compartments.csv"


f1_df = pd.read_csv(f1, header=0, sep="\t")
f2_df = pd.read_csv(f2, header=None, sep="\t", names=["Original_Name"])


corrected_names_df = pd.concat([f2_df, f1_df], axis=1)

output_path = work / "_tmp_optional_missing_vs_original_vs_corrected.csv"
print(f"Output will be written to: {output_path}")
corrected_names_df.to_csv(output_path, index=False, sep="\t", header=False)
print(corrected_names_df.head())
