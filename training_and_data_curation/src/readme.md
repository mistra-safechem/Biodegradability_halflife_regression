# Helper functions folder
This folder contains helper functions for the main scripts, such as data processing, descriptor calculation, and model evaluation.

# Regarding RDKit
The MistraSafeChem code environment uses 2022 RDKit which has fewer descriptors than 2025 RDKit, which was originally used here.

Thus, keeping both descriptors files (rdkitdescriptors*.txt) as reference, but using the 2022 list even in the 2025 environment.

For more information see [`src/feature_selection_in_rdkit_tools.py.md`](src/feature_selection_in_rdkit_tools.py.md).