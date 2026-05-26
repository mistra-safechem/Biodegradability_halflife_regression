#!/bin/bash
# example for inference testing (diclofenac = C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl )
# reminder: there is no air/vega combo
#
# replace uv with Python depending on your environment
uv run ./SVR_interval_inference.py --smiles 'C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl' --model models/SVR_air_hsbd_20260516_145335.json  --json-out test.json
