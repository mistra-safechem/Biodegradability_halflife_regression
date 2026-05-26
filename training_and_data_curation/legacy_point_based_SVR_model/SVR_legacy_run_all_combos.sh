#!/bin/bash
# Script that will run all combinations of compartments vs all combos of data sources
# note: there is no air/vega combo
#
# replace uv run with python depending on your environment
#
uv run ./SVR_legacy_model_and_analysis.py --compartment air --data-source hsbd
uv run ./SVR_legacy_model_and_analysis.py --compartment water --data-source hsbd
uv run ./SVR_legacy_model_and_analysis.py --compartment water --data-source vega
uv run ./SVR_legacy_model_and_analysis.py --compartment soil --data-source hsbd
uv run ./SVR_legacy_model_and_analysis.py --compartment soil --data-source vega
uv run ./SVR_legacy_model_and_analysis.py --compartment sediment --data-source hsbd
uv run ./SVR_legacy_model_and_analysis.py --compartment sediment --data-source vega
uv run ./SVR_legacy_model_and_analysis.py --compartment water --data-source combined
uv run ./SVR_legacy_model_and_analysis.py --compartment soil --data-source combined
uv run ./SVR_legacy_model_and_analysis.py --compartment sediment --data-source combined