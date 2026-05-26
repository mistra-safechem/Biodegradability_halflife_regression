# Missing data Curation

## Workflow
- Incorrect structures from HSBD dataset (see `0a1_1_model_data_extractor_hsbd.ipynb`) saved as  `missing/missing_smiles_{compartment}.txt` containing compound names (unknown/faulty ones).
- merge these four files to single (lots of overlap) with `merge_missing_smiles_to_combofile.py` giving `missing_names_combo_all_compartments.csv`.
- imported this csv manually to EXCEL from where manual curation was/is done (expert understanding or manufacturer lookup esp the "PCBs"). Ambiguous namings or not found compounds remained as N/A. NOTE: the pubchem retrieval was never implemented, this might have been able to resolve more of these compounds, e.g. the PCB ones.
For PCBs: https://www.carlroth.com/
outputfile should be called: `missing_names_combo_manual_corrected.xlsx`
- this excel then exported manually to final file `missing_names_combo_manual_corrected_ready2merge.csv` used second (optional usage) of notebook `0a1_2_optional_missing_hsbd_data.ipynb`.
- Folder `_hsbd_tmp` contains temporary files for optional case of merging with missing. These might still contain duplicate structures!
- note this was not necessary for the vega data which was already curated by source.

**manual, visual check** performed to confirm that duplicate structures contain same t_half values. 

## optional
An optional py script was made mainly for visual comparison/debugging and outputfile isn't used anywhere!
`optional_verify_merge_missing_vs_original_vs_corrected.py`.

## LLM based correction - experimental - not used -
An LLM based version `localollama_check_name.py` was attempted which works in principle, but local models are not powerful enough, even a chemistry centric model. Of note from later testing in other projects, large LLMs (on cloud) do a much better job which would allow for minimal manual curation afterwards. But a question of token costs. For this project here, this line of approach was dismissed.