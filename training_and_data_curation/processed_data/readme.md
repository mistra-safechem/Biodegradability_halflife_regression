# Data curation

## Raw files
Raw files stored under `raw_data/` :
* hsbd_data 
* VEGA_paper_data

check the readme.md in that folder regarding source & data details

---
## Pipeline
Code is found under `0_1_data_clean_descr_calcs_db_build` folder, files are executed in order of filenames.

---
### HSBD data
#### Curation
*Step 1*

notebook: `0a1_1_model_data_extractor_hsbd_data.ipynb`
For HSBD data, due to missing identifiers (smiles, inchi, cas) or incorrect structure names:

- File read with stripping of certain characters.
- Retrieve missing entries (between 5 and 15%)
- check with OPSIN if names can match a smiles
- insert retrieved smiles into datasets (now only 3 - 7% missing)
- remaining missing exported for manual curation
- curation of missing is OPTIONAL
- a first pass smiles cleaning is performed within.
- csv files are created for all compartments

*Step 1b - OPTIONAL*
Since this is a manual expert guided curation step, this is optional. A ready curated file 
is available, no need to rerun those optional steps. Details can be found in `processed_data/missing/readme.md`.

- Once finished, run notebook: `0a1_2_optional_missing_data.ipynb`.
- csv files are created for all compartments


*Step 2*
notebook: `0a2_model_data_extractor_vega_data.ipynb`

- simpler workflow since data already is curated by vega.
- csv files are created for all compartments (not air)

---
### Vega Data
No curation steps are needed, data is already curated by VEGA team. Only step is to extract 
the data and create csv files for all compartments (not air) with notebook: `0a2_model_data_extractor_vega_data.ipynb`.

----
### Database building
Data is handled via SQLite databases. 
- notebooks: `0b1_create_hsbd_database.ipynb` & `0b2_create_vega_database.ipynb` for db and table build
- notebooks: `1a1_descr_hsbd_calc.ipynb` & `1a2_descr_calc_vega.ipynb` calculates features and adds to db
- notebook: `1a3_combined_hsbd_vega.ipynb` combines the data of the two datasets into third db

---
### Files required after pipeline:
- only the database files are required for the next steps, csv files can be deleted if space or file noise is an issue.