# About data extraction, cleaning, feature calculation, and database building
Scripts for data extraction from sources, cleaning of structures, followed by 
feature calculation, and database creation.

# How to run
- Notebooks can (should) be run from within this subfolder directly.
- Start in order of filenames. 
- For explanatory details on data curation, especially the "optional" steps, see the readme.md in "processed_data":
    [..\processed_data\readme.md](..\processed_data\readme.md)

# Interval mapping files (required for DB notebooks)
Two interval mapping CSV files in this folder are used to map continuous half-life values
to interval-based columns stored in the databases.

- `hsbd_interval_mapping.csv`
    - Used by notebook: `0b1_create_hsbd_database.ipynb`

- `vega_interval_mapping.csv`
    - Used by notebook: `0b2_create_vega_database.ipynb`

For details see this [training_and_data_curation/raw_data/readme.md](../raw_data/readme.md).

# New and legacy DB
The created databases are compatible with the legacy code (db_utils.py has 
functions for loading the old and new structure; db really only differs 
with respect to the added interval t_half columns).