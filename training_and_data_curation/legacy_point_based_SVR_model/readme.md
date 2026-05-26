# Legacy point-based SVR model for t_half prediction
This folder contains the code for the legacy point-based SVR model for t_half prediction. 
This was the initial regression development model for the project, and remains available for reference and comparison.

The focus shifted to the interval-based model due to the nature of the underlying data 
as described in the readme.md files in the other subfolders where more details on the data and the modelling approach are provided.

Even older files related to the model development are retained in the `legacy_optimization_scripts` folder, 
such as the RFR, XBGoost approach, etc - see the readme.md file in that folder for details:
[legacy_optimization_scripts/readme.md](legacy_optimization_scripts/readme.md).

## How to run
Run everything from within this subfolder, it will call the necessary files from the src/ folder, and load the data from the data/ folder, 
no need to adjust file paths or copy files from other folders.

