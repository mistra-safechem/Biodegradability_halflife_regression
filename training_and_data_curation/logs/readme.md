# About log review tool scripts
Scripts for automatic log file aggregation and summary are provided in the `logs/` directory.

## Script 1: Aggregates SVR output
Automates concatenation of `logs/combined_logs.txt` from all individual SVR output log files 
in the `logs/` directory, further parsing done with Script 2.

The most basic version of log file concatenation to obtain an overall summary of all SVR runs.

## Script 2: Aggregates all the log files into tables
Automates creation `logs/log_summary.md` from aggregated SVR output in `logs/combined_logs.txt`.  

Script parses each dataset block, extracts key performance, learning-curve, and 
applicability-domain metrics, into markdown table format, for a quick overview of all models' performances,


## Script 3: Generates "presentation-ready" summary
Automates creation of `logs/log_summary_presentation.md` from `logs/log_summary.md`.

It further formates and shortens the summary tables to obtain a quick overview of key performance metrics across all models.


# Date of last run
Final (re)run 2026-05-16 using updated relative paths instead of previous absolute paths;
won't impact the results since they should be reproducible.
An earlier run from 2026-04-30 is available in `.earlier_runs/` for comparison.