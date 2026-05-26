# About Log Review Tools
Script for automatic log file aggregation and summary.

For the legacy point based model and logs this differs from the main intervall based scripts due to, well, legacy, differnt log file formats, differnt metrics etc.

## Script 1: Aggregates SVR output
Automates concatenation of `logs/combined_logs.txt` from all individual SVR output log files 
in the `logs/` directory, further parsing done with Script 2.

The most basic version of log file concatenation to obtain an overall summary of all SVR runs.

## Script 2: Aggregates all the log files
Automates creation `logs/log_review_<timestamp>.md` from aggregated SVR output in `logs/combined_logs.txt`. Script parses each dataset block, extracts key performance, learning-curve, and applicability-domain metrics, into markdown table formats, for a quick overview of all models' performances and some interpretation of the results.

## Script 3: Optional - chemistry analysis 
Automates generation of `chemistry_analysis_interpretation.md` from the chemistry analysis subfolder data. This script parses the chemistry analysis outputs, extracts key insights about chemical space coverage, and generates a narrative interpretation of the chemistry analysis results. 

## Script 4: Generates presentation-ready summary
Automates creation of `logs/log_summary_presentation.md` from `logs/log_review_<timestamp>.md`.


