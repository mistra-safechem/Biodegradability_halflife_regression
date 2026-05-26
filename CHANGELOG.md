# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.2.1] - 2026-05-26 — Post-release updates

### Changed
- typo fixing in readme files and minor adjustment to reflect the current structure.
- best effort to refactor legacy notebooks and py scripts to be runnable from their current location without adjustments, after files had been moved around without reexecuting (to preserve the original state as much as possible for archiving purposes).

---

## [0.2.0] - 2026-05-26 — Full Release (2/2): Modelling Pipeline

### Added
- Entire pipeline for data curation and model creation in `training_and_data_curation/`.
- See `training_and_data_curation/README.md` for usage and documentation

---

## [0.1.0] - 2026-05-24 — Initial Release (1/2): Inference

### Added
- Inference standalone scripts in `inference/` with packaged model assets in `inference/models/`
- See `inference/README.md` for usage and documentation
