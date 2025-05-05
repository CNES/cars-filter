# Changelog

## 0.2.0 Additional statistical filtering parameters (May 2025)

### Added

- New parameters `filtering_constant` and `mean_factor` in cars statistical filtering functions:

`threshold = filtering_constant + mean_factor * mean + dev_factor * std_dev`

## 0.1.2 Statistical filter bug fix (March 2025)

### Fixed

- Robustify percentile utilisation to nan points.

## 0.1.1 Cross-platforms wheels (February 2025)

### Added

 - Add wheels in PyPI package for Windows, MacOS and Ubuntu

### Changed

 - Support for Python 3.8 is over

## 0.1.0 First Official Release (November 2024)

- First open-source release
