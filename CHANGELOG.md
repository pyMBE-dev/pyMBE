# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Switched from `os.makedirs` to `Path().mkdir()` to prevent ocasional failure of the scripts when running them in paralel. (#91)
- `pmb.set_reduced_units()` now redefines the reduced units instead of creating a new instance of `pint.UnitRegistry`. Therefore, the user can do operations between objects defining before and after changing the set of reduced units without getting a `ValueError` (#89)
- `pmb.set_reduced_units()` now checks that the arguments provided by the user have the right dimensionality. (#89)
- The constants stored as attributes in `pyMBE.pymbe_library` are now using their values stabilished in the 2019 SI. The value is taken directly from `scipy.constants` instead of being a hard-coded constant. (#86)
- Switched to Ctest for testing, allowing to run the tests on paralel (#87)

### Added
- Unit testing for reaction methods (#86)

### Fixed
- Removed global state variables, instead they are now created by the constructor of `pyMBE.pymbe_library`. This prevents two instances of the pyMBE library to share the same memory address for their attributes. (#89)
- Required Python dependency versions compatible with ESPResSo 4.2 (#84)
- Fixed several deprecated paths and function names in `tutorials/pyMBE_tutorial.ipynb`. (#77, #78, #79, #80, #81)

## [0.8.0] - 2024-06-18

### Added

* Initial release of pyMBE.
