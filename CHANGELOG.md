# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- `pd.NA` is now enforced as value for empty cells in `pmb.df`, this prevents transformation of variable types from `int` to `float` which was breaking the code when reading `pmb.df` from file (See #102). (#116)
- The sample script `plot_HH.py` has been replaced for specific examples on how to plot data post-processed with pyMBE: `plot_branched_polyampholyte.py`, `plot_peptide.py`, and `plot_peptide_mixture_grxmc_ideal.py`. (#95)
- Sample scripts now take the pH as an argparse input instead of looping over several pH values. This enables paralization of the sample scripts and avoids conflicts with the current post-processing pipeline. (#95)
- Switched from `os.makedirs` to `Path().mkdir()` to prevent ocasional failure of the scripts when running them in paralel. (#91)
- `pmb.set_reduced_units()` now redefines the reduced units instead of creating a new instance of `pint.UnitRegistry`. Therefore, the user can do operations between objects defining before and after changing the set of reduced units without getting a `ValueError` (#89)
- `pmb.set_reduced_units()` now checks that the arguments provided by the user have the right dimensionality. (#89)
- The constants stored as attributes in `pyMBE.pymbe_library` are now using their values established in the 2019 SI. Their value are taken directly from `scipy.constants` instead of being hard-coded constants. (#86)
- Switched to CTest for testing, allowing to run the tests on paralel (#87)

### Added
- Code of conduct of our community `CODE_OF_CONDUCT.md`, adhering to the Contributor Covenant v2.1 (#104) 
- New optional argument `backbone_vector` enabling to build molecules along an input vector using `pmb.create_molecule` and `pmb.create_pmb_object` (#99)
- New boolean flag `--ideal` as argparse argument of `samples/globular_protein.py` enabling to run the script without setting up interactions.
- Unit tests for `pmb.create_protein`, `pmb.enable_motion_of_rigid_object`, `pmb.protein_sequence_parser`, `pmb.define_protein`, `pmb.read_protein_vtf_in_df` (#101)
- Library `lattice.py`, a general builder for crystalline lattices. This library is part of on-going project to support hydrogels in pyMBE. (#93)
- New sample script showing how to use the analysis tools in pyMBE for post-processing time series from the sample scripts `analyze_time_series.py` (#95) 
- A new optional argument `ignore_files`  for `lib.analysis.analyze_time_series`, enabling to provide a list of files to be ignored for post-processing of time series. (#95)
- Functional testing for all sample scripts. (#95)
- Unit testing for reaction methods, bonds, object serialization. (#86, #113)

### Fixed
- Writing and reading `pmb.df` from file does no longer change the variable type from `int` to `float` when there are empty cells in the column. (#116)
- Espresso bond objects stored in `pmb.df` now retain the same value for the  `._bond_id` attribute as the original Espresso objects. (#116)
- Wrong parsing in `pmb.protein_sequence_parser` of input sequences provided as a list of aminoacids using the three letter code. (#101)
- Wrong setup of the rigid object in `pmb.enable_motion_of_rigid_object`, leading to crashes in `samples/globular_protein.py` when enabling the protein motion. (#101)
- The argparse argument `--move_protein` of `samples/globular_protein.py` is now a boolean flag instead of taking arbitrary float values. (#101)
- `lib.analysis.get_dt` now raises a ValueError if the two first two rows of the dataframe have the same values for the time, which break the subsequent code. (#95)
- Removed global state variables, instead they are now created by the constructor of `pyMBE.pymbe_library`. This prevents two instances of the pyMBE library to share the same memory address for their attributes. (#89)
- Required Python dependency versions compatible with ESPResSo 4.2 (#84)
- NumPy 2, Pandas 2 and the development version of ESPResSo are now fully supported. (#106)
- Fixed several deprecated paths and function names in `tutorials/pyMBE_tutorial.ipynb`. (#77, #78, #79, #80, #81)
- Fixed error handling mechanism in the bond search function when the searched bond doesn't exist. (#113)
- Unsupported ESPResSo features now raise a `NotImplementedError` instead of a `ValueError`. (#113)

### Removed
- `pmb.parse_sequence_from_file` has been removed since it is no longer necesary to parse the sequence from pmb.df (#110)

## [0.8.0] - 2024-06-18

### Added

* Initial release of pyMBE.
