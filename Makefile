#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean
.PHONY: tests
.PHONY: docs

ESPResSo_build_path=~/espresso4.2/build

docs:
	pdoc ./pyMBE.py -o ./docs --docformat google 

sample:
	${ESPResSo_build_path}/pypresso sample_scripts/peptide_simulation_example.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	${ESPResSo_build_path}/ipypresso notebook sugar_tutorial.ipynb

tests_peptide:
	 ${ESPResSo_build_path}/pypresso tests/LYS_ASP_peptide.py
	 ${ESPResSo_build_path}/pypresso tests/GLU_HIS_peptide.py
	 ${ESPResSo_build_path}/pypresso tests/histatin5_peptide.py
	 
tests_globular_protein:
	python3 tests/run_test_protein.py --pdb_code 1beb --run_command "${ESPResSo_build_path}/pypresso sample_scripts/globular_protein.py  --pdb 1beb --path_to_cg reference_parameters/coarse_grained_structures/1beb.vtf"
	python3 tests/run_test_protein.py --pdb_code 1f6s --run_command "${ESPResSo_build_path}/pypresso  sample_scripts/globular_protein.py  --pdb 1f6s --metal_ion_name Ca --metal_ion_charge 2 --path_to_cg reference_parameters/coarse_grained_structures/1f6s.vtf"

tests:
	tests_peptide
	tests_globular_protein