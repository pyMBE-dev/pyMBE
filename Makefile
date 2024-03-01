#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean
.PHONY: tests
.PHONY: testsuite
.PHONY: docs

docs:
	pdoc ./pyMBE.py -o ./docs --docformat google 

testsuite:
	python3 testsuite/LYS_ASP_peptide.py

sample:
	python3 sample_scripts/peptide_simulation_example.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tests_peptide:
	python3 tests/LYS_ASP_peptide.py
	python3 tests/GLU_HIS_peptide.py
	python3 tests/histatin5_peptide.py

tests_globular_protein:
	python3 tests/run_test_protein.py --pdb_code 1beb --run_command "python3 sample_scripts/globular_protein.py  --pdb 1beb --path_to_cg reference_parameters/coarse_grained_structures/1beb.vtf"
	python3 tests/run_test_protein.py --pdb_code 1f6s --run_command "python3 sample_scripts/globular_protein.py  --pdb 1f6s --metal_ion_name Ca --metal_ion_charge 2 --path_to_cg reference_parameters/coarse_grained_structures/1f6s.vtf"

tests:
	make tests_peptide
	make tests_globular_protein
