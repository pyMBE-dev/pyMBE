#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean

docs:
	mkdir -p ./documentation
	pdoc ./pyMBE.py -o ./documentation --docformat google

tests:
	python3 testsuite/lj_tests.py
	python3 testsuite/generate_perpendicular_vectors_test.py
	python3 testsuite/create_molecule_position_test.py
	python3 testsuite/read-write-df_test.py
	python3 testsuite/henderson_hasselbalch_tests.py
	python3 testsuite/cph_ideal_tests.py
	python3 testsuite/grxmc_ideal_tests.py
	python3 testsuite/peptide_tests.py
	python3 testsuite/weak_polyelectrolyte_dialysis_test.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	jupyter-lab tutorials/pyMBE_tutorial.ipynb
