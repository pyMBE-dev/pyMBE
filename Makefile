#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean

PYTHON = python3

docs:
	mkdir -p ./documentation
	PDOC_ALLOW_EXEC=0 ${PYTHON} -m pdoc ./pyMBE.py -o ./documentation --docformat google

unit_tests:
	${PYTHON} testsuite/lj_tests.py
	${PYTHON} testsuite/set_particle_acidity_test.py
	${PYTHON} testsuite/bond_tests.py
	${PYTHON} testsuite/generate_perpendicular_vectors_test.py
	${PYTHON} testsuite/create_molecule_position_test.py
	${PYTHON} testsuite/seed_test.py
	${PYTHON} testsuite/read-write-df_test.py
	${PYTHON} testsuite/parameter_test.py
	${PYTHON} testsuite/henderson_hasselbalch_tests.py

functional_tests:
	${PYTHON} testsuite/cph_ideal_tests.py
	${PYTHON} testsuite/grxmc_ideal_tests.py
	${PYTHON} testsuite/peptide_tests.py
	${PYTHON} testsuite/gcmc_tests.py
	${PYTHON} testsuite/weak_polyelectrolyte_dialysis_test.py
	${PYTHON} testsuite/globular_protein_tests.py

tests: unit_tests functional_tests

sample:
	${PYTHON} samples/peptide.py

visual:
	${PYTHON} visualization/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	jupyter-lab tutorials/pyMBE_tutorial.ipynb

pylint:
	${PYTHON} -m pylint pyMBE.py lib/ testsuite/ samples/ maintainer/ visualization/
