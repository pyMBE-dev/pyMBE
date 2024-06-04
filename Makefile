#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean

# whether to run unit tests with code coverage
COVERAGE = 0

# output directory for the code coverage HTML files
COVERAGE_HTML = coverage

# Python executable or launcher, possibly with command line arguments
PYTHON = python3
ifeq ($(COVERAGE),1)
PYTHON := ${PYTHON} -m coverage run --parallel-mode --source=$(CURDIR)
endif

docs:
	mkdir -p ./documentation
	PDOC_ALLOW_EXEC=0 ${PYTHON} -m pdoc ./pyMBE.py -o ./documentation --docformat google

unit_tests:
	${PYTHON} testsuite/serialization_test.py
	${PYTHON} testsuite/lj_tests.py
	${PYTHON} testsuite/set_particle_acidity_test.py
	${PYTHON} testsuite/bond_tests.py
	${PYTHON} testsuite/generate_perpendicular_vectors_test.py
	${PYTHON} testsuite/define_and_create_molecules_unit_tests.py
	${PYTHON} testsuite/create_molecule_position_test.py
	${PYTHON} testsuite/seed_test.py
	${PYTHON} testsuite/read-write-df_test.py
	${PYTHON} testsuite/parameter_test.py
	${PYTHON} testsuite/henderson_hasselbalch_tests.py
	${PYTHON} testsuite/calculate_net_charge_unit_test.py
	${PYTHON} testsuite/setup_salt_ions_unit_tests.py
	${PYTHON} testsuite/globular_protein_unit_tests.py
	${PYTHON} testsuite/analysis_tests.py

functional_tests:
	${PYTHON} testsuite/cph_ideal_tests.py
	${PYTHON} testsuite/grxmc_ideal_tests.py
	${PYTHON} testsuite/peptide_tests.py
	${PYTHON} testsuite/gcmc_tests.py
	${PYTHON} testsuite/weak_polyelectrolyte_dialysis_test.py
	${PYTHON} testsuite/globular_protein_tests.py

tests: unit_tests functional_tests

coverage_xml:
	${PYTHON} -m coverage combine .
	${PYTHON} -m coverage report
	${PYTHON} -m coverage xml

coverage_html:
	${PYTHON} -m coverage combine .
	${PYTHON} -m coverage html --directory="${COVERAGE_HTML}"

sample:
	${PYTHON} samples/peptide.py

visual:
	${PYTHON} visualization/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	jupyter-lab tutorials/pyMBE_tutorial.ipynb

pylint:
	${PYTHON} -m pylint pyMBE.py lib/ testsuite/ samples/ maintainer/ visualization/
