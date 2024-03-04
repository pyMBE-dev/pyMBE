#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean
.PHONY: testsuite
.PHONY: docs

ESPResSo_build_path=~/software/espresso_v4.2/build/

docs:
	pdoc ./pyMBE.py -o ./docs --docformat google 

testsuite:
	${ESPResSo_build_path}/pypresso testsuite/peptide_tests.py

sample:
	${ESPResSo_build_path}/pypresso sample_scripts/peptide_simulation_example.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	${ESPResSo_build_path}/ipypresso notebook sugar_tutorial.ipynb