#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean
.PHONY: testsuite
.PHONY: docs

docs:
	pdoc ./pyMBE.py -o ./docs --docformat google 

testsuite:
	python3 testsuite/peptide_tests.py

sample:
	python3 sample_scripts/peptide_simulation_example.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	jupyter-lab tutorials/pyMBE_tutorial.ipynb
