#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean

docs:
	mkdir -p ./documentation
	pdoc ./pyMBE.py -o ./documentation --docformat google

tests:
	python3 testsuite/peptide_tests.py
	python3 testsuite/set_particle_acidity_test.py

sample:
	python3 sample_scripts/peptide_simulation_example.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	jupyter-lab tutorials/pyMBE_tutorial.ipynb
