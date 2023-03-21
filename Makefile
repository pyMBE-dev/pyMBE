#!/bin/bash
.PHONY: sample
.PHONY: visual 
.PHONY: clean
.PHONY: tests

ESPResSo_build_path=~/espresso4.2/build

sample:
	${ESPResSo_build_path}/pypresso sample_scripts/peptide_simulation_example.py

visual:
	python3 handy_scripts/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	${ESPResSo_build_path}/ipypresso notebook sugar_tutorial.ipynb

tests:
	 ${ESPResSo_build_path}/pypresso tests/LYS_ASP_peptide.py
	 ${ESPResSo_build_path}/pypresso tests/GLU_HIS_peptide.py
	 ${ESPResSo_build_path}/pypresso tests/histatin5_peptide.py
