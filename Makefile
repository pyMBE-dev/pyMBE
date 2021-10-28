#!/bin/bash
.PHONY: run
.PHONY: visual 
.PHONY: clean

ESPResSo_build_path=~/espresso4.1.4/build

run:
	${ESPResSo_build_path}/pypresso  peptide_simulation_example.py

visual:
	python3 vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	${ESPResSo_build_path}/ipypresso notebook sugar_tutorial.ipynb