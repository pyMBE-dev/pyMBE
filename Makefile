#!/bin/bash
.PHONY: run
.PHONY: visual 
.PHONY: clean

run:
	~/program/espresso/build/pypresso minimal_example.py

visual:
	python3 vmd-traj.py
	vmd -e visualization.tcl

clean: 
	rm -r *.txt *.xyz __pycache__
