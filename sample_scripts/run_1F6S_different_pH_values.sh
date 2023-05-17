#!/bin/sh
../pypresso 1F6S_globular_protein.py -pH 2.0
../pypresso 1F6S_globular_protein.py -pH 2.5
../pypresso 1F6S_globular_protein.py -pH 3.0
../pypresso 1F6S_globular_protein.py -pH 3.5
../pypresso 1F6S_globular_protein.py -pH 4.0
../pypresso 1F6S_globular_protein.py -pH 4.5
../pypresso 1F6S_globular_protein.py -pH 5.0
../pypresso 1F6S_globular_protein.py -pH 5.5
../pypresso 1F6S_globular_protein.py -pH 6.0
../pypresso 1F6S_globular_protein.py -pH 6.5
../pypresso 1F6S_globular_protein.py -pH 7.0

mkdir observables_results
mv pH-* observables_results
python3 ../handy_scripts/data_analysis.py observables_results/
python3 plot_analized_observables.py

