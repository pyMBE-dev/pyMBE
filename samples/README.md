# Pipeline of the sample scripts in pyMBE

## Production scripts
Production scripts show examples on how to setup various systems with pyMBE and ESPResSo.
These scripts sample the systems for a specific set of conditions and take as argparse arguments various inputs for the simulations (for example, the pH of the solution).
When run, the production script collect the time series of various quantities and store them in CSV files for later postprocessing.
Such CSV files are systematically named using the input argparse arguments, allowing to backtrace from which specific system are the corresponding time series.
Examples of production scripts are: `branched_polyampholyte.py`, `peptide.py`, and `peptide_mixture_grxmc_ideal.py`. 

## Analysis scripts
Analysis scripts show examples on how to analyze the time series produced with the production scripts of pyMBE.
These scripts read the time series stored from the production scripts and post-process them, calculating the ensemble mean, the error of the mean and the auto-correlation time using the block analysis method. [^1]
These quantities are stored together with their corresponding input conditions, extracted from the filename, into a wrapper CSV file for plotting and further analysis.
Examples of analylis scripts are: `analyze_time_series.py`. 

## Plotting scripts
Plotting scripts show examples on how to plot data post-processed with the analysis scripts of pyMBE and on how to use the toolbox of pyMBE to calculate various analytical solutions.
Examples of plotting scripts are: `plot_branched_polyampholyte.py`, `plot_peptide.py`, and `plot_peptide_mixture_grxmc_ideal.py`.

[^1]: Janke, W. (2002). Statistical analysis of simulations: Data correlations and error estimation. Quantum simulations of complex many-body systems: from theory to algorithms, 10, 423-445. 