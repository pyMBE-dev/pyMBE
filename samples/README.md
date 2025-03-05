# Pipeline of the sample scripts in pyMBE

## Production scripts
Production scripts show examples on how to setup various systems with pyMBE and ESPResSo.
These scripts sample the systems for a specific set of conditions and take as command-line arguments various inputs for the simulations (for example, the pH of the solution).
When run, the production scripts collect the time series of various quantities and store them in CSV files for future postprocessing.
Such CSV files are systematically named using the input arguments, allowing to backtrace from which specific system are the corresponding time series.
Examples of production scripts are: `branched_polyampholyte.py`, `peptide_cpH.py`, `peptide_mixture_grxmc_ideal.py` and `salt_solution_gcmc.py`. 

## Analysis scripts
Analysis scripts show examples on how to analyze the time series produced with the production scripts of pyMBE.
These scripts read the time series stored from the production scripts and post-process them, calculating the ensemble mean, the error of the mean and the auto-correlation time using the block analysis method. [^1]
These quantities are stored together with their corresponding input conditions, extracted from the filename, into a wrapper CSV file for plotting and further analysis.
Examples of analysis scripts are: `analyze_time_series.py`. 

## Plotting scripts
Plotting scripts show examples on how to plot data post-processed with the analysis scripts of pyMBE and on how to use the toolbox of pyMBE to calculate various analytical solutions.
Examples of plotting scripts are: `plot_branched_polyampholyte.py`, `plot_peptide_cpH.py`, and `plot_peptide_mixture_grxmc_ideal.py`.

[^1]: Janke, W. (2002). Statistical analysis of simulations: Data correlations and error estimation. In: Quantum simulations of complex many-body systems: from theory to algorithms, NIC Series volume 10, pp. 423-445. isbn: 3-00-009057-6. URL: <https://www.physik.uni-leipzig.de/~janke/Paper/nic10_423_2002.pdf>

## Example on how to use the pipeline
The sample scripts are designed to be used in the following order: (i) production script, (ii) analysis script and (iii) plotting script. For example:
```bash
python3 branched_polyampholyte.py # By default, stores the time series in `time_series/branched_polyampholyte`
python3 analyze_time_series.py --data_folder time_series/branched_polyampholyte # by default, stores the post-processed data in `time_series/branched_polyampholyte/analyzed_data.csv`
python3 plot_branched_polyampholyte.py # By default, reads the averages data in `time_series/branched_polyampholyte/analyzed_data.csv`
```