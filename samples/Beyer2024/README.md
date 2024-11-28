The scripts in this folder are designed to reproduce the data and plots showcased in our publication [^1].
To reproduce the data, one simply needs to run the following script:
```sh
python3 create_paper_data.py --fig_label 7a --mode long-run --plot
```
where the previous line will run the script to produce Fig. 7a in Ref.[^1] The user can use the argparse argument `--fig_label` to create any of the plots that we presented in that publication as benchmarks: 7a, 7b, 7c, 8a, 8b, 9. The argparse `--mode` controls the statiscal accuracy (i.e. the number of samples) that the script measures. The mode `long-run` should be used to generate data with the same statistical accuracy than in Ref.[^1]. The mode `short-run` can be used for a shorter run for testing or to trying out the scripts for each of our benchmarks:

- peptide.py: for the peptide benchmarks
- globular_protein.py: for the globular protein benchmarks
- weak_polyelectrolyte_dialysis.py: for the weak polyelectrolyte dialysis benchmarks

The optional argparse argument `--plot` controls if these scripts generate the corresponding plot or if the data is simply stored to file. We note that the format of the plots can differ from that of our publication [^1]. Theses scripts are part of the continous integration (CI) pipeline to ensure that future pyMBE releases still reproduce the benchmarks.

[^1]: D. Beyer, P. B. Torres, S. P. Pineda, C. F. Narambuena, J.-N. Grad, P. Košovan, P. M. Blanco. J. Chem. Phys.(2024), 161 (2), 022502. doi: [10.1063/5.0216389](https://doi.org/10.1063/5.0216389).
