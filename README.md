# Sugar library

Sugar is a library that builds on top of the Molecular Dynamics software ESPResSo 
providing tools to simplify the setting up of simulations of peptides and polyelectrolytes.
Internally, Sugar ensures that the units are properly handled before providing the inputs
to ESPResSo using Pint library. 


## Requirements

- ESPResSo v4.2 [link to Espresso](https://espressomd.org/wordpress/download/)
- Pint v0.17 [link to Pint](https://pint.readthedocs.io/en/stable/)
- (for visualization) VMD v1.9.3 [link to VMD](https://www.ks.uiuc.edu/Research/vmd/)

A deprecated version of sugar compatible with ESPResSo v4.1.4 can be found in the branch `sugar_espresso4.1.4`
Note that further developments of sugar will only be developed for ESPResSo v4.2 and thus further support for
that branch is not planned. 

## Contents

- `sugar.py` : source code of Sugar library
- `peptide_simulation_example.py` : example script on how to use Sugar to setup a peptide simulation in ESPResSo
- `sugar_tutorial.ipynb` : jupyter notebook containing a tutorial to Sugar library
- `reference_scripts/`: folder with scripts built to reproduce reference data to ensure the reliability of Sugar
- `handy_scripts/`: folder with various handy scripts complementary with Sugar, contains:
	- `vmd-traj.py` : trajectory processing script. Since trajectories with a varying number of particles cannot be readen by VMD, this script
                processes the output trajectory .vtf files from espresso to make them readable by VMD.
	- `zeta_to_charge.py` : Estimates the charge density of a colloidal object given its zeta potential 

## Usage

### Use Sugar in your simulation scripts 

To use Sugar in your simulations, first clone this repository into your source folder

`git clone git@gitlab.com:blancoapa/sugar_library.git`

then you can load Sugar into your script with the command

`from  sugar_library import sugar`

Please, be aware that Sugar is intended to be a supporting tool to setup simulations with ESPResSo. Thus, for most of its funtionalities ESPResSo must be also loaded to your script

`import espressomd`

and your simulations should be runned using ESPResSo

`{$ESPResSo_build_path}/pypresso your_simulation_script.py`

### Example script of the setting up of a peptide simulation with Sugar

You can run the example script directly using your ESPResSo instalation

`{$ESPResSo_build_path}/pypresso your_simulation_script.py`

or alternatively you can modify the $ESPResSo_build_path variable in `Makefile` to match the path where you have build ESPResSo and run the script with

`make run`

once the simulation finishes, the trajectory can be visualized with the command 

`make visual`

### Run Sugar tutorial

You can run the interactive tutorial of Sugar with the command

`{$ESPResSo_build_path}/ipypresso notebook  sugar_tutorial.ipynb`

or alternatively you can run the command

`make tutorial`

provided that you have modified the $ESPResSo_build_path variable in `Makefile` to match the path where you have build ESPResSo

### Run reference script

You can run the reference script `reference_scripts/histatin5_peptide.py` with the command

`make reference`

or directly with

`{$ESPResSo_build_path}/pypresso reference_scripts/histatin5_peptide.py`

## Contribute

All members of this repository are cordially invited to contribute on the development of this library by providing new functionalities 
or helping fixing any possible bugs. Please, note that any change on Sugar source do should be done in a separate branch of the repository, 
which will be reviewed before being merched to the repository master branch



