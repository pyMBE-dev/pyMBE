# Sugar library

Sugar is a library that builds on top of the Molecular Dynamics software ESPResSo 
providing tools to simplify the setting up of simulations of peptides and polyelectrolytes.

## Requirements

- ESPResSo v4.1 [link to Espresso](https://espressomd.org/wordpress/download/)
- Pint v0.17 [link to Pint](https://pint.readthedocs.io/en/stable/)
- (for visualization) VMD v1.9.3 [link to VMD](https://www.ks.uiuc.edu/Research/vmd/)

## Contents

- sugar.py : source code of Sugar library
- peptide_simulation_example.py : example script on how to use Sugar to setup a peptide simulation in ESPResSo
- sugar_tutorial.ipynb : jupyter notebook containing a tutorial to Sugar library
- vmd-traj.py : trajectory processing script. Since trajectories with a varying number of particles cannot be readen by VMD, this script
                processes the output trajectory .vtf files from espresso to make them readable by VMD.

## Usage

### Use Sugar in your simulation scripts 

To use Sugar in your simulations clone this repository



