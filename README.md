# pyMBE: the python Molecular Brewer for ESPResSo 

pyMBE provides tools to facilitate building up molecules with complex architectures in the Molecular Dynamics software [ESPResSo](https://espressomd.org/wordpress/). Some examples of molecules that can be set up with pyMBE are polyelectrolytes, peptides and proteins. pyMBE bookkeeps all the information about the molecule topology, permiting to link each particle to its corresponding residue and molecule. pyMBE uses the [Pint](https://pint.readthedocs.io/en/stable/) library to permit defining the input parameters in any arbitrary unit system, which is later transformed in the reduced unit system used in  ESPResSo.

## Dependencies

- [ESPResSo](https://espressomd.org/wordpress/) v4.2.1 
- [Pint](https://pint.readthedocs.io/en/stable/) v0.20.01 
- [Pandas](https://pandas.pydata.org/) v1.5.3
- [Numpy](https://numpy.org/)
- [pdoc](https://pdoc.dev/) (for building the docs)

A deprecated version of pyMBE compatible with ESPResSo (under the historical name of pyMBE, Sugar)  can be found in the branch `sugar_espresso4.1.4`. Note that further developments of pyMBE will only be developed for ESPResSo v4.2.1 and its forthcoming realises and no further support for that branch is not planned. 

## Contents

- `pyMBE.py` : source code of pyMBE
- `samples` : folder with various sample scripts showcasing how to use pyMBE to setup different systems.
- `pyMBE_tutorial.ipynb` : Tutorial of pyMBE
- `tests/`: folder with several test scripts to check that new developtments do not break pyMBE
- `handy_scripts/`: folder with various handy scripts and libraries

## Usage

### Use pyMBE in your simulation scripts 

To use pyMBE in your simulations, first clone this repository into your source folder

`git clone git@gitlab.com:blancoapa/sugar_library.git`

then you can load pyMBE into your script with the command

`from  pyMBE import pyMBE`

Please, be aware that pyMBE is intended to be a supporting tool to setup simulations with ESPResSo. Thus, for most of its funtionalities ESPResSo must be also loaded to your script

`import espressomd`

and your simulations should be runned using ESPResSo

`{$ESPResSo_build_path}/pypresso your_simulation_script.py`

### Run the tutorial of pyMBE

You can run the interactive tutorial of pyMBE with the command

`{$ESPResSo_build_path}/ipypresso notebook  pyMBE_tutorial.ipynb`

or alternatively you can run the command

`make tutorial`

provided that you have modified the $ESPResSo_build_path variable in `Makefile` to match the path where you have build ESPResSo.

