# pyMBE: the Python-based Molecule Builder for ESPResSo 

pyMBE provides tools to facilitate building up molecules with complex architectures in the Molecular Dynamics software [ESPResSo](https://espressomd.org/wordpress/). Some examples of molecules that can be set up with pyMBE are polyelectrolytes, peptides and proteins. pyMBE bookkeeps all the information about the molecule topology, permitting to link each particle to its corresponding residue and molecule. pyMBE uses the [Pint](https://pint.readthedocs.io/en/stable/) library to enable input parameters in any arbitrary unit system, which is later transformed in the reduced unit system used in ESPResSo.

## Dependencies

- [ESPResSo](https://espressomd.org/wordpress/) v4.2.1 
- [Pint](https://pint.readthedocs.io/en/stable/) v0.20.01 
- [Pandas](https://pandas.pydata.org/) v1.5.3
- [Pint-Pandas](https://pypi.org/project/Pint-Pandas/) v0.5
- [Numpy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [pdoc](https://pdoc.dev/) (for building the docs)

## Branches

A deprecated version of pyMBE compatible with ESPResSo v4.1.4 (under the historical name of pyMBE, Sugar)  can be found in the branch `sugar_espresso4.1.4`. Note that further development of pyMBE will only be carried out for ESPResSo v4.2.1 and its forthcoming releases, and no further support for that branch is planned.

## Contents

- `docs/`: folder with the API documentation of pyMBE.
- `figs/`: folder with various images used in the tutorials of pyMBE.
- `handy_scripts/`: folder with various handy scripts and libraries.
- `logo/`: folder with the logo of pyMBE.
- `reference_data/`: folder with various reference data set used to validate pyMBE.
- `reference_parameters/`: folder with various sets of parameters from previous works.
- `sample_scripts/`: folder with various sample scripts showcasing how to use pyMBE to setup different systems.
- `tests/`: folder with several test scripts to check that new developments do not break pyMBE.
- `tutorials/`: folder with the available tutorials on pyMBE.
- `visualization/`: folder with helper scripts to aid the visualization of vtf trajectories from constant pH and Grand reaction simulations with [VMD](https://www.ks.uiuc.edu/Research/vmd/).
- `AUTHORS.md`: list of authors and contributors of pyMBE.
- `CONTRIBUTING`: rules on how to contribute to pyMBE.
- `LICENSE.md`: license of pyMBE.
- `pyMBE.py`: source code of pyMBE

## Usage

### Use pyMBE in your simulation scripts 

To use pyMBE in your simulations, first clone this repository into your source folder

```sh
git clone git@gitlab.com:blancoapa/pyMBE.git
```

then you can load pyMBE into your script with the command

```py
from pyMBE import pyMBE
```

Please, be aware that pyMBE is intended to be a supporting tool to setup simulations with ESPResSo. Thus, for most of its functionalities ESPResSo must be also loaded to your script

```py
import espressomd
```

and your simulations should be runned using ESPResSo

```sh
${ESPResSo_build_path}/pypresso your_simulation_script.py
```

### Run the tutorial of pyMBE

You can run the interactive tutorial of pyMBE with the command

```sh
${ESPResSo_build_path}/ipypresso notebook pyMBE_tutorial.ipynb
```

or alternatively you can run the command

```sh
make tutorial
```

provided that you have modified the `$ESPResSo_build_path` variable in `Makefile` to match the path where you have built ESPResSo.

### Run the testsuite

To make sure your code is valid, please run the testsuite before submitting your contribution:

```sh
PYTHONPATH=$(realpath .) make testsuite
```

When contributing new features, consider adding a unit test in the `testsuite/`
folder and a corresponding line in the `testsuite` target of the Makefile.

Every contribution is automatically tested in CI using EESSI (https://www.eessi.io)
and the [EESSI GitHub Action](https://github.com/marketplace/actions/eessi).
