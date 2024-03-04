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

### Set up the pyMBE environment

To use pyMBE in your simulations, first clone this repository locally:

```sh
git clone git@github.com:pm-blanco/pyMBE.git
```

Please, be aware that pyMBE is intended to be a supporting tool to setup simulations with ESPResSo.
Thus, for most of its functionalities ESPResSo must also be available.

Create a virtual environment to install Python dependencies and configure
the path to the ESPResSo build folder:

```sh
python3 -m venv pymbe
source pymbe/bin/activate
python3 maintainer/configure_venv.py --espresso_path=~/espresso/build # adapt path
python3 -m pip install -r requirements.txt
deactivate
```

We highlight that the path `~/espresso/build` is just an example of a possible
path to the ESPResSo build folder. The user should change this path to match
the local path were ESPResSo was installed.

Cluster users who rely on module files to load dependencies should opt for the
following alternative:

```sh
module load ESPResSo/4.2.1-foss-2022a # adapt module name
python3 -m venv --system-site-packages pymbe
source pymbe/bin/activate
python3 maintainer/configure_venv.py
python3 -m pip install -r requirements.txt
deactivate
module purge
```

We highlight that the module files need to be loaded before every activation
of the virtual environment.

Now you can use pyMBE and ESPResSo by activating the virtual environment:

```sh
$ source pymbe/bin/activate
(pymbe) $ python3 -c "import espressomd.version; print(espressomd.version.friendly())"
4.2
(pymbe) $ python3 -c "import pyMBE; print(pyMBE.__file__)"
/home/user/Documents/pyMBE/pyMBE.py
$ deactivate
```

To use pyMBE in JupyterLab, register the virtual environment in a new kernel:

```sh
source pymbe/bin/activate
python3 -m pip install ipykernel "jupyterlab>=4.0.8" "PyOpenGL>=3.1.5"
python3 -m ipykernel install --user --name=pyMBE
deactivate
```

Please be aware the pyMBE kernel will be registered outside the environment,
typically in your home folder. You can later inspect the list of registered
kernels and delete unwanted ones with the following commands:

```sh
jupyter kernelspec list
jupyter kernelspec uninstall pymbe
```

The JupyterLab main menu will now show a new Python kernel called "pyMBE"
that uses the virtual environment.

### Use pyMBE in your simulation scripts

```sh
source pymbe/bin/activate
python3 sample_scripts/peptide.py
deactivate
```

### Run the tutorial of pyMBE

You can run the interactive tutorial of pyMBE with the command:

```sh
source pymbe/bin/activate
jupyter lab tutorials/pyMBE_tutorial.ipynb
deactivate
```

Be sure to use the pyMBE kernel instead of the default Python3 kernel.
The currently active kernel is usually displayed in the top right corner of the notebook.

### Run the testsuite

To make sure your code is valid, please run the testsuite before submitting your contribution:

```sh
source pymbe/bin/activate
make testsuite
deactivate
```

When contributing new features, consider adding a unit test in the `testsuite/`
folder and a corresponding line in the `testsuite` target of the Makefile.

Every contribution is automatically tested in CI using EESSI (https://www.eessi.io)
and the [EESSI GitHub Action](https://github.com/marketplace/actions/eessi).
