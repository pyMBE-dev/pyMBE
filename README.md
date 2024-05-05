# pyMBE: the Python-based Molecule Builder for ESPResSo 

pyMBE provides tools to facilitate building up molecules with complex architectures in the Molecular Dynamics software [ESPResSo](https://espressomd.org/wordpress/). Some examples of molecules that can be set up with pyMBE are polyelectrolytes, peptides and proteins. pyMBE bookkeeps all the information about the molecule topology, permitting to link each particle to its corresponding residue and molecule. pyMBE uses the [Pint](https://pint.readthedocs.io/en/stable/) library to enable input parameters in any arbitrary unit system, which is later transformed in the reduced unit system used in ESPResSo.

## Dependencies

- [ESPResSo](https://espressomd.org/wordpress/) =4.2.1 
- [Pint](https://pint.readthedocs.io/en/stable/) >=0.20.01
- [Pandas](https://pandas.pydata.org/) >=1.5.3
- [Pint-Pandas](https://pypi.org/project/Pint-Pandas/) >=0.3
- [Numpy](https://numpy.org/) >=1.23
- [SciPy](https://scipy.org/) 
- [pdoc](https://pdoc.dev/) (for building the docs)

## Branches

A deprecated version of pyMBE compatible with ESPResSo v4.1.4 (under the historical name of pyMBE, Sugar)  can be found in the branch `sugar_espresso4.1.4`. Note that further development of pyMBE will only be carried out for ESPResSo v4.2.1 and its forthcoming releases, and no further support for that branch is planned.

## Contents

- `docs/`: folder with the API documentation of pyMBE.
- `figs/`: folder with various images used in the tutorials of pyMBE.
- `libs/`: folder with various libraries.
- `logo/`: folder with the logo of pyMBE.
- `maintainer/`: folder with various scripts used by the maintainers.
- `parameters/`: folder with various sets of parameters from previous works.
- `samples/`: folder with various sample scripts showcasing how to use pyMBE to setup different systems.
- `testsuite/`: folder with several test scripts and data for continous integration of the library.
- `tutorials/`: folder with the available tutorials on pyMBE.
- `visualization/`: folder with helper scripts to aid the visualization of vtf trajectories from constant pH and Grand reaction simulations with [VMD](https://www.ks.uiuc.edu/Research/vmd/).
- `AUTHORS.md`: list of authors and contributors of pyMBE.
- `CONTRIBUTING`: rules on how to contribute to pyMBE.
- `LICENSE.md`: license of pyMBE.
- `pyMBE.py`: source code of pyMBE
- `requirements.txt`: list of required libraries to use pyMBE.

## Usage

### Set up the pyMBE virtual environment

To use pyMBE in your simulations, first clone this repository locally:

```sh
git clone git@github.com:pyMBE-dev/pyMBE.git
```

Please, be aware that pyMBE is intended to be a supporting tool to setup simulations with ESPResSo.
Thus, for most of its functionalities ESPResSo must also be available. Following the NEP29 guidelines, we recommend the users of pyMBE to use Python3.10+ when using our module.

The pyMBE module uses its own Python virtual enviroment to avoid incompatibility issues when loading its requierements from other libraries. 
The Python module (`venv`)[https://docs.python.org/3/library/venv.html#module-venv] from the Python Standard Library (starting with Python 3.3)  is needed to set up pyMBE. 
If `venv` is not in the Python distribution of the user, the user will need to first install 'venv' before setting up pyMBE.
For Ubuntu users, this can be done as follows:

```sh
sudo apt install python3-venv
```

To set up pyMBE, the users need to install its virtual environment, install its Python dependencies and configure the path to the ESPResSo build folder as follows:

```sh
python3 -m venv pymbe
source pymbe/bin/activate
python3 maintainer/configure_venv.py --espresso_path=/home/user/espresso/build # adapt path
python3 -m pip install -r requirements.txt
deactivate
```

We highlight that the path `/home/user/espresso/build` is just an example of a possible
path to the ESPResSo build folder. The user should change this path to match
the local absolute path were ESPResSo was installed. 

The pyMBE virtual enviroment can be deactivated at any moment:
```sh
deactivate
```

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
jupyter-lab tutorials/pyMBE_tutorial.ipynb
deactivate
```

Be sure to use the pyMBE kernel instead of the default Python3 kernel.
The currently active kernel is usually displayed in the top right corner of the notebook.

### Run the testsuite

To make sure your code is valid, please run the testsuite before submitting your contribution:

```sh
source pymbe/bin/activate
make tests
deactivate
```

When contributing new features, consider adding a unit test in the `testsuite/`
folder and a corresponding line in the `testsuite` target of the Makefile.

Every contribution is automatically tested in CI using EESSI (https://www.eessi.io)
and the [EESSI GitHub Action](https://github.com/marketplace/actions/eessi).

## References

Check out the corresponding [preprint](https://doi.org/10.48550/arXiv.2401.14954) to learn more about pyMBE.
If you use pyMBE in your research, please cite our preprint:

```bibtex
@article{beyer2024pymbe,
  title={pyMBE: the Python-based Molecule Builder for ESPResSo},
  author={Beyer, David and Torres, Paola B and Pineda, Sebastian P and Narambuena, Claudio F and Grad, Jean-No{\"e}l and Ko{\v{s}}ovan, Peter and Blanco, Pablo M},
  journal={arXiv preprint arXiv:2401.14954},
  year={2024},
  doi={10.48550/arXiv.2401.14954}
}
```

Please also make sure to properly cite the original authors if you use the resources provided in the `parameters/` folder.
The relevant references are provided as comments in the corresponding files.
