#
# Copyright (C) 2024 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pathlib
import pyMBE
import pandas as pd
import numpy as np

pmb = pyMBE.pymbe_library(seed=42)

print("*** Unit test: check that the different pKa sets are correctly formatted ***")

pymbe_root = pathlib.Path(pyMBE.__file__).parent
pka_root = pymbe_root / "parameters" / "pka_sets"

for path in pka_root.glob("*.json"):
    print(f"Checking {path.stem}")
    path_to_pka = pmb.get_resource(path.relative_to(pymbe_root).as_posix())
    assert pathlib.Path(path_to_pka) == path
    pmb.load_pka_set(path_to_pka,
                     verbose=False)

print("*** Test passed ***")

print("*** Unit test: check that the order to execute load_pka_set() and load_interaction_parameters does not change the resulting parameters in pmb.df ***")

path_to_interactions=pmb.get_resource("parameters/peptides/Lunkad2021.json")
path_to_pka=pmb.get_resource("parameters/pka_sets/Hass2015.json")

# First order of loading parameters
pmb.setup_df() # clear the pmb_df
pmb.load_interaction_parameters (filename=path_to_interactions) 
pmb.load_pka_set (filename=path_to_pka)
df_1 = pmb.df.copy()
df_1 = df_1.sort_values(by="name").reset_index(drop=True)
# Drop espresso types (they depend on the order of loading)
df_1 = df_1.drop(labels=('state_one', 'es_type'), axis=1).drop(labels=('state_two', 'es_type'), axis=1)
# Drop bond_object  (assert_frame_equal does not process it well)
df_1 = df_1.sort_index(axis=1).drop(labels="bond_object", axis=1)

# Second order of loading parameters
pmb.setup_df() # clear the pmb_df
pmb.load_pka_set (filename=path_to_pka)
pmb.load_interaction_parameters (filename=path_to_interactions) 
df_2 = pmb.df.copy()
df_2 = df_2.sort_values(by="name").reset_index(drop=True)
# Drop espresso types (they depend on the order of loading)
df_2 = df_2.drop(labels=('state_one', 'es_type'), axis=1).drop(labels=('state_two', 'es_type'), axis=1)
# Drop bond_object  (assert_frame_equal does not process it well)
df_2 = df_2.sort_index(axis=1).drop(labels="bond_object", axis=1)

pd.testing.assert_frame_equal(df_1,df_2)

print("*** Test passed ***")

print("*** Unit test: check that  load_interaction_parameters loads FENE bonds correctly ***")
pmb.setup_df() # clear the pmb_df
path_to_interactions=pmb.get_resource("testsuite/test_parameters/test_FENE.json")
pmb.load_interaction_parameters (filename=path_to_interactions) 

expected_parameters = {'r_0'    : 0.4*pmb.units.nm,
                        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                        'd_r_max': 0.8 * pmb.units.nm}
reduced_units = {'r_0'    : 'reduced_length',
                     'k'      : 'reduced_energy / reduced_length**2',
                     'd_r_max': 'reduced_length'}
parameters_in_df = pmb.df[pmb.df.pmb_type == "bond"].parameters_of_the_potential.values[0]

for key in expected_parameters.keys():
    np.testing.assert_equal(actual=parameters_in_df[key],
                desired=expected_parameters[key].m_as(reduced_units[key]),
                verbose=True)

print("*** Test passed ***")
print("*** Unit test: check that  load_interaction_parameters loads residue, molecule and peptide objects correctly ***")

pmb.setup_df() # clear the pmb_df
path_to_interactions=pmb.get_resource("testsuite/test_parameters/test_molecules.json")
pmb.load_interaction_parameters (filename=path_to_interactions)

expected_residue_parameters={"central_bead":  "A", "side_chains": ["B","C"] }
expected_molecule_parameters={"residue_list":   ["R1","R1", "R1"]}
expected_peptide_parameters= {"sequence":   ['K', 'K', 'K', 'K', 'K', 'D', 'D', 'D', 'D', 'D'], "model": "1beadAA" }

# Check residue
np.testing.assert_equal(actual=pmb.df[pmb.df.name == "R1"].central_bead.values[0],
                desired=expected_residue_parameters["central_bead"],
                verbose=True)

np.testing.assert_equal(actual=frozenset(pmb.df[pmb.df.name == "R1"].side_chains.values[0]),
                desired=frozenset(expected_residue_parameters["side_chains"]),
                verbose=True)
# Check molecule
np.testing.assert_equal(actual=frozenset(pmb.df[pmb.df.name == "M1"].residue_list.values[0]),
                desired=frozenset(expected_molecule_parameters["residue_list"]),
                verbose=True)
# Check peptide
np.testing.assert_equal(actual=pmb.df[pmb.df.name == "P1"].sequence.values[0],
                desired=expected_peptide_parameters["sequence"],
                verbose=True)
np.testing.assert_equal(actual=frozenset(pmb.df[pmb.df.name == "P1"].model.values[0]),
                desired=frozenset(expected_peptide_parameters["model"]),
                verbose=True)
print("*** Test passed ***")
print("*** Unit test: check that  load_interaction_parameters raises a ValueError if one loads a data set with an unknown pmb_type ***")
pmb.setup_df() # clear the pmb_df
path_to_interactions=pmb.get_resource("testsuite/test_parameters/test_non_valid_object.json")
input_parameters={"filename":path_to_interactions}
np.testing.assert_raises(ValueError, pmb.load_interaction_parameters, **input_parameters)
print("*** Test passed ***")
print("*** Unit test: check that  load_interaction_parameters raises a ValueError if one loads a bond not supported by pyMBE ***")
pmb.setup_df() # clear the pmb_df
path_to_interactions=pmb.get_resource("testsuite/test_parameters/test_non_valid_bond.json")
input_parameters={"filename":path_to_interactions}
np.testing.assert_raises(ValueError, pmb.load_interaction_parameters, **input_parameters)
print("*** Test passed ***")
