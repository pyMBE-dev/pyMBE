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

import re 
import ast
import tempfile
import espressomd
import pandas as pd

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(SEED=42)

print ('*** Unit tests: read and write from pyMBE dataframe ***')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm,
                      verbose=False)

#Define particles
pmb.define_particle(
    name = "I",
    sigma = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),
    q = 0,
    acidity = "inert")

pmb.define_particle(
    name = "A",
    acidity = "acidic",
    pka = 4,
    sigma = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),)
    
pmb.define_particle(
    name = "B",
    acidity = "basic",
    pka = 9,
    sigma = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),)

#Define the residues
pmb.define_residue(
    name = "Res_1",
    central_bead = "I",
    side_chains = ["A","B"])

pmb.define_residue(
    name = "Res_2",
    central_bead = "I",
    side_chains = ["Res_1"])   


# Defines a molecule  
molecule_name = "A_molecule"
n_molecules = 1

pmb.define_molecule(
    name = molecule_name,
    residue_list = ["Res_1", "Res_1",
                    "Res_2", "Res_1",
                    "Res_1", "Res_2",
                    "Res_2"])


bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant,
                 }


pmb.define_default_bond(bond_type = bond_type, bond_parameters = harmonic_bond)

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, q=1, sigma=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, sigma=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
volume = n_molecules/(pmb.N_A*molecule_concentration)
L = volume ** (1./3.) # Side of the simulation box
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Setup potential energy

pmb.setup_lj_interactions (espresso_system=espresso_system)
pd.options.display.max_colwidth = 10

# Copy the pmb.df into a new DF for the unit test 
stored_df = pmb.df.copy()

with tempfile.TemporaryDirectory() as tmp_directory:
    # Write the pymbe DF to a csv file
    df_filename = f'{tmp_directory}/df-example_molecule.csv'
    pmb.write_pmb_df (filename = df_filename)
    # Read the same pyMBE df from a csv a load it in pyMBE
    read_df = pmb.read_pmb_df(filename = df_filename)

# Preprocess data for the Unit Test
# The espresso bond object must be converted to a dict in order to compare them using assert_frame_equal
stored_df['bond_object']  = stored_df['bond_object'].apply(lambda x: ast.literal_eval(re.subn('HarmonicBond', '', str(x))[0]) if pd.notnull(x) else x)
read_df['bond_object']  = read_df['bond_object'].apply(lambda x: ast.literal_eval(re.subn('HarmonicBond', '', str(x))[0]) if pd.notnull(x) else x)

print("*** Unit test: check that the dataframe stored by pyMBE to file is the same as the one read from the file (same values and variable types) ***")

pd.testing.assert_frame_equal (stored_df, read_df, check_exact= True)
print("*** Unit test passed***")
