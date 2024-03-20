import espressomd
from espressomd import interactions

import os
import sys
import inspect

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

print ('Read/Write pyMBE dataframe test')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

path_to_interactions=pmb.get_resource("parameters/peptides/Lunkad2021.txt")
path_to_pka=pmb.get_resource("parameters/pka_sets/Hass2015.txt")
pmb.load_interaction_parameters (filename=path_to_interactions) 
pmb.load_pka_set (path_to_pka)

#Define particles

pmb.define_particle(
    name = "I",
    diameter = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),
    q = 0,
    acidity = "inert")

# Acidic particle
pmb.define_particle(
    name = "A",
    acidity = "acidic",
    pka = 4)
    
# Basic particle
pmb.define_particle(
    name = "B",
    acidity = "basic",
    pka = 9)

#Defines the residues

pmb.define_residue(
    name = "Res_1",
    central_bead = "I",
    side_chains = ["A","B"])

pmb.define_residue(
    name = "Res_2",
    central_bead = "I",
    side_chains = ["Res_1"])   


# Defines the example molecule of the paper 
molecule_name = "A_molecule"
n_molecules = 1

pmb.define_molecule(
    name = molecule_name,
    residue_list = ["Res_1", "Res_1",
                    "Res_2", "Res_1",
                    "Res_1", "Res_2",
                    "Res_2"])


generic_bond_lenght=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond = interactions.HarmonicBond(k=generic_harmonic_constant.to('reduced_energy / reduced_length**2').magnitude,
                                 r_0=generic_bond_lenght.to('reduced_length').magnitude)

pmb.define_default_bond(bond_object = generic_bond, bond_type="harmonic")

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

#System parameters
molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
volume = n_molecules/(pmb.N_A*molecule_concentration)
L = volume ** (1./3.) # Side of the simulation box
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

#Setup potential energy

pmb.setup_lj_interactions (espresso_system=espresso_system)

#Write the pymbe DF to a csv file 

df_filename = 'df-example_molecule.csv'

pmb.write_pmb_df (filename=df_filename)

# Read the same pyMBE df from a csv a load it in pyMBE

read_df = pmb.read_pmb_df(df_filename)

#Compare both df

dtype_comparison = pmb.df.dtypes == read_df.dtypes 

print(f"*** Unit test: test that the dtype between the written pmb.df and the file read are equal  ***")

if dtype_comparison.all():
    print(f"*** Unit test passed ***")
else:
    different_dtype_columns = pmb.df.columns[~dtype_comparison]
    for col in different_dtype_columns:
       raise ValueError (f"The following columns have different dtype: {col}: {pmb.df[col].dtype}")

print(f"*** Unit test: test that the rows between the written pmb.df and the file read are equal  ***")

row_comparison = pmb.df.equals(read_df)

if row_comparison:
    print(f"*** Unit test passed ***")
else:
    raise ValueError ('There is some inconsistency between the pmb.df and read_df')