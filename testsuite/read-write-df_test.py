import espressomd
from espressomd import interactions

import os
import sys
import inspect

# Find path to pyMBE
current_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
path_end_index=current_dir.find("pyMBE")
pyMBE_path=current_dir[0:path_end_index]+"pyMBE"
sys.path.insert(0, pyMBE_path) 

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


#Write the pymbe DF to a csv file 

df_filename = 'df-example_molecule.csv'

pmb.write_pmb_df (df=pmb.df, filename=df_filename)

# Read the same pyMBE df from a csv a load it in pyMBE

read_df = pmb.read_pmb_df(df_filename)

#Compare both df



dtype_comparison = pmb.df.dtypes == read_df.dtypes 

print(f"*** Unit test: test that the dtype between the written pmb.df and the file read are equal  ***")

if not dtype_comparison.all():
    print (dtype_comparison)
    raise ValueError ('There are two columns that have different dtype.')
else:
    print(f"*** Unit test passed ***")


print(f"*** Unit test: test that the rows between the written pmb.df and the file read are equal  ***")

row_comparison = pmb.df.equals(read_df)
if not row_comparison:
    raise ValueError ('There is some inconsistency between the pmb.df and read_df')
else:
    print(f"*** Unit test passed ***")