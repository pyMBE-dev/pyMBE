import espressomd
from espressomd import interactions
from espressomd.io.writer import vtf
from espressomd import electrostatics 

import os
import sys
import inspect
import numpy as np

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

#Define the particles 

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

molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L

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


# System parameters
# volume = n_molecules/(pmb.N_A*molecule_concentration)
# L = volume ** (1./3.) # Side of the simulation box
# # calculated_molecule_concentration = n_molecules/(volume*pmb.N_A)

# # Create an instance of an espresso system
# espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# # Add all bonds to espresso system
# pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# # Create your molecules into the espresso system
# pmb.create_pmb_object(name=molecule_name, number_of_objects= n_molecules,espresso_system=espresso_system, use_default_bond=True)

# pmb.create_counterions(object_name=molecule_name,cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system) # Create counterions for the peptide chains

# c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)


#Write the pymbe DF to a csv file 

df_filename = 'df-example_molecule.csv'

pmb.write_pmb_df (df=pmb.df, filename=df_filename)

# Read the same pyMBE df from a csv a load it in pyMBE

read_df = pmb.read_pmb_df(df_filename)

#Compare both df

#comparison of dtype of each column 

dtype_comparison = pmb.df.dtypes == read_df.dtypes 

if not dtype_comparison.all():
    print (dtype_comparison)
    raise ValueError ('There are two columns that have different dtype.')
else:
    print ('All columns have the same format')

#Comparison between each row 
row_comparison = pmb.df.equals(read_df)
if not row_comparison:
    raise ValueError ('There is some inconsistency between the pmb.df and read_df')
