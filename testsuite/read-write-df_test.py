import espressomd
import re 
from ast import literal_eval
from espressomd import interactions
from pandas.testing import assert_frame_equal

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

print ('Read/Write pyMBE dataframe test')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

# path_to_interactions=pmb.get_resource("parameters/peptides/Lunkad2021.txt")
# path_to_pka=pmb.get_resource("parameters/pka_sets/Hass2015.txt")
# pmb.load_interaction_parameters (filename=path_to_interactions) 
# pmb.load_pka_set (path_to_pka)

#Define particles

pmb.define_particle(
    name = "I",
    sigma = 0.3*pmb.units.nm,
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

pmb.define_particle(name=cation_name, q=1, sigma=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, sigma=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

#System parameters
molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
volume = n_molecules/(pmb.N_A*molecule_concentration)
L = volume ** (1./3.) # Side of the simulation box
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

#Setup potential energy

# pmb.setup_lj_interactions (espresso_system=espresso_system)

pmb.pd.options.display.max_colwidth = 10

#Copy the pmb.df into a new DF for the unit test 
stored_df = pmb.df.copy()

#Write the pymbe DF to a csv file 
df_filename = 'df-example_molecule.csv'

pmb.write_pmb_df (filename = df_filename)

# Read the same pyMBE df from a csv a load it in pyMBE
read_df = pmb.read_pmb_df(filename = df_filename)


# Preprocess data for the Unit Test
# The espresso bond object must be converted to a dict in order to compare them using assert_frame_equal

stored_df['parameters_of_the_potential'] = stored_df['parameters_of_the_potential'].apply (lambda x: literal_eval(str(x)) if pmb.pd.notnull(x) else x)
stored_df['bond_object']  = stored_df['bond_object'].apply(lambda x: literal_eval(re.subn('HarmonicBond', '', str(x))[0]) if pmb.pd.notnull(x) else x)

read_df['bond_object']  = read_df['bond_object'].apply(lambda x: literal_eval(re.subn('HarmonicBond', '', str(x))[0]) if pmb.pd.notnull(x) else x)

print(f"*** Unit test: comparison between the written pmb.df and read pmb.df using assert_frame_equal ***")

try:
    assert_frame_equal (stored_df, read_df, check_exact= True)
    print (f"*** Unite test passed***")
except AssertionError as message:
    
    print (f"*** Unite test Failed***")
    print (message)