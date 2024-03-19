# Load espresso, sugar and other necessary libraries
import sys
import os 
import inspect
from matplotlib.style import use
import espressomd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from espressomd.io.writer import vtf
from espressomd import interactions


# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# Load some functions from the handy_scripts library for convinience
from handy_scripts.handy_functions import setup_langevin_dynamics
from handy_scripts.handy_functions import block_analyze

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_range = np.linspace(2, 12, num=20)
Samples_per_pH = 100
MD_steps_per_sample = 0
steps_eq = int(Samples_per_pH/3)
N_samples_print = 10  # Write the trajectory every 100 samples
probability_reaction =0.5 
SEED = 100
dt = 0.001
solvent_permitivity = 78.3

# Define different particles
# Inert particle 
pmb.define_particle(
    name = "I",
    acidity = "inert",
    q = 0,
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))
    
# Acidic particle
pmb.define_particle(
    name = "A",
    acidity = "acidic",
    pka = 4,
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))
    
# Basic particle
pmb.define_particle(
    name = "B",
    acidity = "basic",
    pka = 9,
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

# Define different residues
pmb.define_residue(
    name = "Res_1",
    central_bead = "I",
    side_chains = ["A","B"])
    
pmb.define_residue(
    name = "Res_2",
    central_bead = "I",
    side_chains = ["Res_1"])

# Define the molecule
pmb.define_molecule(
    name = "polyampholyte",
    residue_list = 5*["Res_1"] + 5*["Res_2"])

# Define bonds
generic_bond_length = 0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

generic_bond = interactions.HarmonicBond(
        k=generic_harmonic_constant.to('reduced_energy / reduced_length**2').magnitude,
        r_0=generic_bond_length.to('reduced_length').magnitude)

pmb.define_default_bond(bond_object = generic_bond, 
        bond_type="harmonic")

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, q=1, sigma=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, sigma=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
N_polyampholyte_chains = 5
polyampholyte_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
volume = N_polyampholyte_chains/(pmb.N_A*polyampholyte_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_polyampholyte_concentration = N_polyampholyte_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System(box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_pmb_object(name="polyampholyte", number_of_objects=N_polyampholyte_chains,espresso_system=espresso_system, use_default_bond=True)
pmb.create_counterions(object_name="polyampholyte",cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system)

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)

#List of ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The polyampholyte concentration in your system is ', calculated_polyampholyte_concentration.to('mol/L') , 'with', N_polyampholyte_chains, 'peptides')
print('The ionisable groups in your peptide are ', list_ionisible_groups)

RE, sucessfull_reactions_labels = pmb.setup_cpH(counter_ion=cation_name, constant_pH=2, SEED=SEED)
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

# Setup espresso to track the ionization of the acid/basic groups 
type_map = pmb.get_type_map()
types = list(type_map.values())
espresso_system.setup_type_map(type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

#Setup Langevin
setup_langevin_dynamics(espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = SEED,
                                    time_step=dt,
                                    tune_skin=False)

espresso_system.cell_system.skin=0.4

N_frame=0
Z_pH=[] # List of the average global charge at each pH

#Save the pyMBE dataframe in a CSV file
pmb.df.to_csv('df.csv')

# Main loop for performing simulations at different pH-values

for index in tqdm(range(len(pH_range))):
    
    pH_value=pH_range[index]
    # Sample list inicialization
    Z_sim=[]
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:
            espresso_system.integrator.run(steps=MD_steps_per_sample)        
        else:
            RE.reaction( reaction_steps = total_ionisible_groups)

        if ( step > steps_eq):
            # Get polyampholyte net charge
            charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system, 
                    molecule_name="polyampholyte")      
            Z_sim.append(charge_dict["mean"])

        if (step % N_samples_print == 0) :
            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    Z_pH.append(Z_sim)
    print("pH = {:6.4g} done".format(pH_value))
   
# Estimate the statistical error and the autocorrelation time of the data

print("Net charge analysis")
av_net_charge, err_net_charge, tau_net_charge, block_size_net_charge = block_analyze(input_data=Z_pH)

# Calculate the ideal titration curve of the polyampholyte with Henderson-Hasselbach equation (manually)
pH_range_HH = np.linspace(2, 12, num=1000)
Z_HH_manually = [10 * (1/(1+10**(pH_value-9)) - 1/(1+10**(4-pH_value))) for pH_value in pH_range_HH]

# Calculate the ideal titration curve of the polyampholyte with Henderson-Hasselbach equation (pyMBE)
Z_HH = pmb.calculate_HH(object_name="polyampholyte", 
                        pH_list=pH_range_HH)

fig, ax = plt.subplots(figsize=(10, 7))
plt.errorbar(pH_range, av_net_charge, yerr=err_net_charge, fmt = 'o', capsize=3, label='Simulation')
plt.plot(pH_range_HH, Z_HH_manually, label="Henderson-Hasselbalch (manually)")
ax.plot(pH_range_HH, Z_HH, "-k", label='Henderson-Hasselbach (pyMBE)', linestyle="dashed")
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the polyampholyte / e')

plt.show()
