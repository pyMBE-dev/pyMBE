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
from espressomd import electrostatics

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()
# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'

if not os.path.exists('./frames'):
    os.makedirs('./frames')

#Import functions from handy_functions script 
from handy_scripts.handy_functions import minimize_espresso_system_energy
from handy_scripts.handy_functions import block_analyze

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_range = np.linspace(2, 12, num=20)
Samples_per_pH = 500
MD_steps_per_sample = 0
steps_eq = int(Samples_per_pH)
N_samples_print = 1000  # Write the trajectory every 100 samples
probability_reaction =1
SEED = 42
dt = 0.001
solvent_permitivity = 78.3

# Peptide parameters
sequence1 = 'nGEGHc'
model = '1beadAA'  # Model with 2 beads per each aminoacid
pep1_concentration = 1e-2 *pmb.units.mol/pmb.units.L
N_peptide1_chains = 10

sequence2 = 'nGEEHHc'
pep2_concentration = 1e-2 *pmb.units.mol/pmb.units.L
N_peptide2_chains = 10


# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
path_to_interactions=pmb.get_resource("reference_parameters/interaction_parameters/Lunkad2021.txt")
path_to_pka=pmb.get_resource("reference_parameters/pka_sets/Hass2015.txt")
pmb.load_interaction_parameters (filename=path_to_interactions) 
pmb.load_pka_set (path_to_pka)

# Use a generic parametrization for the aminoacids not parametrized

not_parametrized_neutral_aminoacids = ['A','N','Q','G','I','L','M','F','P','O','S','U','T','W','V','J']
not_parametrized_acidic_aminoacids = ['C','c']
not_parametrized_basic_aminoacids = ['R','n']

already_defined_AA=[]

for aminoacid_key in sequence1+sequence2:
    if aminoacid_key in already_defined_AA:
        continue
    if aminoacid_key in not_parametrized_acidic_aminoacids:
        pmb.define_particle(name=aminoacid_key,
                           acidity='acidic',
                           diameter=0.35*pmb.units.nm, 
                           epsilon=1*pmb.units('reduced_energy'))
    elif aminoacid_key in not_parametrized_basic_aminoacids:
        pmb.define_particle(name=aminoacid_key, acidity='basic',diameter=0.35*pmb.units.nm,epsilon=1*pmb.units('reduced_energy'))
        
    elif aminoacid_key in not_parametrized_neutral_aminoacids:
        pmb.define_particle(name=aminoacid_key,
                           q=0,
                           diameter=0.35*pmb.units.nm, 
                           epsilon=1*pmb.units('reduced_energy'))
    already_defined_AA.append(aminoacid_key)

generic_bond_lenght=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond = interactions.HarmonicBond(k=generic_harmonic_constant.to('reduced_energy / reduced_length**2').magnitude,
                                 r_0=generic_bond_lenght.to('reduced_length').magnitude)

pmb.define_default_bond(bond_object = generic_bond, bond_type="harmonic")

# Defines the peptine in the pyMBE data frame
peptide1 = 'generic_peptide1'
pmb.define_peptide (name=peptide1, sequence=sequence1, model=model)

peptide2 = 'generic_peptide2'
pmb.define_peptide (name=peptide2, sequence=sequence2, model=model)


# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
volume = N_peptide1_chains/(pmb.N_A*pep1_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide1_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_pmb_object(name=peptide1, number_of_objects= N_peptide1_chains,espresso_system=espresso_system, use_default_bond=True)
pmb.create_pmb_object(name=peptide2, number_of_objects= N_peptide2_chains,espresso_system=espresso_system, use_default_bond=True)
pmb.create_counterions(object_name=peptide1, cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system) # Create counterions for the peptide chains
pmb.create_counterions(object_name=peptide2, cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system) # Create counterions for the peptide chains

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)

with open('frames/trajectory0.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

#List of ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))

RE, sucessful_reactions_labels, ionic_strength_res = pmb.setup_grxmc_unified(pH_res=2, c_salt_res=c_salt, cation_name=cation_name, anion_name=anion_name, SEED=SEED)
print('The acid-base reaction has been sucessfully setup for ', sucessful_reactions_labels)

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map(type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Minimize the system energy to avoid huge starting force due to random inicialization of the system
minimize_espresso_system_energy(espresso_system=espresso_system)

espresso_system.time_step = dt

#Save the initial state
with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

# Setup espresso to do langevin dynamics

espresso_system.time_step= dt 
espresso_system.integrator.set_vv()
espresso_system.thermostat.set_langevin(kT=pmb.kT.to('reduced_energy').magnitude, gamma=0.1, seed=SEED)

N_frame=0
Z_pH=[] # List of the average global charge at each pH
xi_plus=[] # List of the average partition coefficient of positive ions

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df (filename='df.csv')

# Main loop for performing simulations at different pH-values

for index in tqdm(range(len(pH_range))):
    
    pH_value=pH_range[index]
    # Sample list inicialization
    Z_sim=[]
    num_plus=[]

    RE, sucessful_reactions_labels, ionic_strength_res = pmb.setup_grxmc_unified(pH_res=pH_value, c_salt_res=c_salt, cation_name=cation_name, anion_name=anion_name, SEED=SEED)

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:
            espresso_system.integrator.run(steps=MD_steps_per_sample)        
        else:
            RE.reaction(reaction_steps = total_ionisible_groups)

        if (step > steps_eq):
            # Get peptide net charge      
            z_one_object=0
            for pid in particle_id_list:
                z_one_object +=espresso_system.part.by_id(pid).q

            Z_sim.append(np.mean((z_one_object)))
            num_plus.append(espresso_system.number_of_particles(type=type_map["Na"]))
            
        if (step % N_samples_print == 0) :
            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    Z_pH.append(Z_sim)
    concentration_plus = (num_plus/(pmb.N_A * L**3)).to('mol/L')
    xi_plus.append(concentration_plus/ionic_strength_res)
    print("pH = {:6.4g} done".format(pH_value))
   
# Estimate the statistical error and the autocorrelation time of the data
print("Net charge analysis")
av_net_charge, err_net_charge, tau_net_charge, block_size_net_charge = block_analyze(input_data=Z_pH)
av_xi_plus, err_xi_plus, tau_xi_plus, block_size_xi_plus = block_analyze(input_data=xi_plus)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
pH_range_HH = np.linspace(2, 12, num=100)
HH_charge_dict = pmb.calculate_HH_Donnan(espresso_system=espresso_system, object_names=[peptide1, peptide2], c_salt=c_salt, pH_list=pH_range_HH)
Z_HH_Donnan = HH_charge_dict["charges_dict"]
pH_sys = HH_charge_dict["pH_system_list"]
xi = HH_charge_dict["partition_coefficients"]

fig, ax = plt.subplots(figsize=(10, 7))
plt.errorbar(pH_range, np.asarray(av_net_charge)/N_peptide1_chains, yerr=err_net_charge/N_peptide1_chains, fmt = 'o', capsize=3, label='Simulation')
ax.plot(pH_range_HH, np.asarray(Z_HH_Donnan[peptide1])+np.asarray(Z_HH_Donnan[peptide2]), "--r", label='HH+Donnan')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the peptide 1 + peptide 2 / e')
plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(10, 7))
plt.errorbar(pH_range, av_xi_plus, yerr=err_xi_plus, fmt = 'o', capsize=3, label='Simulation')
ax.plot(pH_range_HH, np.asarray(xi), "-k", label='HH+Donnan')
plt.legend()
plt.xlabel('pH')
plt.ylabel(r'partition coefficient $\xi_+$')
plt.show()
plt.close()
