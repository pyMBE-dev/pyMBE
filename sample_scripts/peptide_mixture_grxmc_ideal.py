#Load espresso, sugar and other necessary libraries
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

# Find path to pyMBE
current_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
path_end_index=current_dir.find("pyMBE")
pyMBE_path=current_dir[0:path_end_index]+"pyMBE"
sys.path.insert(0, pyMBE_path)

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_range = np.linspace(2, 12, num=20)
Samples_per_pH = 500
MD_steps_per_sample = 0
steps_eq = int(Samples_per_pH)
N_samples_print = 1000  # Write the trajectory every 100 samples
probability_reaction =1
SEED = 100
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
pmb.load_interaction_parameters (filename=pyMBE_path+'/reference_parameters/interaction_parameters/Lunkad2021.txt') 
pmb.load_pka_set (filename=pyMBE_path+'/reference_parameters/pka_sets/Hass2015.txt')

# Use a generic parametrization for the aminoacids not parametrized
not_parametrized_neutral_aminoacids = ['A','N','Q','G','I','L','M','F','P','O','S','U','T','W','V','J']
not_parametrized_acidic_aminoacids = ['C','c']
not_parametrized_basic_aminoacids = ['R','n']

already_defined_AA=[]

for aminoacid_key in sequence1:
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

for aminoacid_key in sequence2:
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

# Defines the peptides in the pyMBE data frame
peptide1 = 'generic_peptide1'
pmb.define_peptide (name=peptide1, sequence=sequence1, model=model)
peptide2 = 'generic_peptide2'
pmb.define_peptide (name=peptide2, sequence=sequence2, model=model)

# Solution parameters
proton_name = 'Hplus'
hydroxide_name = 'OHminus'
sodium_name = 'Na'
chloride_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=proton_name, q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=hydroxide_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=sodium_name, q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=chloride_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))


# System parameters
volume = N_peptide1_chains/(pmb.N_A*pep1_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide1_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_pmb_object_in_espresso(name=peptide1, number_of_objects= N_peptide1_chains,espresso_system=espresso_system, use_default_bond=True)
pmb.create_pmb_object_in_espresso(name=peptide2, number_of_objects= N_peptide2_chains,espresso_system=espresso_system, use_default_bond=True)
pmb.create_counterions_in_espresso(object_name=peptide1,cation_name=proton_name,anion_name=hydroxide_name,espresso_system=espresso_system) # Create counterions for the peptide chains
pmb.create_counterions_in_espresso(object_name=peptide2,cation_name=proton_name,anion_name=hydroxide_name,espresso_system=espresso_system) # Create counterions for the peptide chains

c_salt_calculated = pmb.create_added_salt_in_espresso(espresso_system=espresso_system,cation_name=sodium_name,anion_name=chloride_name,c_salt=c_salt)

with open('frames/trajectory0.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

#List of ionisable groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))

RE, sucessful_reactions_labels, ionic_strength_res = pmb.setup_grxmc_reactions_in_espresso(pH_res=2, c_salt_res=c_salt, proton_name=proton_name, hydroxide_name=hydroxide_name, sodium_name=sodium_name, chloride_name=chloride_name, SEED=SEED)
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
print('\nMinimazing system energy\n')
espresso_system.cell_system.skin = 0.4
espresso_system.time_step = dt 
print('steepest descent')
espresso_system.integrator.set_steepest_descent(f_max=0, gamma=0.1, max_displacement=0.1)
espresso_system.integrator.run(1000)
print('velocity verlet')
espresso_system.integrator.set_vv()  # to switch back to velocity Verlet
espresso_system.integrator.run(1000)
print('\nMinimization finished \n')

espresso_system.time_step = dt

#Save the initial state
with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

# Setup espresso to do langevin dynamics
print (f'Optimizing skin\n')

espresso_system.time_step= dt 
espresso_system.integrator.set_vv()
espresso_system.thermostat.set_langevin(kT=pmb.kT.to('reduced_energy').magnitude, gamma=0.1, seed=SEED)

N_frame=0
Z_pH=[] # List of the average global charge at each pH
xi_plus=[] # List of the average partition coefficient of positive ions

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

#Save the pyMBE dataframe in a CSV file
pmb.df.to_csv('df.csv')

# Main loop for performing simulations at different pH-values

for index in tqdm(range(len(pH_range))):
    
    pH_value=pH_range[index]
    # Sample list inicialization
    Z_sim=[]
    num_plus=[]

    RE, sucessful_reactions_labels, ionic_strength_res = pmb.setup_grxmc_reactions_in_espresso(pH_res=pH_value, c_salt_res=c_salt, proton_name=proton_name, hydroxide_name=hydroxide_name, sodium_name=sodium_name, chloride_name=chloride_name, SEED=SEED)

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
            num_plus.append(espresso_system.number_of_particles(type_map["Na"])+espresso_system.number_of_particles(type_map["Hplus"]))
            
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

def block_analyze (input_data,n_blocks=16):

    data = np.asarray(input_data)
    block = 0
    # this number of blocks is recommended by Janke as a reasonable compromise
    # between the conflicting requirements on block size and number of blocks
    block_size = int(data.shape[1] / n_blocks)
    print(f"block_size: {block_size}")
    # initialize the array of per-block averages
    block_average = np.zeros((n_blocks, data.shape[0]))
    # calculate averages per each block
    for block in range(n_blocks):
        block_average[block] = np.average(data[:, block * block_size: (block + 1) * block_size], axis=1)
    # calculate the average and average of the square
    av_data = np.average(data, axis=1)
    av2_data = np.average(data * data, axis=1)
    # calculate the variance of the block averages
    block_var = np.var(block_average, axis=0)
    # calculate standard error of the mean
    err_data = np.sqrt(block_var / (n_blocks - 1))
    # estimate autocorrelation time using the formula given by Janke
    # this assumes that the errors have been correctly estimated
    tau_data = np.zeros(av_data.shape)

    for val in range(av_data.shape[0]):
        if av_data[val] == 0:
            # unphysical value marks a failure to compute tau
            tau_data[val] = -1.0
        else:
            tau_data[val] = 0.5 * block_size * n_blocks / (n_blocks - 1) * block_var[val] \
                / (av2_data[val] - av_data[val] * av_data[val])

    # check if the blocks contain enough data for reliable error estimates
    print("uncorrelated samples per block:\nblock_size/tau = ", block_size/tau_data)
    threshold = 10.  # block size should be much greater than the correlation time
    if np.any(block_size / tau_data < threshold):
        print("\nWarning: some blocks may contain less than ", threshold, "uncorrelated samples."
        "\nYour error estimated may be unreliable."
        "\nPlease, check them using a more sophisticated method or run a longer simulation.")
        print("? block_size/tau > threshold ? :", block_size/tau_data > threshold)
    else:
        print("\nAll blocks seem to contain more than ", threshold, "uncorrelated samples.\
        Error estimates should be OK.")
    return av_data, err_data,tau_data,block_size

print("Net charge analysis")
av_net_charge, err_net_charge, tau_net_charge, block_size_net_charge = block_analyze(input_data=Z_pH)
av_xi_plus, err_xi_plus, tau_xi_plus, block_size_xi_plus = block_analyze(input_data=xi_plus)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
pH_range_HH = np.linspace(2, 12, num=100)
Z_HH1 = pmb.calculate_HH(object_name=peptide1, pH_list=pH_range_HH)
Z_HH2 = pmb.calculate_HH(object_name=peptide2, pH_list=pH_range_HH)

HH_charge_dict = pmb.calculate_HH_Donnan(espresso_system=espresso_system, object_names=[peptide1, peptide2], c_salt=c_salt, pH_list=pH_range_HH)
Z_HH_Donnan = HH_charge_dict["charges_dict"]
pH_sys = HH_charge_dict["pH_system_list"]
xi = HH_charge_dict["partition_coefficients"]

fig, ax = plt.subplots(figsize=(10, 7))
plt.errorbar(pH_range, np.asarray(av_net_charge)/N_peptide1_chains, yerr=err_net_charge/N_peptide1_chains, fmt = 'o', capsize=3, label='Simulation')
plt.errorbar(np.asarray(pH_range)-np.log10(np.asarray(av_xi_plus)), np.asarray(av_net_charge)/N_peptide1_chains, yerr=err_net_charge/N_peptide1_chains, fmt = 'x', capsize=3, label='Simulation (corrected)')
ax.plot(pH_range_HH, np.asarray(Z_HH1)+np.asarray(Z_HH2), "-k", label='Henderson-Hasselbalch')
ax.plot(pH_range_HH, np.asarray(Z_HH_Donnan[peptide1])+np.asarray(Z_HH_Donnan[peptide2]), "--r", label='HH+Donnan')
ax.plot(pH_sys, np.asarray(Z_HH_Donnan[peptide1])+np.asarray(Z_HH_Donnan[peptide2]), "-.y", label='HH+Donnan (corrected)')
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
