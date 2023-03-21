"""
Script that simulates a dilute solution of GlU5-HIS5 peptide 
It calculates the peptide average charge and radius of gyration (using a constant pH simulation) 
All parameters are taken from Ref. 1, whose results are here reproduced as reference.
This script is part of pyMBE library and it is meant to serve as reference for it. 

Authors: Dr. Pablo M. Blanco (Charles University)

[1] Lunkad, R., Murmiliuk, A., Hebbeker, P., Boublík, M., Tošner, Z., Štěpánek, M., & Košovan, P.  
Quantitative prediction of charge regulation in oligopeptides. 
Molecular Systems Design & Engineering, 2021, 6(2), 122-131.

"""
# Load espresso, pyMBE and other necessary libraries
import os
import sys
import inspect
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.style import use
import espressomd
from espressomd import interactions
from espressomd.io.writer import vtf
from espressomd import electrostatics 

# For loading pyMBE from parent folder
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters

pH_range = np.linspace(2, 12, num=21)
Samples_per_pH = 100
MD_steps_per_sample =100
steps_eq =int(Samples_per_pH/3)
N_samples_print = 10 # Write the trajectory every 100 samples
probability_reaction = 0.5 
SEED = 1
dt = 0.01
solvent_permitivity = 78.3 

L = 25.513*pmb.units.nm

# Peptide parameters
N_aminoacids = 5
sequence = "E"*N_aminoacids+"H"*N_aminoacids
model ='2beadAA'  # Model with 2 beads per each aminoacid
pep_concentration = 1e-4 *pmb.units.mol/pmb.units.L
residue_positions = [0,3,5,len(sequence)-1] # Residue positions to calculate its average charge

# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
pmb.load_interaction_parameters (filename='reference_parameters/interaction_parameters/Lunkad2021.txt') 
pmb.load_pka_set (filename='reference_parameters/pka_sets/CRC1991.txt')

# Create a molecule entry on the pyMBE dataframe 
pmb.define_peptide (name=sequence, sequence=sequence, model=model)

# Salt parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=1e-2 * pmb.units.mol/ pmb.units.L

pmb.define_particle( name=cation_name,  q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle( name=anion_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
volume = L**3
N_peptide_chains = int ( volume * pmb.N_A * pep_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system = espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_pmb_object_in_espresso (name=sequence, number_of_objects= N_peptide_chains,espresso_system=espresso_system, use_default_bond=True)

# Create counterions for the peptide chains
pmb.create_counterions_in_espresso(pmb_object='particle',cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system) 
c_salt_calculated=pmb.create_added_salt_in_espresso(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)

#List of ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide are ', list_ionisible_groups)

# Setup the acid-base reactions of the peptide using the constant pH ensemble
RE, sucessfull_reactions_labels=pmb.setup_constantpH_reactions_in_espresso(counter_ion=cation_name, constant_pH=2, SEED = SEED)
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

# Setup espresso to track the each type defined in type_map
type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy
pmb.setup_lj_interactions_in_espresso(espresso_system=espresso_system)

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

#Electrostatic energy setup

print ('Electrostatics setup')
bjerrum_length = pmb.e.to('reduced_charge')**2 / (4 *pmb.units.pi*pmb.units.eps0* solvent_permitivity * pmb.kT.to('reduced_energy'))
coulomb_prefactor = bjerrum_length.to('reduced_length') * pmb.kT.to('reduced_energy')
coulomb = espressomd.electrostatics.P3M ( prefactor = coulomb_prefactor.magnitude, accuracy=1e-3)

print('\nBjerrum length ', bjerrum_length.to('nm'), '=', bjerrum_length.to('reduced_length'))

espresso_system.time_step = dt
espresso_system.actors.add(coulomb)

# save the optimal parameters and add them by hand
p3m_params = coulomb.get_params()
espresso_system.actors.remove(coulomb)

coulomb = espressomd.electrostatics.P3M (
                            prefactor = coulomb_prefactor.magnitude,
                            accuracy = 1e-3,
                            mesh = p3m_params['mesh'],
                            alpha = p3m_params['alpha'] ,
                            cao = p3m_params['cao'],
                            r_cut = p3m_params['r_cut'],
                            tune = False,
                                )

espresso_system.actors.add(coulomb)

print('\nElectrostatics successfully added to the system \n')

# Write the initial state
with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

print (f'Optimizing skin\n')
espresso_system.time_step = dt 
espresso_system.integrator.set_vv()
espresso_system.thermostat.set_langevin(kT=pmb.kT.to('reduced_energy').magnitude, gamma=0.1, seed=SEED)

espresso_system.cell_system.tune_skin ( min_skin = 1, 
                                        max_skin = espresso_system.box_l[0]/2, tol=1e-3, 
                                        int_steps=1000, adjust_max_skin=True)

print('Optimized skin value: ', espresso_system.cell_system.skin, '\n')

N_frame=0
Z_pH=[] # List of the average global charge at each pH
Rg_pH=[] 

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()
first_peptide_id = min(particle_id_list)
print (particle_id_list,'chain_length',len(particle_id_list))

#Save the pyMBE dataframe in a CSV file
pmb.df.to_csv('df.csv')

for index in tqdm(range(len(pH_range))):
    # Sample list inicialization
    pH_value=pH_range[index]
    print (pH_value)
    Z_sim=[]
    Rg_sim=[]
    Z_groups_time_series=[]       
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value
    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:
            espresso_system.integrator.run(steps=MD_steps_per_sample)
        else:
            RE.reaction(reaction_steps=total_ionisible_groups)

        if ( step > steps_eq):
            # Get peptide net charge        
            z_one_object=0
            for pid in particle_id_list:
                z_one_object +=espresso_system.part.by_id(pid).q
            
            Z_sim.append(z_one_object)
            # Get peptide radius of gyration
            Rg = espresso_system.analysis.calc_rg(chain_start=first_peptide_id, number_of_chains=N_peptide_chains, chain_length=len(particle_id_list))
            Rg_value = pmb.units.Quantity(Rg[0], 'reduced_length')
            Rg_nm = Rg_value.to('nm').magnitude
            Rg_sim.append(Rg_nm)

        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    
    Z_pH.append(np.array(Z_sim))
    Rg_pH.append(Rg_sim)

    print("pH = {:6.4g} done".format(pH_value))

#  Function to Estimate the statistical error and the autocorrelation time of the data
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

    return av_data, err_data, tau_data, block_size

print("Net charge analysis")
av_charge, err_charge, tau_charge, block_size = block_analyze(input_data=pmb.np.array(Z_pH))

print("Rg analysis")
av_rg, err_rg, tau_rg, block_size = block_analyze(input_data=Rg_pH)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
Z_HH = pmb.calculate_HH(sequence=sequence, pH=pH_range)

# Load the reference data 

reference_file_Path = str(parentdir)+"/reference_data/Glu-HisMSDE.csv"
reference_data = pd.read_csv(reference_file_Path)

Z_ref = N_aminoacids*-1*reference_data['aaa']+N_aminoacids*reference_data['aab']         
Rg_ref = reference_data['arg']*0.37

# Plot the results

plt.figure(figsize=[11, 9])
plt.subplot(1, 2, 1)
plt.suptitle('Peptide sequence: '+''.join(pmb.protein_sequence_parser(sequence=sequence)))
plt.errorbar(pH_range, av_charge, yerr=err_charge, fmt = '-o', capsize=3, label='Simulation')
plt.plot(pH_range, Z_ref, '-o', color = 'b', label='Lunkad2021')
plt.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.axhline(y=0.0, color="gray", linestyle="--")
plt.xlabel('pH')
plt.ylabel('Net charge / e')
plt.legend()

plt.subplot(1, 2, 2)
plt.errorbar(pH_range, av_rg, yerr=err_rg, fmt = '-o', capsize=3, label='Simulation')
plt.plot(pH_range, Rg_ref, '-o', color='b',  label='Lunkad2021')
plt.xlabel('pH')
plt.ylabel('Radius of gyration / nm')
plt.legend()

plt.show()