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

# Find path to pyMBE
current_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
path_end_index=current_dir.find("pyMBE")
pyMBE_path=current_dir[0:path_end_index]+"pyMBE"
sys.path.insert(0, pyMBE_path)

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# Load some functions from the handy_scripts library for convinience
from handy_scripts.handy_functions import setup_electrostatic_interactions_in_espresso
from handy_scripts.handy_functions import minimize_espresso_system_energy
from handy_scripts.handy_functions import setup_langevin_dynamics_in_espresso
from handy_scripts.handy_functions import block_analyze

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_range = np.linspace(2, 12, num=21)
Samples_per_pH = 1000
MD_steps_per_sample = 1000
steps_eq =int(Samples_per_pH/3)
N_samples_print = 10 # Write the trajectory every 100 samples
probability_reaction = 0.5 
SEED = 1
dt = 0.01
solvent_permitivity = 78.3 

L = 25.513*pmb.units.nm

# Peptide parameters
N_aminoacids = 5
sequence = 'K'*N_aminoacids+'D'*N_aminoacids
model = '2beadAA'  # Model with 2 beads per each aminoacid
pep_concentration = 1e-4 *pmb.units.mol/pmb.units.L

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt = 1e-2 * pmb.units.mol/ pmb.units.L

# Define salt parameters

pmb.define_particle( name=cation_name,  q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle( name=anion_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

pmb.load_interaction_parameters (filename=pyMBE_path+'/reference_parameters/interaction_parameters/Lunkad2021.txt') 
pmb.load_pka_set (filename=pyMBE_path+'/reference_parameters/pka_sets/CRC1991.txt')

# Define the peptide on the pyMBE dataframe 
pmb.define_peptide( name=sequence, sequence=sequence, model=model)

# System parameters
volume = L**3
N_peptide_chains = int ( volume * pmb.N_A * pep_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system = espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso (espresso_system=espresso_system)

# Create your molecules into the espresso system

pmb.create_pmb_object_in_espresso (name=sequence, number_of_objects= N_peptide_chains,espresso_system=espresso_system, use_default_bond=True)


# Create counterions for the peptide chains
pmb.create_counterions_in_espresso (object_name=sequence,cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system) 
c_salt_calculated = pmb.create_added_salt_in_espresso (espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)


#List of ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide are ', list_ionisible_groups)

# Setup the acid-base reactions of the peptide using the constant pH ensemble
RE, sucessfull_reactions_labels = pmb.setup_constantpH_reactions_in_espresso(counter_ion=cation_name, constant_pH=2, SEED = SEED)
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

# Setup espresso to track the each type defined in type_map
type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy
pmb.setup_lj_interactions_in_espresso (espresso_system=espresso_system)
minimize_espresso_system_energy (espresso_system=espresso_system)
setup_electrostatic_interactions_in_espresso(units=pmb.units,
                                            espresso_system=espresso_system,
                                            kT=pmb.kT)
minimize_espresso_system_energy (espresso_system=espresso_system)


setup_langevin_dynamics_in_espresso (espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = SEED,
                                    time_step=dt,
                                    tune_skin=False)

espresso_system.cell_system.skin=0.4

# Save the initial state
with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

N_frame=0
Z_pH=[] # List of the average global charge at each pH
Rg_pH=[] 

particle_id_list = pmb.get_particle_id_map(object_name=sequence)["all"]
first_peptide_id = min(particle_id_list)

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
            charge_dict=pmb.calculate_net_charge_in_molecules(espresso_system=espresso_system, 
                                                                object_name=sequence)      
            Z_sim.append(charge_dict["mean"])
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

print("Net charge analysis")
av_charge, err_charge, tau_charge, block_size = block_analyze(input_data=pmb.np.array(Z_pH))

print("Rg analysis")
av_rg, err_rg, tau_rg, block_size = block_analyze(input_data=Rg_pH)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
Z_HH = pmb.calculate_HH(object_name=sequence,
                         pH_list=pH_range)

# Load the reference data 
reference_file_Path = pyMBE_path+"/reference_data/Lys-AspMSDE.csv"
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
