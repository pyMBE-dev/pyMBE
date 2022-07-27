"""
Script that simulates a dilute solution of GlU5-HIS5 peptide 
It calculates the peptide average charge and radius of gyration (using a constant pH simulation) 
All parameters are taken from Ref. 1, whose results are here reproduced as reference.
This script is part of Sugar library and it is meant to serve as reference for it. 

Authors: Dr. Pablo M. Blanco (Charles University)
[1] Lunkad, R., Murmiliuk, A., Hebbeker, P., Boublík, M., Tošner, Z., Štěpánek, M., & Košovan, P.  
Quantitative prediction of charge regulation in oligopeptides. 
Molecular Systems Design & Engineering, 2021, 6(2), 122-131.
"""

# Load espresso, sugar and other necessary libraries

from matplotlib.style import use
import espressomd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import inspect
from espressomd import interactions

# For loading sugar from parent folder

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
import sugar 
# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
from espressomd.io.writer import vtf
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Create an instance of sugar library

sg=sugar.sugar_library()
# Simulation parameters

pH_range = np.linspace(2, 12, num=21)

Samples_per_pH= 1000
MD_steps_per_sample=1000
steps_eq=int(Samples_per_pH/3)
N_samples_print= 100  # Write the trajectory every 100 samples
probability_reaction=0.5 
L=25.513*sg.units.nm

# Peptide parameters
N_aminoacids=5
sequence="K"*N_aminoacids+"D"*N_aminoacids
model='2beadAA'  # Model with 2 beads per each aminoacid
pep_concentration=1e-4 *sg.units.mol/sg.units.L
residue_positions=[0,3,5,len(sequence)-1] # Residue positions to calculate its average charge

    # Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

sg.load_interaction_parameters(filename='reference_parameters/interaction_parameters/Lunkad2021.txt') 
sg.load_pka_set(filename='reference_parameters/pka_sets/CRC1991.txt')

    # Create an instance of a sugar molecule object for the peptide

peptide = sg.peptide(name=sequence, sequence=sequence, model=model)

# Salt parameters

c_salt=1e-2 * sg.units.mol/ sg.units.L
cation=sg.particle(name='Na',  q=1, diameter=0.35*sg.units.nm, epsilon=1*sg.units('reduced_energy'))
anion=sg.particle(name='Cl',  q=-1, diameter=0.35*sg.units.nm,  epsilon=1*sg.units('reduced_energy'))

# System parameters

volume=L**3
N_peptide_chains=int(volume*sg.N_A*pep_concentration)
L=volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration=N_peptide_chains/(volume*sg.N_A)
dict_titrable_groups=sg.count_titrable_particles(sugar_object=peptide)
total_ionisible_groups=sum(dict_titrable_groups.values())
print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide is ', dict_titrable_groups)

    # Create an instance of an espresso system

espresso_system=espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system

sg.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system

for _ in range(N_peptide_chains):
    sg.create_sugar_object_in_espresso(sugar_object=peptide, espresso_system=espresso_system, use_default_bond=True)

sg.create_counterions_in_espresso(sugar_object=peptide,cation=cation,anion=anion,espresso_system=espresso_system) # Create counterions for the peptide chains
c_salt_calculated=sg.create_added_salt_in_espresso(espresso_system=espresso_system,cation=cation,anion=anion,c_salt=c_salt)

# Setup the acid-base reactions of the peptide using the constant pH ensemble

RE, sucessfull_reactions_labels=sg.setup_constantpH_reactions_in_espresso(counter_ion=cation, constant_pH=2)
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

# Setup espresso to track the ionization of the acid/basic groups in peptide

sg.setup_espresso_to_track_ionization(espresso_system=espresso_system)

# Setup the non-interacting type for speeding up the sampling of the reactions

type_dict=sg.get_all_stored_types()
non_interacting_type=max(type_dict.keys())+1
RE.set_non_interacting_type(type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

sg.setup_lj_interactions_in_espresso(espresso_system=espresso_system)
sg.minimize_espresso_system_energy(espresso_system=espresso_system)
sg.setup_electrostatic_interactions_in_espresso(espresso_system=espresso_system, c_salt=c_salt)
sg.minimize_espresso_system_energy(espresso_system=espresso_system)

# Setup espresso to do langevin dynamics

sg.setup_langevin_dynamics_in_espresso(espresso_system=espresso_system)

# Write the initial state

with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

N_frame=0
Z_pH=[] # List of the average global charge at each pH
Rg_pH=[] 
first_peptide_id=sg.get_ids_from_sugar(sugar_object=peptide)[0][0]
# Main loop for performing simulations at different pH-values

for pH_value in pH_range:

    # Sample list inicialization

    Z_sim=[]
    Rg_sim=[]
    Z_groups_time_series=[]       

    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:

            espresso_system.integrator.run(steps=MD_steps_per_sample)
        
        else:
        
            RE.reaction(steps=total_ionisible_groups)

        if ( step > steps_eq):

            # Get peptide net charge

            Z_net_list = sg.get_net_charge_from_espresso(espresso_system=espresso_system, sugar_object=peptide)
            Z_sim.append(np.mean(np.array(Z_net_list)))

            Rg = espresso_system.analysis.calc_rg(chain_start=first_peptide_id, number_of_chains=N_peptide_chains, chain_length=len(sequence)*2)
            Rg_value = sg.units.Quantity(Rg[0], 'reduced_length')
            Rg_nm = Rg_value.to('nm').magnitude
            Rg_sim.append(Rg_nm)

        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    sg.write_progress(step=list(pH_range).index(pH_value), total_steps=len(pH_range))
    Z_pH.append(np.array(Z_sim))
    Rg_pH.append(Rg_sim)

    print("pH = {:6.4g} done".format(pH_value))


# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = sg.calculate_HH(sequence=peptide.sequence, pH=pH_range)

# Plot the results

# Estimate the statistical error and the autocorrelation time of the data
print("Charge analysis")
av_charge, err_charge, tau_charge, block_size = sg.block_analyze(input_data=sg.np.array(Z_pH))
print("Rg analysis")
av_rg, err_rg, tau_rg, block_size = sg.block_analyze(input_data=Rg_pH)

# Load the reference data 

reference_file_Path=str(parentdir)+"/reference_data/Lys-AspMSDE.csv"
import pandas
reference_data=pandas.read_csv(reference_file_Path)

Z_ref=N_aminoacids*-1*reference_data['aaa']+N_aminoacids*reference_data['aab']         
Rg_ref=reference_data['arg']*0.37

# Plot the results

plt.figure(figsize=[11, 9])
plt.subplot(1, 2, 1)
plt.suptitle('Peptide sequence: '+''.join(sg.protein_sequence_parser(sequence=sequence)))
plt.errorbar(pH_range, av_charge, yerr=err_charge, fmt = '-o', capsize=3, label='Simulation')
plt.plot(pH_range, Z_ref, '-o', color = 'b', label='Lunkad2021')
plt.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.axhline(y=0.0, color="gray", linestyle="--")
plt.xlabel('pH')
plt.ylabel('Net charge / e')
plt.legend()

plt.subplot(1, 2, 2)
plt.errorbar(pH_range, av_rg, yerr=err_rg, fmt = '-o', capsize=3, label='Simulation')
plt.plot(pH_range, Rg_ref, "-o", color='b', label='Lunkad2021')
plt.xlabel('pH')
plt.ylabel('Radius of gyration / nm')
plt.legend()

plt.show()
