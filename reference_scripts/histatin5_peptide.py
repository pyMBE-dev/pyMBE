"""
Script that simulates a dilute solution of histatin-5 peptide at very low salt concentration.
It calculates the peptide average charge and radius of gyration (using a constant pH simulation) 
All parameters are taken from Ref. 1, whose results are here reproduced as reference.
This script is part of Sugar library and it is meant to serve as reference for it. 

Authors: Dr. Pablo M. Blanco (Charles University) and MsC. Albert Martinez (Royal College of Surgeons in Ireland)

[1] Blanco, P. M., Madurga, S., Garcés, J. L., Mas, F., & Dias, R. S. 
 Influence of macromolecular crowding on the charge regulation of intrinsically disordered proteins. 
 Soft Matter, 17(3), 655-669,2021.

"""

# Load espresso, sugar and other necessary libraries

import espressomd

import os
import sys
import inspect
import numpy as np
# For loading sugar from parent folder

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import sugar

import matplotlib.pyplot as plt
import espressomd.reaction_ensemble
from espressomd import interactions

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
from espressomd.io.writer import vtf
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Create an instance of sugar library

sg=sugar.sugar_library()

# Simulation parameters
sg.set_reduced_units(unit_length=0.4*sg.units.nm)
pH_range = np.linspace(2, 12, num=21)
Samples_per_pH= 1000
MD_steps_per_sample=1000
steps_eq=int(Samples_per_pH/3)
N_samples_print= 100  # Write the trajectory every 100 samples
probability_reaction=0.5 

# Peptide parameters
sequence='NH2-ASP-SER-HIS-ALA-LYS-ARG-HIS-HIS-GLY-TYR-LYS-ARG-LYS-PHE-HIS-GLU-LYS-HIS-HIS-SER-HIS-ARG-GLY-TYR-COOH'
model='1beadAA'
pep_concentration=1.56e-4 *sg.units.mol/sg.units.L
bead_size=0.4*sg.units.nm


# Solution parameters 

c_salt=1e-3 * sg.units.mol/ sg.units.L
cation=sg.particle(name='Na', type=sg.propose_unused_type(), q=1, diameter=0.2*sg.units.nm, epsilon=1*sg.units('reduced_energy'))
anion=sg.particle(name='Cl', type=sg.propose_unused_type(), q=-1, diameter=0.36*sg.units.nm,  epsilon=1*sg.units('reduced_energy'))

    # Load pKa set
sg.load_pka_set(filename='reference_parameters/pka_sets/Nozaki1967.txt')
    # use generic paramenters for the peptide

acidic_aminoacids=['c','E','D','Y','C']
basic_aminoacids=['R','n','K','H']
N_aminoacids=len(sg.protein_sequence_parser(sequence=sequence))

for aminoacid_key in sg.protein_sequence_parser(sequence=sequence):
    if aminoacid_key not in sg.stored_sugar_objects['particle'].keys():
        if aminoacid_key in acidic_aminoacids:
            sg.particle(name=aminoacid_key, acidity='acidic', diameter=bead_size, epsilon=1*sg.units('reduced_energy'))
        elif aminoacid_key in basic_aminoacids:
            sg.particle(name=aminoacid_key, acidity='basic', diameter=bead_size, epsilon=1*sg.units('reduced_energy'))
        else:
            sg.particle(name=aminoacid_key, q=0, diameter=bead_size, epsilon=1*sg.units('reduced_energy'))

generic_bond_lenght=0.4 * sg.units.nm
generic_harmonic_constant=0.41 * sg.units.N / sg.units.m
generic_bond = interactions.HarmonicBond(k=generic_harmonic_constant.to('reduced_energy / reduced_length**2').magnitude,
                                 r_0=generic_bond_lenght.to('reduced_length').magnitude)
sg.define_default_bond(bond=generic_bond)

    # Create an instance of a sugar molecule object for the peptide

histidin5 = sg.peptide(name='histidin5', sequence=sequence, model=model)

# System parameters from Ref. 1

L=22 * sg.units.nm # Side of the simulation box

    # Create an instance of an espresso system

espresso_system=espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system

sg.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system

N_peptide_chains=int(L**3*pep_concentration*sg.N_A)

for _ in range(N_peptide_chains):
    sg.create_sugar_object_in_espresso(sugar_object=histidin5, espresso_system=espresso_system, use_default_bond=True)

dict_titrable_groups=sg.count_titrable_particles(sugar_object=histidin5)
total_ionisible_groups=sum(dict_titrable_groups.values())
sg.create_counterions_in_espresso(sugar_object=histidin5,cation=cation,anion=anion,espresso_system=espresso_system) # Create counterions for the peptide chains
c_salt_calculated=sg.create_added_salt_in_espresso(espresso_system=espresso_system,cation=cation,anion=anion,c_salt=c_salt)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', pep_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide are ', dict_titrable_groups)

# Setup the acid-base reactions of the peptide using the constant pH ensemble

RE=sg.setup_constantpH_reactions_in_espresso(counter_ion=cation)

# Setup espresso to track the ionization of the acid/basic groups in peptide

sg.setup_espresso_to_track_ionization(espresso_system=espresso_system)

# Setup the non-interacting type for speeding up the sampling of the reactions

type_dict=sg.get_all_stored_types()
non_interacting_type=max(type_dict.keys())+1
RE.set_non_interacting_type(non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

sg.setup_lj_interactions_in_espresso(espresso_system=espresso_system)
sg.setup_electrostatic_interactions_in_espresso(espresso_system=espresso_system, c_salt=c_salt)

# Minimize the system energy to avoid huge starting force due to random inicialization of the system

sg.minimize_espresso_system_energy(espresso_system=espresso_system)

# Setup espresso to do langevin dynamics

sg.setup_langevin_dynamics_in_espresso(espresso_system=espresso_system)

# Write the initial state

with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates) 

Rg_pH=[]
Z_pH=[]
N_frame=0
first_peptide_id=sg.get_ids_from_sugar(sugar_object=histidin5)[0][0]

# Main loop for performing simulations at different pH-values

for pH_value in pH_range:

    Z_sim=[]
    Rg_sim=[]
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if sg.np.random.random() > probability_reaction:

            espresso_system.integrator.run(steps=MD_steps_per_sample)
        
        else:
        
            RE.reaction(total_ionisible_groups)

        if ( step > steps_eq):

            Z_net_list = sg.get_net_charge_from_espresso(espresso_system=espresso_system, sugar_object=histidin5)
            Z_sim.append(np.mean(np.array(Z_net_list)))
            Rg = espresso_system.analysis.calc_rg(chain_start=first_peptide_id, number_of_chains=N_peptide_chains, chain_length=N_aminoacids)
            Rg_value = sg.units.Quantity(Rg[0], 'reduced_length')
            Rg_nm = Rg_value.to('nm').magnitude
            Rg_sim.append(Rg_nm)

        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    sg.write_progress(step=list(pH_range).index(pH_value), total_steps=len(pH_range))

    Z_pH.append(Z_sim)
    Rg_pH.append(Rg_sim)
    print("pH = {:6.4g} done".format(pH_value))
    
# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = sg.calculate_HH(sequence=histidin5.sequence, pH=pH_range)

# Estimate the statistical error and the autocorrelation time of the data
print("Charge analysis")
av_charge, err_charge, tau_charge, block_size = sg.block_analyze(input_data=sg.np.array(Z_pH))
print("Rg analysis")
av_rg, err_rg, tau_rg, block_size = sg.block_analyze(input_data=Rg_pH)

# Load the reference data 

reference_file_Path=str(parentdir)+"/reference_data/histatin5_SoftMatter.txt"
reference_data=sg.np.loadtxt(reference_file_Path, delimiter=",")

Z_ref=reference_data[:,1]         
Z_err_ref=reference_data[:,2]

Rg_ref=reference_data[:,4]/10
Rg_err_ref=reference_data[:,5]/10

# Plot the results

plt.figure(figsize=[11, 9])
plt.subplot(1, 2, 1)
plt.suptitle('Peptide sequence: '+''.join(sg.protein_sequence_parser(sequence=sequence)))
plt.errorbar(pH_range, av_charge, yerr=err_charge, fmt = '-o', capsize=3, label='Simulation')
plt.errorbar(pH_range, Z_ref, yerr=Z_err_ref, color='b', fmt = '-o', ecolor = 'b', capsize=3, label='Blanco2021')
plt.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.axhline(y=0.0, color="gray", linestyle="--")
plt.xlabel('pH')
plt.ylabel('Charge of Histatin-5 / e')
plt.legend()

plt.subplot(1, 2, 2)
plt.errorbar(pH_range, av_rg, yerr=err_rg, fmt = '-o', capsize=3, label='Simulation')
plt.errorbar(pH_range, Rg_ref, yerr=Rg_err_ref, color='b', fmt = '-o', ecolor = 'b', capsize=3, label='Blanco2021')
plt.xlabel('pH')
plt.ylabel('Radius of gyration of Histatin-5 / nm')
plt.legend()

plt.show()
