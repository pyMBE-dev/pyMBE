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

# For loading sugar from parent folder

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import sugar

import matplotlib.pyplot as plt
import os
import espressomd.reaction_ensemble

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
from espressomd.io.writer import vtf
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Create an instance of sugar library

sg=sugar.sugar_library()

# Set the reduced length to match the bead size in Ref. 1
sg.set_reduced_units(unit_length=0.4*sg.units.nm)

# Peptide parameters
sequence='NH2-ASP-SER-HIS-ALA-LYS-ARG-HIS-HIS-GLY-TYR-LYS-ARG-LYS-PHE-HIS-GLU-LYS-HIS-HIS-SER-HIS-ARG-GLY-TYR-COOH'

# Model used in Ref. 1

model='1beadpeptide'
pKa_set='nozaki' 
sg.param.default.default_bond.k=0.41 * sg.units.N / sg.units.m

# Create a sugar molecule object for histatin-5 peptide

peptide=sg.molecule(sequence=sequence, model=model, pKa_set=pKa_set)

# Solution parameters from Ref. 1

c_salt=1e-3 * sg.units.mol/ sg.units.L
pep_concentration=1.56e-4 *sg.units.mol/sg.units.L

    # Create an instance of a sugar particle object for the added salt cations and anions

added_salt_cation=sg.particle()
added_salt_cation.q=1
added_salt_cation.type=0
added_salt_cation.radius=0.1*sg.units.nm

added_salt_anion=sg.particle()
added_salt_anion.q=-1
added_salt_anion.type=1
added_salt_anion.radius=0.18*sg.units.nm

# System parameters from Ref. 1

L=22 * sg.units.nm # Side of the simulation box

    # Create an instance of an espresso system

system=espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Simulation parameters

pH_range = sg.np.linspace(2, 12, num=21)
Samples_per_pH= 100
MD_steps_per_sample=1000
steps_eq=int(Samples_per_pH/3)
N_samples_print= 100  # Write the trajectory every 100 samples
probability_reaction=0.5 

# Add peptides to your simulation box

volume=system.volume()*sg.units('reduced_length**3')
peptide.N=int(volume*pep_concentration*sg.N_A)
sg.create_molecule(peptide, system)
calculated_peptide_concentration=peptide.N/(volume*sg.N_A)
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', peptide.N, 'molecules')

# Count the number of titrable groups in your peptide

N_titrable_groups=sg.count_titrable_groups(mol=peptide)
print('The number of ionisable groups in your peptide is ', N_titrable_groups)

# Add added salt ions to your simulation box

c_salt_calculated=sg.create_added_salt(system=system, cation=added_salt_cation, anion=added_salt_anion, c_salt=c_salt)

# Add counter-ions to neutralize the peptide charge

positive_counterion, negative_counterion = sg.create_counterions(system=system, mol=peptide)

# Set the radius of the counter-ions to the reference values

positive_counterion.radius=0.1*sg.units.nm
negative_counterion.radius=0.18*sg.units.nm

# Setup the acid-base reactions of the peptide (in the constant pH ensemble)

RE=sg.setup_acidbase_reactions(mol=peptide, counter_ion=added_salt_cation)

# Setup espresso to track the ionization of the acid/basic groups in peptide

sg.track_ionization(system=system, mol=peptide)

# Setup the non-interacting type for speeding up the sampling of the reactions

non_interacting_type=max(system.part[:].type)+1
RE.set_non_interacting_type(non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

    # Setup Lennard-Jones interactions (By default Weeks-Chandler-Andersen, i.e. only steric repulsion)

mol_list=[peptide,positive_counterion,negative_counterion,added_salt_cation,added_salt_anion] # list of sugar molecules/particles for setting up the LJ interactions

sg.setup_lj_interactions(mol_list=mol_list, system=system)

    # Setup the electrostatic potential (By default, it is used the p3m method but the Debye-huckel potential can also be set up by method='DH')

sg.setup_electrostatic_interactions(system=system, c_salt=c_salt)

# Minimize the system energy to avoid huge starting force due to random inicialization of the system

sg.minimize_system_energy(system=system)

# Setup espresso to do langevin dynamics

sg.setup_langevin_dynamics(system=system)

# Write the initial state

with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(system, coordinates)
    vtf.writevcf(system, coordinates)

N_frame=0

# Lists for the observables

Z_pH=[] 
Rg_pH=[]

# Main loop for performing simulations at different pH-values

for pH_value in pH_range:

    Z_sim=[]
    Rg_sim=[]
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if sg.np.random.random() > probability_reaction:

            system.integrator.run(steps=MD_steps_per_sample)
        
        else:
        
            RE.reaction(N_titrable_groups)

        if ( step > steps_eq):

            Z, Z2 = sg.calculate_molecule_charge(system=system, mol=peptide)
            Z_sim.append(Z)
            Rg = system.analysis.calc_rg(chain_start=0, number_of_chains=1, chain_length=26)
            Rg_value = sg.units.Quantity(Rg[0], 'reduced_length')
            Rg_nm = Rg_value.to('nm').magnitude
            Rg_sim.append(Rg_nm)

        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(system, coordinates)
                vtf.writevcf(system, coordinates)

    sg.write_progress(step=list(pH_range).index(pH_value), total_steps=len(pH_range))

    Z_pH.append(Z_sim)
    Rg_pH.append(Rg_sim)
    print("pH = {:6.4g} done".format(pH_value))
    
# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = sg.calculate_HH(mol=peptide, pH=list(pH_range))

# Estimate the statistical error and the autocorrelation time of the data
print("Charge analysis")
av_charge, err_charge, tau_charge, block_size = sg.block_analyze(input_data=sg.np.array(Z_pH))
print("Rg analysis")
av_rg, err_rg, tau_rg, block_size = sg.block_analyze(input_data=Rg_pH)

# Load the reference data 

reference_file_Path=str(parentdir)+"/reference_scripts/histatin5_reference_data.txt"
reference_data=sg.np.loadtxt(reference_file_Path, delimiter=",")

Z_ref=reference_data[:,1]         
Z_err_ref=reference_data[:,2]

Rg_ref=reference_data[:,4]/10
Rg_err_ref=reference_data[:,5]/10

# Plot the results
peptide_seq = ''.join([str(elem) for elem in peptide.sequence])

plt.figure(figsize=[11, 9])
plt.subplot(1, 2, 1)
plt.suptitle('Peptide sequence: '+ peptide_seq)
plt.errorbar(pH_range, av_charge, yerr=err_charge, fmt = '-o', capsize=3, label='Simulation')
plt.errorbar(pH_range, Z_ref, yerr=Z_err_ref, color='b', fmt = '-o', ecolor = 'b', capsize=3, label='Reference values')
plt.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.axhline(y=0.0, color="gray", linestyle="--")
plt.xlabel('pH')
plt.ylabel('Charge of Histatin-5 / e')
plt.legend()

plt.subplot(1, 2, 2)
plt.errorbar(pH_range, av_rg, yerr=err_rg, fmt = '-o', capsize=3, label='Simulation')
plt.errorbar(pH_range, Rg_ref, yerr=Rg_err_ref, color='b', fmt = '-o', ecolor = 'b', capsize=3, label='Reference values')
plt.xlabel('pH')
plt.ylabel('Radius of gyration of Histatin-5 / nm')
plt.legend()

plt.show()


exit()
