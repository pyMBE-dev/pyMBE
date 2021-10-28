# Load espresso, sugar and other necessary libraries

import espressomd
import sugar 
import numpy as np
import matplotlib.pyplot as plt
import os
import espressomd.reaction_ensemble

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
from espressomd.io.writer import vtf
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Create an instance of sugar library

sg=sugar.sugar_library()

# Peptide parameters

sequence="nHHHEEEc"
model='1beadpeptide'
pKa_set='nozaki'
pep_concentration=1.56e-4 *sg.units.mol/sg.units.L

    # Create an instance of a sugar molecule object for the peptide

peptide = sg.molecule(sequence=sequence, model=model,  pKa_set=pKa_set)

# Salt parameters

c_salt=1e-3 * sg.units.mol/ sg.units.L

    # Create an instance of a sugar particle object for the added salt cations and anions

added_salt_cation=sg.particle()
added_salt_cation.q=1
added_salt_cation.type=0

added_salt_anion=sg.particle()
added_salt_anion.q=-1
added_salt_anion.type=1

# System parameters

L=22 * sg.units.nm # Side of the simulation box

    # Create an instance of an espresso system

system=espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Simulation parameters

pH_range = np.linspace(2, 12, num=20)
Samples_per_pH= 100
MD_steps_per_sample=1000
steps_eq=int(Samples_per_pH/3)
N_samples_print= 100  # Write the trajectory every 100 samples
probability_reaction=0.5 

# Add peptides to your simulation box

volume=system.volume()*sg.units('reduced_length**3')
peptide.N=int(volume*pep_concentration*sg.units.N_A)
sg.create_molecule(peptide, system)
calculated_peptide_concentration=peptide.N/(volume*sg.units.N_A)
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', peptide.N, 'molecules')

# Count the number of titrable groups in your peptide

N_titrable_groups=sg.count_titrable_groups(mol=peptide)
print('The number of ionisable groups in your peptide is ', N_titrable_groups)

# Add added salt ions to your simulation box

c_salt_calculated=sg.create_added_salt(system=system, cation=added_salt_cation, anion=added_salt_anion, c_salt=c_salt)

# Add counter-ions to neutralize the peptide charge

positive_counterion, negative_counterion = sg.create_counterions(system=system, mol=peptide)

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
Z_pH=[] # Average charge list

# Main loop for performing simulations at different pH-values

for pH_value in pH_range:

    Z_sim=[]
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:

            system.integrator.run(steps=MD_steps_per_sample)
        
        else:
        
            RE.reaction(N_titrable_groups)

        if ( step > steps_eq):

            Z, Z2 = sg.calculate_molecule_charge(system=system, mol=peptide)
            Z_sim.append(Z)

        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(system, coordinates)
                vtf.writevcf(system, coordinates)

    Z_sim=np.array(Z_sim)
    Z_pH.append(Z_sim.mean())
    print("pH = {:6.4g} done".format(pH_value))
    
# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = sg.calculate_HH(mol=peptide, pH=list(pH_range))

# Plot the results

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(pH_range, Z_pH, "ro", label='Simulation')
ax.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.title('Peptide sequence: '+ sequence)
plt.show()

exit()
