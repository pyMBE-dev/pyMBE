# Load espresso, sugar and other necessary libraries

from matplotlib.style import use
import espressomd
import sugar 
import numpy as np
import matplotlib.pyplot as plt
import os
import espressomd.reaction_ensemble
from espressomd import interactions

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
from espressomd.io.writer import vtf
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters

pH_range = np.linspace(2, 12, num=20)

Samples_per_pH= 100
MD_steps_per_sample=1000
steps_eq=int(Samples_per_pH/3)
N_samples_print= 100  # Write the trajectory every 100 samples
probability_reaction=0.5 

# Create an instance of sugar library

sg=sugar.sugar_library()

# Peptide parameters

sequence="nHHHEEEc"
model='2beadAA'  # Model with 2 beads per each aminoacid
pep_concentration=5.56e-4 *sg.units.mol/sg.units.L
N_peptide_chains=5
residue_positions=[0,3,5,len(sequence)-1] # Residue positions to calculate its average charge

    # Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

sg.load_parameters(filename='reference_parameters/Lunkad2021.txt') 
sg.load_parameters(filename='reference_parameters/Hass2015.txt')

    # Use a generic parametrization for the aminoacids not parametrized

not_parametrized_neutral_aminoacids=['A','N','Q','G','I','L','M','F','P','O','S','U','T','W','V','J']
not_parametrized_acidic_aminoacids=['C','c']
not_parametrized_basic_aminoacids=['R','n']

for aminoacid_key in sequence:
    if aminoacid_key in not_parametrized_acidic_aminoacids:
        sg.particle(name=aminoacid_key, acidity='acidic', diameter=0.35*sg.units.nm, epsilon=1*sg.units('reduced_energy'))
    elif aminoacid_key in not_parametrized_basic_aminoacids:
        sg.particle(name=aminoacid_key, acidity='basic', diameter=0.35*sg.units.nm, epsilon=1*sg.units('reduced_energy'))
    elif aminoacid_key in not_parametrized_neutral_aminoacids:
        sg.particle(name=aminoacid_key, q=0, diameter=0.35*sg.units.nm, epsilon=1*sg.units('reduced_energy'))

generic_bond_lenght=0.4*sg.units.nm
generic_harmonic_constant=400*sg.units('reduced_energy / nm**2')
generic_bond = interactions.HarmonicBond(k=generic_harmonic_constant.to('reduced_energy / nm**2').magnitude,
                                 r_0=generic_bond_lenght.to('reduced_length').magnitude)
sg.define_default_bond(bond=generic_bond)

    # Create an instance of a sugar molecule object for the peptide

peptide = sg.peptide(name='generic_peptide', sequence=sequence, model=model)

# Salt parameters

c_salt=5e-3 * sg.units.mol/ sg.units.L
cation=sg.particle(name='Na', type=sg.propose_unused_type(), q=1, diameter=0.35*sg.units.nm, epsilon=1*sg.units('reduced_energy'))
anion=sg.particle(name='Cl', type=sg.propose_unused_type(), q=-1, diameter=0.35*sg.units.nm,  epsilon=1*sg.units('reduced_energy'))

# System parameters

volume=N_peptide_chains/(sg.N_A*pep_concentration)
L=volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration=N_peptide_chains/(volume*sg.N_A)
dict_titrable_groups=sg.count_titrable_particles(object=peptide)
total_ionisible_groups=sum(dict_titrable_groups.values())
print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide is ', dict_titrable_groups)

    # Create an instance of an espresso system

system=espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system

sg.add_bonds_to_system(system=system)

# Create your molecules into the espresso system

for _ in range(N_peptide_chains):
    sg.create_object_in_system(object=peptide, system=system, use_default_bond=True)

sg.create_counterions_in_system(object=peptide,cation=cation,anion=anion,system=system) # Create counterions for the peptide chains
c_salt_calculated=sg.create_added_salt_in_system(system=system,cation=cation,anion=anion,c_salt=c_salt)

# Setup the acid-base reactions of the peptide using the constant pH ensemble

RE=sg.setup_constantpH_reactions(counter_ion=cation)

# Setup espresso to track the ionization of the acid/basic groups in peptide

sg.track_ionization(system=system)

# Setup the non-interacting type for speeding up the sampling of the reactions

type_dict=sg.get_all_stored_types()
non_interacting_type=max(type_dict.keys())+1
RE.set_non_interacting_type(non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

sg.setup_lj_interactions(system=system)
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
Z_pH=[] # List of the average global charge at each pH
Z_groups_pH=[] # List of the average charge of groups in residue_positions at each pH

# Main loop for performing simulations at different pH-values

for pH_value in pH_range:

    # Sample list inicialization

    Z_sim=[]
    Z_groups_time_series=[]       

    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:

            system.integrator.run(steps=MD_steps_per_sample)
        
        else:
        
            RE.reaction(total_ionisible_groups)

        if ( step > steps_eq):

            # Get peptide net charge

            Z_net_list = sg.get_net_charge(system=system, object=peptide)
            Z_sim.append(np.mean(np.array(Z_net_list)))
            
            # Get the charge of the residues in residue_position

            Z_charge_res_dict = sg.get_charge_in_residues(system=system, molecule=peptide)
            Z_groups_av=[]
            for residue_position in residue_positions:
                z_pos_av=[]
                for molecule_dict in Z_charge_res_dict:
                    z_pos_av.append(molecule_dict[residue_position][sequence[residue_position]])
                Z_groups_av.append(np.mean(np.array(z_pos_av)))
            Z_groups_time_series.append(Z_groups_av)
        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(system, coordinates)
                vtf.writevcf(system, coordinates)

    sg.write_progress(step=list(pH_range).index(pH_value), total_steps=len(pH_range))
    Z_pH.append(np.array(Z_sim))
    Z_groups_pH.append(np.array(Z_groups_time_series))

    print("pH = {:6.4g} done".format(pH_value))


# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = sg.calculate_HH(object=peptide, pH=list(pH_range))

# Estimate the statistical error and the autocorrelation time of the data
av_net_charge, err_net_charge, tau_net_charge, block_size_net_charge = sg.block_analyze(input_data=Z_pH)

group_averages=[]
for time_serie_pH in Z_groups_pH:
    group_averages.append(sg.block_analyze(input_data=np.transpose(time_serie_pH)))

# Plot the results

fig, ax = plt.subplots(figsize=(10, 7))
plt.errorbar(pH_range, av_net_charge, yerr=err_net_charge, fmt = '-o', capsize=3, label='Simulation')
ax.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.title('Peptide sequence: '+ sequence)
plt.show()

exit()
