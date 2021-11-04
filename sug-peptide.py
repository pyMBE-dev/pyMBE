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

sg.set_reduced_units(unit_length=0.4*sg.units.nm)
sg.param.default.default_bond.k=0.41 * sg.units.N / sg.units.m

# Peptide parameters
sequence='NH2-ASP-SER-HIS-ALA-LYS-ARG-HIS-HIS-GLY-TYR-LYS-ARG-LYS-PHE-HIS-GLU-LYS-HIS-HIS-SER-HIS-ARG-GLY-TYR-COOH'
model='1beadpeptide'
pKa_set='nozaki' 
pep_concentration=1.56e-4 *sg.units.mol/sg.units.L

peptide=sg.molecule(sequence=sequence, model=model, pKa_set=pKa_set)

#sg.write_parameters(peptide)

# Salt parameters

c_salt=0e-3 * sg.units.mol/ sg.units.L

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

pH_range = np.linspace(2, 12, num=21)
Samples_per_pH= 10
MD_steps_per_sample=1000
steps_eq=int(Samples_per_pH/3)
N_samples_print= 100  # Write the trajectory every 100 samples
probability_reaction=0.5 

# Add peptides to your simulation box

volume=system.volume()*sg.units('reduced_length**3')
#print(volume,(L.to('reduced_length')**3))
peptide.N=int(volume*pep_concentration*sg.N_A)
sg.create_molecule(peptide, system)
sg.write_parameters(peptide)
calculated_peptide_concentration=peptide.N/(volume*sg.N_A)
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
Rg_pH=[]

# Main loop for performing simulations at different pH-values

for pH_value in pH_range:

    Z_sim=[]
    Rg_sim=[]
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
            Rg = system.analysis.calc_rg(chain_start=0, number_of_chains=1, chain_length=26)
            Rg_value = sg.units.Quantity(Rg[0], 'reduced_length')
            Rg_nm = Rg_value.to('nm').magnitude
            Rg_sim.append(Rg_nm)

        if (step % N_samples_print == 0) :

            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(system, coordinates)
                vtf.writevcf(system, coordinates)

    Z_pH.append(Z_sim)
    Rg_pH.append(Rg_sim)
    print("pH = {:6.4g} done".format(pH_value))
    
# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = sg.calculate_HH(mol=peptide, pH=list(pH_range))



# statistical analysis of the results
def block_analyze(input_data, n_blocks=16):
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
    return av_data, err_data, tau_data, block_size


charge_file=open("prot_charge.dat","w")
rg_file=open("prot_rg.dat","w")


# estimate the statistical error and the autocorrelation time using the formula given by Janke
av_charge, err_charge, tau_charge, block_size = block_analyze(Z_pH)
print("average charge = ", av_charge)
print("err = ", err_charge)
print("tau = ", tau_charge)
for i in range(len(av_charge)):
    charge_file.write(str(pH_range[i]) + "  " + str(av_charge[i]) + "  " + str(err_charge[i]) + " \n")
charge_file.close()


av_rg, err_rg, tau_rg, block_size = block_analyze(Rg_pH)
print("average radius of gyration = ", av_rg)
print("err = ", err_rg)
print("tau = ", tau_rg)
for i in range(len(av_rg)):
    rg_file.write(str(pH_range[i]) + "  " + str(av_rg[i]) + "  " + str(err_rg[i]) + " \n")
rg_file.close()


# check if the blocks contain enough data for reliable error estimates
print("uncorrelated samples per block:\nblock_size/tau = ",
      block_size/tau_charge)
threshold = 10.  # block size should be much greater than the correlation time
if np.any(block_size / tau_charge < threshold):
    print("\nWarning: some blocks may contain less than ", threshold, "uncorrelated samples."
          "\nYour error estimated may be unreliable."
          "\nPlease, check them using a more sophisticated method or run a longer simulation.")
    print("? block_size/tau > threshold ? :", block_size/tau_charge > threshold)
else:
    print("\nAll blocks seem to contain more than ", threshold, "uncorrelated samples.\
    Error estimates should be OK.")


# Reference data
Z_ref=[13.674, 12.891, 12.266, 11.788, 11.169, 10.167, 8.853, 7.447, 6.252, 5.416, 4.837, 4.278, 3.69, 3.112, 2.501, 1.729, 0.826, -0.115, -0.946 , -1.654,  -2.275]
Z_err_ref=[0.00692820323028, 0.00640312423743, 0.00469041575982, 0.00761577310586, 0.00721110255093, 0.00721110255093, 0.00714142842854, 0.00707106781187, 0.0636631761696, 0.0637024332345, 0.0637024332345, 0.0895321171424, 0.0895823643358, 0.089554452709, 0.089554452709, 0.0895377015564, 0.0895879456177, 0.0634428877022, 0.0634428877022, 0.0634428877022, 0.00721110255093]
Rg_ref=[2.063699, 2.013557, 1.956493, 1.920424, 1.890579, 1.832828, 1.752297, 1.664951, 1.597367, 1.557794, 1.533217, 1.505957, 1.487000, 1.462986, 1.456899, 1.447990, 1.438085, 1.447151, 1.449776, 1.472057, 1.486251]
Rg_err_ref=[0.04091, 0.04239, 0.04242, 0.04136, 0.04389, 0.04479, 0.04384, 0.04380, 0.04232, 0.04066, 0.04011, 0.04106, 0.04113, 0.03784, 0.04007, 0.03859, 0.03748, 0.03768, 0.03693, 0.03761, 0.03809]


# Plot the results
peptide_seq = ''.join([str(elem) for elem in peptide.sequence])
err_charge=np.array(err_charge)
Z_err_ref=np.array(Z_err_ref)
fig, ax = plt.subplots(figsize=(10, 7))
ax.errorbar(pH_range, av_charge, yerr=err_charge, fmt = '-o', capsize=3, label='Simulation')
ax.errorbar(pH_range, Z_ref, yerr=Z_err_ref, color='b', fmt = '-o', ecolor = 'b', capsize=3, label='Reference values')
ax.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
ax.axhline(y=0.0, color="gray", linestyle="--")
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.title('Peptide sequence: '+ peptide_seq)
plt.show()

err_rg=np.array(err_rg)
Rg_ref=np.array(Rg_ref)
Rg_err_ref=np.array(Rg_err_ref)
fig, ax = plt.subplots(figsize=(10, 7))
ax.errorbar(pH_range, av_rg, yerr=err_rg, fmt = '-o', capsize=3, label='Simulation')
ax.errorbar(pH_range, Rg_ref, yerr=Rg_err_ref, color='b', fmt = '-o', ecolor = 'b', capsize=3, label='Reference values')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Radius of gyration of Histatin-5 / e')
plt.title('Peptide sequence: '+ peptide_seq)
plt.show()


exit()
