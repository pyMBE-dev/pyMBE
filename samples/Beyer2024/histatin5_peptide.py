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

# Load some functions from the handy_scripts library for convinience
from handy_scripts.handy_functions import setup_electrostatic_interactions
from handy_scripts.handy_functions import minimize_espresso_system_energy
from handy_scripts.handy_functions import setup_langevin_dynamics
from handy_scripts.handy_functions import block_analyze

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_range = np.linspace(2, 12, num=21)
Samples_per_pH = 1000
MD_steps_per_sample = 1000
steps_eq = int(Samples_per_pH/3)
N_samples_print = 10  # Write the trajectory every 100 samples
probability_reaction = 0.5 
SEED = 100
dt = 0.001
solvent_permitivity = 78.3 

L = 22 * pmb.units.nm # Side of the simulation box

# Peptide parameters
peptide_name = 'histidin5'
sequence='NH2-ASP-SER-HIS-ALA-LYS-ARG-HIS-HIS-GLY-TYR-LYS-ARG-LYS-PHE-HIS-GLU-LYS-HIS-HIS-SER-HIS-ARG-GLY-TYR-COOH'
model = '1beadAA'
calculated_peptide_concentration = 1.56e-4 *pmb.units.mol/pmb.units.L
bead_size = 0.4*pmb.units.nm
epsilon = 1*pmb.units('reduced_energy')

# Solution parameters 
cation_name = 'Na'
anion_name = 'Cl'
c_salt = 5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name,  q=1, diameter=0.2*pmb.units.nm, epsilon=epsilon)
pmb.define_particle(name=anion_name,  q=-1, diameter=0.36*pmb.units.nm,  epsilon=epsilon)

# use generic parameters for all beads in the peptide
acidic_aminoacids = ['c','E','D','Y','C']
basic_aminoacids  = ['R','n','K','H']
N_aminoacids = len (pmb.protein_sequence_parser(sequence=sequence))

# Load pKa set
pmb.load_pka_set (filename=pyMBE_path+'/parameters/pka_sets/Nozaki1967.txt')

already_defined_AA=[]
for aminoacid_key in pmb.protein_sequence_parser(sequence=sequence):
    if aminoacid_key in already_defined_AA:
        continue
    if aminoacid_key in acidic_aminoacids:
        pmb.define_particle (name=aminoacid_key, acidity='acidic', diameter=bead_size, epsilon=epsilon)
    elif aminoacid_key in basic_aminoacids:
        pmb.define_particle (name=aminoacid_key, acidity='basic', diameter=bead_size, epsilon=epsilon)
    else:
        pmb.define_particle (name=aminoacid_key, q=0, diameter=bead_size, epsilon=epsilon)
    already_defined_AA.append(aminoacid_key)

generic_bond_lenght=0.4 * pmb.units.nm
generic_harmonic_constant=0.41 * pmb.units.N / pmb.units.m
generic_bond = interactions.HarmonicBond(k=generic_harmonic_constant.to('reduced_energy / reduced_length**2').magnitude,
                                 r_0=generic_bond_lenght.to('reduced_length').magnitude)

pmb.define_default_bond(bond_object=generic_bond, bond_type="harmonic")

# Define the peptide in the pyMBE data frame 

pmb.define_peptide (name=peptide_name, sequence=sequence, model=model)

# Create an instance of an espresso system

espresso_system = espressomd.System(box_l=[L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system

N_peptide_chains = int(L**3*calculated_peptide_concentration*pmb.N_A)
pmb.create_pmb_object (name=peptide_name, number_of_objects= N_peptide_chains, espresso_system=espresso_system, use_default_bond=True, position = [[L.to('reduced_length').magnitude/2]*3])

with open('frames/trajectory0.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates) 

# Create counterions for the peptide chains
pmb.create_counterions(object_name=peptide_name,cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system)
c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)

#List of ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.drop_duplicates().to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.drop_duplicates().to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print('The box length of your system is', L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide are ', list_ionisible_groups)

# Setup the acid-base reactions of the peptide using the constant pH ensemble
RE, sucessfull_reactions_labels = pmb.setup_cpH(counter_ion=cation_name, constant_pH=2, SEED = SEED )
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

#Setup the potential energy
pmb.setup_lj_interactions(espresso_system=espresso_system)
minimize_espresso_system_energy (espresso_system=espresso_system)
setup_electrostatic_interactions(units=pmb.units,
                                            espresso_system=espresso_system,
                                            kT=pmb.kT)
minimize_espresso_system_energy (espresso_system=espresso_system)


setup_langevin_dynamics(espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = SEED,
                                    time_step=dt,
                                    tune_skin=False)

espresso_system.cell_system.skin=0.4

#Save the initial state 
with open('frames/trajectory1.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates) 

Rg_pH=[]
Z_pH=[]
N_frame=0

particle_id_list = pmb.get_particle_id_map(object_name=peptide_name)["all"]
first_peptide_id = min(particle_id_list)

#Save thepyMBE dataframe in a CSV file
pmb.df.to_csv('df.csv')

# Main loop for performing simulations at different pH-values

for index in tqdm(range(len(pH_range))):

    # Sample list inicialization
    pH_value=pH_range[index]
    print (pH_value)
    Z_sim=[]
    Rg_sim=[]
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value
    for step in range(Samples_per_pH+steps_eq):
        if np.random.random() > probability_reaction:
            espresso_system.integrator.run(steps=MD_steps_per_sample)
        else:
            RE.reaction(reaction_steps=total_ionisible_groups)

        if ( step > steps_eq):        
            # Get peptide net charge
            charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system, 
                                                                molecule_name=peptide_name)      
            Z_sim.append(charge_dict["mean"])
            
            Rg = espresso_system.analysis.calc_rg(chain_start=first_peptide_id, number_of_chains=N_peptide_chains, chain_length=len(particle_id_list))
            Rg_value = pmb.units.Quantity(Rg[0], 'reduced_length')
            Rg_nm = Rg_value.to('nm').magnitude
            Rg_sim.append(Rg_nm)

        if (step % N_samples_print == 0) :
            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    Z_pH.append(Z_sim)
    Rg_pH.append(Rg_sim)
    print("pH = {:6.4g} done".format(pH_value))

# Estimate the statistical error and the autocorrelation time of the data

print("Net charge analysis")
av_charge, err_charge, tau_charge, block_size = block_analyze(input_data=pmb.np.array(Z_pH))

print("Rg analysis")
av_rg, err_rg, tau_rg, block_size = block_analyze(input_data=Rg_pH)
    
# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation

Z_HH = pmb.calculate_HH(object_name=peptide_name, 
                        pH_list=pH_range)

# Load the reference data 
reference_file_Path = pyMBE_path+"/reference_data/histatin5_SoftMatter.txt"
reference_data = np.loadtxt(reference_file_Path, delimiter=",")

Z_ref=reference_data[:,1]         
Z_err_ref=reference_data[:,2]

Rg_ref=reference_data[:,4]/10
Rg_err_ref=reference_data[:,5]/10

# Plot the results

plt.figure(figsize=[11, 9])
plt.subplot(1, 2, 1)
plt.suptitle('Peptide sequence: '+''.join(pmb.protein_sequence_parser(sequence=sequence)))
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
