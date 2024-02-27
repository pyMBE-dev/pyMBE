# Load espresso, sugar and other necessary libraries
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

# Load some functions from the handy_scripts library for convinience
from handy_scripts.handy_functions import setup_electrostatic_interactions
from handy_scripts.handy_functions import minimize_espresso_system_energy
from handy_scripts.handy_functions import setup_langevin_dynamics
from handy_scripts.handy_functions import block_analyze

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_range = np.linspace(2, 12, num=20)
Samples_per_pH = 1000
MD_steps_per_sample = 1000
steps_eq = int(Samples_per_pH/3)
N_samples_print = 10  # Write the trajectory every 100 samples
probability_reaction =0.5 
SEED = 100
dt = 0.001
solvent_permitivity = 78.3

# Peptide parameters

sequence = 'nGEGGHc'
model = '2beadAA'  # Model with 2 beads per each aminoacid
pep_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
N_peptide_chains = 4

# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
pmb.load_interaction_parameters (filename=pyMBE_path+'/reference_parameters/interaction_parameters/Lunkad2021.txt') 
pmb.load_pka_set (filename=pyMBE_path+'/reference_parameters/pka_sets/Hass2015.txt')

# Use a generic parametrization for the aminoacids not parametrized

not_parametrized_neutral_aminoacids = ['A','N','Q','G','I','L','M','F','P','O','S','U','T','W','V','J']
not_parametrized_acidic_aminoacids = ['C','c']
not_parametrized_basic_aminoacids = ['R','n']

already_defined_AA=[]

for aminoacid_key in sequence:
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

# Defines the peptine in the pyMBE data frame
peptide_name = 'generic_peptide'
pmb.define_peptide (name=peptide_name, sequence=sequence, model=model)

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, q=1, diameter=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, diameter=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
volume = N_peptide_chains/(pmb.N_A*pep_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_pmb_object(name=peptide_name, number_of_objects= N_peptide_chains,espresso_system=espresso_system, use_default_bond=True)
pmb.create_counterions(object_name=peptide_name,cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system) # Create counterions for the peptide chains

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)

with open('frames/trajectory0.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

#List of ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
print('The peptide concentration in your system is ', calculated_peptide_concentration.to('mol/L') , 'with', N_peptide_chains, 'peptides')
print('The ionisable groups in your peptide are ', list_ionisible_groups)

RE, sucessfull_reactions_labels = pmb.setup_cpH(counter_ion=cation_name, constant_pH=2, SEED=SEED)
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

#Setup the potential energy
pmb.setup_lj_interactions (espresso_system=espresso_system)
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

N_frame=0
Z_pH=[] # List of the average global charge at each pH

#Save the pyMBE dataframe in a CSV file
pmb.df.to_csv('df.csv')

# Main loop for performing simulations at different pH-values

for index in tqdm(range(len(pH_range))):
    
    pH_value=pH_range[index]
    # Sample list inicialization
    Z_sim=[]
    RE.constant_pH = pH_value

    # Inner loop for sampling each pH value

    for step in range(Samples_per_pH+steps_eq):
        
        if np.random.random() > probability_reaction:
            espresso_system.integrator.run(steps=MD_steps_per_sample)        
        else:
            RE.reaction( reaction_steps = total_ionisible_groups)

        if ( step > steps_eq):
            # Get peptide net charge
            charge_dict=pmb.calculate_net_charge (  espresso_system=espresso_system, 
                                                    molecule_name=peptide_name)      
            Z_sim.append(charge_dict["mean"])

        if (step % N_samples_print == 0) :
            N_frame+=1
            with open('frames/trajectory'+str(N_frame)+'.vtf', mode='w+t') as coordinates:
                vtf.writevsf(espresso_system, coordinates)
                vtf.writevcf(espresso_system, coordinates)

    Z_pH.append(Z_sim)
    print("pH = {:6.4g} done".format(pH_value))
   
# Estimate the statistical error and the autocorrelation time of the data

print("Net charge analysis")
av_net_charge, err_net_charge, tau_net_charge, block_size_net_charge = block_analyze(input_data=Z_pH)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
Z_HH = pmb.calculate_HH(object_name=peptide_name, 
                        pH_list=pH_range)

fig, ax = plt.subplots(figsize=(10, 7))
plt.errorbar(pH_range, av_net_charge, yerr=err_net_charge, fmt = '-o', capsize=3, label='Simulation')
ax.plot(pH_range, Z_HH, "-k", label='Henderson-Hasselbach')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.title('Peptide sequence: '+ sequence)

plt.show()
