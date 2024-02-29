# Load espresso, sugar and other necessary libraries
import sys
import os 
import inspect
from matplotlib.style import use
import espressomd
import numpy as np
import matplotlib.pyplot as plt
import argparse 
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

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

parser = argparse.ArgumentParser(description='Script to run the peptide test cases for pyMBE')
parser.add_argument('--sequence', type=str, required= True,  help='sequence of the peptide')
parser.add_argument('--pH', type=float, required= True,  help='pH of the solution')
parser.add_argument('--verbose', type=bool, default= False,  help='switch to activate prints')
args = parser.parse_args()

# Inputs
sequence=args.sequence
pH=args.pH

# Simulation parameters
Nsamples = 1000
MD_steps_per_sample = 1000
SEED = 100
dt = 0.01
solvent_permitivity = 78.3
pep_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L

# Sanity check

Lunkad_test_sequences=["E"*5+"H"*5,"K"*5+"D"*5]
Blanco_test_sequence=["nDSHAKRHHGYKRKFHHSHRGYc"]

valid_sequences=Lunkad_test_sequences+Blanco_test_sequence

if sequence not in valid_sequences:
    raise ValueError(f"ERROR: the only valid peptide sequence for this test script are {valid_sequences}")

if sequence in Lunkad_test_sequences:
    pmb.load_interaction_parameters (filename=pyMBE_path+'/parameters/interaction_parameters/Lunkad2021.txt') 
    pmb.load_pka_set (filename=pyMBE_path+'/parameters/pka_sets/CRC1991.txt')
    model = '2beadAA'  # Model with 2 beads per each aminoacid
    N_peptide_chains = 4
    diameter_Na=0.35*pmb.units.nm
    diameter_Cl=0.35*pmb.units.nm
    c_salt=1e-2 * pmb.units.mol/ pmb.units.L

elif sequence in Blanco_test_sequence:
    pmb.load_interaction_parameters (filename=pyMBE_path+'/parameters/interaction_parameters/Blanco2020.txt') 
    pmb.load_pka_set (filename=pyMBE_path+'/parameters/pka_sets/Nozaki1967.txt')
    model = '1beadAA'
    N_peptide_chains = 1
    c_salt = 5e-3 * pmb.units.mol/ pmb.units.L
    diameter_Na=0.2*pmb.units.nm
    diameter_Cl=0.36*pmb.units.nm

pmb.define_peptide (name=sequence, sequence=sequence, model=model)

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, 
                    q=1, 
                    diameter=0.35*pmb.units.nm, 
                    epsilon=1*pmb.units('reduced_energy'))

pmb.define_particle(name=anion_name,  
                    q=-1, 
                    diameter=0.35*pmb.units.nm,  
                    epsilon=1*pmb.units('reduced_energy'))

# System parameters
volume = N_peptide_chains/(pmb.N_A*pep_concentration)
L = volume ** (1./3.) # Side of the simulation box

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_pmb_object(name=sequence, 
                    number_of_objects= N_peptide_chains,
                    espresso_system=espresso_system)

# Create counterions for the peptide chains
pmb.create_counterions(object_name=sequence,
                    cation_name=cation_name,
                    anion_name=anion_name,
                    espresso_system=espresso_system) 

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                        cation_name=cation_name,
                                        anion_name=anion_name,
                                        c_salt=c_salt)


RE, sucessfull_reactions_labels = pmb.setup_cpH(counter_ion=cation_name, 
                                                constant_pH=pH, 
                                                SEED=SEED)

if args.verbose:

    print("The box length of your system is", L.to('reduced_length'), L.to('nm'))
    print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map =pmb.get_type_map()
espresso_system.setup_type_map(type_list = list(type_map.values()))

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

# Main loop 
  
Z_samples=[]
for sample in (pbar := tqdm(range(Nsamples))):
    
    espresso_system.integrator.run(steps=MD_steps_per_sample)        
    RE.reaction( reaction_steps = len(sequence))

    charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system, 
                                            molecule_name=sequence)      
    Z_samples.append(charge_dict["mean"])
    pbar.set_description(f"Sample = {sample}/{Nsamples}")

if args.verbose:
    print("*** DONE ***")
   
