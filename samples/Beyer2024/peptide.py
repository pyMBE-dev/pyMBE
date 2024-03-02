# Load espresso, sugar and other necessary libraries
import sys
import os 
import inspect
import espressomd
import numpy as np
import pandas as pd
import argparse 
from tqdm import tqdm
from espressomd.io.writer import vtf
from espressomd import interactions
from espressomd import electrostatics

# Import pyMBE
import pyMBE
from lib import analysis
from lib import handy_functions as hf

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

parser = argparse.ArgumentParser(description='Script to run the peptide test cases for pyMBE')
parser.add_argument('--sequence', 
                    type=str, 
                    required= True,  
                    help='sequence of the peptide')
parser.add_argument('--pH', 
                    type=float, 
                    required= True,  
                    help='pH of the solution')
parser.add_argument('--verbose', 
                    type=bool, 
                    default= False,  
                    help='switch to activate prints')
parser.add_argument('--mode', 
                    type=str, 
                    default= "short-run",  
                    help='sets for how long the simulation runs, valid modes are "short-run" and "long-run"')
args = parser.parse_args()

# Inputs
sequence=args.sequence
pH=args.pH
inputs={"pH": args.pH,
        "sequence": args.sequence}

mode=args.mode

valid_modes=["short-run","long-run"]
if mode not in valid_modes:
    raise ValueError(f"Mode {mode} is not currently supported, valid modes are {valid_modes}")

# Simulation parameters
if mode == "short-run":
    Nsamples = 1000
    MD_steps_per_sample = 1000

if mode == "long-run":
    Nsamples = 5000
    MD_steps_per_sample = 5000

SEED = 100
dt = 0.01
solvent_permitivity = 78.3
pep_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L

# Sanity check
Lunkad_test_sequences=["E"*5+"H"*5,"K"*5+"D"*5]
Blanco_test_sequence=["nDSHAKRHHGYKRKFHEKHHSHRGYc"]

valid_sequences=Lunkad_test_sequences+Blanco_test_sequence

if sequence not in valid_sequences:
    raise ValueError(f"ERROR: the only valid peptide sequence for this test script are {valid_sequences}")

if sequence in Lunkad_test_sequences:
    pmb.load_interaction_parameters (filename='parameters/peptides/Lunkad2021.txt') 
    pmb.load_pka_set (filename='parameters/pka_sets/CRC1991.txt')
    model = '2beadAA'  # Model with 2 beads per each aminoacid
    N_peptide_chains = 4
    diameter_Na=0.35*pmb.units.nm
    diameter_Cl=0.35*pmb.units.nm
    c_salt=1e-2 * pmb.units.mol/ pmb.units.L
    chain_length=len(sequence)*2

elif sequence in Blanco_test_sequence:
    pmb.load_interaction_parameters (pmb.get_resource(path='parameters/peptides/Blanco2020.txt')) 
    pmb.load_pka_set (pmb.get_resource(path='parameters/pka_sets/Nozaki1967.txt'))
    model = '1beadAA'
    N_peptide_chains = 1
    c_salt = 5e-3 * pmb.units.mol/ pmb.units.L
    diameter_Na=0.2*pmb.units.nm
    diameter_Cl=0.36*pmb.units.nm
    chain_length=len(sequence)

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
hf.minimize_espresso_system_energy (espresso_system=espresso_system)
hf.setup_electrostatic_interactions(units=pmb.units,
                                        espresso_system=espresso_system,
                                        kT=pmb.kT)
hf.minimize_espresso_system_energy (espresso_system=espresso_system)


hf.setup_langevin_dynamics(espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = SEED,
                                    time_step=dt,
                                    tune_skin=False)

espresso_system.cell_system.skin=0.4

# Main loop 
  
labels_obs=["time","charge","rg"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

for sample in (pbar := tqdm(range(Nsamples))):
    # Run LD
    espresso_system.integrator.run(steps=MD_steps_per_sample)
    # Run MC        
    RE.reaction(reaction_steps=len(sequence))
    # Sample observables
    charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system, 
                                            molecule_name=sequence)      
    
    Rg = espresso_system.analysis.calc_rg(chain_start=0, 
                                        number_of_chains=N_peptide_chains, 
                                        chain_length=chain_length)
    # Store observables
    time_series["time"].append(espresso_system.time)
    time_series["charge"].append(charge_dict["mean"])
    time_series["rg"].append(Rg[0])

data_path=pmb.get_resource(path="samples/Beyer2024")+"/time_series"

if not os.path.exists(data_path):
    os.makedirs(data_path)

time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(f"{data_path}/{filename}_time_series.csv", index=False)

if args.verbose:
    print("*** DONE ***")
   
