#
# Copyright (C) 2024 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
from tqdm import tqdm
import espressomd
import argparse
import numpy as np
import pandas as pd 
from espressomd.io.writer import vtf
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(SEED=42)

#Import functions from handy_functions script 
from lib.handy_functions import setup_electrostatic_interactions
from lib.handy_functions import minimize_espresso_system_energy
from lib.handy_functions import setup_langevin_dynamics
from lib import analysis
# Here you can adjust the width of the panda columns displayed when running the code 
pd.options.display.max_colwidth = 10

#This line allows you to see the complete amount of rows in the dataframe
pd.set_option('display.max_rows', None)

valid_modes=["short-run","long-run", "test"]

parser = argparse.ArgumentParser(description='Script to run globular protein simulation in espressomd')

parser.add_argument('--pdb',
                    type=str, 
                    required= True,  
                    help='PDB code of the protein')
parser.add_argument('--pH', 
                    type=float, 
                    required= True,  
                    help='pH value')
parser.add_argument('--path_to_cg', 
                    type=str, 
                    required= True,  
                    help='Path to the CG structure of the protein')
parser.add_argument('--move_protein', 
                    type=float, 
                    required= False, 
                    default=False,  
                    help='Activates the motion of the protein')

parser.add_argument('--mode',
                    type=str,
                    default= "short-run",
                    help='sets for how long the simulation runs, valid modes are {valid_modes}')

parser.add_argument('--output',
                    type=str,
                    required= False,
                    help='output directory')

parser.add_argument('--no_verbose', action='store_false', help="Switch to deactivate verbose",default=True)

args = parser.parse_args ()
mode=args.mode
verbose=args.no_verbose
protein_name = args.pdb
pH_value = args.pH 

inputs={"pH": args.pH,
        "pdb": args.pdb}

#System Parameters 
LANGEVIN_SEED = 77 

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 
solvent_permitivity = 78.3
epsilon = 1*pmb.units('reduced_energy')
sigma = 1*pmb.units("reduced_length")
ion_size = 0.4*pmb.units.nm

#Simulation Parameters
 #  in LJ units of time
dt = 0.01
stride_traj = 100 # in LJ units of time

if mode == 'short-run':
    stride_obs = 10 #  in LJ units of time
    integ_steps = int (stride_obs/dt)
    t_max = 5e3
    N_samples = int (t_max / stride_obs)
    
elif mode == 'long-run':
    stride_obs = 10 #  in LJ units of time
    integ_steps = int (stride_obs/dt)
    t_max = 1e5
    N_samples = int (t_max / stride_obs)
    
elif mode == 'test':
    t_max = 1e3
    stride_obs = 10 #  in LJ units of time
    integ_steps = int (stride_obs/dt)
    N_samples = int (t_max / stride_obs)
else: 
    raise RuntimeError()    

#Switch for Electrostatics and WCA interactions 
WCA = True
Electrostatics = True

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

espresso_system = espressomd.System(box_l=[Box_L.to('reduced_length').magnitude] * 3)
espresso_system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

#Reads the VTF file of the protein model
path_to_cg=pmb.get_resource(args.path_to_cg)
topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)
#Defines the protein in the pmb.df
pmb.define_protein (name=protein_name, 
                    topology_dict=topology_dict, 
                    model = '2beadAA',
                    lj_setup_mode = "wca")

# Here we define the solution particles in the pmb.df 
cation_name = 'Na'
anion_name = 'Cl'

pmb.define_particle(name = cation_name, 
                    q = 1, 
                    sigma=sigma, 
                    epsilon=epsilon,
                    offset=ion_size-sigma)

pmb.define_particle(name = anion_name,  
                    q =-1, 
                    sigma=sigma, 
                    epsilon=epsilon,
                    offset=ion_size-sigma)

# Here we upload the pka set from the reference_parameters folder
path_to_pka=pmb.get_resource('parameters/pka_sets/Nozaki1967.json') 
pmb.load_pka_set(filename=path_to_pka)

#We create the protein in espresso 
pmb.create_protein(name=protein_name,
                    number_of_proteins=1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

#Here we activate the motion of the protein 
if args.move_protein:
    pmb.enable_motion_of_rigid_object(espresso_system=espresso_system,
                                        name=protein_name)

# Here we put the protein on the center of the simulation box
protein_id = pmb.df.loc[pmb.df['name']==protein_name].molecule_id.values[0]
pmb.center_molecule_in_simulation_box (molecule_id=protein_id,
                                    espresso_system=espresso_system)

# Creates counterions and added salt 
pmb.create_counterions (object_name=protein_name,
                        cation_name=cation_name,
                        anion_name=anion_name,
                        espresso_system=espresso_system)

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                          cation_name=cation_name,
                                          anion_name=anion_name,
                                          c_salt=c_salt)

#Here we calculated the ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

if verbose:
    print('The box length of the system is', Box_L.to('reduced_length'), Box_L.to('nm'))
    print('The ionisable groups in the protein are ', list_ionisible_groups)
    print ('The total amount of ionizable groups are:',total_ionisible_groups)

#Setup of the reactions in espresso 
cpH, labels = pmb.setup_cpH(counter_ion=cation_name, 
                            constant_pH= pH_value)
if verbose:
    print('The acid-base reaction has been sucessfully setup for ', labels)

type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
cpH.set_non_interacting_type (type=non_interacting_type)
if verbose:
    print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

if WCA:
    pmb.setup_lj_interactions (espresso_system=espresso_system)
    minimize_espresso_system_energy (espresso_system=espresso_system)
    if Electrostatics:
        setup_electrostatic_interactions (units=pmb.units,
                                        espresso_system=espresso_system,
                                        kT=pmb.kT)

#Save the initial state 
n_frame = 0
with open('frames/trajectory'+str(n_frame)+'.vtf', mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

setup_langevin_dynamics (espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = LANGEVIN_SEED)

observables_df = pd.DataFrame()
time_step = []
net_charge_list = []

Z_sim=[]
particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df (filename='df.csv')

#Here we start the main loop over the Nsamples 

labels_obs=["time","charge"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

charge_dict=pmb.calculate_net_charge (espresso_system=espresso_system, 
                                            molecule_name=protein_name)
    
net_charge_residues = charge_dict ['residues']
net_charge_amino_save = {}
AA_label_list=[]    
for amino in net_charge_residues.keys():
    amino_part_row=pmb.df[(pmb.df['residue_id']== amino) & (pmb.df['acidity'] != "inert")]
    if not amino_part_row.empty:
        label = f'charge_{amino_part_row["name"].values[0]}'
        if label not in AA_label_list:
            AA_label_list.append(label)
            net_charge_amino_save[label] = []
            time_series[label] = []

for step in tqdm(range(N_samples),disable=not verbose):      
    espresso_system.integrator.run (steps = integ_steps)
    cpH.reaction(reaction_steps = total_ionisible_groups)
    charge_dict=pmb.calculate_net_charge (espresso_system=espresso_system, 
                                            molecule_name=protein_name)
    charge_residues = charge_dict['residues']
    charge_residues_per_type={}

    for label in AA_label_list:
        charge_residues_per_type[label]=[]

    for amino in charge_residues.keys():
        amino_part_row=pmb.df[(pmb.df['residue_id']== amino) & (pmb.df['acidity'] != "inert")]
        if not amino_part_row.empty:
            label = f'charge_{amino_part_row["name"].values[0]}'
            if label in AA_label_list:
                charge_residues_per_type[label].append(charge_residues[amino])

    if step % stride_traj == 0 :
        n_frame +=1
        with open('frames/trajectory'+str(n_frame)+'.vtf', mode='w+t') as coordinates:
            vtf.writevsf(espresso_system, coordinates)
            vtf.writevcf(espresso_system, coordinates)

    # Store observables
    time_series["time"].append(espresso_system.time)
    time_series["charge"].append(charge_dict["mean"])
    
    for label in AA_label_list:
        charge_amino = np.mean(charge_residues_per_type[label])
        time_series[label].append(charge_amino)

data_path = args.output
if data_path is None:
    data_path=pmb.get_resource(path="samples/Beyer2024/")+"/time_series/globular_protein"

if not os.path.exists(data_path):
    os.makedirs(data_path)
    
time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(f"{data_path}/{filename}_time_series.csv", index=False)

if verbose:
    print("*** DONE ***")
