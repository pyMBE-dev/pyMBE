import os
import sys
import inspect
from tqdm import tqdm
import espressomd
import argparse

from espressomd import interactions
from espressomd.io.writer import vtf
from espressomd import electrostatics 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

#Import functions from handy_functions script 
from lib.handy_functions import setup_electrostatic_interactions
from lib.handy_functions import minimize_espresso_system_energy
from lib.handy_functions import setup_langevin_dynamics

# Here you can adjust the width of the panda columns displayed when running the code 
pmb.pd.options.display.max_colwidth = 10

#This line allows you to see the complete amount of rows in the dataframe
pmb.pd.set_option('display.max_rows', None)

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

parser = argparse.ArgumentParser(description='Script to run globular protein simulation in espressomd')

parser.add_argument('--pdb', type=str, required= True,  help='PDB code of the protein')
parser.add_argument('--pH', type=float, required= True,  help='pH value')
parser.add_argument('--path_to_cg', type=str, required= True,  help='Path to the CG structure of the protein')
parser.add_argument('--move_protein', type=float, required= False, default=False,  help='Activates the motion of the protein')
parser.add_argument('--metal_ion_name', type=str, required= False, default=None,  help='Name of the metal ion in the protein')
parser.add_argument('--metal_ion_charge', type=int, required= False, default=None,  help='Charge of the metal ion in the protein')


args = parser.parse_args ()
protein_name = args.pdb

#System Parameters 
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

SEED = 77 
pH_value = args.pH 
c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 
solvent_permitivity = 78.3
epsilon = 1*pmb.units('reduced_energy')

#Simulation Parameters
t_max = 1e3 #  in LJ units of time
stride_obs = 10 #  in LJ units of time
stride_traj = 100 # in LJ units of time
dt = 0.01
N_samples = int (t_max / stride_obs)
integ_steps = int (stride_obs/dt)

#Switch for Electrostatics and WCA interactions 
WCA = True
Electrostatics = True

espresso_system = espressomd.System(box_l=[Box_L.to('reduced_length').magnitude] * 3)
espresso_system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

#Reads the VTF file of the protein model
path_to_cg=pmb.get_resource(args.path_to_cg)
topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)
#Defines the protein in the pmb.df
pmb.define_protein (name=protein_name, topology_dict=topology_dict, model = '2beadAA')

#Create dictionary with the value of epsilon and sigma for each residue
clean_sequence = pmb.df.loc[pmb.df['name']== protein_name].sequence.values[0]

epsilon_dict = {}
sigma_dict = {}

for residue in clean_sequence:
    if residue not in epsilon_dict.keys():
        epsilon_dict [residue] = epsilon
        sigma_dict [residue] = 0.355*pmb.units.nm
    epsilon_dict  ['CA'] = epsilon
    sigma_dict ['CA'] = 0.355*pmb.units.nm

#Define epsilon and sigma for each particle into pmb.df

pmb.define_particles_parameter_from_dict (param_dict = epsilon_dict,
                                            param_name ='epsilon')
pmb.define_particles_parameter_from_dict (param_dict = sigma_dict,
                                            param_name ='sigma')

#Defines the metal ion present in the protein 
if args.metal_ion_name is not None:
    pmb.define_particle(name = args.metal_ion_name, 
                        q=args.metal_ion_charge, 
                        sigma=0.355*pmb.units.nm, 
                        epsilon=epsilon)

# Here we define the solution particles in the pmb.df 
cation_name = 'Na'
anion_name = 'Cl'

pmb.define_particle(name = cation_name, q = 1, sigma=0.2*pmb.units.nm, epsilon=epsilon)
pmb.define_particle(name = anion_name,  q =-1, sigma=0.2*pmb.units.nm, epsilon=epsilon)

# Here we upload the pka set from the reference_parameters folder
path_to_pka=pmb.get_resource('parameters/pka_sets/Nozaki1967.txt') 
pmb.load_pka_set (filename=path_to_pka)

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

c_salt_calculated = pmb.create_added_salt (espresso_system=espresso_system,
                                                    cation_name=cation_name,
                                                    anion_name=anion_name,
                                                    c_salt=c_salt)

#Here we calculated the ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print('The box length of the system is', Box_L.to('reduced_length'), Box_L.to('nm'))
print('The ionisable groups in the protein are ', list_ionisible_groups)
print ('The total amount of ionizable groups are:',total_ionisible_groups)

#Setup of the reactions in espresso 
RE, sucessfull_reactions_labels = pmb.setup_cpH(counter_ion=cation_name, 
                                                constant_pH= pH_value, 
                                                SEED = SEED )
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

if (WCA):

    print ('Setup of LJ interactions.. ')

    pmb.setup_lj_interactions (espresso_system=espresso_system)
    minimize_espresso_system_energy (espresso_system=espresso_system)

    if (Electrostatics):

        setup_electrostatic_interactions (units=pmb.units,
                                        espresso_system=espresso_system,
                                        kT=pmb.kT)

#Save the initial state 
n_frame = 0
pmb.write_output_vtf_file(espresso_system=espresso_system,
                        filename=f"frames/trajectory{n_frame}.vtf")

setup_langevin_dynamics (espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = SEED)

observables_df = pmb.pd.DataFrame()
time_step = []
net_charge_list = []
net_charge_amino_save = {}

Z_sim=[]
particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df (filename='df.csv')

#Here we start the main loop over the Nsamples 

for step in tqdm(range(N_samples)):
        
        espresso_system.integrator.run (steps = integ_steps)
        RE.reaction( reaction_steps = total_ionisible_groups)

        charge_dict=pmb.calculate_net_charge (espresso_system=espresso_system, 
                                                molecule_name=protein_name)
        net_charge = charge_dict['mean']
        net_charge_residues = charge_dict ['residues']

        time_step.append (str(espresso_system.time))
        net_charge_list.append (net_charge)

        if len(net_charge_amino_save.keys()) == 0:
            for amino in net_charge_residues.keys():
                net_charge_amino_save [amino] = []
        for amino in net_charge_residues.keys():            
            net_charge_amino_save [amino].append (net_charge_residues[amino])

        if (step % stride_traj == 0  ):
            n_frame +=1
            pmb.write_output_vtf_file(espresso_system=espresso_system,
                                        filename=f"frames/trajectory{n_frame}.vtf")

#We save the calculated observables into a new df 
observables_df['time'] = time_step 
observables_df['Znet'] = net_charge_list

observables_df=pmb.pd.concat([observables_df,pmb.pd.DataFrame(net_charge_amino_save)], axis=1)

observables_df.to_csv(f'pH-{pH_value}_observables.csv',index=False)
