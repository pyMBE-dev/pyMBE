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
import sys
import numpy as np 
import espressomd
import pyMBE
import re
import json
from pint import UnitRegistry, Quantity
from lib.handy_functions import setup_electrostatic_interactions
from lib.handy_functions import minimize_espresso_system_energy
from lib.handy_functions import setup_langevin_dynamics

ureg = UnitRegistry()

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

def custom_serializer(obj):
    if isinstance(obj, Quantity):
        return {"value": obj.magnitude, "unit": str(obj.units)}  
    raise TypeError(f"Type {type(obj)} not serializable")

def custom_deserializer(dct):
    if "value" in dct and "unit" in dct:
        return ureg.Quantity(dct["value"], dct["unit"])  
    return dct  

mode = "test"
protein_pdb = '1beb'
valid_modes = ["save","test"]

if mode not in valid_modes:
    raise ValueError(f"mode {mode} not supported, valid modes are {valid_modes}")

print("*** Unit test: check that read_protein_vtf_in_df() loads the protein topology correctly ***")

filename = "testsuite/tests_data/protein_topology_dict.json"
path_to_cg=pmb.get_resource(f'parameters/globular_proteins/{protein_pdb}.vtf')
topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)
path_to_parfile=pmb.get_resource(filename)

if mode == "save":
    with open (path_to_parfile, "w") as output:
        json.dump(topology_dict, output,default=custom_serializer)
    sys.exit()
elif mode == "test":
    with open (path_to_parfile, "r") as file:
        load_json = json.load(file,object_hook=custom_deserializer)

np.testing.assert_equal(actual= topology_dict, 
                        desired= load_json,
                        verbose = True)

print("*** Unit test passed ***")


print("*** Unit test: check that define_protein() defines the aminoacids in the protein correctly ***")

protein_model = '2beadAA'

pmb.define_protein (name=protein_pdb, 
                    topology_dict=topology_dict, 
                    model = protein_model,
                    lj_setup_mode = "wca")
sequence = []
clean_sequence= []

for aminoacid in topology_dict.keys():
    
    input_parameters=topology_dict[aminoacid]
    residue_name = re.split(r'\d+', aminoacid)[0]
    sequence.append(residue_name)
    
    if residue_name not in ['CA', 'n', 'c']:
        clean_sequence.append(residue_name)

    
    for index in pmb.df[pmb.df['name']==residue_name].index:
        if residue_name not in sequence:           
            np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="particle", 
                                verbose=True)

residue_list = pmb.define_AA_residues(sequence= clean_sequence,
                                      model = protein_model)

for residue in residue_list:
    for index in pmb.df[pmb.df['name']==residue].index:
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                        desired="residue", 
                        verbose=True)

protein_index = pmb.df[pmb.df['name']==protein_pdb].index

np.testing.assert_equal(actual=str(pmb.df.loc[protein_index, "name"].values[0]), 
                            desired=protein_pdb, 
                            verbose=True)        

np.testing.assert_equal(actual=pmb.df.loc[protein_index, ('sequence','')].values[0], 
                    desired=clean_sequence, 
                    verbose=True)


np.testing.assert_equal(actual=pmb.df.loc[protein_index, ('residue_list','')].values[0], 
                    desired=residue_list, 
                    verbose=True)  

print("*** Unit test passed ***")

print("*** Unit test: check that define_protein() raises a ValueError if the user provides a wrong model")

input_parameters={"name": protein_pdb,
                 "topology_dict": topology_dict,
                 "model" : "3beadAA",
                "lj_setup_mode": "wca"}

np.testing.assert_raises(ValueError, pmb.define_protein, **input_parameters)

input_parameters={"name": protein_pdb,
                 "topology_dict": topology_dict,
                 "model" : protein_model,
                "lj_setup_mode": "awc"}

np.testing.assert_raises(ValueError, pmb.define_protein, **input_parameters)

print("*** Unit test passed ***")


print("*** Unit test: check that create_protein() creates all the particles in the protein into the espresso_system with the properties defined in pmb.df  ***")

espresso_system=espressomd.System(box_l = [10]*3)
espresso_system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

# Here we upload the pka set from the reference_parameters folder
path_to_pka=pmb.get_resource('parameters/pka_sets/Nozaki1967.json') 
pmb.load_pka_set(filename=path_to_pka)

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

residue_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].residue_id.dropna().to_list()

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

molecule_id = pmb.df.loc[pmb.df['name']==protein_pdb].molecule_id.values[0]

center_of_mass_es = pmb.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)

center_of_mass = np.zeros(3)
axis_list = [0,1,2]

for aminoacid in topology_dict.keys():
    initial_pos = topology_dict[aminoacid]['initial_pos']

    for axis in axis_list:
        center_of_mass[axis] +=  initial_pos[axis]
center_of_mass = center_of_mass/ len(topology_dict.keys())

distance_es = np.zeros(3)
distance_topology = np.zeros(3)

for id in particle_id_list:

    initial_pos_es = espresso_system.part.by_id(id).pos
    charge = espresso_system.part.by_id(id).q
    es_type = espresso_system.part.by_id(id).type

    residue_id = pmb.df.loc[pmb.df['particle_id']==id].residue_id.values[0]
    residue_name = pmb.df.loc[pmb.df['particle_id']==id].name.values[0]

    initial_pos = topology_dict[residue_name+residue_id]['initial_pos']
    index = pmb.df.loc[pmb.df['particle_id']==id].index

    for axis in axis_list:
        distance_es[axis] = (initial_pos_es[axis] - center_of_mass_es[axis])**2
        distance_topology[axis] = (initial_pos[axis] - center_of_mass[axis])**2

    relative_distance_es = np.sqrt(np.sum(distance_es))
    relative_distance = np.sqrt(np.sum(distance_es))
    
    np.testing.assert_equal(actual=relative_distance_es, 
                        desired=relative_distance, 
                        verbose=True)

    np.testing.assert_equal(actual=charge, 
                        desired=pmb.df.loc[index, ("state_one","z")].values[0], 
                        verbose=True)
      
    np.testing.assert_equal(actual=es_type, 
                    desired=pmb.df.loc[index, ("state_one","es_type")].values[0], 
                    verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_protein() does not create any protein for number_of_proteins <= 0  ***")

starting_number_of_particles=len(espresso_system.part.all())

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=0,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=-1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that enable_motion_of_rigid_object() moves the protein correctly ***")

positions = []
for pid in particle_id_list:
    positions.append(espresso_system.part.by_id(pid).pos)

pmb.enable_motion_of_rigid_object(espresso_system=espresso_system,
                                  name=protein_pdb)

dt = 0.01
t_max = 1e3
stride_obs = 10 #  in LJ units of time
integ_steps = int (stride_obs/dt)
N_samples = int (t_max / stride_obs)
epsilon = 1*pmb.units('reduced_energy')
sigma = 1*pmb.units("reduced_length")
ion_size = 0.4*pmb.units.nm

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 

pH_value= 3

# Here we define the solution particles in the pmb.df 
cation_name = 'Na'
anion_name = 'Cl'

pmb.define_particle(name = cation_name, 
                    z = 1, 
                    sigma=sigma, 
                    epsilon=epsilon,
                    offset=ion_size-sigma)

pmb.define_particle(name = anion_name,  
                    z =-1, 
                    sigma=sigma, 
                    epsilon=epsilon,
                    offset=ion_size-sigma)

# Estimate the radius of the protein
protein_ids = pmb.get_particle_id_map(object_name=protein_pdb)["all"]
protein_radius=0
protein_center=espresso_system.box_l/2
for pid in protein_ids:
    protein_part = espresso_system.part.by_id(pid)
    dist = protein_part.pos - protein_center
    dist = np.linalg.norm(dist)
    if dist > protein_radius:
        protein_radius = dist

# Create counter-ions 
protein_net_charge = pmb.calculate_net_charge(espresso_system=espresso_system,
                                              molecule_name=protein_pdb,
                                              dimensionless=True)["mean"]

## Get coordinates outside the volume occupied by the protein
counter_ion_coords=pmb.generate_coordinates_outside_sphere(center=protein_center,
                                                            radius=protein_radius,
                                                            max_dist=Box_L.m_as("reduced_length"),
                                                            n_samples=abs(protein_net_charge))
if protein_net_charge > 0:
    counter_ion_name=anion_name
elif protein_net_charge < 0:
    counter_ion_name=cation_name

pmb.create_particle(name=counter_ion_name,
                    espresso_system=espresso_system,
                    number_of_particles=int(abs(protein_net_charge)),
                    position=counter_ion_coords)

# Create added salt ions
N_ions= int((Box_V.to("reduced_length**3")*c_salt.to('mol/reduced_length**3')*pmb.N_A).magnitude)

## Get coordinates outside the volume occupied by the protein
added_salt_ions_coords=pmb.generate_coordinates_outside_sphere(center=protein_center,
                                                            radius=protein_radius,
                                                            max_dist=Box_L.m_as("reduced_length"),
                                                            n_samples=N_ions*2)
## Create cations
pmb.create_particle(name=cation_name,
                    espresso_system=espresso_system,
                    number_of_particles=N_ions,
                    position=added_salt_ions_coords[:N_ions])
## Create anions
pmb.create_particle(name=anion_name,
                    espresso_system=espresso_system,
                    number_of_particles=N_ions,
                    position=added_salt_ions_coords[N_ions:])

#Here we calculated the ionisible groups 
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

#Setup of the reactions in espresso 
cpH, labels = pmb.setup_cpH(counter_ion=cation_name, 
                            constant_pH= pH_value)

type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
cpH.set_non_interacting_type (type=non_interacting_type)

pmb.setup_lj_interactions (espresso_system=espresso_system)
minimize_espresso_system_energy (espresso_system=espresso_system)
setup_electrostatic_interactions (units=pmb.units,
                                        espresso_system=espresso_system,
                                        kT=pmb.kT)

setup_langevin_dynamics (espresso_system=espresso_system, 
                        kT = pmb.kT, 
                        SEED = 77)

print("*** Running simulation ***")
for step in range(N_samples):      
    espresso_system.integrator.run (steps = integ_steps)
    cpH.reaction(reaction_steps = total_ionisible_groups)
    

positions_enable_motion = []
for pid in particle_id_list:
    positions_enable_motion.append(espresso_system.part.by_id(pid).pos)

np.testing.assert_raises(AssertionError, np.testing.assert_array_equal, positions, positions_enable_motion)

print("*** Unit test passed ***")


print("*** Unit test: check that enable_motion_of_rigid_object() raises a ValueError if a wrong pmb_type is provided***")

input_parameters = {"espresso_system":espresso_system,
                    "name": "CA"}

np.testing.assert_raises(ValueError, pmb.enable_motion_of_rigid_object, **input_parameters)

print("*** Unit test passed ***")

print("*** Unit test: check that protein_sequence_parser() correctly returns que protein sequence ***")

output = ["R", "E", "C", "H"]
sequence = "ARG-GLU-CYS-HIS"

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

output = ["R", "E", "C", "H"]
sequence = "R-E-C-H"

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

output = ["R", "E", "C", "H"]
sequence = ["R","E", "C", "H"]

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

output = ["R", "E", "C", "H"]
sequence = ["ARG","GLU", "CYS", "HIS"]

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

print("*** Unit test: check that protein_sequence_parser() raises a ValueError if a wrong residue key is provided***")

input_parameters = {"sequence":"ARG-GLU-TUR-HIS"}

np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":"A-E-E-X"}

np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":["A", "E","X"]}

np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

print("*** Unit test passed ***")