#
# Copyright (C) 2024-2025 pyMBE-dev team
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

import numpy as np
import random
import unittest as ut
import pyMBE
from lib.lattice import DiamondLattice
import espressomd


pmb = pyMBE.pymbe_library(seed=42)

# Define node particle
NodeType = "node_type"
pmb.define_particle(name=NodeType, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

CounterIon = "counter_ion"
pmb.define_particle(name=CounterIon, sigma=0.5*pmb.units.nm, epsilon=1.5*pmb.units("reduced_energy"))

# define monomers
BeadType1 = "C"
pmb.define_particle(name=BeadType1, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
BeadType2 = "M"
pmb.define_particle(name=BeadType2, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

Res1 = "res_1"
pmb.define_residue(
    name=Res1,  # Name of the residue
    central_bead=BeadType1,  # Define the central bead name
    side_chains=[]  # Assuming no side chains for the monomer
)

Res2 = "res_2"
pmb.define_residue(
    name=Res2,  # Name of the residue
    central_bead=BeadType2,  # Define the central bead name
    side_chains=[]  # Assuming no side chains for the monomer
)

molecule_name = 'alternating_residue'
pmb.define_molecule(name=molecule_name, residue_list = [Res1, Res2, Res1, Res2, Res1, Res2])

# define bond parameters
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_length = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_length,
                      'k'      : generic_harmonic_constant}

pmb.define_bond(bond_type = 'harmonic',
                        bond_parameters = HARMONIC_parameters, particle_pairs = [[BeadType1, BeadType1],
                                                                                 [BeadType1, BeadType2],
                                                                                 [BeadType2, BeadType2]])
pmb.define_bond(bond_type = 'harmonic',
                        bond_parameters = HARMONIC_parameters, particle_pairs = [[NodeType, BeadType1],
                                                                                 [NodeType, BeadType2]])

print("*** Unit Test: check that only non-negative values of monomers per chain are allowed ***")
MPC=0
np.testing.assert_raises(ValueError, DiamondLattice, MPC, generic_bond_length)
MPC="invalid"
np.testing.assert_raises(ValueError, DiamondLattice, MPC, generic_bond_length)
MPC=-5
np.testing.assert_raises(ValueError, DiamondLattice, MPC, generic_bond_length)
print("*** Unit Test passed ***")

print("*** Unit test: check that any objects are other than DiamondLattice passed to initialize_lattice_builder raises a TypeError ***")
np.testing.assert_raises(TypeError, pmb.initialize_lattice_builder, None)
print("*** Unit test passed ***")


MPC=8
diamond_lattice = DiamondLattice(MPC, generic_bond_length)
box_l = diamond_lattice.BOXL
espresso_system = espressomd.System(box_l = [box_l]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

pmb.create_particle(name=CounterIon,
                    espresso_system=espresso_system,
                    number_of_particles=1,
                    position=[[np.random.uniform(0,box_l)]*3])

pmb.create_molecule(name=molecule_name,
                        number_of_molecules=1,
                        espresso_system=espresso_system,
                        use_default_bond=False,
                        list_of_first_residue_positions = [[np.random.uniform(0,box_l)]*3])

# Setting up node topology
indices = diamond_lattice.indices
node_topology = {}

for index in range(len(indices)):
    node_topology[index]={"particle_name": NodeType,
                          "lattice_index": indices[index]}

# Setting up chain topology
#connectivity = diamond_lattice.connectivity
node_labels = lattice_builder.node_labels
chain_labels = lattice_builder.chain_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
#connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j]) for i, j in connectivity}
chain_topology = {}
residue_list = [Res1]*(MPC//2) + [Res2]*(MPC//2)
#chain_id = 0
for chain_data in chain_labels.items():
    chain_label = chain_data[1]
    node_label_pair = chain_data[0]
    node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
    chain_topology[chain_label]={'node_start':reverse_node_labels[node_label_s],
                              'node_end': reverse_node_labels[node_label_e],
                              'residue_list':residue_list}

last_chain_id = max(chain_topology.keys())
chain_topology[last_chain_id]['residue_list'] = chain_topology[last_chain_id]['residue_list'][::-1]

#last_chain_id = max(chain_topology.keys())
#chain_topology[last_chain_id]['residue_list'] = ["res_1" if i % 2 == 0 else "res_2" for i in range(len(residue_list))]
################################################
incomplete_node_map = {
        0: {"particle_name": NodeType, "lattice_index": [0, 0, 0]},
        1: {"particle_name": NodeType, "lattice_index": [1, 1, 1]},
    }
incomplete_chain_map = {0: {"node_start": "[0 0 0]", "node_end":"[1 1 1]" , "residue_list": residue_list}}

np.testing.assert_raises(ValueError, pmb.define_hydrogel, "test_hydrogel", incomplete_node_map, chain_topology)

np.testing.assert_raises(ValueError, pmb.define_hydrogel, "test_hydrogel", node_topology, incomplete_chain_map)

def compare_node_maps(map1, map2):
    # Ensure keys are comparable
    map1 = {int(k): v for k, v in map1.items()}
    map2 = {int(k): v for k, v in map2.items()}

    # Compare key sets
    if set(map1.keys()) != set(map2.keys()):
        return False

    # Compare each node's details
    for key in map1:
        if map1[key]["particle_name"] != map2[key]["particle_name"]:
            return False
        if not np.array_equal(np.array(map1[key]["lattice_index"]), np.array(map2[key]["lattice_index"])):
            return False
    return True

def parse_string_to_array(string):
    """
    Convert a string representation of a list (e.g., '[3 1 3]') into a numpy array.
    """
    string = string.strip("[]")  # Remove brackets
    elements = map(int, string.split())  # Split by spaces and convert to integers
    return np.array(list(elements))

def compare_chain_maps(map1, map2):
    # Normalize keys to integers
    map1 = {int(k): v for k, v in map1.items()}
    map2 = {int(k): v for k, v in map2.items()}

    # Compare key sets
    if set(map1.keys()) != set(map2.keys()):
        return False

    # Compare each chain's details
    for key in map1:
        # Compare node_start and node_end
        node_start_1 = parse_string_to_array(map1[key]["node_start"])
        node_start_2 = parse_string_to_array(map2[key]["node_start"])
        if not np.array_equal(node_start_1, node_start_2):
            return False

        node_end_1 = parse_string_to_array(map1[key]["node_end"])
        node_end_2 = parse_string_to_array(map2[key]["node_end"])
        if not np.array_equal(node_end_1, node_end_2):
            return False

        # Compare residue_list
        if map1[key]["residue_list"] != map2[key]["residue_list"]:
            return False
    return True

#######################################################
pmb.define_hydrogel("my_hydrogel",node_topology, chain_topology)
existing_hydrogel_name = "my_hydrogel"
pmb.define_hydrogel(existing_hydrogel_name,node_topology, chain_topology)
hydrogel_count = len(pmb.df[pmb.df["name"] == existing_hydrogel_name])
assert hydrogel_count == 1, f"Hydrogel '{existing_hydrogel_name}' should not be redefined."
assert existing_hydrogel_name in pmb.df["name"].values
assert pmb.df.loc[pmb.df["name"] == existing_hydrogel_name, "pmb_type"].values[0] == "hydrogel"

# Verify node_map and chain_map are correctly added
assert compare_node_maps(pmb.df.loc[pmb.df["name"] == existing_hydrogel_name, "node_map"].values[0], node_topology)
assert compare_chain_maps(pmb.df.loc[pmb.df["name"] == existing_hydrogel_name, "chain_map"].values[0], chain_topology)
for chain_id in chain_topology:
    molecule_name = f"chain_{chain_topology[chain_id]['node_start']}_{chain_topology[chain_id]['node_end']}"
    assert molecule_name in pmb.df["name"].values

# Creating hydrogel
hydrogel_info = pmb.create_hydrogel("my_hydrogel", espresso_system)

print("*** Hydrogel created: Unit test to verify their name and positions ***")
############################

class Test(ut.TestCase):
 
    def test_format_node(self):
        assert pmb.format_node([1, 2, 3]) == "[1 2 3]"
        assert pmb.format_node([4, 5, 6]) == "[4 5 6]"

    def test_hydrogel_info(self):
        assert hydrogel_info["name"] == "my_hydrogel"

    def test_node_positions(self):
        for _, node_id in hydrogel_info["nodes"].items():
            node_pos = espresso_system.part.by_id(int(node_id[0])).pos
            node_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "particle") & (pmb.df["particle_id"] == node_id[0])]["name"].values[0]
            node_label = node_labels[pmb.format_node(list((node_pos*(4/lattice_builder.BOXL)).astype(int)))]
            node_data = node_topology[node_label]
            node_name = node_data["particle_name"] 
            # Assert node's name and position are correctly set
            np.testing.assert_equal(node_name_in_espresso, node_name)
            np.testing.assert_allclose(node_pos, np.array(node_data["lattice_index"]) * 0.25 * diamond_lattice.BOXL, atol=1e-7)

    def test_chain_placement_and_connectivity(self):
        for molecule_id, molecule_data in hydrogel_info["chains"].items():
            # Ensure that chain's node_start and node_end are correctly set
            node_start = molecule_data["node_start"]
            node_end = molecule_data["node_end"]
            chain_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "molecule") & (pmb.df["molecule_id"] == molecule_id)]["name"].values[0]
            # Assert chain's node_start and node_end
            np.testing.assert_equal(chain_name_in_espresso, f"chain_{node_start}_{node_end}")
            # Check if chain is connected in the espresso system (e.g., check bond or distance between node_start and node_end)
            node_start_id = hydrogel_info["nodes"][node_start][0]
            node_end_id = hydrogel_info["nodes"][node_end][0]
            start_pos = espresso_system.part.by_id(int(node_start_id)).pos
            end_pos = espresso_system.part.by_id(int(node_end_id)).pos
            vec_between_nodes = end_pos - start_pos
            # Ensure that the chain is connected (check distance, should be within acceptable bond length range)
            vec_between_nodes = vec_between_nodes -  diamond_lattice.BOXL * np.round(vec_between_nodes / diamond_lattice.BOXL)
            distance_between_nodes = np.linalg.norm(vec_between_nodes)
            np.testing.assert_allclose(distance_between_nodes, (diamond_lattice.MPC+1)*generic_bond_length.magnitude, atol=0.0000001)
    
    def test_all_residue_placement(self):
        for _, chain_data in hydrogel_info["chains"].items():
            node_start = chain_data["node_start"]
            node_end = chain_data["node_end"]
            node_start_label = lattice_builder.node_labels[node_start]
            node_end_label = lattice_builder.node_labels[node_end]
            node_start_pos = np.array([float(x) for x in node_start.strip('[]').split()]) * 0.25 * lattice_builder.BOXL
            node_end_pos = np.array([float(x) for x in node_end.strip('[]').split()]) * 0.25 * lattice_builder.BOXL
            node_start_id = espresso_system.part.select(lambda p: (p.pos == node_start_pos).all()).id[0]
            node_start_name = pmb.df[(pmb.df["particle_id"]==node_start_id) & (pmb.df["pmb_type"]=="particle")]["name"].values[0]
            prev_particle_name = node_start_name 
            vec_between_nodes = (np.array([float(x) for x in node_end.strip('[]').split()]) -
                                     np.array([float(x) for x in node_start.strip('[]').split()])) * 0.25 * lattice_builder.BOXL
            vec_between_nodes = vec_between_nodes - lattice_builder.BOXL * np.round(vec_between_nodes / lattice_builder.BOXL)
            backbone_vector = vec_between_nodes / (diamond_lattice.MPC + 1)
   
            # Get the list of residues (keys in chain_data, excluding node_start/node_end)
            residues_in_chain = {k: v for k, v in chain_data.items() if isinstance(k, int)}
             
            # Loop through each residue in the chain
            for (res_id, res_data) in residues_in_chain.items():
                central_bead_id = res_data["central_bead_id"]
                
                # Get the position of the central bead from the espresso system
                central_bead_pos = espresso_system.part.by_id(central_bead_id).pos

                # Calculate the expected position of the residue's central bead
                residue_index = list(residues_in_chain.keys()) .index(res_id)
                current_part_name = pmb.df[(pmb.df["particle_id"]==central_bead_id) & (pmb.df["pmb_type"]=="particle")]["name"].values[0] 
                bond_length = pmb.get_bond_length(prev_particle_name, current_part_name)
                expected_position = np.array([float(x) for x in node_start.strip('[]').split()]) * 0.25 * diamond_lattice.BOXL + (residue_index + 1) * backbone_vector * bond_length
                
                # Validate that the central bead's position matches the expected position
                np.testing.assert_allclose(central_bead_pos, expected_position, atol=1e-7)
                chain_label_key = "("+str(node_start_label)+", "+str(node_end_label)+")"
                chain_label = chain_labels[chain_label_key]
                expected_res_name = chain_topology[chain_label]["residue_list"][residue_index]
                expected_node_start = chain_topology[chain_label]["node_start"]
                expected_node_end = chain_topology[chain_label]["node_end"]
                residue_name = pmb.df[(pmb.df["pmb_type"]=="residue") & (pmb.df["residue_id"]==res_id)]["name"].values[0]
                np.testing.assert_equal(node_start, expected_node_start)
                np.testing.assert_equal(node_end, expected_node_end)
                np.testing.assert_equal(residue_name, expected_res_name)
                prev_particle_name = current_part_name

if __name__ == "__main__":
    ut.main(exit=False)

#####-- Invalid hydrogel name --#####
# Test if create_hydrogel raises an exception when provided with invalid data
print("*** Unit Test: Check invalid inputs for create_hydrogel ***")
np.testing.assert_raises(ValueError, pmb.create_hydrogel, "invalid_hydrogel", espresso_system)
print("*** Invalid Input Test Passed ***")

# Check if the molecules (chains) are correctly stored in the hydrogel data
for ((molecule_id, molecule_data),chain_id_in_input_dict) in zip(hydrogel_info["chains"].items(),chain_topology.items()):
    molecule_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "molecule") & (pmb.df["molecule_id"] == molecule_id)]["name"].values[0]
    np.testing.assert_equal(molecule_name_in_espresso, f"chain_{molecule_data['node_start']}_{molecule_data['node_end']}")
    Residue_list_in_espresso = pmb.df[(pmb.df["pmb_type"]=="molecule") & (pmb.df["molecule_id"]==molecule_id)]["residue_list"].values[0]
    residue_list = chain_id_in_input_dict[1]["residue_list"]
    np.testing.assert_equal(Residue_list_in_espresso, residue_list)

print("*** Checking if the ends of the randomly chosen chain is connected to node_start and node_end ***")
random_chain_id = random.choice(list(chain_topology.keys())) # Choose a random chain of 0-15
molecule_random = hydrogel_info["chains"][random_chain_id]
numeric_keys = {key: value for key, value in molecule_random.items() if isinstance(key, int)}

# Extract the first and last elements
keys = list(numeric_keys.keys())
first_key = keys[0]
last_key = keys[-1]

Res_node_start = numeric_keys[first_key]
Res_node_end = numeric_keys[last_key]

central_bead_near_node_start = Res_node_start["central_bead_id"]
central_bead_near_node_end = Res_node_end["central_bead_id"]

#node_ids = [0,1,2,3,4,5,6,7]
node_ids = []
for indice in node_labels.keys():
    index_pos = np.array(list(int(x) for x in indice.strip('[]').split()))*0.25*lattice_builder.BOXL
    node_id = espresso_system.part.select(lambda p: (p.pos == index_pos).all()).id[0]
    node_ids.append(node_id)

bead_ids_in_random_molecule = [i for i in range(central_bead_near_node_start, central_bead_near_node_end+1)]
#print(pmb.df[pmb.df["molecule_id"]==random_chain_id])
filtered_df = pmb.df[
    pmb.df["particle_id"].isin(node_ids) & 
    pmb.df["particle_id2"].isin(bead_ids_in_random_molecule)
    ]

# Extract scalar values for central_bead_node_start and central_bead_node_end
central_bead_node_start = filtered_df[filtered_df["particle_id2"] == central_bead_near_node_start]["particle_id"].iloc[0]
central_bead_node_end = filtered_df[filtered_df["particle_id2"] == central_bead_near_node_end]["particle_id"].iloc[0]

bond_name_node_start = filtered_df[
    (filtered_df["particle_id"] == central_bead_node_start) & 
    (filtered_df["particle_id2"] == central_bead_near_node_start)
]["name"].iloc[0]

bond_name_node_end = filtered_df[
    (filtered_df["particle_id"] == central_bead_node_end) & 
    (filtered_df["particle_id2"] == central_bead_near_node_end)
]["name"].iloc[0]

for _, row in filtered_df.iterrows():
    bond_object = row["bond_object"]
    if bond_object is None:
        raise ValueError("Bond object is not defined near nodes")

central_bead_name_near_node_start = pmb.df[pmb.df["particle_id"]==central_bead_near_node_start]["name"].values[0]
central_bead_name_near_node_end = pmb.df[pmb.df["particle_id"]==central_bead_near_node_end]["name"].values[0]

if central_bead_name_near_node_start == BeadType1:
    possible_bond_names = [NodeType+"-"+BeadType1, BeadType1+"-"+NodeType]
    assert bond_name_node_start in possible_bond_names

elif central_bead_name_near_node_start == BeadType2:
    possible_bond_names = [NodeType+"-"+BeadType2, BeadType2+"-"+NodeType]
    assert bond_name_node_start in possible_bond_names

if central_bead_name_near_node_end == BeadType1:
    possible_bond_names = [NodeType+"-"+BeadType1, BeadType1+"-"+NodeType]
    assert bond_name_node_end in possible_bond_names

elif central_bead_name_near_node_end == BeadType2:
    possible_bond_names = [NodeType+"-"+BeadType2, BeadType2+"-"+NodeType]
    assert bond_name_node_end in possible_bond_names

print("*** Unit Test passed ***")
