import numpy as np
import random
import unittest as ut
import pyMBE
from lib.lattice import DiamondLattice
import espressomd
from unittest.mock import patch
import io

pmb = pyMBE.pymbe_library(seed=42)
MPC = 0 
BOND_LENGTH = 0.355 * pmb.units.nm

print("*** Unit test: check that any objects are other than DiamondLattice passed to initialize_lattice_builder raises a NameError ***")
np.testing.assert_raises(TypeError, pmb.initialize_lattice_builder, None)
print("*** Unit test passed ***")
# Define node particle
NodeType = "node_type"
pmb.define_particle(name=NodeType, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
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

with ut.TestCase().assertRaisesRegex(ValueError, "MPC must be a non-zero positive integer."):
    diamond_lattice = DiamondLattice(MPC, generic_bond_length)

with ut.TestCase().assertRaisesRegex(ValueError, "MPC must be a non-zero positive integer."):
    diamond_lattice = DiamondLattice("invalid", generic_bond_length)  # MPC is a string

with ut.TestCase().assertRaisesRegex(ValueError, "MPC must be a non-zero positive integer."):
    diamond_lattice = DiamondLattice(-5, generic_bond_length)  # MPC is negative

print("*** Unit Test passed ***")

MPC=8
diamond_lattice = DiamondLattice(MPC, generic_bond_length)
espresso_system = espressomd.System(box_l = [diamond_lattice.BOXL]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

# Setting up node topology
indices = diamond_lattice.indices
node_topology = {}

for index in range(len(indices)):
    node_topology[index]={"particle_name": NodeType,
                          "lattice_index": indices[index]}

# Setting up chain topology
connectivity = diamond_lattice.connectivity
node_labels = lattice_builder.node_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j]) for i, j in connectivity}
chain_topology = {}
residue_list = [Res1]*(MPC//2) + [Res2]*(MPC//2)
chain_id = 0
for node_s, node_e in connectivity_with_labels:
    chain_topology[chain_id]={'node_start':node_s,
                              'node_end': node_e,
                              'residue_list':residue_list}
    chain_id+=1

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

def test_incomplete_node_map():
    #node_map = {0: {"particle_name": "node0", "lattice_index": [0, 0, 0]}}
    #chain_map = {i: {"node_start": str(i), "node_end": str(i+1), "residue_list": ["A", "B", "C"]} for i in range(16)}

    with patch('sys.stdout', new=io.StringIO()) as mock_stdout:
        pmb.define_hydrogel(name="test_hydrogel", node_map = incomplete_node_map, chain_map = chain_topology)
        output = mock_stdout.getvalue().strip()

    expected_output = (
        "Incomplete hydrogel: node_map must contain exactly 8 node_labels (0 through 7). "
        "Found: [0, 1]"
    )
    assert output == expected_output, f"Unexpected print output: {output}"


def test_incomplete_chain_map():
    #node_map = {i: {"particle_name": f"node{i}", "lattice_index": [i, i, i]} for i in range(8)}
    #chain_map = {i: {"node_start": str(i), "node_end": str(i+1), "residue_list": ["A", "B", "C"]} for i in range(10)}

    with patch('sys.stdout', new=io.StringIO()) as mock_stdout:
        pmb.define_hydrogel(name="test_hydrogel", node_map = node_topology, chain_map = incomplete_chain_map)
        output = mock_stdout.getvalue().strip()

    expected_output = (
        "Incomplete hydrogel: chain_map must contain exactly 16 chain_ids (0 through 15). "
        "Found: [0]"
    )
    assert output == expected_output, f"Unexpected print output: {output}"
test_incomplete_node_map()
test_incomplete_chain_map()
##################################################################
#pmb.define_hydrogel("test_hydrogel", incomplete_node_map, chain_map)
#pmb.create_hydrogel("test_hydrogel", espresso_system)
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
#def test_hydrogel_already_defined():
existing_hydrogel_name = "my_hydrogel"
pmb.define_hydrogel(existing_hydrogel_name,node_topology, chain_topology)
hydrogel_count = len(pmb.df[pmb.df["name"] == existing_hydrogel_name])
assert hydrogel_count == 1, f"Hydrogel '{existing_hydrogel_name}' should not be redefined."
#test_hydrogel_already_defined()
assert existing_hydrogel_name in pmb.df["name"].values
assert pmb.df.loc[pmb.df["name"] == existing_hydrogel_name, "pmb_type"].values[0] == "hydrogel"

# Verify node_map and chain_map are correctly added
assert compare_node_maps(pmb.df.loc[pmb.df["name"] == existing_hydrogel_name, "node_map"].values[0], node_topology)
assert compare_chain_maps(pmb.df.loc[pmb.df["name"] == existing_hydrogel_name, "chain_map"].values[0], chain_topology)
for chain_id in chain_topology:
    molecule_name = f"chain_{chain_topology[chain_id]['node_start']}_{chain_topology[chain_id]['node_end']}"
    assert molecule_name in pmb.df["name"].values

print("*** Unit Test: check node map and chain map parameters ***")

#random_chain_id = random.choice(list(chain_topology.keys())) # Choose a random chain of 0-15
# extract node_start, node_end and residue_list of random chain
#chain_data = chain_topology[random_chain_id] 
#node_start = chain_data["node_start"] 
#node_end = chain_data["node_end"]
#residue_list = chain_data['residue_list']
# choosing random residue in the residue_list
#random_res_in_res_list = random.choice(list(residue_list))
#res_indices = [i for i, residue in enumerate(residue_list) if residue == random_res_in_res_list]
#chosen_index = random.choice(res_indices)
# Expected res_id of the randomly chosen residue
#random_res_id = random_chain_id*(len(residue_list)) + chosen_index

#hydrogel_data = pmb.df[pmb.df["name"] == "my_hydrogel"]
#stored_chain_map = hydrogel_data["chain_map"].values[0]

#ut.TestCase().assertEqual(stored_chain_map[str(random_chain_id)]["node_start"], node_start)
#ut.TestCase().assertEqual(stored_chain_map[str(random_chain_id)]["node_end"], node_end)
#ut.TestCase().assertEqual(stored_chain_map[str(random_chain_id)]["residue_list"],residue_list)
#molecule_name = f"chain_{node_start}_{node_end}"
#molecule_data = pmb.df[
#    (pmb.df["pmb_type"] == "molecule") & (pmb.df["name"] == molecule_name)
#]

# Ensure that exactly one molecule matches
#ut.TestCase().assertEqual(len(molecule_data), 1, f"Molecule {molecule_name} not found in pmb.df")

#random_node = random.choice(list(node_topology.keys()))
#node_data = node_topology[random_node]
#node_name = node_data['particle_name']
#node_index = node_data['lattice_index']

#stored_node_map = hydrogel_data["node_map"].values[0]
#ut.TestCase().assertEqual(stored_node_map[str(random_node)]["particle_name"], node_name)
#ut.TestCase().assertEqual(all(stored_node_map[str(random_node)]["lattice_index"]), all(node_index))

print("*** Unit Test passed ***")

# Creating hydrogel
hydrogel_info = pmb.create_hydrogel("my_hydrogel", espresso_system)
print("*** Hydrogel created: Unit test to verify their name and positions ***")

############################
def test_format_node():
    assert pmb.format_node([1, 2, 3]) == "[1 2 3]"
    assert pmb.format_node([4, 5, 6]) == "[4 5 6]"
test_format_node()
############################

def test_hydrogel_info():
    assert hydrogel_info["name"] == "my_hydrogel"

test_hydrogel_info()
###########################
def test_node_positions():
    for node_index, node_id in hydrogel_info["nodes"].items():
        # Ensure that the node is added to the espresso system and has a position
        node_pos = espresso_system.part.by_id(int(node_id[0])).pos
        node_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "particle") & (pmb.df["particle_id"] == node_id[0])]["name"].values[0]
        node_data = node_topology[node_id[0]]
        node_name = node_data["particle_name"]
    
        # Assert node's name and position are correctly set
        ut.TestCase().assertEqual(node_name_in_espresso, node_name)
        ut.TestCase().assertTrue(np.allclose(node_pos, np.array(node_data["lattice_index"]) * 0.25 * diamond_lattice.BOXL, atol=1e-7))
test_node_positions()
#############################
def test_chain_placement_and_connectivity():
    for molecule_id, molecule_data in hydrogel_info["chains"].items():
        # Ensure that chain's node_start and node_end are correctly set
        node_start = molecule_data["node_start"]
        node_end = molecule_data["node_end"]
        chain_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "molecule") & (pmb.df["molecule_id"] == molecule_id)]["name"].values[0]

        # Assert chain's node_start and node_end
        ut.TestCase().assertEqual(chain_name_in_espresso, f"chain_{node_start}_{node_end}")

        # Check if chain is connected in the espresso system (e.g., check bond or distance between node_start and node_end)
        node_start_id = hydrogel_info["nodes"][node_start][0]
        node_end_id = hydrogel_info["nodes"][node_end][0]
        start_pos = espresso_system.part.by_id(int(node_start_id)).pos
        end_pos = espresso_system.part.by_id(int(node_end_id)).pos
        vec_between_nodes = end_pos - start_pos
        # Ensure that the chain is connected (check distance, should be within acceptable bond length range)
        vec_between_nodes = vec_between_nodes -  diamond_lattice.BOXL * np.round(vec_between_nodes / diamond_lattice.BOXL)
        distance_between_nodes = np.linalg.norm(vec_between_nodes)
        ut.TestCase().assertAlmostEqual(distance_between_nodes, (diamond_lattice.MPC+1)*generic_bond_length.magnitude, delta=0.0000001)
test_chain_placement_and_connectivity()
########################################
def test_all_residue_placement():
    for chain_id, chain_data in hydrogel_info["chains"].items():
        node_start = chain_data["node_start"]
        node_end = chain_data["node_end"]

        # Get the list of residues (keys in chain_data, excluding node_start/node_end)
        residues_in_chain = {k: v for k, v in chain_data.items() if isinstance(k, int)}

        # Loop through each residue in the chain
        for (res_id, res_data) in residues_in_chain.items():
            central_bead_id = res_data["central_bead_id"]

            # Get the position of the central bead from the espresso system
            central_bead_pos = espresso_system.part.by_id(central_bead_id).pos

            # Calculate the expected position of the residue's central bead
            residue_index = list(residues_in_chain.keys()).index(res_id)

            vec_between_nodes = (np.array([float(x) for x in node_end.strip('[]').split()]) - 
                                 np.array([float(x) for x in node_start.strip('[]').split()])) * 0.25 * diamond_lattice.BOXL
            vec_between_nodes = vec_between_nodes - diamond_lattice.BOXL * np.round(vec_between_nodes / diamond_lattice.BOXL)
            backbone_vector = vec_between_nodes / (diamond_lattice.MPC + 1)
            expected_position = (np.array([float(x) for x in node_start.strip('[]').split()]) * 0.25 * diamond_lattice.BOXL + 
                                 (residue_index + 1) * backbone_vector)
            
            # Validate that the central bead's position matches the expected position
            np.allclose(central_bead_pos, expected_position, atol=1e-2)
            expected_res_name = chain_topology[chain_id]["residue_list"][residue_index]
            expected_node_start = chain_topology[chain_id]["node_start"]
            expected_node_end = chain_topology[chain_id]["node_end"]
            residue_name = pmb.df[(pmb.df["pmb_type"]=="residue") & (pmb.df["residue_id"]==res_id)]["name"].values[0]
            ut.TestCase().assertEqual(node_start, expected_node_start)
            ut.TestCase().assertEqual(node_end, expected_node_end)
            ut.TestCase().assertEqual(residue_name, expected_res_name)
test_all_residue_placement()
###############################################################################################################################

#####-- Invalid hydrogel name --#####
# Test if create_hydrogel raises an exception when provided with invalid data
print("*** Unit Test: Check invalid inputs for create_hydrogel ***")

with ut.TestCase().assertRaises(ValueError):
    pmb.create_hydrogel("invalid_hydrogel", espresso_system)  # No such hydrogel in df

print("*** Invalid Input Test Passed ***")
#######################################################################

# Check if the molecules (chains) are correctly stored in the hydrogel data
for ((molecule_id, molecule_data),chain_id_in_input_dict) in zip(hydrogel_info["chains"].items(),chain_topology.items()):
    molecule_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "molecule") & (pmb.df["molecule_id"] == molecule_id)]["name"].values[0]
    ut.TestCase().assertEqual(molecule_name_in_espresso, f"chain_{molecule_data['node_start']}_{molecule_data['node_end']}")
    Residue_list_in_espresso = pmb.df[(pmb.df["pmb_type"]=="molecule") & (pmb.df["molecule_id"]==molecule_id)]["residue_list"].values[0]
    residue_list = chain_id_in_input_dict[1]["residue_list"]
    ut.TestCase().assertEqual(Residue_list_in_espresso, residue_list)

########################################################################


#Node_name_in_espresso = pmb.df[(pmb.df["pmb_type"]=="particle") & (pmb.df["particle_id"]==random_node)]["name"].values[0]
#ut.TestCase().assertEqual(Node_name_in_espresso, node_name)
#ut.TestCase().assertEqual(all(espresso_system.part.by_id(random_node).pos), all(node_index*0.25*diamond_lattice.BOXL))

#Chain_name_in_espresso = pmb.df[(pmb.df["pmb_type"]=="molecule") & (pmb.df["molecule_id"]==random_chain_id)]["name"].values[0] 
#Residue_list_in_espresso = pmb.df[(pmb.df["pmb_type"]=="molecule") & (pmb.df["molecule_id"]==random_chain_id)]["residue_list"].values[0]
#Residue_random_name_in_espresso = pmb.df[(pmb.df["pmb_type"]=="residue") & (pmb.df["residue_id"]==random_res_id)]["name"].values[0]
#Residue_random_mol_id = pmb.df[(pmb.df["name"]==Residue_random_name_in_espresso) & (pmb.df["residue_id"]==random_res_id)]["molecule_id"].values[0]

#ut.TestCase().assertEqual(Chain_name_in_espresso, "chain_"+node_start+"_"+node_end)
#ut.TestCase().assertEqual(Residue_list_in_espresso, residue_list)
#ut.TestCase().assertEqual(Residue_random_name_in_espresso,  random_res_in_res_list)

#vec_between_nodes = (np.array(list(int(x) for x in node_end.strip('[]').split())) -  np.array(list(int(x) for x in node_start.strip('[]').split())))*0.25*diamond_lattice.BOXL
#vec_between_nodes = vec_between_nodes - diamond_lattice.BOXL * np.round(vec_between_nodes/diamond_lattice.BOXL)
#backbone_vector = np.array(list(vec_between_nodes/(diamond_lattice.MPC + 1)))
#random_res_pos_central_bead =  np.array(list(int(x) for x in node_start.strip('[]').split()))*0.25*diamond_lattice.BOXL + (chosen_index + 1)*backbone_vector
#random_res_central_bead_id = hydrogel_info["chains"][Residue_random_mol_id][random_res_id]["central_bead_id"]

#np.allclose(espresso_system.part.by_id(random_res_central_bead_id).pos, random_res_pos_central_bead, atol=1e-2)
#print("*** Unit Test passed ***")

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

node_ids = [0,1,2,3,4,5,6,7]
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
