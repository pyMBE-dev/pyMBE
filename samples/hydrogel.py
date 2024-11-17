import pyMBE
import espressomd
import numpy as np
from tqdm import tqdm
import re
from pathlib import Path
import matplotlib.pyplot as plt
import json
import pandas as pd
from lib.handy_functions import minimize_espresso_system_energy
from lib.handy_functions import setup_langevin_dynamics
from lib.analysis import block_analyze
from lib.lattice import DiamondLattice

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)
# Creating instance from diamond lattice 
reduced_unit_set = pmb.get_reduced_units()
# Monomers per chain
MPC = 8
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

residue_list = [Res1]*(MPC//2) + [Res2]*(MPC//2)

# Defining bonds in the hydrogel for all different pairs
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

# Provide MPC and BOND_LENGTH to Diamond Lattice
diamond_lattice = DiamondLattice(MPC,generic_bond_length)
# Get the box length from the object of DiamondLattice 
# and use it to create espresso system
espresso_system = espressomd.System(box_l = [diamond_lattice.BOXL]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)
# Currently pass the pyMBE object into lattice builder
# so you have only one instance of the pyMBE by setting 
# self.pmb = pmb in constructor of the LatticeBuilder
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice) # Dont need pmb when integrated to pyMBE.py

######---creating a hydrogel---###########
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
# reverse_node_labels looks like {0:"[0 0 0]",1:"[1 1 1]", ..}
connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j]) for i, j in connectivity}
chain_topology = {}

def parse_node(node_str):
    return list(map(int, node_str.strip("[]").split()))

chain_id = 0
for node_s, node_e in connectivity_with_labels:
    chain_topology[chain_id]={'node_start':node_s,
                              'node_end': node_e,
                              'residue_list':residue_list}
    chain_id+=1

pmb.define_hydrogel("my_hydrogel",node_topology, chain_topology)
node_positions = pmb.create_hydrogel("my_hydrogel", espresso_system)
   
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
lattice_builder.draw_lattice(ax)
lattice_builder.draw_simulation_box(ax)
plt.legend(fontsize=12)
plt.show()




