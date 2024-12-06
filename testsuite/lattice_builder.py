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
#
import numpy as np
import lib.lattice
import unittest as ut
import matplotlib
import matplotlib.pyplot as plt
import pyMBE
import espressomd

matplotlib.use("Agg") # use a non-graphic backend

pmb = pyMBE.pymbe_library(seed=42)
MPC = 4
BOND_LENGTH = 0.355 * pmb.units.nm

# Define node particle
NodeType1 = "node_type1"
NodeType2 = "node_type2"


# define monomers
BeadType1 = "carbon"
BeadType2 = "oxygen"
BeadType3 = "nitrogen"

Res1 = "res_1"
Res2 = "res_2"
Res3 = "res_3"

# Defining bonds in the hydrogel for all different pairs
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_length = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_length,
                       'k'      : generic_harmonic_constant}

class Test(ut.TestCase):
    colormap = {
        "default_linker":"green",
        "default_monomer":"blue",
        Res3: "red",
        NodeType2: "orange",
        NodeType1: "cyan",
        Res1: "yellow",
        Res2: "magenta"
        }
    
    @classmethod
    def setUpClass(cls): 
        pmb.define_particle(name=NodeType1, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name=NodeType2, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name=BeadType1, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name=BeadType2, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name=BeadType3, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_residue(
            name=Res1,
            central_bead=BeadType1,
            side_chains=[]
            )
        pmb.define_residue(
            name=Res2,
            central_bead=BeadType2,  
            side_chains=[]
            )
        pmb.define_residue(
            name=Res3,
            central_bead=BeadType3,
            side_chains=[]
            )
        pmb.define_bond(bond_type = 'harmonic',
                        bond_parameters = HARMONIC_parameters, particle_pairs = [[BeadType1, BeadType1],
                                                                                 [BeadType1, BeadType2],
                                                                                 [BeadType1, BeadType3],
                                                                                 [BeadType2, BeadType2],
                                                                                 [BeadType2, BeadType3],
                                                                                 [BeadType3, BeadType3],
                                                                                 [BeadType1, NodeType1],
                                                                                 [BeadType1, NodeType2],
                                                                                 [BeadType2, NodeType1],
                                                                                 [BeadType2, NodeType2],
                                                                                 [BeadType3, NodeType1],
                                                                                 [BeadType3, NodeType2]])

    def test_lattice_setup(self):
        
        diamond = lib.lattice.DiamondLattice(MPC, BOND_LENGTH)
        espresso_system = espressomd.System(box_l = [diamond.BOXL]*3)
        pmb.add_bonds_to_espresso(espresso_system = espresso_system)
        lattice = pmb.initialize_lattice_builder(diamond)
        sequence = [Res3, Res1, Res2, Res1]
        # build default structure
        assert len(lattice.nodes) == len(diamond.indices)
        assert len(lattice.chains) ==  0

        # this function need some work
        lattice.add_default_chains(mpc=2)
        assert len(lattice.chains) ==  len(diamond.connectivity)

        # define custom nodes
        assert lattice.get_node("[1 1 1]") == "default_linker"
        assert lattice.get_node("[0 0 0]") == "default_linker"

        pos_node1 = pmb.set_node("[1 1 1]", NodeType1, espresso_system=espresso_system)
        np.testing.assert_equal(actual = lattice.get_node("[1 1 1]"), desired = NodeType1, verbose=True)
        pos_node2 = pmb.set_node("[0 0 0]", NodeType2, espresso_system=espresso_system)
        np.testing.assert_equal(actual = lattice.get_node("[0 0 0]"), desired = NodeType2, verbose=True)
        pos_node3 = pmb.set_node("[2 2 0]", NodeType2, espresso_system=espresso_system)
        np.testing.assert_equal(actual = lattice.get_node("[2 2 0]"), desired = NodeType2, verbose=True)

        node_positions={}
        node1_label = lattice.node_labels["[1 1 1]"]
        node_positions[node1_label]=pos_node1
        node2_label = lattice.node_labels["[0 0 0]"]
        node_positions[node2_label]=pos_node2
        node3_label = lattice.node_labels["[2 2 0]"]
        node_positions[node3_label]=pos_node3

        # define molecule in forward direction
        molecule_name = "chain_[1 1 1]_[0 0 0]"
        pmb.define_molecule(name=molecule_name, residue_list=sequence)
        pmb.set_chain("[1 1 1]", "[0 0 0]", node_positions, espresso_system=espresso_system)
        np.testing.assert_equal(actual = lattice.get_chain("[1 1 1]", "[0 0 0]"), desired = sequence, verbose=True)
        np.testing.assert_equal(actual = lattice.get_chain("[0 0 0]", "[1 1 1]"), desired = sequence[::-1], verbose=True)

        # define custom chain in reverse direction
        molecule_name = "chain_[0 0 0]_[1 1 1]"
        pmb.define_molecule(name=molecule_name, residue_list=sequence)
        pmb.set_chain("[0 0 0]", "[1 1 1]", node_positions, espresso_system=espresso_system)
        np.testing.assert_equal(lattice.get_chain("[1 1 1]", "[0 0 0]"), sequence[::-1])
        np.testing.assert_equal(lattice.get_chain("[0 0 0]", "[1 1 1]"), sequence)

        ####---Raise Exceptions---####
        # define custom chain between normally unconnected nodes
        molecule_name = "chain_[0 0 0]_[2 2 0]"
        pmb.define_molecule(name=molecule_name, residue_list=sequence)
        np.testing.assert_raises(AssertionError, 
                         pmb.set_chain, 
                         "[0 0 0]", "[2 2 0]", node_positions, espresso_system=espresso_system)

        # define custom chain that loops
        molecule_name = "chain_[0 0 0]_[0 0 0]"
        pmb.define_molecule(name=molecule_name, residue_list=sequence)
        np.testing.assert_raises(AssertionError,
                         pmb.set_chain,
                         "[0 0 0]", "[0 0 0]", node_positions, espresso_system=espresso_system)

        lattice.set_colormap(self.colormap)
        for index, (label, color) in enumerate(self.colormap.items()):
            np.testing.assert_equal(actual = lattice.get_monomer_color(label),desired = color, verbose=True)
            np.testing.assert_equal(actual = lattice.get_monomer_color_index(label),desired = index, verbose=True)


        # Test invalid operations
        with self.assertRaisesRegex(RuntimeError, "monomer 'unknown' has no associated color in the colormap"):
            lattice.get_monomer_color("unknown")
        with self.assertRaises(AssertionError):
            lattice.set_colormap("red")

        # Test node operations
        with self.assertRaisesRegex(AssertionError, r"node '\[0 5 13\]' doesn't exist in a diamond lattice"):
            lattice.get_node("[0 5 13]")

        # Test plot
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection="3d", computed_zorder=False)
        lattice.draw_lattice(ax)
        lattice.draw_simulation_box(ax)
        plt.close(fig)

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection="3d", computed_zorder=False)
        lattice.set_colormap(self.colormap)
        lattice.draw_lattice(ax)
        lattice.draw_simulation_box(ax)
        ax.legend()
        plt.close(fig)

if __name__ == "__main__":
    ut.main()

