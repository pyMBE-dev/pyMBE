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

import lib.lattice
import unittest as ut
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg") # use a non-graphic backend


class Test(ut.TestCase):
    colormap = {
        "default_linker": "C0",
        "default_monomer": "C1",
        "silicon": "C2",
        "carbon": "C3",
        "oxygen": "C4",
        "nitrogen": "C5",
    }

    def test_builder(self):
        diamond = lib.lattice.DiamondLattice
        lattice = lib.lattice.LatticeBuilder(diamond, strict=False)
        sequence = ["nitrogen", "carbon", "oxygen", "carbon"]
        # build default structure
        self.assertEqual(len(lattice.nodes), len(diamond.indices))
        self.assertEqual(len(lattice.chains), 0)
        lattice.add_default_chains(mpc=2)
        self.assertEqual(len(lattice.chains), len(diamond.connectivity))
        # define custom nodes
        self.assertEqual(lattice.get_node("[1 1 1]"), "default_linker")
        lattice.set_node(node="[1 1 1]", residue="silicon")
        self.assertEqual(lattice.get_node("[1 1 1]"), "silicon")
        # define custom chain in forward direction
        lattice.set_chain("[0 0 0]", "[1 1 1]", sequence=sequence)
        self.assertEqual(lattice.get_chain("[0 0 0]", "[1 1 1]"), sequence)
        self.assertEqual(lattice.get_chain("[1 1 1]", "[0 0 0]"), sequence[::-1])
        # define custom chain in reverse direction
        lattice.set_chain("[1 1 1]", "[0 0 0]", sequence=sequence)
        self.assertEqual(lattice.get_chain("[0 0 0]", "[1 1 1]"), sequence[::-1])
        self.assertEqual(lattice.get_chain("[1 1 1]", "[0 0 0]"), sequence)
        # define custom chain between normally unconnected nodes
        lattice.set_chain("[0 0 0]", "[2 2 0]", sequence=sequence)
        self.assertEqual(lattice.get_chain("[0 0 0]", "[2 2 0]"), sequence)
        self.assertEqual(lattice.get_chain("[2 2 0]", "[0 0 0]"), sequence[::-1])
        # define custom chain that loops
        lattice.set_chain("[0 0 0]", "[0 0 0]", sequence=["carbon"])
        self.assertEqual(lattice.get_chain("[0 0 0]", "[0 0 0]"), ["carbon"])
        # colors
        lattice.set_colormap(self.colormap)
        for index, (label, color) in enumerate(self.colormap.items()):
            self.assertEqual(lattice.get_monomer_color(label), color)
            self.assertEqual(lattice.get_monomer_color_index(label), index)

    def test_exceptions(self):
        diamond = lib.lattice.DiamondLattice
        lattice = lib.lattice.LatticeBuilder(diamond, strict=True)
        with self.assertRaisesRegex(RuntimeError, "monomer 'unknown' has no associated color in the colormap"):
            lattice.get_monomer_color("unknown")
        with self.assertRaises(AssertionError):
            lattice.set_colormap("red")
        with self.assertRaises(AssertionError):
            lattice.set_colormap([("oxygen", "red")])
        with self.assertRaisesRegex(AssertionError, r"node '\[0 5 13\]' doesn't exist in a diamond lattice"):
            lattice.get_node("[0 5 13]")
        with self.assertRaisesRegex(AssertionError, r"node '\[-1 0 0\]' doesn't exist in a diamond lattice"):
            lattice.set_node("[-1 0 0]", "carbon")
        with self.assertRaisesRegex(RuntimeError, r"no chain has been defined between '\[0 0 0\]' and '\[1 1 1\]' yet"):
            lattice.get_chain("[0 0 0]", "[1 1 1]")
        with self.assertRaisesRegex(AssertionError, r"node '\[0 5 13\]' doesn't exist in a diamond lattice"):
            lattice.get_chain("[0 5 13]", "[1 1 1]")
        with self.assertRaisesRegex(AssertionError, r"node '\[0 5 13\]' doesn't exist in a diamond lattice"):
            lattice.get_chain("[0 0 0]", "[0 5 13]")
        with self.assertRaisesRegex(AssertionError, r"there is no chain between '\[0 0 0\]' and '\[2 2 0\]' in a diamond lattice \(strict mode is enabled\)"):
            lattice.get_chain("[0 0 0]", "[2 2 0]")
        with self.assertRaisesRegex(AssertionError, r"there is no chain between '\[0 0 0\]' and '\[0 0 0\]' in a diamond lattice \(strict mode is enabled\)"):
            lattice.get_chain("[0 0 0]", "[0 0 0]")
        with self.assertRaises(AssertionError):
            lattice.set_chain("[0 0 0]", "[1 1 1]", sequence=[])
        with self.assertRaises(AssertionError):
            lattice.set_chain("[0 0 0]", "[1 1 1]", sequence="oxygen")
        with self.assertRaisesRegex(AssertionError, r"chain cannot be defined between '\[0 0 0\]' and '\[0 0 0\]' since it would form a loop with a non-symmetric sequence \(under-defined stereocenter\)"):
            lattice.strict = False
            lattice.set_chain("[0 0 0]", "[0 0 0]", sequence=["carbon", "oxygen"])

    def test_plot(self):
        # smoke test: check for runtime errors when drawing on a 3D canvas
        sequence = ["nitrogen", "carbon", "oxygen", "carbon"]
        diamond = lib.lattice.DiamondLattice
        lattice = lib.lattice.LatticeBuilder(diamond)
        lattice.add_default_chains(mpc=2)
        lattice.set_node(node="[1 1 1]", residue="silicon")
        lattice.set_chain("[0 0 0]", "[1 1 1]", sequence=sequence)
        # with default colormap
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection="3d", computed_zorder=False)
        lattice.draw_lattice(ax)
        lattice.draw_simulation_box(ax)
        ax.legend()
        plt.close(fig)
        # with custom colormap
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection="3d", computed_zorder=False)
        lattice.set_colormap(self.colormap)
        lattice.draw_lattice(ax)
        lattice.draw_simulation_box(ax)
        ax.legend()
        plt.close(fig)


if __name__ == "__main__":
    ut.main()
