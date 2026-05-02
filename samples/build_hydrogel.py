#
# Copyright (C) 2024-2026 pyMBE-dev team
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
# Demonstrates:
#   1. Original hydrogel build (mpc=40, two bead types, no angles) with a 3D plot.
#   2. Angle-potential hydrogel (mpc=4, three bead types).
#      Run with --include to also add crosslinker-adjacent angle templates.

import argparse
import numpy as np
import espressomd
import matplotlib.pyplot as plt

import pyMBE
from pyMBE.lib.lattice import DiamondLattice


def run_original_hydrogel(espresso_system):
    """Build and visualize the original hydrogel (mpc=40, no angles)."""
    pmb = pyMBE.pymbe_library(seed=42)
    mpc = 40
    NodeType = "node_type"
    pmb.define_particle(name=NodeType,
                        sigma=0.355*pmb.units.nm,
                        epsilon=1*pmb.units('reduced_energy'))
    BeadType1 = "C"
    pmb.define_particle(name=BeadType1,
                        sigma=0.355*pmb.units.nm,
                        epsilon=1*pmb.units('reduced_energy'))
    BeadType2 = "M"
    pmb.define_particle(name=BeadType2,
                        sigma=0.355*pmb.units.nm,
                        epsilon=1*pmb.units('reduced_energy'))

    Res1 = "res_1"
    pmb.define_residue(name=Res1, central_bead=BeadType1, side_chains=[])
    Res2 = "res_2"
    pmb.define_residue(name=Res2, central_bead=BeadType2, side_chains=[])

    residue_list = [Res1]*(mpc//2) + [Res2]*(mpc//2)
    pmb.define_molecule(name="hydrogel_chain", residue_list=residue_list)

    generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
    generic_bond_l = 0.355*pmb.units.nm
    HARMONIC_parameters = {'r_0': generic_bond_l, 'k': generic_harmonic_constant}
    pmb.define_bond(bond_type='harmonic',
                    bond_parameters=HARMONIC_parameters,
                    particle_pairs=[[BeadType1, BeadType1],
                                    [BeadType1, BeadType2],
                                    [BeadType2, BeadType2],
                                    [NodeType, BeadType1],
                                    [NodeType, BeadType2]])

    diamond_lattice = DiamondLattice(mpc, generic_bond_l)
    espresso_system.box_l = [diamond_lattice.box_l] * 3
    lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

    indices = diamond_lattice.indices
    node_topology = []
    for index in range(len(indices)):
        node_topology.append({"particle_name": NodeType,
                               "lattice_index": indices[index]})

    connectivity = diamond_lattice.connectivity
    node_labels = lattice_builder.node_labels
    reverse_node_labels = {v: k for k, v in node_labels.items()}
    connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j])
                                 for i, j in connectivity}
    chain_topology = []
    for node_s, node_e in connectivity_with_labels:
        chain_topology.append({'node_start': node_s,
                                'node_end': node_e,
                                'molecule_name': "hydrogel_chain"})

    lattice_builder.chains = chain_topology
    pmb.define_hydrogel("my_hydrogel", node_topology, chain_topology)
    pmb.create_hydrogel("my_hydrogel", espresso_system)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    lattice_builder.draw_lattice(ax=ax, pmb=pmb)
    lattice_builder.draw_simulation_box(ax)
    plt.legend(fontsize=12)
    plt.show()


def define_angle_hydrogel_system(espresso_system, include_crosslinker_angles):
    """Build the angle-potential hydrogel (mpc=4)."""
    pmb = pyMBE.pymbe_library(seed=42)

    node_type = "D"
    internal_bead_type = "A"
    chain_end_bead_type = "C"
    terminal_residue_name = "TerminalRes"
    internal_residue_name = "InternalRes"
    molecule_name = "hydrogel_chain"

    bond_length = 0.355 * pmb.units.nm
    mpc = 4

    pmb.define_particle(name=node_type,
                        sigma=bond_length,
                        epsilon=1*pmb.units("reduced_energy"))
    pmb.define_particle(name=internal_bead_type,
                        sigma=bond_length,
                        epsilon=1*pmb.units("reduced_energy"))
    pmb.define_particle(name=chain_end_bead_type,
                        sigma=bond_length,
                        epsilon=1*pmb.units("reduced_energy"))

    pmb.define_residue(name=terminal_residue_name,
                       central_bead=chain_end_bead_type,
                       side_chains=[])
    pmb.define_residue(name=internal_residue_name,
                       central_bead=internal_bead_type,
                       side_chains=[])
    pmb.define_molecule(name=molecule_name,
                        residue_list=[
                            terminal_residue_name,
                            internal_residue_name,
                            internal_residue_name,
                            terminal_residue_name,
                        ])

    harmonic_bond_parameters = {
        "r_0": bond_length,
        "k": 400 * pmb.units("reduced_energy / reduced_length**2"),
    }
    pmb.define_bond(bond_type="harmonic",
                    bond_parameters=harmonic_bond_parameters,
                    particle_pairs=[
                        [internal_bead_type, internal_bead_type],
                        [internal_bead_type, chain_end_bead_type],
                        [node_type, chain_end_bead_type],
                    ])

    angle_parameters = {
        "k": 25 * pmb.units("reduced_energy"),
        "phi_0": np.pi * pmb.units(""),
    }
    pmb.define_angle(angle_type="harmonic",
                     angle_parameters=angle_parameters,
                     particle_triplets=[(chain_end_bead_type, internal_bead_type, internal_bead_type)])

    if include_crosslinker_angles:
        pmb.define_angle(angle_type="harmonic",
                         angle_parameters=angle_parameters,
                         particle_triplets=[
                             (node_type, chain_end_bead_type, internal_bead_type),  # canonical: A-C-D
                             (chain_end_bead_type, node_type, chain_end_bead_type),  # canonical: C-D-C
                         ])

    diamond_lattice = DiamondLattice(mpc, bond_length)
    espresso_system.box_l = [diamond_lattice.box_l] * 3
    pmb.initialize_lattice_builder(diamond_lattice)

    node_topology = [
        {"particle_name": node_type, "lattice_index": index}
        for index in diamond_lattice.indices
    ]
    chain_topology = []
    for node_connectivity in diamond_lattice.connectivity:
        node_start = str(diamond_lattice.indices[node_connectivity[0]])
        node_end = str(diamond_lattice.indices[node_connectivity[1]])
        chain_topology.append({
            "node_start": node_start,
            "node_end": node_end,
            "molecule_name": molecule_name,
        })

    pmb.define_hydrogel(name="simple_hydrogel",
                        node_map=node_topology,
                        chain_map=chain_topology)
    hydrogel_id = pmb.create_hydrogel(name="simple_hydrogel",
                                      espresso_system=espresso_system,
                                      gen_angle=True)
    return pmb, hydrogel_id


def summarize_angles(pmb):
    """Return angle templates, instances, and counts by canonical angle name."""
    angle_templates_df = pmb.get_templates_df(pmb_type="angle")
    angle_instances_df = pmb.get_instances_df(pmb_type="angle")
    angle_counts = angle_instances_df["name"].value_counts().sort_index()
    return angle_templates_df, angle_instances_df, angle_counts


def run_angle_case(label, include_crosslinker_angles, expected_angle_names, espresso_system):
    """Clear the system, build one angle-hydrogel case, and assert expected angle names."""
    espresso_system.part.clear()
    espresso_system.bonded_inter.clear()
    espresso_system.non_bonded_inter.reset()

    print(f"\n############ {label} ############")
    pmb, hydrogel_id = define_angle_hydrogel_system(espresso_system, include_crosslinker_angles)
    angle_templates_df, angle_instances_df, angle_counts = summarize_angles(pmb)
    print(f"Hydrogel assembly id: {hydrogel_id}")
    print("\nAngle templates:")
    print(angle_templates_df)
    print("\nAngle instance dataframe:")
    print(angle_instances_df)
    print("Angle instance counts:")
    print(angle_counts)

    angle_names = set(angle_counts.index)
    if angle_names != expected_angle_names:
        raise AssertionError(
            f"{label}: expected angle names {sorted(expected_angle_names)}, "
            f"got {sorted(angle_names)}"
        )
    print(f"{label}: OK")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build hydrogels: original (mpc=40) then angle-potential (mpc=4)."
    )
    parser.add_argument(
        "--include",
        action="store_true",
        help="Include crosslinker-adjacent angles (D-C-A and C-D-C) in the angle demo.",
    )
    args = parser.parse_args()

    espresso_system = espressomd.System(box_l=[1.0] * 3)

    # Part 1: original hydrogel (mpc=40, no angles)
    run_original_hydrogel(espresso_system)

    # Part 2: angle-potential hydrogel (mpc=4); clearing is done inside run_angle_case
    if args.include:
        run_angle_case(label="Chain + crosslinker angles",
                       include_crosslinker_angles=True,
                       expected_angle_names={"A-A-C", "A-C-D", "C-D-C"},
                       espresso_system=espresso_system)
    else:
        run_angle_case(label="Chain angles only",
                       include_crosslinker_angles=False,
                       expected_angle_names={"A-A-C"},
                       espresso_system=espresso_system)
