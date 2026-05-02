#
# Example: demonstrating hydrogel angle generation in pyMBE.
#
# This script builds a simple hydrogel of:
#   - N: crosslinker/node particles
#   - C: residue/chain particles
#
# By default, it defines only C-C-C and generates only chain-internal angles.
# If run with `--include`, it also defines the crosslinker-adjacent
# angle templates N-C-C and C-N-C.
#
# Note:
#   - The canonical angle name for N-C-C is C-C-N, because pyMBE sorts the
#     two side particles alphabetically when building angle template names.
#

import argparse
import numpy as np
import espressomd

import pyMBE
from pyMBE.lib.lattice import DiamondLattice


def define_simple_hydrogel_system(include_crosslinker_angles):
    """
    Build a simple hydrogel with C residues and N crosslinkers.
    """
    pmb = pyMBE.pymbe_library(seed=42)

    node_type = "N"
    bead_type = "C"
    residue_name = "Res"
    molecule_name = "hydrogel_chain"

    bond_length = 0.355 * pmb.units.nm
    mpc = 4

    pmb.define_particle(name=node_type,
                        sigma=bond_length,
                        epsilon=1 * pmb.units("reduced_energy"))
    pmb.define_particle(name=bead_type,
                        sigma=bond_length,
                        epsilon=1 * pmb.units("reduced_energy"))

    pmb.define_residue(name=residue_name,
                       central_bead=bead_type,
                       side_chains=[])
    pmb.define_molecule(name=molecule_name,
                        residue_list=[residue_name] * mpc)

    harmonic_bond_parameters = {
        "r_0": bond_length,
        "k": 400 * pmb.units("reduced_energy / reduced_length**2"),
    }
    pmb.define_bond(bond_type="harmonic",
                    bond_parameters=harmonic_bond_parameters,
                    particle_pairs=[
                        [bead_type, bead_type],
                        [node_type, bead_type],
                    ])

    angle_parameters = {
        "k": 25 * pmb.units("reduced_energy"),
        "phi_0": np.pi * pmb.units(""),
    }
    pmb.define_angle(angle_type="harmonic",
                     angle_parameters=angle_parameters,
                     particle_triplets=[(bead_type, bead_type, bead_type)])

    if include_crosslinker_angles:
        pmb.define_angle(angle_type="harmonic",
                         angle_parameters=angle_parameters,
                         particle_triplets=[
                             (node_type, bead_type, bead_type),  # canonical: C-C-N
                             (bead_type, node_type, bead_type),  # canonical: C-N-C
                         ])

    diamond_lattice = DiamondLattice(mpc, bond_length)
    espresso_system = espressomd.System(box_l=[diamond_lattice.box_l] * 3)
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
    return pmb, espresso_system, hydrogel_id


def summarize_angles(pmb):
    """
    Return angle templates, angle instances, and counts by canonical angle name.
    """
    angle_templates_df = pmb.get_templates_df(pmb_type="angle")
    angle_instances_df = pmb.get_instances_df(pmb_type="angle")
    angle_counts = angle_instances_df["name"].value_counts().sort_index()
    return angle_templates_df, angle_instances_df, angle_counts


def run_case(label, include_crosslinker_angles, expected_angle_names):
    """
    Build one hydrogel case, print angle counts, and assert expectations.
    """
    print(f"\n############ {label} ############")
    pmb, _, hydrogel_id = define_simple_hydrogel_system(
        include_crosslinker_angles=include_crosslinker_angles
    )
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
        description="Build a simple hydrogel and optionally include crosslinker-adjacent angles."
    )
    parser.add_argument(
        "--include",
        action="store_true",
        help="Include crosslinker-adjacent angles (N-C-C and C-N-C).",
    )
    args = parser.parse_args()

    if args.include:
        run_case(label="Chain + crosslinker angles",
                 include_crosslinker_angles=True,
                 expected_angle_names={"C-C-C", "C-C-N", "C-N-C"})
    else:
        run_case(label="Chain angles only",
                 include_crosslinker_angles=False,
                 expected_angle_names={"C-C-C"})
