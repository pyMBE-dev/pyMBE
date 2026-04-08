#
# Example: demonstrating angle potential support in pyMBE.
#
# This script shows two complementary usages of `gen_angle=True`:
#   1. create_residue with gen_angle=True
#        - builds a single residue with auto-generated angles
#   2. create_molecule with gen_angle=True
#        - builds a chain of residues; angles are auto-generated across the
#          full molecule, including angles that span residue boundaries
#
# Note: angle auto-generation is currently supported only in `create_residue`
# and `create_molecule`. Other higher-level builders (e.g. `create_hydrogel`)
# do not yet propagate `gen_angle`.
#

import pyMBE
import numpy as np
import espressomd

# ----------------------------------------------------------------------
# Set up the ESPResSo system and pyMBE
# ----------------------------------------------------------------------
espresso_system = espressomd.System(box_l=[20] * 3)
pmb = pyMBE.pymbe_library(seed=42)

# ----------------------------------------------------------------------
# Particle definitions
# ----------------------------------------------------------------------
pmb.define_particle(name='A',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

pmb.define_particle(name='B',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

pmb.define_particle(name='C',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

# ----------------------------------------------------------------------
# Default bond used everywhere
# ----------------------------------------------------------------------
HARMONIC_bond_parameters = {
    'r_0': 0.4 * pmb.units.nm,
    'k': 40 * pmb.units('reduced_energy / reduced_length**2'),
}
pmb.define_default_bond(bond_type='harmonic',
                        bond_parameters=HARMONIC_bond_parameters)

# ----------------------------------------------------------------------
# Default angle used by both demos
# ----------------------------------------------------------------------
default_angle_params = {
    'k': 50 * pmb.units('reduced_energy'),
    'phi_0': np.pi * pmb.units(''),
}
pmb.define_default_angle(angle_type="harmonic",
                         angle_parameters=default_angle_params)


def show_database(label):
    print(f"\n############ {label} ############")
    for pmb_type in ["particle", "residue", "molecule", "bond", "angle"]:
        print(f"\n=== {pmb_type} templates ===")
        print(pmb.get_templates_df(pmb_type=pmb_type))
        print(f"\n=== {pmb_type} instances ===")
        print(pmb.get_instances_df(pmb_type=pmb_type))


# ======================================================================
# Demo 1: create_residue with gen_angle=True
# ======================================================================
# Topology: central A bonded to side particles A, B, C.
# Auto-generated angles (central A has 3 neighbors): A-A-B, A-A-C, B-A-C
pmb.define_residue(name="Res_1",
                   central_bead="A",
                   side_chains=["A", "B", "C"])

residue_id = pmb.create_residue(name="Res_1",
                                espresso_system=espresso_system,
                                use_default_bond=True,
                                gen_angle=True)

show_database("Demo 1: single residue with gen_angle=True")

# ----------------------------------------------------------------------
# Clean up: remove the residue (and its particles, bonds, angles) from
# both ESPResSo and the pyMBE database before running Demo 2.
# ----------------------------------------------------------------------
pmb.delete_instances_in_system(instance_id=residue_id,
                               pmb_type="residue",
                               espresso_system=espresso_system)

# ======================================================================
# Demo 2: create_molecule with gen_angle=True
# ======================================================================
# Topology: a chain of 4 Res_2 residues, each with central A and one
# side chain B. After auto-generation, angles include both intra-residue
# (A-A-B) and inter-residue (A-A-A) triplets.
pmb.define_residue(name="Res_2",
                   central_bead="A",
                   side_chains=["B"])

pmb.define_molecule(name="Mol_1",
                    residue_list=["Res_2"] * 4)

pmb.create_molecule(name="Mol_1",
                    number_of_molecules=1,
                    espresso_system=espresso_system,
                    backbone_vector=[1.0, 0.0, 0.0],
                    use_default_bond=True,
                    gen_angle=True)

show_database("Demo 2: linear molecule with gen_angle=True")
