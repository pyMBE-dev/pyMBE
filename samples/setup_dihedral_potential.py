#
# Copyright (C) 2026 pyMBE-dev team
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
# Example: demonstrating dihedral (torsion) potential support in pyMBE.
#
# This script shows:
#   1. create_molecule with gen_dihedral=True
#        - builds a linear chain of residues; dihedrals are auto-generated
#          across the full molecule from the bond topology
#   2. a multi-term dihedral on a single quadruplet
#        - GROMACS-style force fields stack several [dihedrals] lines (with
#          different multiplicities) on the same four atoms; pyMBE registers one
#          ESPResSo Dihedral object per term
#   3. validation of the ESPResSo <-> GROMACS proper-dihedral (type 1) mapping
#        - ESPResSo Dihedral: V(phi) = k*[1 - cos(n*phi - phi_0)]
#        - GROMACS type 1:    V(phi) = k*[1 + cos(n*phi - phi_s)]
#        - identical for phi_0 = phi_s + pi (phi_s converted from degrees)


import numpy as np
import espressomd
import pyMBE

# ----------------------------------------------------------------------
# Set up the ESPResSo system and pyMBE
# ----------------------------------------------------------------------
box_l = [20] * 3
espresso_system = espressomd.System(box_l=box_l)
pmb = pyMBE.pymbe_library(seed=42)
pmb.set_simulation_engine(espresso_system)

# ----------------------------------------------------------------------
# Particle and default bond definitions
# ----------------------------------------------------------------------
for particle_name in ("A", "B", "C", "D"):
    pmb.define_particle(name=particle_name,
                        z=0,
                        sigma=0.4 * pmb.units.nm,
                        epsilon=1 * pmb.units('reduced_energy'))

pmb.define_default_bond(bond_type='harmonic',
                        bond_parameters={'r_0': 0.4 * pmb.units.nm,
                                         'k': 40 * pmb.units('reduced_energy / reduced_length**2')})


def dihedral_terms_on_particle(particle_id):
    """Return the ESPResSo Dihedral interaction objects bonded to `particle_id`."""
    terms = []
    for bond_tuple in espresso_system.part.by_id(particle_id).bonds:
        if len(bond_tuple) == 4 and type(bond_tuple[0]).__name__ == "Dihedral":
            terms.append(bond_tuple[0])
    return terms


# ======================================================================
# Demo 1: create_molecule with gen_dihedral=True
# ======================================================================
# Topology: a linear chain A-B-C-D. The single bonded quadruplet
# (A, B, C, D) picks up the periodic dihedral defined below.
pmb.define_dihedral_potential(dihedral_type="periodic",
                              dihedral_parameters={'k': 5 * pmb.units('reduced_energy'),
                                                   'phi_0': np.pi * pmb.units(''),
                                                   'mult': 1},
                              particle_quadruplets=[('A', 'B', 'C', 'D')])

pmb.define_residue(name="Res_A", central_bead="A", side_chains=[])
pmb.define_residue(name="Res_B", central_bead="B", side_chains=[])
pmb.define_residue(name="Res_C", central_bead="C", side_chains=[])
pmb.define_residue(name="Res_D", central_bead="D", side_chains=[])
pmb.define_molecule(name="Mol_1", residue_list=["Res_A", "Res_B", "Res_C", "Res_D"])

molecule_ids = pmb.create_molecule(name="Mol_1",
                                   number_of_molecules=1,
                                   box_l=box_l,
                                   backbone_vector=[1.0, 0.0, 0.0],
                                   use_default_bond=True,
                                   gen_dihedral=True)
pmb.add_instances_to_engine()

molecule_pids = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                       attribute="molecule_id",
                                                       value=molecule_ids[0])
print("\n############ Demo 1: linear molecule with gen_dihedral=True ############")
print(pmb.get_instances_df(pmb_type="dihedral"))
assert len(pmb.get_instances_df(pmb_type="dihedral")) == 1, "Demo 1: expected one dihedral"
assert len(dihedral_terms_on_particle(min(molecule_pids))) == 1
print("Demo 1: one backbone dihedral generated OK")

pmb.delete_instances_in_system(instance_id=molecule_ids[0], pmb_type="molecule")
pmb.db.delete_templates(pmb_type="dihedral")


# ======================================================================
# Demo 2: multi-term dihedral on one quadruplet
# ======================================================================
multi_term_parameters = [
    {'k': 4 * pmb.units('reduced_energy'), 'phi_0': 0.0 * pmb.units(''), 'mult': 1},
    {'k': 2 * pmb.units('reduced_energy'), 'phi_0': np.pi * pmb.units(''), 'mult': 3},
]
pmb.define_dihedral_potential(dihedral_type="periodic",
                              dihedral_parameters=multi_term_parameters,
                              particle_quadruplets=[('A', 'B', 'C', 'D')])

pids = [pmb.create_particle(name=name, box_l=box_l, number_of_particles=1)[0]
        for name in ("A", "B", "C", "D")]
pmb.create_bond(particle_id1=pids[0], particle_id2=pids[1], use_default_bond=True)
pmb.create_bond(particle_id1=pids[1], particle_id2=pids[2], use_default_bond=True)
pmb.create_bond(particle_id1=pids[2], particle_id2=pids[3], use_default_bond=True)
pmb.create_dihedral_potential(particle_id1=pids[0],
                              particle_id2=pids[1],
                              particle_id3=pids[2],
                              particle_id4=pids[3])
pmb.add_instances_to_engine()

print("\n############ Demo 2: multi-term dihedral ############")
print(f"registered {len(dihedral_terms_on_particle(pids[0]))} ESPResSo Dihedral objects "
      f"for a {len(multi_term_parameters)}-term dihedral")
assert len(dihedral_terms_on_particle(pids[0])) == len(multi_term_parameters)
print("Demo 2: one Dihedral object per term OK")

for inst_id in pids:
    pmb.delete_instances_in_system(instance_id=inst_id, pmb_type="particle")
pmb.db.delete_templates(pmb_type="dihedral")


# ======================================================================
# Demo 3: ESPResSo <-> GROMACS proper dihedral (type 1) mapping
# ======================================================================
# A GROMACS proper dihedral line "phi_s k n" (phi_s in degrees) maps onto a
# pyMBE "periodic" dihedral with the same k, mult = n and phi_0 = radians(phi_s) + pi.
gromacs_phi_s_deg = 180.0
gromacs_k = 3.5          # reduced units here; kJ/mol in a real GROMACS topology
gromacs_n = 2
phi_0 = np.deg2rad(gromacs_phi_s_deg) + np.pi

pmb.define_dihedral_potential(dihedral_type="periodic",
                              dihedral_parameters={'k': gromacs_k * pmb.units('reduced_energy'),
                                                   'phi_0': phi_0 * pmb.units(''),
                                                   'mult': gromacs_n},
                              particle_quadruplets=[('A', 'B', 'C', 'D')])

pids = [pmb.create_particle(name=name, box_l=box_l, number_of_particles=1)[0]
        for name in ("A", "B", "C", "D")]
pmb.create_bond(particle_id1=pids[0], particle_id2=pids[1], use_default_bond=True)
pmb.create_bond(particle_id1=pids[1], particle_id2=pids[2], use_default_bond=True)
pmb.create_bond(particle_id1=pids[2], particle_id2=pids[3], use_default_bond=True)
pmb.create_dihedral_potential(particle_id1=pids[0],
                              particle_id2=pids[1],
                              particle_id3=pids[2],
                              particle_id4=pids[3])
pmb.add_instances_to_engine()

params = dihedral_terms_on_particle(pids[0])[0].get_params()
print("\n############ Demo 3: GROMACS proper dihedral (type 1) mapping ############")
print(f"  bend  = {params['bend']}   (GROMACS k = {gromacs_k})")
print(f"  mult  = {params['mult']}   (GROMACS n = {gromacs_n})")
print(f"  phase = {params['phase']:.6f}   (radians(phi_s) + pi = {phi_0:.6f})")
assert np.isclose(params['bend'], gromacs_k)
assert params['mult'] == gromacs_n
assert np.isclose(params['phase'], phi_0)
print("Demo 3: GROMACS proper dihedral (type 1) mapping verified OK")

for inst_id in pids:
    pmb.delete_instances_in_system(instance_id=inst_id, pmb_type="particle")
pmb.db.delete_templates(pmb_type="dihedral")
