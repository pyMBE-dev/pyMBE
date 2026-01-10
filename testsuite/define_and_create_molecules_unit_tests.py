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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import pandas as pd
import espressomd
import logging
import io
# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)])

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# The unit tests for define_particle are in lj_tests.py and set_particle_acidity


particle_parameters={"S1":{"name":"S1",
                           "sigma":1*pmb.units.nm,
                           "offset":0.5*pmb.units.nm,
                           "epsilon":1.0*pmb.units.reduced_energy,
                           "z":1},
                     "S2":{"name":"S2",
                           "sigma":2*pmb.units.nm,
                           "epsilon":1.0*pmb.units.reduced_energy,
                           "offset":1.5*pmb.units.nm,
                           "z": 1},
                     "S3":{"name":"S3",
                           "sigma":3*pmb.units.nm,
                           "epsilon":1.0*pmb.units.reduced_energy,
                           "offset":2.5*pmb.units.nm,
                           "z":2}}
for particle_set in particle_parameters.values():
    pmb.define_particle(**particle_set)

print("*** Unit test: check that define_residue() stores the parameters correctly in the pyMBE database ***")

residue_parameters={"R1":{"name": "R1",
                        "central_bead": "S1",
                        "side_chains": []},
                    "R2":{"name": "R2",
                        "central_bead": "S1",
                        "side_chains": ["S2","S3"]},
                    "R3":{"name": "R3",
                        "central_bead": "S2",
                        "side_chains": ["R2"]}}

for parameter_set in residue_parameters.values():
    pmb.define_residue(**parameter_set)

for residue_name in residue_parameters.keys():
    input_parameters=residue_parameters[residue_name]
    res_tpl = pmb.db.get_template(pmb_type="residue", 
                                  name=residue_name)
    np.testing.assert_equal(actual=res_tpl.pmb_type, 
                            desired="residue", 
                            verbose=True)
    np.testing.assert_equal(actual=res_tpl.central_bead, 
                            desired=input_parameters["central_bead"], 
                            verbose=True)
    np.testing.assert_equal(actual=res_tpl.side_chains, 
                            desired=input_parameters["side_chains"], 
                            verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that define_molecule() stores the parameters correctly in the pyMBE database ***")

molecule_parameters={"M1":{"name": "M1",
                     "residue_list": []},
                    "M2":{"name": "M2",
                    "residue_list": ["R1","R2","R3"]}}


for parameter_set in molecule_parameters.values():
    pmb.define_molecule(**parameter_set)

for molecule_name in molecule_parameters.keys():
    input_parameters=molecule_parameters[molecule_name]
    mol_tpl = pmb.db.get_template(pmb_type="molecule", 
                                  name=molecule_name)
    np.testing.assert_equal(actual=mol_tpl.pmb_type, 
                            desired="molecule", 
                            verbose=True)
    np.testing.assert_equal(actual=mol_tpl.residue_list, 
                            desired=input_parameters["residue_list"], 
                            verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_particle() creates particles into the espresso_system with the properties defined in the pyMBE database ***")
# Create an instance of an espresso system
espresso_system=espressomd.System(box_l = [10]*3)
particle_positions=[[0,0,0],[1,1,1]]
retval = pmb.create_particle(name="S1",
                             espresso_system=espresso_system,
                             number_of_particles=2,
                             fix=True,
                             position=particle_positions)
np.testing.assert_array_equal(retval, [0, 1])

particle_ids=pmb.get_particle_id_map(object_name="S1")["all"]
type_map=pmb.get_type_map()
for pid in particle_ids:
    particle=espresso_system.part.by_id(pid)
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map["S1"], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters["S1"]["z"], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.fix, 
                                desired=[True]*3, 
                                verbose=True)
    np.testing.assert_equal(actual=particle.pos, 
                                desired=particle_positions[pid], 
                                verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_particle() does not create any particle for number_of_particles <= 0  ***")
starting_number_of_particles=len(espresso_system.part.all())

for number_of_particles in [0, -1]:
    retval = pmb.create_particle(name="S1",
                                 espresso_system=espresso_system,
                                 number_of_particles=number_of_particles)
    np.testing.assert_equal(len(retval), 0)
# If no particles have been created, only two particles should be in the system (from the previous test)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)
print("*** Unit test passed ***")

print("*** Unit test: check that create_particle() does not create any particle if one provides an undefined name  ***")
pmb.create_particle(name="S23",
                    espresso_system=espresso_system,
                    number_of_particles=1)
    
# If no particles have been created, only two particles should be in the system (from the previous test)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)
print("*** Unit test passed ***")

# Unit tests for delete particle
print("*** Unit test: check that delete_particle deletes the particle and cleans the pyMBE database  ***")
starting_number_of_particles=len(espresso_system.part.all())
starting_number_of_rows=len(pmb.get_instances_df(pmb_type="particle"))
# This should delete one particle instance
pmb.delete_instances_in_system(instance_id=0,
                               pmb_type="particle",
                               espresso_system=espresso_system)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles-1, 
                        verbose=True)
particle_df = pmb.get_instances_df(pmb_type="particle")
np.testing.assert_equal(actual=len(particle_df), 
                        desired=starting_number_of_rows-1, 
                        verbose=True)

print(pmb.get_instances_df(pmb_type="particle"))

# Delete the other particle instance to simplify the rest of the tests
pmb.delete_instances_in_system(instance_id=1,
                               pmb_type="particle",
                               espresso_system=espresso_system)

print("*** Unit test: check that create_residue() creates a simple residue into the espresso_system with the properties defined in pmb.df  ***")

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

pmb.define_default_bond(bond_type = bond_type,
                        bond_parameters = bond)

central_bead_position=[[0,0,0]]
backbone_vector=np.array([1.,2.,3.])

pmb.create_residue(name="R2",
                   espresso_system=espresso_system,
                   central_bead_position=central_bead_position,
                   backbone_vector=backbone_vector,
                   use_default_bond=True)

particle_ids=pmb.get_particle_id_map(object_name="R2")["all"]

# Check that the particle properties are correct
for pid in particle_ids:
    particle=espresso_system.part.by_id(pid)
    particle_tpl = pmb.db.get_instance(pmb_type="particle",
                                       instance_id=pid)
    particle_name = particle_tpl.name
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map[particle_name], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters[particle_name]["z"], 
                                verbose=True)
    # Check that the position are correct
    # Central bead
    if particle_name == "S1":
        np.testing.assert_equal(actual=particle.pos, 
                                desired=central_bead_position[0], 
                                verbose=True)
    else: # Side chains should be in positions perpendicular to the backbone vector
        np.testing.assert_almost_equal(actual=np.dot(particle.pos,backbone_vector), 
                                desired=0, 
                                verbose=True)
    # Check that particles have the correct residue id
    residue_id = particle_tpl.residue_id
    np.testing.assert_equal(actual=residue_id, 
                                desired=0, 
                                verbose=True)

# Check that particles are correctly bonded
# Central bead S1 (id 0) should be bonded to S2 (id 1) and S3 (id 2)
bonded_pairs=[]
bond_df = pmb.get_instances_df(pmb_type="bond")
for bond_index in bond_df.index:
    particle_id1= bond_df.loc[bond_index,"particle_id1"]
    particle_id2= bond_df.loc[bond_index,"particle_id2"]
    bonded_pair=frozenset([particle_id1,particle_id2])
    bonded_pairs.append(bonded_pair)
    bonded_in_espresso = False
    for pid in bonded_pair:
        for bond in espresso_system.part.by_id(pid).bonds[:]:
            bond_object = bond[0]
            partner_id  = bond[1]
            if partner_id in bonded_pair:
                bonded_in_espresso=True
    np.testing.assert_equal(actual=bonded_in_espresso, 
                            desired=True, 
                            verbose=True)

np.testing.assert_equal(actual=frozenset(bonded_pairs), 
                        desired=frozenset([frozenset([0,1]),frozenset([0,2])]), 
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_residue() creates a nested residue into the espresso_system with the properties defined in pmb.df  ***")

pmb.create_residue(name="R3",
                    espresso_system=espresso_system,
                    use_default_bond=True)

particle_ids=pmb.get_particle_id_map(object_name="R3")["all"]

# Check that the particle properties are correct
for pid in particle_ids:
    particle=espresso_system.part.by_id(pid)
    particle_tpl = pmb.db.get_instance(pmb_type="particle",
                                        instance_id=pid)
    particle_name = particle_tpl.name
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map[particle_name], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters[particle_name]["z"], 
                                verbose=True)
    # Check that particles have the correct residue id
    residue_id = particle_tpl.residue_id
    np.testing.assert_equal(actual=residue_id, 
                                desired=1, 
                                verbose=True)

# Check that particles are correctly bonded, new bonds are:
# Central bead S2 (id 3) should be bonded to R2 central bead S1 (id 4)
# Central bead S1 (id 4) should be bonded to side chains S2 (id 5) and S3 (id 6)

bonded_pairs=[]
bond_df = pmb.get_instances_df(pmb_type="bond")

for bond_index in bond_df.index:
    particle_id1= bond_df.loc[bond_index,"particle_id1"]
    particle_id2= bond_df.loc[bond_index,"particle_id2"]
    bonded_pair=frozenset([particle_id1,particle_id2])
    bonded_pairs.append(bonded_pair)
    bonded_in_espresso = False
    for pid in bonded_pair:
        for bond in espresso_system.part.by_id(pid).bonds[:]:
            bond_object = bond[0]
            partner_id  = bond[1]
            if partner_id in bonded_pair:
                bonded_in_espresso=True
                
    np.testing.assert_equal(actual=bonded_in_espresso, 
                            desired=True, 
                            verbose=True)

np.testing.assert_equal(actual=frozenset(bonded_pairs), 
                        desired=frozenset([frozenset([0,1]),
                                           frozenset([0,2]),
                                            frozenset([3,4]),
                                            frozenset([4,5]),
                                            frozenset([4,6])]), 
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_residue() does not create any residue if name is not defined in the pyMBE database ***")
starting_number_of_particles=len(espresso_system.part.all())
pmb.create_residue(name="R51",
                    espresso_system=espresso_system,
                    use_default_bond=True)
# If no particles have been created, the number of particles should be the same as before
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)

# Tests for delete_residue
print("*** Unit test: check that delete_residue deletes the particle and cleans the pyMBE database ***")

# This should delete 3 particles (residue 0 is a R2 residue)
starting_number_of_particles=len(espresso_system.part.all())
pmb.delete_instances_in_system(instance_id=0,
                               pmb_type="residue",
                               espresso_system=espresso_system)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles-3, 
                        verbose=True)
# There should be only one residue instance now in the pyMBE database
np.testing.assert_equal(actual=len(pmb.get_instances_df(pmb_type="residue")), 
                        desired=1, 
                        verbose=True)
# And there should be only 4 particles (central bead + 2 side chains + central bead of R3)
np.testing.assert_equal(actual=len(pmb.get_instances_df(pmb_type="particle")), 
                        desired=4, 
                        verbose=True)
# Delete the other residue instance to simplify the rest of the tests
pmb.delete_instances_in_system(instance_id=1,
                               pmb_type="residue",
                               espresso_system=espresso_system)
print("*** Unit test passed ***")
# Additional unit tests for define_molecule are in create_molecule_position_test
print("*** Unit test: check that create_molecule() creates a simple molecule into the espresso_system with the properties defined in pmb.df  ***")

backbone_vector = np.array([1,3,-4])
magnitude = np.linalg.norm(backbone_vector)
backbone_vector = backbone_vector/magnitude
molecule_info_M2 = pmb.create_molecule(name="M2",
                    number_of_molecules=1,
                    espresso_system=espresso_system,
                    backbone_vector = backbone_vector,
                    use_default_bond=True)

particle_ids=pmb.get_particle_id_map(object_name="M2")["all"]

# Residue and molecule IDs expected
# For the  M2 molecule created, the residue and molecule IDs should be as follows:
# R1 (residue_id=0, molecule_id=0), R2 (residue_id=1, molecule_id=0), R3 (residue_id=2, molecule_id=0)

residue_ids={0: 0, 1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 2, 7: 2}
molecule_ids={0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}

# Check that the particle properties are correct
for pid in particle_ids:
    particle=espresso_system.part.by_id(pid)
    particle_tpl = pmb.db.get_instance(pmb_type="particle",
                                        instance_id=pid)
    particle_name = particle_tpl.name
    np.testing.assert_equal(actual=particle.type, 
                            desired=type_map[particle_name], 
                            verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters[particle_name]["z"], 
                                verbose=True)
    # Check that particles have the correct residue id
    residue_id = particle_tpl.residue_id
    np.testing.assert_equal(actual=residue_id, 
                                desired=residue_ids[pid], 
                                verbose=True)
    # Check that particles have the correct molecule id
    molecule_id = particle_tpl.molecule_id
    np.testing.assert_equal(actual=molecule_id, 
                                desired=molecule_ids[pid], 
                                verbose=True)

# Check that the molecule has the right residues
residue_list=[]
residue_df = pmb.get_instances_df(pmb_type="residue")
for res_index in residue_df[residue_df['molecule_id']==0].index:
    resname = residue_df.loc[res_index,"name"]
    residue_list.append(resname)
np.testing.assert_equal(actual=frozenset(residue_list), 
                        desired=frozenset(molecule_parameters["M2"]["residue_list"]), 
                        verbose=True)

# Expected bonded pairs for the molecule
# Molecule 0:
# S1(0)-S1(1) (R1-R2)
# S1(1)-S2(2) (R2)
# S1(1)-S3(3) (R2)
# S2(1)-S2(4) (R2-R3)
# S2(4)-S1(5) (R3)
# S1(5)-S2(6) (R3)
# S1(5)-S3(7) (R3)

bonded_pairs_ref=[frozenset([0,1]),
                  frozenset([1,2]),
                  frozenset([1,3]),
                  frozenset([1,4]),
                  frozenset([4,5]),
                  frozenset([5,6]),
                  frozenset([5,7])]

# Check that particles are correctly bonded
bonded_pairs=[]
bond_df = pmb.get_instances_df(pmb_type="bond")
for bond_index in bond_df.index:
    particle_id1= bond_df.loc[bond_index,"particle_id1"]
    particle_id2= bond_df.loc[bond_index,"particle_id2"]
    bonded_pair=frozenset([particle_id1,particle_id2])
    bonded_pairs.append(bonded_pair)
    bonded_in_espresso = False
    for pid in bonded_pair:
        for bond in espresso_system.part.by_id(pid).bonds[:]:
            bond_object = bond[0]
            partner_id  = bond[1]
            if partner_id in bonded_pair:
                bonded_in_espresso=True
    np.testing.assert_equal(actual=bonded_in_espresso, 
                            desired=True, 
                            verbose=True)
np.testing.assert_equal(actual  = frozenset(bonded_pairs), 
                        desired = frozenset(bonded_pairs_ref), 
                        verbose = True)

print("*** Unit test passed ***")
print("*** Unit test: check the backbone vector of the molecule in espresso and the given input backbone vector are same ***")

central_bead_positions = []

residue_map=pmb.get_particle_id_map(object_name="M2")["residue_map"]
for res_id in residue_map.keys():
    central_bead_id = min(residue_map[res_id])
    central_bead_pos = espresso_system.part.by_id(central_bead_id).pos
    central_bead_positions.append(central_bead_pos)

# Here one expects 3 central bead positions for residues R1, R2, and R3
np.testing.assert_equal(len(central_bead_positions),len(molecule_parameters["M2"]["residue_list"]))
   
backbone_direction_1 = central_bead_positions[1] - central_bead_positions[0]
backbone_direction_2 = central_bead_positions[2] - central_bead_positions[1]
backbone_direction_1 /= np.linalg.norm(backbone_direction_1)
backbone_direction_2 /= np.linalg.norm(backbone_direction_2)
np.testing.assert_almost_equal(
    actual = backbone_direction_1,
    desired = backbone_vector,
    verbose = True)
np.testing.assert_almost_equal(
    actual = backbone_direction_2,
    desired = backbone_vector,
    verbose = True)


print("*** Unit test passed ***")

print("*** Unit test: check that create_molecule() does not create any molecule for number_of_molecules <= 0  ***")

starting_number_of_particles=len(espresso_system.part.all())
pmb.create_molecule(name="M2",
                    number_of_molecules=0,
                    espresso_system=espresso_system,
                    use_default_bond=True)
pmb.create_molecule(name="M2",
                    number_of_molecules=-1,
                    espresso_system=espresso_system,
                    use_default_bond=True)
# If no particles have been created, only two particles should be in the system (from the previous test)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)
print("*** Unit test passed ***")

print("*** Unit test: check that create_molecule() does not create any molecule if one provides an undefined name  ***")

starting_number_of_particles=len(espresso_system.part.all())
pmb.create_molecule(name="M23",
                    number_of_molecules=1,
                    espresso_system=espresso_system,
                    use_default_bond=True)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)
print("*** Unit test passed ***")

# Unit tests for get_radius_map
print("*** Unit test: check that get_radius_map() provides the right amount of radii corresponding to the number of different particles in the simulation box ***")

np.testing.assert_equal(actual=len(pmb.get_radius_map()),
                        desired=len(particle_parameters.values()),
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that get_radius_map() provides the right values of the radii of the particles, which corresponds to (sigma+offset)/2 ***")

desired_radii=[]
for particle in particle_parameters.values():
    desired_radii.append((particle['sigma'].magnitude+particle['offset'].magnitude)/2)

actual_radii=[pmb.get_radius_map()[0],
               pmb.get_radius_map()[1],
               pmb.get_radius_map()[2],]

np.testing.assert_equal(actual=actual_radii,
                        desired=desired_radii,
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that the default value for the argument 'dimensionless' in get_radius_map() is True ***")

np.testing.assert_equal(actual=isinstance(pmb.get_radius_map()[0],float),
                        desired=True,
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that if the argument 'dimensionless' is False in get_radius_map() then we obtain the corresponding units ***")


np.testing.assert_equal(actual=pmb.get_radius_map(dimensionless=False)[0].dimensionality,
                        desired=pmb.units.nm.dimensionality,
                        verbose=True)

print("*** Unit test passed ***")

# Tests for delete_residue
print("*** Unit test: check that delete_molecule deletes the particle and cleans pmb.df  ***")
# create another molecule just to have two molecules in the system
pmb.create_molecule(name="M2",
                    number_of_molecules=1,
                    espresso_system=espresso_system,
                    backbone_vector = backbone_vector,
                    use_default_bond=True)

# This should delete 8 particles (molecule 0 is a M2 molecule)

starting_number_of_particles=len(espresso_system.part.all())
pmb.delete_instances_in_system(instance_id=0,
                               pmb_type="molecule",
                               espresso_system=espresso_system)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles-8, 
                        verbose=True)

# There should only one molecule instance now in the pyMBE database
np.testing.assert_equal(actual=len(pmb.get_instances_df(pmb_type="molecule")), 
                        desired=1, 
                        verbose=True)
# There should be only 3 residues (from the remaining M2 molecule)
np.testing.assert_equal(actual=len(pmb.get_instances_df(pmb_type="residue")), 
                        desired=3, 
                        verbose=True)
# There should be only 8 particles (from the remaining M2 molecule)
np.testing.assert_equal(actual=len(pmb.get_instances_df(pmb_type="particle")), 
                        desired=8,
                        verbose=True)

print("*** Unit test passed ***")