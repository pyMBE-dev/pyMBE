
# Import pyMBE and other libraries
import pyMBE
import numpy as np
import espressomd

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

# The unit tests for define_particle are in lj_tests.py and set_particle_acidity

print("*** Unit test: check that define_particles() does not setup any particle if no parameters are provided ***")
output = pmb.define_particles(parameters={})
np.testing.assert_equal(actual=output, 
                        desired=0, 
                        verbose=True)

print("*** Unit test: check that define_particles() defines a set of particles correctly ***")
particle_parameters={"S1":{"name": "S1",
                            "sigma":1*pmb.units.nm,
                            "q":0},
                    "S2":{"name": "S2",
                            "sigma":2*pmb.units.nm,
                            "q": 1},
                    "S3":{"name": "S3",
                            "sigma":3*pmb.units.nm,
                            "q":2}}
pmb.define_particles(parameters=particle_parameters)

for particle_name in particle_parameters.keys():
    input_parameters=particle_parameters[particle_name]
    for index in pmb.df[pmb.df['name']==particle_name].index:
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="particle", 
                                verbose=True)
        np.testing.assert_equal(actual=pmb.df.loc[index, "sigma"].values[0].to("nm").magnitude, 
                                desired=input_parameters["sigma"].to("nm").magnitude, 
                                verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that define_residue() stores the parameters correctly in pmb.df ***")

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
    for index in pmb.df[pmb.df['name']==residue_name].index:
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="residue", 
                                verbose=True)
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "central_bead"].values[0]), 
                                desired=input_parameters["central_bead"], 
                                verbose=True)
        np.testing.assert_equal(actual=pmb.df.loc[index, "side_chains"].values[0], 
                                desired=input_parameters["side_chains"], 
                                verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that define_residue() raises a ValueError if the user provides an already defined name  ***")
input_parameters={"name": "S3",
                 "central_bead": "S2",
                 "side_chains": ["R1"]}
np.testing.assert_raises(ValueError, pmb.define_residue, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that define_molecule() stores the parameters correctly in pmb.df ***")

molecule_parameters={"M1":{"name": "M1",
                     "residue_list": []},
                    "M2":{"name": "M2",
                    "residue_list": ["R1","R2","R3"]}}

for parameter_set in molecule_parameters.values():
    pmb.define_molecule(**parameter_set)

for molecule_name in molecule_parameters.keys():
    input_parameters=molecule_parameters[molecule_name]
    for index in pmb.df[pmb.df['name']==molecule_name].index:
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="molecule", 
                                verbose=True)
        np.testing.assert_equal(actual=pmb.df.loc[index, "residue_list"].values[0], 
                                desired=input_parameters["residue_list"], 
                                verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that define_molecule() raises a ValueError if the user provides an already defined name  ***")
input_parameters={"name": "S3",
                 "residue_list": ["R1"]}
np.testing.assert_raises(ValueError, pmb.define_molecule, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that create_particle() creates particles into the espresso_system with the properties defined in pmb.df  ***")
# Create an instance of an espresso system
espresso_system=espressomd.System(box_l = [10]*3)
particle_positions=[[0,0,0],[1,1,1]]
pmb.create_particle(name="S1",
                    espresso_system=espresso_system,
                    number_of_particles=2,
                    fix=True,
                    position=particle_positions)

particle_ids=pmb.get_particle_id_map(object_name="S1")["all"]
type_map=pmb.get_type_map()

for pid in particle_ids:
    particle=espresso_system.part.by_id(pid)
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map["S1"], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters["S1"]["q"], 
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
pmb.create_particle(name="S1",
                    espresso_system=espresso_system,
                    number_of_particles=0)
pmb.create_particle(name="S1",
                    espresso_system=espresso_system,
                    number_of_particles=-1)
# If no particles have been created, only two particles should be in the system (from the previous test)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)
print("*** Unit test passed ***")

print("*** Unit test: check that create_particle() raises a ValueError if the user provides the name of an object that is not a particle ***")
input_parameters={"name": "R2",
                 "espresso_system": espresso_system,
                 "number_of_particles": 1}
np.testing.assert_raises(ValueError, pmb.create_particle, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that create_residue() creates a simple residue into the espresso_system with the properties defined in pmb.df  ***")

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

pmb.define_default_bond(bond_type = bond_type,
                        bond_parameters = bond)

pmb.add_bonds_to_espresso(espresso_system=espresso_system)

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
    particle_name = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["name"].values[0]
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map[particle_name], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters[particle_name]["q"], 
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
    residue_id = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["residue_id"].values[0]
    np.testing.assert_equal(actual=residue_id, 
                                desired=0, 
                                verbose=True)

# Check that particles are correctly bonded
bonded_pairs=[]
for bond_index in pmb.df[pmb.df['pmb_type']=="bond"].index:
    particle_id1= pmb.df.loc[bond_index,"particle_id"].values[0]
    particle_id2= pmb.df.loc[bond_index,"particle_id2"].values[0]
    bonded_pair=frozenset([particle_id1,particle_id2])
    bonded_pairs.append(bonded_pair)
    bonded_in_espresso = False
    for pid in bonded_pair:
        for bond in espresso_system.part.by_id(pid).bonds[:]:
            bond_object = bond[0]
            partner_id  = bond[1]
            if partner_id in bonded_pair:
                bonded_in_espresso=True
                # Test that the bond object is correctly stored in pyMBE
                np.testing.assert_equal(actual=pmb.df.loc[bond_index,"bond_object"].values[0], 
                            desired=bond_object, 
                            verbose=True)
                np.testing.assert_equal(actual=pmb.df.loc[bond_index,"residue_id"].values[0], 
                            desired=0, 
                            verbose=True)
    np.testing.assert_equal(actual=bonded_in_espresso, 
                            desired=True, 
                            verbose=True)

np.testing.assert_equal(actual=frozenset(bonded_pairs), 
                        desired=frozenset([frozenset([2,3]),frozenset([2,4])]), 
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
    particle_name = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["name"].values[0]
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map[particle_name], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters[particle_name]["q"], 
                                verbose=True)
    # Check that particles have the correct residue id
    residue_id = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["residue_id"].values[0]
    np.testing.assert_equal(actual=residue_id, 
                                desired=1, 
                                verbose=True)

# Check that particles are correctly bonded

bonded_pairs=[]
for bond_index in pmb.df[(pmb.df['pmb_type']=="bond") & (pmb.df['residue_id']==1)].index:
    particle_id1= pmb.df.loc[bond_index,"particle_id"].values[0]
    particle_id2= pmb.df.loc[bond_index,"particle_id2"].values[0]
    bonded_pair=frozenset([particle_id1,particle_id2])
    bonded_pairs.append(bonded_pair)
    bonded_in_espresso = False
    for pid in bonded_pair:
        for bond in espresso_system.part.by_id(pid).bonds[:]:
            bond_object = bond[0]
            partner_id  = bond[1]
            if partner_id in bonded_pair:
                bonded_in_espresso=True
                # Test that the bond object is correctly stored in pyMBE
                np.testing.assert_equal(actual=pmb.df.loc[bond_index,"bond_object"].values[0], 
                            desired=bond_object, 
                            verbose=True)
                np.testing.assert_equal(actual=pmb.df.loc[bond_index,"residue_id"].values[0], 
                            desired=1, 
                            verbose=True)

    np.testing.assert_equal(actual=bonded_in_espresso, 
                            desired=True, 
                            verbose=True)

np.testing.assert_equal(actual=frozenset(bonded_pairs), 
                        desired=frozenset([frozenset([5,6]),
                                            frozenset([6,7]),
                                            frozenset([6,8])]), 
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_residue() raises a ValueError if the user provides the name of an object that is not a residue ***")
input_parameters={"name": "S2",
                 "espresso_system": espresso_system}
np.testing.assert_raises(ValueError, pmb.create_residue, **input_parameters)
print("*** Unit test passed ***")
print("*** Unit test: check that create_residue() raises a ValueError if the any of the names in side_chains does not correspond to a previously defined particle ***")
pmb.define_residue(name="test",
               central_bead="S1",
               side_chains=["test1"])
input_parameters={"name": "test",
                 "espresso_system": espresso_system}
np.testing.assert_raises(ValueError, pmb.create_residue, **input_parameters)
print("*** Unit test passed ***")
# Additional unit tests for define_molecule are in create_molecule_position_test
print("*** Unit test: check that create_molecule() creates a simple molecule into the espresso_system with the properties defined in pmb.df  ***")

pmb.create_molecule(name="M2",
                    number_of_molecules=2,
                    espresso_system=espresso_system,
                    use_default_bond=True)

particle_ids=pmb.get_particle_id_map(object_name="M2")["all"]

residue_ids={9: 2, 10: 3, 11: 3, 12: 3, 13: 4, 14: 4, 15: 4, 16: 4, # First molecule
            17: 5, 18: 6, 19: 6, 20: 6, 21: 7, 22: 7, 23: 7, 24: 7} # Second molecule

molecule_ids={9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, # First molecule
              17: 1, 18: 1, 19: 1, 20: 1, 21: 1, 22: 1, 23: 1, 24: 1} # Second molecule

# Check that the particle properties are correct
for pid in particle_ids:
    particle=espresso_system.part.by_id(pid)
    particle_name = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["name"].values[0]
    np.testing.assert_equal(actual=particle.type, 
                                desired=type_map[particle_name], 
                                verbose=True)
    np.testing.assert_equal(actual=particle.q, 
                                desired=particle_parameters[particle_name]["q"], 
                                verbose=True)
    # Check that particles have the correct residue id
    residue_id = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["residue_id"].values[0]
    np.testing.assert_equal(actual=residue_id, 
                                desired=residue_ids[pid], 
                                verbose=True)
    # Check that particles have the correct molecule id
    molecule_id = pmb.df[(pmb.df['particle_id']==pid) & (pmb.df['pmb_type']=="particle")]["molecule_id"].values[0]
    np.testing.assert_equal(actual=molecule_id, 
                                desired=molecule_ids[pid], 
                                verbose=True)

# Check that the molecules have the right residues
for mol_id in [0,1]:
    residue_list=[]
    for res_index in pmb.df[(pmb.df['pmb_type']=="residue") & (pmb.df['molecule_id']==mol_id)].index:
        resname=pmb.df.loc[res_index,"name"].values[0]
        residue_list.append(resname)
    np.testing.assert_equal(actual=frozenset(residue_list), 
                            desired=frozenset(molecule_parameters["M2"]["residue_list"]), 
                            verbose=True)

bonded_pairs_ref={0: [frozenset([9,10]),
                      frozenset([10,11]),
                      frozenset([10,12]),
                      frozenset([10,13]),
                      frozenset([13,14]),
                      frozenset([14,15]),
                      frozenset([14,16])],
                  1: [frozenset([17,18]),
                      frozenset([18,19]),
                      frozenset([18,20]),
                      frozenset([18,21]),
                      frozenset([21,22]),
                      frozenset([22,23]),
                      frozenset([22,24])]}

# Check that particles are correctly bonded
bonded_pairs={}
for mol_id in [0,1]:
    bonded_pairs[mol_id]=[]
    for bond_index in pmb.df[(pmb.df['pmb_type']=="bond") & (pmb.df['molecule_id']==mol_id)].index:
        particle_id1= pmb.df.loc[bond_index,"particle_id"].values[0]
        particle_id2= pmb.df.loc[bond_index,"particle_id2"].values[0]
        bonded_pair=frozenset([particle_id1,particle_id2])
        bonded_pairs[mol_id].append(bonded_pair)
        bonded_in_espresso = False
        for pid in bonded_pair:
            for bond in espresso_system.part.by_id(pid).bonds[:]:
                bond_object = bond[0]
                partner_id  = bond[1]
                if partner_id in bonded_pair:
                    bonded_in_espresso=True
                    # Test that the bond object is correctly stored in pyMBE
                    np.testing.assert_equal(actual=pmb.df.loc[bond_index,"bond_object"].values[0], 
                                desired=bond_object, 
                                verbose=True)
        np.testing.assert_equal(actual=bonded_in_espresso, 
                                desired=True, 
                                verbose=True)
    np.testing.assert_equal(actual=frozenset(bonded_pairs[mol_id]), 
                            desired=frozenset(bonded_pairs_ref[mol_id]), 
                            verbose=True)

print("*** Unit test passed ***")
print("*** Unit test: check that create_molecule() does not create any molcule for number_of_molecules <= 0  ***")

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