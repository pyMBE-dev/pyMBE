import numpy as np
import espressomd
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(SEED=42)

# Define a set of ions
pmb.define_particle(name="Na", 
                    q=1)
pmb.define_particle(name="Ca", 
                    q=2)
pmb.define_particle(name="Cl", 
                    q=-1)
pmb.define_particle(name="SO4", 
                    q=-2)

type_map=pmb.get_type_map()
# System parameters
c_salt_input = 0.01 * pmb.units.mol/ pmb.units.L
N_SALT_ION_PAIRS = 50
volume = N_SALT_ION_PAIRS/(pmb.N_A*c_salt_input)
L = volume ** (1./3.) # Side of the simulation box

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)
espresso_system.setup_type_map(type_map.values())

#### Unit tests for the added salt

def check_salt_concentration(espresso_system,cation_name,anion_name,c_salt,N_SALT_ION_PAIRS, verbose=False):
    charge_map=pmb.get_charge_map()
    type_map=pmb.get_type_map()
    espresso_system.setup_type_map(type_map.values())
    c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                            cation_name=cation_name,
                                            anion_name=anion_name,
                                            c_salt=c_salt,
                                            verbose=verbose)

    np.testing.assert_equal(espresso_system.number_of_particles(type_map[cation_name]),N_SALT_ION_PAIRS*abs(charge_map[type_map[anion_name]]))
    np.testing.assert_equal(espresso_system.number_of_particles(type_map[anion_name]),N_SALT_ION_PAIRS*abs(charge_map[type_map[cation_name]]))
    np.testing.assert_almost_equal(c_salt_calculated.m_as("mol/L"), c_salt.m_as("mol/L"))
    espresso_system.part.clear()
print("*** Unit test: test that create_added_salt works for a 1:1 salt (NaCl-like). Should print the added salt concentration and number of ions ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Na",
                        anion_name="Cl",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS,
                        verbose=True)
print("*** Unit test passed***")
print("*** Unit test: test that create_added_salt works for a 2:1 salt (CaCl_2-like) ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Ca",
                        anion_name="Cl",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: test that create_added_salt works for a 1:2 salt (Na_2SO_4-like) ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Na",
                        anion_name="SO4",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: test that create_added_salt works for a 2:2 salt (CaSO_4-like) ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Ca",
                        anion_name="SO4",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: check that create_added_salt works for an input c_salt in [particle/lenght**3]. Should print the concentration and number of ions")
c_salt_part=c_salt_input*pmb.N_A
espresso_system.setup_type_map(type_map.values())
c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                            cation_name="Na",
                                            anion_name="Cl",
                                            c_salt=c_salt_part)
np.testing.assert_equal(espresso_system.number_of_particles(type_map["Na"]),N_SALT_ION_PAIRS)
np.testing.assert_equal(espresso_system.number_of_particles(type_map["Cl"]),N_SALT_ION_PAIRS)
np.testing.assert_almost_equal(c_salt_calculated.m_as("reduced_length**-3"), c_salt_part.m_as("reduced_length**-3"))
espresso_system.part.clear()

print("*** Unit test: check that create_added_salt raises a ValueError if one provides a cation_name of an object that has been defined with a non-positive charge ***")
input_parameters={"cation_name":"Cl",
                    "anion_name":"SO4",
                    "c_salt":c_salt_input,
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_added_salt, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that create_added_salt raises a ValueError if one provides a anion_name of an object that has been defined with a non-negative charge ***")
input_parameters={"cation_name":"Na",
                    "anion_name":"Ca",
                    "c_salt":c_salt_input,
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_added_salt, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that create_added_salt raises a ValueError if one provides a c_salt with the wrong dimensionality ***")
input_parameters={"cation_name":"Na",
                    "anion_name":"Cl",
                    "c_salt":1*pmb.units.nm,
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_added_salt, **input_parameters)
print("*** Unit test passed ***")

### Unit tests for the counter ions:


def setup_molecules():
    pmb.define_particle(name='0P',
                            q=0)
    pmb.define_particle(name='+1P',
                        q=+1)
    pmb.define_particle(name='-1P',
                        q=-1)
    pmb.define_residue(
        name = 'R1',
        central_bead = '0P',
        side_chains = ['+1P']
        )

    pmb.define_residue(
        name = 'R2',
        central_bead = '0P',
        side_chains = ['-1P']
        )

    bond_type = 'harmonic'
    generic_bond_length=0.4 * pmb.units.nm
    generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

    harmonic_bond = {'r_0'    : generic_bond_length,
                    'k'      : generic_harmonic_constant,
                    }

    pmb.define_default_bond(bond_type = bond_type, 
                            bond_parameters = harmonic_bond)
    # Add all bonds to espresso system
    pmb.add_bonds_to_espresso(espresso_system=espresso_system)
    molecule_name = 'positive_polyampholyte'
    pmb.define_molecule(name=molecule_name, 
                        residue_list = ['R1']*3+['R2']*2)

    molecule_name = 'isoelectric_polyampholyte'
    pmb.define_molecule(name=molecule_name, 
                        residue_list = ['R1']*3+['R2']*3)

    molecule_name = 'negative_polyampholyte'
    pmb.define_molecule(name=molecule_name, 
                        residue_list = ['R1']*2+['R2']*3)

def test_counterions(molecule_name, cation_name, anion_name, espresso_system, expected_numbers,verbose=False):
    setup_molecules()
    pmb.create_molecule(name=molecule_name,
                        number_of_molecules= 2,
                        espresso_system=espresso_system,
                        use_default_bond=True)
    pmb.create_counterions(object_name=molecule_name,
                            cation_name=cation_name,
                            anion_name=anion_name,
                            espresso_system=espresso_system,
                            verbose=verbose)
    espresso_system.setup_type_map(type_map.values())
    np.testing.assert_equal(espresso_system.number_of_particles(type_map[cation_name]),expected_numbers[cation_name])
    np.testing.assert_equal(espresso_system.number_of_particles(type_map[anion_name]),expected_numbers[anion_name])
    pmb.destroy_pmb_object_in_system(espresso_system=espresso_system,
                                    name=molecule_name)
    espresso_system.part.clear()

print("*** Unit test: check that create_counterions creates the right number of monovalent counter ions for a polyampholyte with positive net charge. Should print the number of ions. ***")

test_counterions(molecule_name='positive_polyampholyte', 
                cation_name="Na", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Na":4,
                                  "Cl":6},
                verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions creates the right number of monovalent counter ions for a polyampholyte at its isoelectric point ***")

test_counterions(molecule_name='isoelectric_polyampholyte', 
                cation_name="Na", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Na":6,
                                  "Cl":6})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions creates the right number of monovalent counter ions for a polyampholyte with a negative net charge ***")

test_counterions(molecule_name='negative_polyampholyte', 
                cation_name="Na", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Na":6,
                                  "Cl":4})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions creates the right number of multivalent counter ions for a polyampholyte ***")

test_counterions(molecule_name='negative_polyampholyte', 
                cation_name="Ca", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Ca":3,
                                  "Cl":4})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions raises a ValueError if the charge number of the cation is not divisible by the negative charge of the polyampholyte ***")
setup_molecules()
pmb.create_molecule(name='isoelectric_polyampholyte',
                        number_of_molecules= 1,
                        espresso_system=espresso_system,
                        use_default_bond=True)
input_parameters={"cation_name":"Ca",
                    "anion_name":"Cl",
                    "object_name":'isoelectric_polyampholyte',
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_counterions, **input_parameters)
print("*** Unit test passed ***")
print("*** Unit test: check that create_counterions raises a ValueError if the charge number of the anion is not divisible by the positive charge of the polyampholyte ***")
setup_molecules()
input_parameters={"cation_name":"Na",
                    "anion_name":"SO4",
                    "object_name":'isoelectric_polyampholyte',
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_counterions, **input_parameters)
pmb.destroy_pmb_object_in_system(espresso_system=espresso_system,
                                    name='isoelectric_polyampholyte')

print("*** Unit test passed ***")
print("*** Unit test: check that no create_counterions does not create counterions for molecules with no charge")
pmb.define_particle(name='0P',
                        q=0)
pmb.define_residue(
    name = 'R0',
    central_bead = '0P',
    side_chains = []
    )

pmb.define_molecule(name='neutral_molecule', 
                    residue_list = ['R0'])

pmb.create_molecule(name='neutral_molecule',
                    number_of_molecules= 1,
                    espresso_system=espresso_system)
pmb.create_counterions(object_name='neutral_molecule',
                            cation_name="Na",
                            anion_name="Cl",
                            espresso_system=espresso_system,
                            verbose=False)
espresso_system.setup_type_map(type_map.values())
np.testing.assert_equal(espresso_system.number_of_particles(type_map["Na"]),0)
np.testing.assert_equal(espresso_system.number_of_particles(type_map["Cl"]),0)
print("*** Unit test passed ***")