# Import pyMBE and other libraries
import pyMBE
import numpy as np
import warnings
        
# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

#### DEFINE_BOND TESTS ###

print(f"*** Unit test: check that define_bond sets up a harmonic bond correctly ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

pmb.define_bond(bond_type = bond_type,
                bond_parameters = bond,
                particle_pairs = [['A', 'A']])

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    np.testing.assert_equal(actual=str(np.array(pmb.filter_df(pmb_type = 'bond')['bond_object'])[0])[:8], 
                        desired='Harmonic', 
                        verbose=True)

print("*** Unit test passed ***")

# Clean pmb.df

pmb.setup_df()

print(f"*** Unit test: check that define_bond sets up a FENE bond correctly ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'FENE'

bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
        'd_r_max': 0.8 * pmb.units.nm}

pmb.define_bond(bond_type = bond_type,
                bond_parameters = bond,
                particle_pairs = [['A', 'A']])

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    np.testing.assert_equal(actual=str(np.array(pmb.filter_df(pmb_type = 'bond')['bond_object'])[0])[:4], 
            desired='Fene', 
            verbose=True)

print("*** Unit test passed ***")

# Clean pmb.df

pmb.setup_df()

print(f"*** Unit test: check that define_bond sets up a harmonic and a FENE bonds correctly in the same script ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name='B', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type_1 = 'harmonic'
bond_1 = {'r_0'    : 0.4*pmb.units.nm,
          'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

pmb.define_bond(bond_type = bond_type_1,
                bond_parameters = bond_1,
                particle_pairs = [['A', 'A']])

bond_type_2 = 'FENE'

bond_2 = {'r_0'    : 0.4*pmb.units.nm,
          'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
          'd_r_max': 0.8 * pmb.units.nm}

pmb.define_bond(bond_type = bond_type_2,
                bond_parameters = bond_2,
                particle_pairs = [['B', 'B']])

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    np.testing.assert_equal(actual=[str(np.array(pmb.filter_df(pmb_type = 'bond')['bond_object'])[0])[:8], str(np.array(pmb.filter_df(pmb_type = 'bond')['bond_object'])[1])[:4]],
                        desired=['Harmonic','Fene'],
                        verbose=True)

print("*** Unit test passed ***")

# Clean pmb.df

pmb.setup_df()

print(f"*** Unit test: check that define_bond raises a ValueError if the provided bond_type differs from the two currently implemented (harmonic and FENE) ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'Quartic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  "particle_pairs" : [['A', 'A']]
                  }

np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_bond raises a ValueError if the provided dictionary does not have the Equilibrium length (r_0) for setting up a harmonic bond ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'harmonic'
bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  "particle_pairs" : [['A', 'A']]
                  }

np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_bond raises a ValueError if the provided dictionary does not have the Magnitud of the potential (k) for setting up a harmonic bond ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  "particle_pairs" : [['A', 'A']]
                  }

np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_bond raises a ValueError if the provided dictionary does not have the Magnitud of the potential (k) for setting up a FENE bond ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'FENE'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'd_r_max': 0.8 * pmb.units.nm}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  "particle_pairs" : [['A', 'A']]
                  }

np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_bond raises a ValueError if the provided dictionary does not have the Maximal stretching length (d_r_max) for setting up a FENE bond ***")

# Define dummy particle

pmb.define_particle(name='A', q=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

# Define dummy bond

bond_type = 'FENE'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  "particle_pairs" : [['A', 'A']]
                  }

np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

print(f"*** Unit test passed ***")

#### DEFINE_DEFAULT_BOND TESTS ###

print(f"*** Unit test: check that define_default_bond sets up a harmonic bond correctly ***")

# Define dummy bond

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

pmb.define_default_bond(bond_type = bond_type,
                        bond_parameters = bond)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    np.testing.assert_equal(actual=str(np.array(pmb.filter_df(pmb_type = 'bond')['bond_object'])[0])[:8],
                        desired='Harmonic',
                        verbose=True)

print("*** Unit test passed ***")

# Clean pmb.df

pmb.setup_df()

print(f"*** Unit test: check that define_dafult_bond sets up a FENE bond correctly ***")

# Define dummy bond

bond_type = 'FENE'

bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
        'd_r_max': 0.8 * pmb.units.nm}

pmb.define_default_bond(bond_type = bond_type,
                        bond_parameters = bond)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    np.testing.assert_equal(actual=str(np.array(pmb.filter_df(pmb_type = 'bond')['bond_object'])[0])[:4],
            desired='Fene',
            verbose=True)

print("*** Unit test passed ***")

# Clean pmb.df

pmb.setup_df()

print(f"*** Unit test: check that define_default_bond raises a ValueError if the provided bond_type differs from the two currently implemented (harmonic and FENE) ***")

# Define dummy bond

bond_type = 'Quartic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  }

np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_default_bond raises a ValueError if the provided dictionary does not have the Equilibrium length (r_0) for setting up a harmonic bond ***")

# Define dummy bond

bond_type = 'harmonic'
bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  }

np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_default_bond raises a ValueError if the provided dictionary does not have the Magnitud of the potential (k) for setting up a harmonic bond ***")

# Define dummy bond

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  }

np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_default_bond raises a ValueError if the provided dictionary does not have the Magnitud of the potential (k) for setting up a FENE bond ***")

# Define dummy bond

bond_type = 'FENE'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'd_r_max': 0.8 * pmb.units.nm}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  }

np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_default_bond raises a ValueError if the provided dictionary does not have the Maximal stretching length (d_r_max) for setting up a FENE bond ***")

# Define dummy bond

bond_type = 'FENE'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

input_parameters={"bond_type": bond_type,
                  "bond_parameters" : bond,
                  }

np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

print(f"*** Unit test passed ***")

