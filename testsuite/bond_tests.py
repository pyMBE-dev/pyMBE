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
import unittest as ut

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

class Test(ut.TestCase):

    def setUp(self):
        pmb.setup_df()

    def check_bond_setup(self, bond_object, input_parameters, bond_type):
        """
        Checks that pyMBE sets up a bond object correctly.

        Args:
                bond_object(`espressomd.interactions`): instance of a espresso bond object.
                input_parameters(`dict`): dictionary with the parameters for the bond potential.
                bond_type(`str`): label identifying the bond type
        """
        # Check that pyMBE stores the correct type of bond
        self.assertIn(bond_type.lower(), str(bond_object).lower(),
                      msg="pyMBE does not store the correct type of bond")
        # Check that pyMBE defines the bond with the correct parameters
        bond_params = bond_object.get_params()
        reduced_units = {'r_0'    : 'reduced_length',
                         'k'      : 'reduced_energy / reduced_length**2',
                         'd_r_max': 'reduced_length'}
        for key in input_parameters.keys():
            np.testing.assert_equal(actual=bond_params[key],
                    desired=input_parameters[key].m_as(reduced_units[key]),
                    verbose=True)

    def test_custom_bond_harmonic(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        pmb.define_bond(bond_type = bond_type,
                        bond_parameters = bond,
                        particle_pairs = [['A', 'A']])

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[0]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

    def test_custom_bond_fene(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        bond_type = 'FENE'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}

        pmb.define_bond(bond_type = bond_type,
                        bond_parameters = bond,
                        particle_pairs = [['A', 'A']])

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[0]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

    def test_custom_bond_harmonic_and_fene(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name='B', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        bond_type_1 = 'harmonic'
        bond_1 = {'r_0'    : 0.4 * pmb.units.nm,
                  'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}
        pmb.define_bond(bond_type = bond_type_1,
                        bond_parameters = bond_1,
                        particle_pairs = [['A', 'A']])

        bond_type_2 = 'FENE'
        bond_2 = {'r_0'    : 0.4 * pmb.units.nm,
                  'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                  'd_r_max': 0.8 * pmb.units.nm}

        pmb.define_bond(bond_type = bond_type_2,
                        bond_parameters = bond_2,
                        particle_pairs = [['B', 'B']])

        bond_object_1 = pmb.filter_df(pmb_type='bond')['bond_object'][2]
        bond_object_2 = pmb.filter_df(pmb_type='bond')['bond_object'][3]

        self.check_bond_setup(bond_object=bond_object_1,
                              input_parameters=bond_1,
                              bond_type=bond_type_1)
        self.check_bond_setup(bond_object=bond_object_2,
                              input_parameters=bond_2,
                              bond_type=bond_type_2)

    def test_custom_bond_raised_exceptions(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        # check exceptions for unknown bond types
        bond_type = 'Quartic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          "particle_pairs" : [['A', 'A']]
                          }

        np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

        # check exceptions for missing bond equilibrium length
        bond_type = 'harmonic'
        bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          "particle_pairs" : [['A', 'A']]
                          }

        np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4 * pmb.units.nm}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          "particle_pairs" : [['A', 'A']]
                          }

        np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'd_r_max': 0.8 * pmb.units.nm}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          "particle_pairs" : [['A', 'A']]
                          }

        np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

        # check exceptions for missing bond maximal length
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          "particle_pairs" : [['A', 'A']]
                          }

        np.testing.assert_raises(ValueError, pmb.define_bond, **input_parameters)

    def test_default_bond_harmonic(self):
        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        pmb.define_default_bond(bond_type = bond_type,
                                bond_parameters = bond)

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[0]
        self.check_bond_setup(bond_object=bond_object,
                        input_parameters=bond,
                        bond_type=bond_type)

    def test_default_bond_fene(self):
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}

        pmb.define_default_bond(bond_type = bond_type,
                                bond_parameters = bond)

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[0]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

    def test_default_bond_raised_exceptions(self):
        # check exceptions for unknown bond types
        bond_type = 'Quartic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          }

        np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

        # check exceptions for missing bond equilibrium length
        bond_type = 'harmonic'
        bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          }

        np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4*pmb.units.nm}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          }

        np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'd_r_max': 0.8 * pmb.units.nm}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          }

        np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)

        # check exceptions for missing bond maximal length
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type,
                          "bond_parameters" : bond,
                          }

        np.testing.assert_raises(ValueError, pmb.define_default_bond, **input_parameters)


if __name__ == '__main__':
    ut.main()
