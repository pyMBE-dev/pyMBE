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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import unittest as ut
import io
import logging
import espressomd

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)] )

# Create an instance of pyMBE library
espresso_system=espressomd.System (box_l = [10]*3)

pmb = pyMBE.pymbe_library(seed=42)
pmb.define_particle(name='A', 
            z=0, 
            sigma=0.4*pmb.units.nm, 
            epsilon=1*pmb.units('reduced_energy'))

pmb.define_particle(name='B', 
                    z=0, 
                    sigma=0.4*pmb.units.nm, 
                    epsilon=1*pmb.units('reduced_energy'))

harmonic_params = {'r_0'    : 0.4 * pmb.units.nm,
                    'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

FENE_params = {'r_0'    : 0.4 * pmb.units.nm,
             'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
             'd_r_max': 0.8 * pmb.units.nm}

class Test(ut.TestCase):

    def get_bond_object(self, particle_id_pair):
        """
        Returns the bond object stored in espresso betwen a given pair of bonded particle ids.
        """
        for pid in particle_id_pair:
            if espresso_system.part.by_id(pid).bonds:
                return espresso_system.part.by_id(pid).bonds[0][0]

    def check_bond_setup(self, bond_object, input_parameters, bond_type):
        """
        Checks that pyMBE sets up a harmonic bond object correctly.

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

    def test_bond_setup(self):
        """
        Unit test to check the setup of bonds in pyMBE
        """
        #Define bond
        # check particle bond
        pmb.define_bond(bond_type = "harmonic",
                        bond_parameters = harmonic_params,
                        particle_pairs = [['A', 'A']])
        # Create two particles
        pids = pmb.create_particle(name="A",
                                   espresso_system=espresso_system,
                                   number_of_particles=2)

        pmb.create_bond(particle_id1=pids[0],
                        particle_id2=pids[1],
                        espresso_system=espresso_system,
                        use_default_bond=False)
        
        bond_object = self.get_bond_object(particle_id_pair=pids)
        
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=harmonic_params,
                              bond_type="harmonic")
        # Clean-up database
        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle",
                                           espresso_system=espresso_system)
        pid_A = pmb.create_particle(name="A",
                                   espresso_system=espresso_system,
                                   number_of_particles=1)

        harmonic_params_test = {'r_0'    : 0.5 * pmb.units.nm,
                                'k'      : 500 * pmb.units('reduced_energy / reduced_length**2')}
        pmb.define_bond(bond_type = "harmonic",
                bond_parameters = harmonic_params_test,
                particle_pairs = [['A', 'B']])
        
        pid_B = pmb.create_particle(name="B",
                                    espresso_system=espresso_system,
                                    number_of_particles=1)
        
        # Test that the bond is properly setup when there is a default bond
        pmb.define_default_bond(bond_type = "harmonic",
                                bond_parameters = harmonic_params)

        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_A[0],
                        espresso_system=espresso_system,
                        use_default_bond=True)
        
        bond_object = self.get_bond_object(particle_id_pair=[pid_B[0],pid_A[0]])
        
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=harmonic_params_test,
                              bond_type="harmonic")
        # Clean-up database
        for inst_id in pid_B+pid_A:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle",
                                           espresso_system=espresso_system)
        pmb.db.delete_templates(pmb_type="bond")
        
        # Test setup of FENE bonds
        pmb.define_bond(bond_type = "FENE",
                        bond_parameters = FENE_params,
                        particle_pairs = [['A', 'A']])
        # Create two particles
        pids = pmb.create_particle(name="A",
                                   espresso_system=espresso_system,
                                   number_of_particles=2)

        pmb.create_bond(particle_id1=pids[0],
                        particle_id2=pids[1],
                        espresso_system=espresso_system,
                        use_default_bond=False)
        
        bond_object = self.get_bond_object(particle_id_pair=pids)
        
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=FENE_params,
                              bond_type="FENE")
        # Clean-up database
        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle",
                                           espresso_system=espresso_system)
        pid_A = pmb.create_particle(name="A",
                                   espresso_system=espresso_system,
                                   number_of_particles=1)

        FENE_params_test =  {'r_0'    : 0.5 * pmb.units.nm,
                             'k'      : 500 * pmb.units('reduced_energy / reduced_length**2'),
                             'd_r_max': 0.5 * pmb.units.nm}
        pmb.define_bond(bond_type = "FENE",
                bond_parameters = FENE_params_test,
                particle_pairs = [['A', 'B']])
        
        pid_B = pmb.create_particle(name="B",
                                    espresso_system=espresso_system,
                                    number_of_particles=1)
        
        # Test that the FENE bond is properly setup when there is a default bond
        pmb.define_default_bond(bond_type = "harmonic",
                                bond_parameters = harmonic_params)

        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_A[0],
                        espresso_system=espresso_system,
                        use_default_bond=True)
        
        bond_object = self.get_bond_object(particle_id_pair=[pid_B[0],pid_A[0]])
        
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=FENE_params_test,
                              bond_type="FENE")
        # Clean-up database
        for inst_id in pid_B+pid_A:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle",
                                           espresso_system=espresso_system)
        
        pmb.db.delete_templates(pmb_type="bond")
        
        # Test setup of the default bond
        pmb.define_default_bond(bond_type = "harmonic",
                                bond_parameters = harmonic_params)

        pids = pmb.create_particle(name="A",
                                   espresso_system=espresso_system,
                                   number_of_particles=2)

        pmb.create_bond(particle_id1=pids[0],
                        particle_id2=pids[1],
                        espresso_system=espresso_system,
                        use_default_bond=True)
        
        bond_object = self.get_bond_object(particle_id_pair=pids)
        
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=harmonic_params,
                              bond_type="harmonic")
        # Clean-up database
        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle",
                                           espresso_system=espresso_system)
        pmb.db.delete_templates(pmb_type="bond")
        
        # Test setup of default bond when there are other bonds defined
        pid_A = pmb.create_particle(name="A",
                                   espresso_system=espresso_system,
                                   number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    espresso_system=espresso_system,
                                    number_of_particles=1)

        pmb.define_default_bond(bond_type = "FENE",
                                bond_parameters = FENE_params_test)
        pmb.define_bond(bond_type = "harmonic",
                        bond_parameters = harmonic_params_test,
                        particle_pairs = [['A', 'A'], ['B','B']])

        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_A[0],
                        espresso_system=espresso_system,
                        use_default_bond=True)
        
        bond_object = self.get_bond_object(particle_id_pair=[pid_B[0],pid_A[0]])
        
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=FENE_params_test,
                              bond_type="FENE")
        # Clean-up database
        for inst_id in pid_B+pid_A:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle",
                                           espresso_system=espresso_system)
        pmb.db.delete_templates(pmb_type="bond")

    def test_bond_raised_exceptions(self):
        for callback in [pmb.define_bond, pmb.define_default_bond]:
            with self.subTest(msg=f'using method {callback.__qualname__}()'):
                self.check_bond_exceptions(callback)

    def check_bond_exceptions(self, callback):
        # check exceptions for unknown bond types
        bond_type = 'Quartic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(NotImplementedError, callback, **input_parameters)

        # check exceptions for missing bond equilibrium length
        bond_type = 'harmonic'
        bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4 * pmb.units.nm}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'd_r_max': 0.8 * pmb.units.nm}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # check exceptions for missing bond maximal length
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)


if __name__ == '__main__':
    ut.main()
