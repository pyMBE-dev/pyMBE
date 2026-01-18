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
import numpy as np
import pandas as pd
import pyMBE
import unittest as ut

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)


class Test(ut.TestCase):
    def test_inert_particles_setup(self):
        """
        Test that an inert particle is correctly set up in the pyMBE database.
        """
        input_parameters={"name":"I", 
                  "acidity": pd.NA,
                  "pka": pd.NA,
                  "z":2,
                  "sigma": 1.0*pmb.units.reduced_length,
                  "epsilon": 1.0*pmb.units.reduced_energy}
        pmb.define_particle(**input_parameters)
        part_tpl = pmb.db.get_template(name="I", 
                                    pmb_type="particle")
        self.assertTrue(hasattr(part_tpl, "states"))
        self.assertEqual(len(part_tpl.states), 1)
        state_one = part_tpl.states["I"]
        self.assertEqual(state_one.name, "I")
        self.assertEqual(state_one.z, 2)
        pmb.db.delete_template(name="I", pmb_type="particle")   

    def test_acidic_particles_setup(self):
        """
        Test that an acidic particle is correctly set up in the pyMBE database.
        """
        input_parameters={"name":"A", 
                  "acidity": "acidic",
                  "pka":4,
                  "sigma": 1.0*pmb.units.reduced_length,
                  "epsilon": 1.0*pmb.units.reduced_energy}
        pmb.define_particle(**input_parameters)
        part_tpl = pmb.db.get_template(name="A", 
                                    pmb_type="particle")
        self.assertTrue(hasattr(part_tpl, "states"))
        self.assertEqual(len(part_tpl.states), 2)
        state_one = part_tpl.states["AH"]
        self.assertEqual(state_one.name, "AH")
        self.assertEqual(state_one.z, 0)
        state_two = part_tpl.states["A"]
        self.assertEqual(state_two.name, "A")
        self.assertEqual(state_two.z, -1)
        self.assertNotEqual(state_one.es_type, state_two.es_type)
        pmb.db.delete_template(name="A", pmb_type="particle")

    def test_basic_particles_setup(self):
        """
        Test that a basic particle is correctly set up in the pyMBE database.
        """
        input_parameters={"name":"B", 
                  "acidity": "basic",
                  "pka":9,
                  "sigma": 1.0*pmb.units.reduced_length,
                  "epsilon": 1.0*pmb.units.reduced_energy}
        pmb.define_particle(**input_parameters)
        part_tpl = pmb.db.get_template(name="B", 
                                    pmb_type="particle")
        self.assertTrue(hasattr(part_tpl, "states"))
        self.assertEqual(len(part_tpl.states), 2)
        state_one = part_tpl.states["BH"]
        self.assertEqual(state_one.name, "BH")
        self.assertEqual(state_one.z, 1)
        state_two = part_tpl.states["B"]
        self.assertEqual(state_two.name, "B")
        self.assertEqual(state_two.z, 0)
        self.assertNotEqual(state_one.es_type, state_two.es_type)
        pmb.db.delete_template(name="B", pmb_type="particle")

    def sanity_tests(self):
        """
        Unit tests to check that set_particle_acidity raises ValueErrors when expected.
        """
        # Check that set_particle_acidity raises a ValueError if pKa is not provided and pKa is acidic or basic
        input_parametersA={"name":"A", 
                        "acidity": "acidic" }

        input_parametersB= {"name": "B",
                        "acidity": "basic"}
        self.assertRaises(ValueError, pmb.set_particle_acidity,**input_parametersA)
        self.assertRaises(ValueError, pmb.set_particle_acidity, **input_parametersB)
        # Check that set_particle_acidity raises a ValueError if a non-supported acidity is provided
        input_parametersA={"name":"A", 
                        "acidity": "random" }
        self.assertRaises(ValueError, pmb.set_particle_acidity,**input_parametersA)

if __name__ == "__main__":
    ut.main()

