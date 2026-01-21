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
        state_tpl = pmb.db.get_template(name="I", 
                                    pmb_type="particle_state")

        self.assertEqual(state_tpl.name, "I")
        self.assertEqual(state_tpl.z, 2)
        pmb.db.delete_template(name="I", pmb_type="particle")   
        pmb.db.delete_template(name="I", pmb_type="particle_state")   

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
        protonated_state = pmb.db.get_template(name="AH", 
                                    pmb_type="particle_state")
        deprotonated_state = pmb.db.get_template(name="A", 
                                    pmb_type="particle_state")
        self.assertEqual(protonated_state.name, "AH")
        self.assertEqual(protonated_state.z, 0)
        self.assertEqual(deprotonated_state.name, "A")
        self.assertEqual(deprotonated_state.z, -1)
        self.assertNotEqual(protonated_state.es_type, deprotonated_state.es_type)
        pmb.db.delete_template(name="A", pmb_type="particle")
        pmb.db.delete_template(name="AH", pmb_type="particle_state")
        pmb.db.delete_template(name="A", pmb_type="particle_state")

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

        protonated_state = pmb.db.get_template(name="BH", 
                                    pmb_type="particle_state")
        deprotonated_state = pmb.db.get_template(name="B", 
                                    pmb_type="particle_state")

        self.assertEqual(protonated_state.name, "BH")
        self.assertEqual(protonated_state.z, 1)
        self.assertEqual(deprotonated_state.name, "B")
        self.assertEqual(deprotonated_state.z, 0)
        self.assertNotEqual(protonated_state.es_type, deprotonated_state.es_type)
        pmb.db.delete_template(name="B", pmb_type="particle")
        pmb.db.delete_template(name="BH", pmb_type="particle_state")
        pmb.db.delete_template(name="B", pmb_type="particle_state")

    def test_sanity_acidity(self):
        """
        Unit tests to check that define_monoprototic_acidbase_reaction raises ValueErrors when expected.
        """
        # Check that define_monoprototic_acidbase_reaction raises a ValueError if a non-supported acidity is provided
        input_parametersA={"particle_name":"A", 
                           "acidity": "random",
                           "pka":4,}
        self.assertRaises(ValueError, pmb.define_monoprototic_acidbase_reaction,**input_parametersA)

if __name__ == "__main__":
    ut.main()

