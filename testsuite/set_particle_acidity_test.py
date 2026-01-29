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
import numpy as np
import pandas as pd
import pyMBE
import unittest as ut
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
# Create an instance of pyMBE library

participants = [ReactionParticipant(particle_name="A",
                                    state_name="HA",
                                    coefficient=1),
                ReactionParticipant(particle_name="A",
                                    state_name="A",
                                    coefficient=1)]


class Test(ut.TestCase):
    def test_inert_particles_setup(self):
        """
        Test that an inert particle is correctly set up in the pyMBE database.
        """
        pmb = pyMBE.pymbe_library(seed=42)
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
        pmb = pyMBE.pymbe_library(seed=42)
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
        pmb = pyMBE.pymbe_library(seed=42)
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

        self.assertEqual(protonated_state.name, 
                         "BH")
        self.assertEqual(protonated_state.z, 
                         1)
        self.assertEqual(deprotonated_state.name, 
                         "B")
        self.assertEqual(deprotonated_state.z, 
                         0)
        self.assertNotEqual(protonated_state.es_type, 
                            deprotonated_state.es_type)
        pmb.db.delete_template(name="B", 
                               pmb_type="particle")
        pmb.db.delete_template(name="BH", 
                               pmb_type="particle_state")
        pmb.db.delete_template(name="B", 
                               pmb_type="particle_state")

    def test_sanity_acidity(self):
        """
        Unit tests to check that define_monoprototic_acidbase_reaction raises ValueErrors when expected.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        # Check that define_monoprototic_acidbase_reaction raises a ValueError if a non-supported acidity is provided
        input_parametersA={"particle_name":"A", 
                           "acidity": "random",
                           "pka":4,}
        self.assertRaises(ValueError, 
                          pmb.define_monoprototic_acidbase_reaction,
                          **input_parametersA)
        # Check that define_monoprototic_particle_states raises a ValueError if a non-supported acidity is provided
        input_parametersA={"particle_name":"A", 
                           "acidity": "random",}
        self.assertRaises(ValueError, 
                          pmb.define_monoprototic_particle_states,
                          **input_parametersA)

    def test_get_pka_set_empty(self):
        """
        Unit test to check that get_pka_set() returns an empty dict if no reactions have been defined
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_set = pmb.get_pka_set()
        self.assertEqual(pka_set, 
                         {})

    def test_get_pka_set_monoprotic_acid(self):
        """
        Unit test to check that get_pka_set() returns the right output for a monoprotic acid
        """
        pmb = pyMBE.pymbe_library(seed=42)
        reaction = Reaction(reaction_type="monoprotic_acid",
                            pK=4.5,
                            particle_name="A",
                            participants=participants)
        pmb.db._reactions["r1"] = reaction
        pka_set = pmb.get_pka_set()
        expected = {"A": {"pka_value": 4.5,
                           "acidity": "acidic"}}
        self.assertEqual(pka_set, expected)

    def test_get_pka_set_monoprotic_base(self):
        """
        Unit test to check that get_pka_set() returns the right output for a monoprotic base
        """
        pmb = pyMBE.pymbe_library(seed=42)
        reaction = Reaction(reaction_type="monoprotic_base",
                            pK=9.2,
                            particle_name="A",
                            participants=participants)
        pmb.db._reactions["r1"] = reaction
        pka_set = pmb.get_pka_set()
        expected = {"A": {"pka_value": 9.2,
                          "acidity": "basic"}}
        self.assertEqual(pka_set, expected)

    def test_get_pka_set_unsupported_reaction_skipped(self):
        """
        Unit test to check that get_pka_set() ignores unsupported reactions
        """
        pmb = pyMBE.pymbe_library(seed=42)
        supported = Reaction(reaction_type="monoprotic_acid",
                            pK=5.0,
                            particle_name="A",
                            participants=participants)
        unsupported = Reaction(reaction_type="redox",
                                pK=1.0,
                                particle_name="X",
                                participants=participants)

        pmb.db._reactions["r1"] = supported
        pmb.db._reactions["r2"] = unsupported

        pka_set = pmb.get_pka_set()

        self.assertEqual(len(pka_set), 1)
        self.assertIn("A", pka_set)
        self.assertNotIn("X", pka_set)

    
    def test_get_pka_set_duplicate_particle_raises(self):
        """
        Checks the sanity test for particles involved in multiple reactions
        """
        pmb = pyMBE.pymbe_library(seed=42)
        r1 = Reaction(reaction_type="monoprotic_acid",
                    pK=4.0,
                    particle_name="A",
                    participants=participants)
        r2 = Reaction(reaction_type="monoprotic_base",
                    pK=9.0,
                    particle_name="A",
                    participants=participants)

        pmb.db._reactions["r1"] = r1
        pmb.db._reactions["r2"] = r2

        with self.assertRaisesRegex(ValueError, "Multiple acid/base reactions found for particle 'A'"):
            pmb.get_pka_set()

if __name__ == "__main__":
    ut.main()

