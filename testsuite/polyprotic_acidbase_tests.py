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

import json
import os
import tempfile
import pyMBE
import unittest as ut


class TestPolyprotic(ut.TestCase):

    @staticmethod
    def _expected_polyprotic_charge(pka_list, pH, acidity):
        cumulative_pka = 0.0
        terms = []
        for index, pka in enumerate(pka_list, start=1):
            cumulative_pka += pka
            terms.append(10.0 ** (-cumulative_pka + index * pH))
        numerator = sum(index * term for index, term in enumerate(terms, start=1))
        denominator = 1.0 + sum(terms)
        if acidity == "acidic":
            return -numerator / denominator
        return len(pka_list) - numerator / denominator

    def _assert_charge_profile_matches_reference(self, pmb, molecule_name, pka_list, acidity, pH_values):
        expected = [
            self._expected_polyprotic_charge(pka_list=pka_list, pH=pH, acidity=acidity)
            for pH in pH_values
        ]
        observed = pmb.calculate_HH(template_name=molecule_name, pH_list=pH_values)
        for pH, observed_charge, expected_charge in zip(pH_values, observed, expected):
            self.assertAlmostEqual(
                observed_charge,
                expected_charge,
                places=10,
                msg=f"Charge mismatch for {acidity} particle at pH={pH}",
            )

    def test_triprotic_acid_states(self):
        """
        Test that a triprotic acid creates 4 states with correct names and charges.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.15, 7.20, 12.35])
        states = pmb.db.get_particle_states_templates("PO4")
        self.assertEqual(len(states), 4)
        expected = {"H3PO4": 0, "H2PO4": -1, "HPO4": -2, "PO4": -3}
        for name, z in expected.items():
            self.assertIn(name, states)
            self.assertEqual(states[name].z, z)
        # All es_types must be unique
        es_types = [s.es_type for s in states.values()]
        self.assertEqual(len(es_types), len(set(es_types)))

    def test_triprotic_acid_initial_state(self):
        """
        Test that the initial state is the most protonated one.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.15, 7.20, 12.35])
        tpl = pmb.db.get_template("particle", "PO4")
        self.assertEqual(tpl.initial_state, "H3PO4")

    def test_triprotic_acid_reactions(self):
        """
        Test that a triprotic acid creates 3 stepwise reactions with correct pKa values.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.15, 7.20, 12.35])
        reactions = list(pmb.db.get_reactions())
        self.assertEqual(len(reactions), 3)
        for r in reactions:
            self.assertEqual(r.reaction_type, "polyprotic_acid")
        pks = sorted([r.pK for r in reactions])
        self.assertAlmostEqual(pks[0], 2.15)
        self.assertAlmostEqual(pks[1], 7.20)
        self.assertAlmostEqual(pks[2], 12.35)

    def test_triprotic_acid_reaction_participants(self):
        """
        Test that each reaction has the correct reactant/product pair.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="A",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.0, 7.0, 12.0])
        reactions = sorted(pmb.db.get_reactions(), key=lambda r: r.pK)
        expected_pairs = [("H3A", "H2A"), ("H2A", "HA"), ("HA", "A")]
        for r, (reactant, product) in zip(reactions, expected_pairs):
            reactants = [p for p in r.participants if p.coefficient < 0]
            products = [p for p in r.participants if p.coefficient > 0]
            self.assertEqual(len(reactants), 1)
            self.assertEqual(len(products), 1)
            self.assertEqual(reactants[0].state_name, reactant)
            self.assertEqual(products[0].state_name, product)

    def test_triprotic_acid_reaction_names(self):
        """
        Test that reaction names follow the expected adjacent-state convention.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="A",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.0, 7.0, 12.0])
        reactions = sorted(pmb.db.get_reactions(), key=lambda r: r.pK)
        expected_names = ["H3A <-> H2A", "H2A <-> HA", "HA <-> A"]
        self.assertEqual([reaction.name for reaction in reactions], expected_names)

    def test_diprotic_base_states(self):
        """
        Test that a diprotic base creates 3 states with correct charges.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="B",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=2,
                                       acidity="basic",
                                       pka_list=[6.0, 10.0])
        states = pmb.db.get_particle_states_templates("B")
        self.assertEqual(len(states), 3)
        expected = {"H2B": 2, "HB": 1, "B": 0}
        for name, z in expected.items():
            self.assertIn(name, states)
            self.assertEqual(states[name].z, z)

    def test_diprotic_base_reactions(self):
        """
        Test that a diprotic base creates 2 reactions with type polyprotic_base.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="B",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=2,
                                       acidity="basic",
                                       pka_list=[6.0, 10.0])
        reactions = list(pmb.db.get_reactions())
        self.assertEqual(len(reactions), 2)
        for r in reactions:
            self.assertEqual(r.reaction_type, "polyprotic_base")

    def test_monoprotic_backward_compat(self):
        """
        Test that a single pka float still works (monoprotic path).
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_particle(name="A",
                            sigma=0.35*pmb.units.nm,
                            epsilon=1*pmb.units("reduced_energy"),
                            acidity="acidic",
                            pka=4.0)
        states = pmb.db.get_particle_states_templates("A")
        self.assertEqual(len(states), 2)
        self.assertIn("AH", states)
        self.assertIn("A", states)
        reactions = list(pmb.db.get_reactions())
        self.assertEqual(len(reactions), 1)
        self.assertEqual(reactions[0].reaction_type, "monoprotic_acid")

    def test_get_pka_set_polyprotic(self):
        """
        Test that get_pka_set returns pka_values list for polyprotic particles.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.15, 7.20, 12.35])
        pka_set = pmb.get_pka_set()
        self.assertIn("PO4", pka_set)
        self.assertIn("pka_values", pka_set["PO4"])
        self.assertEqual(pka_set["PO4"]["acidity"], "acidic")
        self.assertEqual(pka_set["PO4"]["pka_values"], [2.15, 7.20, 12.35])

    def test_get_pka_set_polyprotic_preserves_pka_order(self):
        """
        Test that get_pka_set preserves the user-provided pKa ordering.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_values = [7.20, 2.15, 12.35]
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=pka_values)
        pka_set = pmb.get_pka_set()
        self.assertEqual(pka_set["PO4"]["pka_values"], pka_values)

    def test_get_pka_set_mixed(self):
        """
        Test get_pka_set with both monoprotic and polyprotic particles.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_particle(name="A",
                            sigma=0.35*pmb.units.nm,
                            epsilon=1*pmb.units("reduced_energy"),
                            acidity="acidic",
                            pka=4.0)
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=[2.15, 7.20, 12.35])
        pka_set = pmb.get_pka_set()
        self.assertIn("pka_value", pka_set["A"])
        self.assertIn("pka_values", pka_set["PO4"])

    def test_load_pka_set_polyprotic(self):
        """
        Test that load_pka_set creates reactions only (not states), matching monoprotic behavior.
        User creates states separately via define_polyprotic_particle_states.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_data = {
            "metadata": {"source": "test"},
            "data": {
                "PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]}
            }
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(pka_data, f)
            tmpfile = f.name
        try:
            metadata = pmb.load_pka_set(tmpfile)
            self.assertEqual(metadata, {"source": "test"})
            # Only reactions created, no states yet
            reactions = list(pmb.db.get_reactions())
            self.assertEqual(len(reactions), 3)
            for r in reactions:
                self.assertEqual(r.reaction_type, "polyprotic_acid")
            # Now create states separately (as user would do after load_pka_set)
            pmb.define_polyprotic_particle_states("PO4", n=3, acidity="acidic")
            states = pmb.db.get_particle_states_templates("PO4")
            self.assertEqual(len(states), 4)
        finally:
            os.unlink(tmpfile)

    def test_load_pka_set_mixed(self):
        """
        Test that load_pka_set handles a mix of monoprotic and polyprotic in one file.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_data = {
            "metadata": {},
            "data": {
                "A": {"acidity": "acidic", "pka_value": 4.5},
                "PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]}
            }
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(pka_data, f)
            tmpfile = f.name
        try:
            pmb.load_pka_set(tmpfile)
            reactions = list(pmb.db.get_reactions())
            mono_rxns = [r for r in reactions if r.reaction_type == "monoprotic_acid"]
            poly_rxns = [r for r in reactions if r.reaction_type == "polyprotic_acid"]
            self.assertEqual(len(mono_rxns), 1)
            self.assertEqual(len(poly_rxns), 3)
        finally:
            os.unlink(tmpfile)

    def test_check_pka_set_validation(self):
        """
        Test that _check_pka_set rejects invalid formats.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        # Missing both pka_value and pka_values
        with self.assertRaises(ValueError):
            pmb._check_pka_set({"X": {"acidity": "acidic"}})
        # Missing acidity
        with self.assertRaises(ValueError):
            pmb._check_pka_set({"X": {"pka_value": 4.0}})
        # Both pka_value and pka_values
        with self.assertRaises(ValueError):
            pmb._check_pka_set({"X": {"acidity": "acidic", "pka_value": 4.0, "pka_values": [1, 2]}})

    def test_define_polyprotic_particle_states_errors(self):
        """
        Test error handling in define_polyprotic_particle_states.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        # n < 2
        with self.assertRaises(ValueError):
            pmb.define_polyprotic_particle_states("X", n=1, acidity="acidic")
        # Invalid acidity
        with self.assertRaises(ValueError):
            pmb.define_polyprotic_particle_states("X", n=2, acidity="neutral")

    def test_define_polyprotic_reactions_pka_count_mismatch(self):
        """
        Test that mismatched pka_list length raises ValueError.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        state_names = ["H3A", "H2A", "HA", "A"]
        # 4 states need 3 pKa values, not 2
        with self.assertRaises(ValueError):
            pmb.define_polyprotic_acidbase_reactions(
                particle_name="A",
                state_names=state_names,
                pka_list=[2.0, 7.0],
                acidity="acidic")

    def test_define_polyprotic_particle_pka_count_mismatch(self):
        """
        Test that define_polyprotic_particle raises ValueError when len(pka_list) != n.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        with self.assertRaises(ValueError):
            pmb.define_polyprotic_particle(name="Bad",
                                           sigma=0.35*pmb.units.nm,
                                           epsilon=1*pmb.units("reduced_energy"),
                                           n=3,
                                           acidity="acidic",
                                           pka_list=[2.0, 7.0])

    def test_define_polyprotic_acidbase_reactions_invalid_acidity(self):
        """
        Test that define_polyprotic_acidbase_reactions raises ValueError for invalid acidity.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        with self.assertRaises(ValueError):
            pmb.define_polyprotic_acidbase_reactions(
                particle_name="A",
                state_names=["H2A", "HA", "A"],
                pka_list=[2.0, 7.0],
                acidity="neutral")

    def test_auto_naming_convention(self):
        """
        Test the H{k}{name} naming convention for various n values.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        # n=2 (diprotic): H2X, HX, X
        names = pmb.define_polyprotic_particle_states("X", n=2, acidity="acidic")
        self.assertEqual(names, ["H2X", "HX", "X"])
        # n=4 (tetraprotic): H4Y, H3Y, H2Y, HY, Y
        names = pmb.define_polyprotic_particle_states("Y", n=4, acidity="acidic")
        self.assertEqual(names, ["H4Y", "H3Y", "H2Y", "HY", "Y"])


    def test_calculate_HH_triprotic_acid(self):
        """
        Test triprotic-acid Henderson-Hasselbalch charge across multiple pH values.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_values = [2.15, 7.20, 12.35]
        pmb.define_polyprotic_particle(name="PO4",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=3,
                                       acidity="acidic",
                                       pka_list=pka_values)
        pmb.define_residue(name="res", central_bead="PO4", side_chains=[])
        pmb.define_molecule(name="test_mol", residue_list=["res"])
        pH_values = [0.0, 2.15, 5.0, 7.20, 10.0, 12.35, 14.0]
        self._assert_charge_profile_matches_reference(
            pmb=pmb,
            molecule_name="test_mol",
            pka_list=pka_values,
            acidity="acidic",
            pH_values=pH_values,
        )
        Z = pmb.calculate_HH(template_name="test_mol", pH_list=[0, 14])
        # At pH 0: fully protonated, charge ≈ 0
        self.assertAlmostEqual(Z[0], 0.0, places=1)
        # At pH 14: fully deprotonated, charge ≈ -3
        self.assertAlmostEqual(Z[1], -3.0, places=1)

    def test_calculate_HH_diprotic_base(self):
        """
        Test diprotic-base Henderson-Hasselbalch charge across multiple pH values.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_values = [6.0, 10.0]
        pmb.define_polyprotic_particle(name="B",
                                       sigma=0.35*pmb.units.nm,
                                       epsilon=1*pmb.units("reduced_energy"),
                                       n=2,
                                       acidity="basic",
                                       pka_list=pka_values)
        pmb.define_residue(name="res", central_bead="B", side_chains=[])
        pmb.define_molecule(name="test_mol", residue_list=["res"])
        pH_values = [0.0, 3.0, 6.0, 8.0, 10.0, 12.0, 14.0]
        self._assert_charge_profile_matches_reference(
            pmb=pmb,
            molecule_name="test_mol",
            pka_list=pka_values,
            acidity="basic",
            pH_values=pH_values,
        )
        Z = pmb.calculate_HH(template_name="test_mol", pH_list=[0, 14])
        # At pH 0: fully protonated, charge ≈ +2
        self.assertAlmostEqual(Z[0], 2.0, places=2)
        # At pH 14: fully deprotonated, charge ≈ 0
        self.assertAlmostEqual(Z[1], 0.0, places=2)


if __name__ == "__main__":
    ut.main()
