#
# Copyright (C) 2026 pyMBE-dev team
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

import unittest as ut
import pyMBE
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.instances.hydrogel import HydrogelInstance
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.pint_quantity import PintQuantity
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
import pint

class Test(ut.TestCase):
    def test_sanity_db(self):
        """
        Sanity tests for the pyMBE database
        """
        pmb = pyMBE.pymbe_library(23)
        pmb.define_molecule(name ="test",
                            residue_list=[])
        pmb.define_peptide(name="test",
                           sequence="",
                           model="1beadAA")
        inputs = {"name": "test",
                  "allowed_types": {"molecule", "peptide"}}
        self.assertRaises(ValueError,
                          pmb._get_template_type,
                          **inputs)
        
    def test_instance_id_validators(self):
        """
        Tests that negative values of instances raise a ValueError in the pyMBE database
        """
        inputs = {"name":"A",
                  "particle_id":-1,
                   "initial_state":"A"}
        self.assertRaises(ValueError,
                          ParticleInstance,
                          **inputs)
        inputs = {"name":"A",
                  "residue_id":-1}
        self.assertRaises(ValueError,
                          ResidueInstance,
                          **inputs)
        inputs = {"name":"A",
                  "molecule_id":-1}
        self.assertRaises(ValueError,
                          MoleculeInstance,
                          **inputs)
        inputs = {"name":"A",
                  "molecule_id":-1}
        self.assertRaises(ValueError,
                          PeptideInstance,
                          **inputs)
        inputs = {"name":"A",
                  "molecule_id":-1}
        self.assertRaises(ValueError,
                          ProteinInstance,
                          **inputs)
        inputs = {"name":"A",
                  "assembly_id":-1}
        self.assertRaises(ValueError,
                          HydrogelInstance,
                          **inputs)
        inputs = {"name":"A",
                  "bond_id":-1,
                  "particle_id1":1,
                  "particle_id2":2}
        self.assertRaises(ValueError,
                          BondInstance,
                          **inputs)
        inputs = {"name":"A",
                  "bond_id":1,
                  "particle_id1":-1,
                  "particle_id2":2}
        self.assertRaises(ValueError,
                          BondInstance,
                          **inputs)
        inputs = {"name":"A",
                  "bond_id":1,
                  "particle_id1":1,
                  "particle_id2":-2}
        self.assertRaises(ValueError,
                          BondInstance,
                          **inputs)
        
    def test_make_name_bond_template(self):
        inputs = {"bond_type": "harmonic",
                  "parameters": {"r": PintQuantity(magnitude=1,
                                                   units="nm",
                                                   dimension="length")}}
        bond_tpl = BondTemplate(**inputs)
        self.assertRaises(RuntimeError,
                          bond_tpl._make_name)
        
    def test_exceptions_pint_quantity(self):
        units = pint.UnitRegistry()
        inputs = {"q":1,
                  "expected_dimension": "length",
                  "ureg": units}
        self.assertRaises(TypeError,
                          PintQuantity.from_quantity,
                          **inputs)
        inputs = {"q":1*units.nm,
                  "expected_dimension": "unknown",
                  "ureg": units}
        self.assertRaises(ValueError,
                          PintQuantity.from_quantity,
                          **inputs)
        inputs = {"q":1*units.nm**2,
                  "expected_dimension": "length",
                  "ureg": units}
        self.assertRaises(ValueError,
                          PintQuantity.from_quantity,
                          **inputs)
        
    def test_exceptions_reaction_template(self):
        """
        Tests sanity of the Reaction template
        """
        # Reactions with less than 2 participants trigger a value error
        inputs = {"participants":[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=1)],
                 "pK":1,
                 "reaction_type":"test"}
        
        self.assertRaises(ValueError,
                          Reaction,
                          **inputs)
        # Reactions with a participant with a 0 stechiometric coeff. trigger a value error
        inputs = {"participants":[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=0),
                                  ReactionParticipant(particle_name="B",
                                                   state_name="B",
                                                   coefficient=1)],
                 "pK":1,
                 "reaction_type":"test"}
        
        self.assertRaises(ValueError,
                          Reaction,
                          **inputs)
        # Reactions with a participant with a 0 stechiometric coeff. triggers a ValueError
        inputs = {"participants":[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=0),
                                  ReactionParticipant(particle_name="B",
                                                   state_name="B",
                                                   coefficient=1)],
                 "pK":1,
                 "reaction_type":"test"}
        
        self.assertRaises(ValueError,
                          Reaction,
                          **inputs)
        # Adding a new participant with a 0 stechiometric coeff. triggers a ValueError
        react_tpl =Reaction(participants=[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=-1),
                                            ReactionParticipant(particle_name="B",
                                                            state_name="B",
                                                            coefficient=1)],
                            pK=1,
                            reaction_type="test")
        inputs={"particle_name": "C",
                "state_name":"C",
                "coefficient":0}
        self.assertRaises(ValueError,
                          react_tpl.add_participant,
                          **inputs)

if __name__ == '__main__':
    ut.main()