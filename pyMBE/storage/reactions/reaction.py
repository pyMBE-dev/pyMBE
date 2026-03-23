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
#

from typing import List, Dict, Optional
from pydantic import BaseModel, validator, root_validator

class ReactionParticipant(BaseModel):
    """
    Represents one participant in a chemical reaction.

    A reaction participant is defined by a particle name, a specific
    state of that particle, and an integer stoichiometric coefficient.
    Negative coefficients indicate reactants, whereas positive
    coefficients indicate products.

    Attributes:
        particle_name (str):
            The name of the particle template participating in the reaction.
        state_name (str):
            The name of the particle state.
        coefficient (int):
            Stoichiometric coefficient of the participant:
            - ``coefficient ∈ {..., -3, -2, -1}`` → reactant
            - ``coefficient ∈ {1, 2, 3, ...}`` → product

    Notes:
        - Coefficients of zero are forbidden.
        - Together, ``particle_name`` and ``state_name`` identify a unique
          chemical species in the simulation framework.
    """
    particle_name: str
    state_name: str
    coefficient: int

    @validator("coefficient")
    def no_zero_coeff(cls, coefficient):
        """
        Ensures that the stoichiometric coefficient is non-zero.

        Args:
            coefficient ('int'):
                Stoichiometric coefficient of the participant.

        Returns:
            ('int'):
                The validated non-zero coefficient.
        """
        if coefficient == 0:
            raise ValueError("Stoichiometric coefficient cannot be zero.")
        return coefficient

class Reaction(BaseModel):
    """
    Defines a chemical reaction between particle states.

    Attributes:
        name ('str'):
            Unique identifier for the reaction.

        participants ('List[ReactionParticipant]'):
            List of reactants and products with stoichiometric coefficients.
            Must include at least two participants.

        reaction_type ('str'):
            A categorical descriptor of the reaction, such as ``"acid_base"``

        pK ('float'):
            Reaction equilibrium parameter (e.g., pKa, -log K). The meaning
            depends on ``reaction_type``.

        simulation_method ('str', optional):
            Simulation method used to study the reaction.

        metadata ('dict', optional):
            Optional free-form metadata for additional reaction details,
            notes, or model-specific configuration.
    """
    participants: List[ReactionParticipant]
    reaction_type: str 
    pK: float 
    metadata: Optional[Dict] = None
    simulation_method: Optional[str] = None
    name: Optional[str] = None 

    @validator("participants")
    def at_least_two_participants(cls, participants):
        """
        Ensures that the reaction contains at least two participants.

        Args:
            participants ('List[ReactionParticipant]'):
                List of reaction participants.

        Returns:
            ('List[ReactionParticipant]'):
                The validated list of participants.
        """
        if len(participants) < 2:
            raise ValueError("A reaction must have at least 2 participants.")
        return participants

    @classmethod
    def _build_name_from_participants(cls, participants):
        """
        Builds a reaction name from a list of reaction participants.

        Args:
            participants ('List[ReactionParticipant]'):
                List of participants to split into reactants and products.

        Returns:
            ('str'):
                Reaction name in the format ``A + B <-> C + D``.
        """
        reactants = []
        products = []
        for participant in participants:
            if participant.coefficient < 0:
                reactants.append(participant.state_name)
            else:
                products.append(participant.state_name)
        reactants.sort()
        products.sort()
        return f"{' + '.join(reactants)} <-> {' + '.join(products)}"

    @root_validator
    def generate_name(cls, values):
        """
        Automatically generates a reaction name from the participants.

        The name is constructed by separating reactants and products
        based on the sign of their stoichiometric coefficients and
        joining them with a reversible reaction symbol.

        Returns:
            ('dict'):
                Updated model values including the generated reaction name.
        """
        participants = values.get("participants", [])
        values["name"] = cls._build_name_from_participants(participants)
        return values

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add_participant(self, particle_name, state_name, coefficient):
        """
        Adds a new participant to the reaction.

        Args:
            particle_name ('str'):
                Name of the particle template.

            state_name ('str'):
                Name of the particle state.

            coefficient ('int'):
                Stoichiometric coefficient of the participant.
                Must be non-zero.
        """
        new_participant = ReactionParticipant(
            particle_name=particle_name,
            state_name=state_name,
            coefficient=coefficient,
        )
        self.participants.append(new_participant)

        # Explicitly regenerate name after mutation
        self.name = self._build_name_from_participants(self.participants)

    def add_simulation_method(self, simulation_method):
        """
        Sets the simulation method used to model the reaction.

        Args:
            simulation_method ('str'):
                Label identifying the simulation method.
        """
        self.simulation_method = simulation_method
