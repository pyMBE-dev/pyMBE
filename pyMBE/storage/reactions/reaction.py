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
from pydantic import BaseModel, Field, field_validator, model_validator


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
            The state of the particle (e.g., protonation state, charge state).
        coefficient (int):
            Stoichiometric coefficient of the participant:
            - ``coefficient < 0`` → reactant
            - ``coefficient > 0`` → product

    Notes:
        - Coefficients of zero are forbidden.
        - Together, ``particle_name`` and ``state_name`` identify a unique
          chemical species in the simulation framework.
    """
    particle_name: str
    state_name: str
    coefficient: int

class Reaction(BaseModel):
    """
    Defines a chemical reaction between particle states.

    Attributes:
        name ('str'):
            Unique identifier for the reaction.

        participants ('List[ReactionParticipant]'):
            List of reactants and products with stoichiometric coefficients.
            Must include at least two participants.

        pK ('float'):
            Reaction equilibrium parameter (e.g., pKa, log K). The meaning
            depends on ``reaction_type``.

        reaction_type ('str'):
            A categorical descriptor of the reaction, such as ``"acid_base"``

        simulation_method ('str', optional):
            Simulation method used to study the reaction.

        metadata ('dict', optional):
            Optional free-form metadata for additional reaction details,
            notes, or model-specific configuration.
    """
    participants: List[ReactionParticipant]
    pK: float 
    reaction_type: str 
    metadata: Optional[Dict] = None
    simulation_method: Optional[str] = None
    name: Optional[str] = None 

    @model_validator(mode="after")
    def generate_name(self):
        """Automatically generate reaction name from participants."""
        reactants = []
        products = []

        for p in self.participants:
            species = f"{p.state_name}"
            if p.coefficient < 0:
                reactants.append(species)
            else:
                products.append(species)

        reactants = sorted(reactants)
        products = sorted(products)

        left = " + ".join(reactants)
        right = " + ".join(products)

        # reversible reaction symbol
        self.name = f"{left} <-> {right}"
        return self

    @field_validator("participants")
    def at_least_two_participants(cls, v):
        if len(v) < 2:
            raise ValueError("A reaction must have at least 2 participants.")
        return v

    @field_validator("participants")
    def no_zero_coeff(cls, v):
        for p in v:
            if p.coefficient == 0:
                raise ValueError(f"Participant {p.state_name} has coefficient 0.")
        return v

    def add_participant(self, particle_name, state_name, coefficient):
        """
        Add a new reaction participant to the reaction.

        Args:
            particle_name ('str'):
                Name of the particle participating in the reaction.
            state_name ('str'):
                Specific state of the particle.
            coefficient ('int'):
                Stoichiometric coefficient for the participant.

        """
        if coefficient == 0:
            raise ValueError("Stoichiometric coefficient cannot be zero.")
        new_participant = ReactionParticipant(particle_name=particle_name,
                                            state_name=state_name,
                                            coefficient=coefficient)
        self.participants.append(new_participant)
        self.generate_name()
        
    
    def add_simulation_method(self, simulation_method):
        """
        Adds which simulation is used to simulate the reaction

        Args:
            simulation_method ('str'): label of the simulation method
        """
        self.simulation_method = simulation_method