#
# Copyright (C) 2025 pyMBE-dev team
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
from pydantic import BaseModel, Field, field_validator


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

    A ``Reaction`` object captures the stoichiometry and thermodynamic
    properties of a chemical equilibrium. 
    This can represent phenomena such as acid–base reactionsor any multi-species reaction scheme 
    supported by the simulation engine.

    Attributes:
        name (str):
            Unique identifier for the reaction.
        participants (List[ReactionParticipant]):
            List of reactants and products with stoichiometric coefficients.
            Must include at least two participants.
        pK (float):
            Reaction equilibrium parameter (e.g., pKa, log K). The meaning
            depends on ``reaction_type``.
        reaction_type (str):
            A categorical descriptor of the reaction, such as ``"acid_base"``
        metadata (Optional[Dict]):
            Optional free-form metadata for additional reaction details,
            notes, or model-specific configuration.

    Validation:
        - At least one participant are required.
        - All participants must have non-zero stoichiometric coefficients.

    Examples:
        Acid dissociation of HA:
            HA ↔ H⁺ + A⁻

        Represented as:
            Reaction(
                name="acid_dissociation",
                participants=[
                    ReactionParticipant("A", "HA", -1),
                    ReactionParticipant("A", "A-", 1),
                    ReactionParticipant("H", "H+", 1),
                ],
                pK=4.75,
                reaction_type="acid_base",
            )
    """
    name: str
    participants: List[ReactionParticipant]
    pK: float = Field(..., description="pKa, logK, eq constant, etc.")
    reaction_type: str = Field(..., description="acid_base, binding, redox, ...")
    metadata: Optional[Dict] = None

    @field_validator("participants")
    def at_least_two_participants(cls, v):
        if len(v) < 2:
            raise ValueError("A reaction must have at least 1 participant.")
        return v

    @field_validator("participants")
    def no_zero_coeff(cls, v):
        for p in v:
            if p.coefficient == 0:
                raise ValueError(f"Participant {p.name} has coefficient 0.")
        return v
