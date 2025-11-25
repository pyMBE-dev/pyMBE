from typing import List, Dict, Optional
from pydantic import BaseModel, Field, field_validator


class ReactionParticipant(BaseModel):
    """
    One participant in a reaction.
    coefficient < 0  -> reactant
    coefficient > 0  -> product
    """
    particle_name: str
    state_name: str
    coefficient: int

class Reaction(BaseModel):
    name: str
    participants: List[ReactionParticipant]
    pK: float = Field(..., description="pKa, logK, eq constant, etc.")
    reaction_type: str = Field(..., description="acid_base, binding, redox, ...")
    metadata: Optional[Dict] = None

    @field_validator("participants")
    def at_least_two_participants(cls, v):
        if len(v) < 2:
            raise ValueError("A reaction must have at least 2 participants.")
        return v

    @field_validator("participants")
    def no_zero_coeff(cls, v):
        for p in v:
            if p.coefficient == 0:
                raise ValueError(f"Participant {p.name} has coefficient 0.")
        return v
