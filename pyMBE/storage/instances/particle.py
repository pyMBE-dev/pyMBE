from pydantic import field_validator
from ..base_type import PMBBaseModel


class ParticleInstance(PMBBaseModel):
    """
    A placed particle within the simulation.
    """
    pmb_type: str = "particle"
    particle_id: int
    initial_state: str
    residue_id: int | None = None
    molecule_id: int | None = None

    @field_validator("particle_id")
    def validate_particle_id(cls, pid):
        if pid < 0:
            raise ValueError("particle_id must be a non-negative integer.")
        return pid
