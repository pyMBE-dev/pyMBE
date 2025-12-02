from pyMBE.storage.base_type import PMBBaseModel
from pydantic import field_validator


class ResidueInstance(PMBBaseModel):
    pmb_type: str = "residue"
    name: str            # residue template name
    residue_id: int
    molecule_id: int | None = None
    
    @field_validator("residue_id")
    def validate_residue_id(cls, rid):
        if rid < 0:
            raise ValueError("residue_id must be a non-negative integer.")
        return rid
