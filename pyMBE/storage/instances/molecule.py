from pyMBE.storage.base_type import PMBBaseModel
from pydantic import field_validator


class MoleculeInstance(PMBBaseModel):
    pmb_type: str = "molecule"
    name: str            # molecule template name
    molecule_id: int 
    
    @field_validator("molecule_id")
    def validate_residue_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid
