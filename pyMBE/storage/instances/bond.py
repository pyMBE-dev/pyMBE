from pyMBE.storage.base_type import PMBBaseModel
from pydantic import field_validator

class BondInstance(PMBBaseModel):
    pmb_type: str = "bond"
    bond_id: int
    name : str            # bond template name
    particle_id1: int
    particle_id2: int 
    

    @field_validator("bond_id")
    def validate_bond_id(cls, bid):
        if bid < 0:
            raise ValueError("bond_id must be a non-negative integer.")
        return bid