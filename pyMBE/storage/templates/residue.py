from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class ResidueTemplate(PMBBaseModel):
    pmb_type: str = Field(default="residue", frozen=True)
    name: str
    central_bead: str 
    side_chains: list[str] = []
                