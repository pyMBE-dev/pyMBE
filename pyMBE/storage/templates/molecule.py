from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class MoleculeTemplate(PMBBaseModel):
    pmb_type: str = Field(default="molecule", frozen=True)
    name: str
    residue_list: list[str] 


