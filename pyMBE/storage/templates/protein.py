from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class ProteinTemplate(PMBBaseModel):
    pmb_type: str = Field(default="protein", frozen=True)
    name: str
    model: str
    residue_list: list[str] 
    sequence: list[str] 