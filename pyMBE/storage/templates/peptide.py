from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class PeptideTemplate(PMBBaseModel):
    pmb_type: str = Field(default="peptide", frozen=True)
    name: str
    model: str
    residue_list: list[str] 
    sequence: list[str] 