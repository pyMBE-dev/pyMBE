from typing import List
from pydantic import Field
from ..base_type import PMBBaseModel

class HydrogelInstance(PMBBaseModel):
    pmb_type: str = Field(default="hydrogel", frozen=True)
    hydrogel_id: int
    name: str
    molecule_ids: List[str] = Field(default_factory=list)
