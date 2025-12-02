from typing import List
from pydantic import Field, BaseModel
from ..base_type import PMBBaseModel


class HydrogelNode(BaseModel):
    particle_name: str
    lattice_index: List[int]  # must be length 3


class HydrogelChain(BaseModel):
    node_start: str
    node_end: str
    residue_list: List[str]   # list of residue names


class HydrogelTemplate(PMBBaseModel):
    """
    A hydrogel definition consists of:
    - node_map: list of nodes with particle names and lattice positions
    - chain_map: list of node-node polymer chains with residue lists
    """
    pmb_type: str = Field(default="hydrogel", frozen=True)
    name: str

    node_map: List[HydrogelNode] = Field(default_factory=list)
    chain_map: List[HydrogelChain] = Field(default_factory=list)

