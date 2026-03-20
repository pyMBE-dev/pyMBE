#
# Copyright (C) 2026 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from typing import List
from pydantic import Field, BaseModel, validator
from ..base_type import PMBBaseModel

class HydrogelNode(BaseModel):
    """
    Represents a node in a hydrogel network.

    Attributes:
        particle_name ('str'): 
            Name of the particle at this node.

        lattice_index ('List[int]'): 
            3D lattice position of the node. Must be a list of length 3.
    """
    particle_name: str
    lattice_index: List[int]  # must be length 3
    @validator("lattice_index", pre=True)
    def coerce_lattice_index(cls, v):
        # Accept tuple, list, numpy array, etc.
        try:
            v = list(v)
        except TypeError:
            raise ValueError("lattice_index must be an iterable of 3 integers") 
        return v

class HydrogelChain(BaseModel):
    """
    Represents a polymer chain between two hydrogel nodes.

    Attributes:
        molecule_name ('str'): 
            Name of the molecule representing the polymer chain.

        node_start ('str'): 
            Name of the starting node.

        node_end ('str'): 
            Name of the ending node.
    """
    molecule_name: str
    node_start: str
    node_end: str   
    
class HydrogelTemplate(PMBBaseModel):
    """
    Template defining a hydrogel network in the pyMBE database.

    Attributes:
        pmb_type ('str'): 
            Fixed type identifier for this template. Always "hydrogel".

        name ('str'): 
            Unique name of the hydrogel template.

        node_map ('List[HydrogelNode]'): 
            List of nodes defining the hydrogel lattice.

        chain_map ('List[HydrogelChain]'): 
            List of polymer chains connecting nodes.
    """
    pmb_type: str = Field(default="hydrogel", frozen=True)
    name: str
    node_map: List[HydrogelNode] = Field(default_factory=list)
    chain_map: List[HydrogelChain] = Field(default_factory=list)