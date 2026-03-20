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

from typing import List, Literal
from pyMBE.storage.base_type import PMBBaseModel

class ResidueTemplate(PMBBaseModel):
    """
    Template defining a residue in the pyMBE database.

    Attributes:
        pmb_type (Literal["residue"]): Fixed type identifier. Always "residue".
        name (str): Unique name of the residue template.
        central_bead (str): Name of the central bead representing the residue.
        side_chains (List[str]): List of side-chain names attached to the central bead.
            Defaults to an empty list if no side chains are present.
    """
    pmb_type: Literal["residue"] = "residue"
    name: str
    central_bead: str 
    side_chains: List[str] = []
                
