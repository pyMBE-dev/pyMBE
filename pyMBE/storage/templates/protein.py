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

class ProteinTemplate(PMBBaseModel):
    """
    Template defining a protein in the pyMBE database.

    Attributes:
        pmb_type ('Literal["protein"]'):
            Fixed type identifier. Always "protein".

        name ('str'): 
            Unique name of the protein template.

        model ('str'): 
            Name or type of the model used for this protein.

        residue_list ('List[str]'): 
            Ordered list of residue names that compose the protein.

        sequence ('List[str]'): 
            Ordered sequence of residues representing the protein's structure.
    """
    pmb_type: Literal["protein"] = "protein"
    name: str
    model: str
    residue_list: List[str]
    sequence: str 
