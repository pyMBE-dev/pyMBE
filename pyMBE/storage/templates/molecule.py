#
# Copyright (C) 2025 pyMBE-dev team
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

from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class MoleculeTemplate(PMBBaseModel):
    """
    Template defining a molecule in pyMBE.

    Attributes:
        pmb_type (str): Fixed type identifier for this template. Always "molecule".
        name (str): Unique name of the molecule template.
        residue_list (List[str]): Ordered list of residue names that make up the molecule.
    """
    pmb_type: str = Field(default="molecule", frozen=True)
    name: str
    residue_list: list[str] 


