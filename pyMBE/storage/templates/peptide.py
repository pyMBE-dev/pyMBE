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

from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class PeptideTemplate(PMBBaseModel):
    """
    Template defining a peptide in a pyMBE simulation.

    Attributes:
        pmb_type (str): Fixed type identifier. Always "peptide".
        name (str): Unique name of the peptide template.
        model (str): Name or type of the model used for this peptide.
        residue_list (List[str]): Ordered list of residue names that make up the peptide.
        sequence (List[str]): Ordered sequence of residues representing the peptide's structure.
    """
    pmb_type: str = Field(default="peptide", frozen=True)
    name: str
    model: str
    residue_list: list[str] 
    sequence: str 