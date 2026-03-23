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
from pydantic import validator
from typing import Literal, Optional


class ProteinInstance(PMBBaseModel):
    """
    Instance of a protein molecule placed in the simulation.

    Attributes:
        pmb_type ('Literal["protein"]'):
            Fixed string identifying this object as a protein instance. Always ``"protein"``.
        
        name ('str'):
            Name of the protein template from which this instance was created. 

        molecule_id ('int'):
            Unique non-negative integer identifying this protein within the database. 

        assembly_id ('Optional[int]'):
            Identifier of the super-parent assembly (e.g. hydrogel) to which this residue belongs. ``None`` indicates that the residue is not assigned to any assembly.

    Notes:
        - A ``ProteinInstance`` only records the identity of the protein and its template association.
        - Residues and particles that belong to the protein reference this instance through their ``molecule_id`` values.
        - The structural connectivity (residue sequence, domains) is  handled at the template level or by the builder modules.
    """
    pmb_type: Literal["protein"] = "protein"
    name: str            # molecule template name
    molecule_id: int 
    assembly_id: Optional[int] = None
    
    @validator("molecule_id")
    def validate_residue_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid
