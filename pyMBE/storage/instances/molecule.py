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

class MoleculeInstance(PMBBaseModel):
    """
    Persistent instance representation of a molecule.

    Attributes:
        pmb_type ('Literal["molecule"]'):
            Fixed string identifying this object as a molecule instance.  Always ``"molecule"``.
        
        name ('str'):
            Name of the molecule **template** from which this instance was created. This must correspond to an existing ``MoleculeTemplate`` in the database.

        molecule_id ('int'):
            Unique non-negative integer identifying this molecule instance within the database.

        assembly_id (Optional[int]):
            Identifier of the super-parent assembly (e.g. hydrogel) to which this molecule belongs. ``None`` indicates that the molecule is not assigned to any assembly.

    Notes:
        - Validation of whether ``name`` corresponds to a registered  molecule template is performed at the database level.
        - Structural or connectivity information (e.g., residue ordering) is maintained outside this class in the instance registry.
    """

    pmb_type: Literal["molecule"] = "molecule"
    name: str            # molecule template name
    molecule_id: int 
    assembly_id: Optional[int] = None
    
    @validator("molecule_id")
    def validate_residue_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid
