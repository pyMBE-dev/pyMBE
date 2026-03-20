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
from typing import Optional


class PeptideInstance(PMBBaseModel):
    """
    Instance of a peptide molecule placed in the simulation.

    Attributes:
        pmb_type ('str'):
            Fixed string identifying this object as a peptide instance. Always ``"peptide"``.
        
        name ('str'):
            Name of the peptide template from which this instance was created. 

        molecule_id ('int'):
            Unique non-negative integer identifying this peptide within the database. 

        assembly_id ('int' | 'None'):
            Identifier of the super-parent assembly (e.g. hydrogel) to which this residue belongs. ``None`` indicates that the residue is not assigned to any assembly.

    Notes:
        - This class only tracks the identity of the peptide instance. Residues and particles belonging to the peptide reference this instance through their ``molecule_id`` fields.
        - Connectivity (ordering of residues), spatial arrangement, and bonding interactions are managed separately by the database or simulation engine.
    """
    pmb_type: str = "peptide"
    name: str            # molecule template name
    molecule_id: int 
    assembly_id: Optional[int] = None

    @validator("molecule_id")
    def validate_residue_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid
