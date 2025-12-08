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
from pydantic import field_validator


class MoleculeInstance(PMBBaseModel):
    """
    Persistent instance representation of a molecule.

    A ``MoleculeInstance`` links a concrete molecule in the system to a
    molecule template (through its ``name``) and assigns it a unique
    integer identifier. Molecule instances typically serve as containers
    for ordered lists of residue instances, which are managed in the
    database layer outside of this class.

    This class is intentionally minimal and fully serializable. It stores
    no engine-specific data or structural objects.

    Attributes:
        pmb_type (str):
            Fixed string identifying this object as a molecule instance.
            Always ``"molecule"``.
        name (str):
            Name of the molecule **template** from which this instance
            was created. This must correspond to an existing
            ``MoleculeTemplate`` in the database.
        molecule_id (int):
            Unique non-negative integer identifying this molecule
            instance within the database.
        assembly_id (int | None):
            Identifier of the super-parent assembly (e.g. hydrogel) to which this residue belongs.
            ``None`` indicates that the residue is not assigned to any assembly.

    Notes:
        - Validation of whether ``name`` corresponds to a registered
          molecule template is performed at the database level.
        - Structural or connectivity information (e.g., residue ordering)
          is maintained outside this class in the instance registry.
    """

    pmb_type: str = "molecule"
    name: str            # molecule template name
    molecule_id: int 
    assembly_id: int | None = None
    
    @field_validator("molecule_id")
    def validate_residue_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid
