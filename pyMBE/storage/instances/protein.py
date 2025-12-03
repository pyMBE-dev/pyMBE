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


class ProteinInstance(PMBBaseModel):
    """
    Instance of a protein molecule placed in the simulation.

    ``ProteinInstance`` represents a concrete protein object created
    from a protein template defined in the database. Each instance
    corresponds to one full protein chain and is uniquely identified
    by its ``molecule_id``.

    Attributes:
        pmb_type (str):
            Fixed string identifying this object as a protein instance.
            Always ``"protein"``.
        name (str):
            Name of the protein template from which this instance was
            created. This usually corresponds to a user-defined or
            imported protein type or sequence identifier.
        molecule_id (int):
            Unique non-negative integer identifying this protein within
            the database. Assigned by the database manager upon creation.

    Notes:
        - A ``ProteinInstance`` only records the identity of the protein
          and its template association.
        - Residues and particles that belong to the protein reference
          this instance through their ``molecule_id`` values.
        - The structural connectivity (residue sequence, domains) is
          handled at the template level or by the builder modules.
    """
    pmb_type: str = "protein"
    name: str            # molecule template name
    molecule_id: int 
    
    @field_validator("molecule_id")
    def validate_residue_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid