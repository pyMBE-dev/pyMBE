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

from typing import List
from pydantic import Field
from ..base_type import PMBBaseModel

class HydrogelInstance(PMBBaseModel):
    """
    Persistent instance representation of a hydrogel object.

    A ``HydrogelInstance`` stores the high-level composition of a
    hydrogel in terms of the constituent polymer chain molecules.  
    Each hydrogel is assigned a unique integer ID and has a human-readable
    name, along with a list of molecule identifiers referencing previously
    registered molecule instances.

    This class is intentionally lightweight and fully serializable.  
    It does **not** store simulation-engine internal objects
    (such as lattice builders, Espresso handles, network topologies, etc.).
    These are expected to be constructed externally at run time.

    Attributes:
        pmb_type (str):
            Fixed string identifier for this instance type. Always
            ``"hydrogel"``.
        assembly_id (int):
            Unique non-negative integer identifying this hydrogel instance.
        name (str):
            Human-readable name for the hydrogel (e.g., ``"HG_001"``).

    Notes:
        - This class represents the *instance* level (what specific
          hydrogel exists in the system), not a template describing generic
          hydrogel types.
        - The integrity of ``molecule_ids`` (e.g., references to existing
          molecule instances) should be validated in the database layer
          during creation or update and not inside this class.
    """
    pmb_type: str = Field(default="hydrogel", frozen=True)
    assembly_id: int
    name: str
