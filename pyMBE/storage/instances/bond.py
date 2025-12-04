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

class BondInstance(PMBBaseModel):
    """
    Instance representation of a bond between two particles.

    A ``BondInstance`` links two particle instances using a specified
    bond template. This class stores only lightweight, serializable
    identifiers (not Espresso objects or interaction handles), ensuring
    that the object can be safely persisted, exported, and reloaded from
    CSV or other storage formats.

    Attributes:
        pmb_type (str):
            Fixed identifier set to ``"bond"`` for all bond instances.
        bond_id (int):
            Unique non-negative integer identifying this bond instance.
        name (str):
            Name of the bond template from which this instance was created.
        particle_id1 (int):
            ID of the first particle involved in the bond.
        particle_id2 (int):
            ID of the second particle involved in the bond.

    Validators:
        validate_bond_id:
            Ensures that ``bond_id`` is a non-negative integer.

    Notes:
        - ``particle_id1`` and ``particle_id2`` must correspond to
          particle instance IDs already registered in the database.
        - This class does **not** store simulation engine–specific
          objects (e.g., Espresso bond handles). Those should be created
          by a runtime builder separate from the persistent database.
    """
    pmb_type: str = "bond"
    bond_id: int
    name : str            # bond template name
    particle_id1: int
    particle_id2: int 
    

    @field_validator("bond_id")
    def validate_bond_id(cls, bid):
        if bid < 0:
            raise ValueError("bond_id must be a non-negative integer.")
        return bid