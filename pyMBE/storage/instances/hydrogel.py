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

from typing import Literal
from ..base_type import PMBBaseModel
from pydantic import validator

class HydrogelInstance(PMBBaseModel):
    """
    Persistent instance representation of a hydrogel object.

    Attributes:
        pmb_type ('Literal["hydrogel"]'):
            Fixed string identifier for this instance type. Always ``"hydrogel"``.

        assembly_id ('int'):
            Unique non-negative integer identifying this hydrogel instance.

        name ('str'):
            Human-readable name for the hydrogel (e.g., ``"HG_001"``).

    Notes:
        - This class represents the *instance* level (what specific
          hydrogel exists in the system), not a template describing generic
          hydrogel types.
    """
    pmb_type: Literal["hydrogel"] = "hydrogel"
    assembly_id: int
    name: str
    @validator("assembly_id")
    def validate_assembly_id(cls, aid):
        if aid < 0:
            raise ValueError("assembly_id must be a non-negative integer.")
        return aid
