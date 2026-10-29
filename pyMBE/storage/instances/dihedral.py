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
from pyMBE.storage.base_type import PMBBaseModel
from pydantic import validator

class DihedralInstance(PMBBaseModel):
    """
    Instance representation of a dihedral (torsion) between four particles.

    Attributes:
        pmb_type ('Literal["dihedral"]'):
            Fixed identifier set to ``"dihedral"`` for all dihedral instances.

        dihedral_id ('int'):
            Unique non-negative integer identifying this dihedral instance.

        name ('str'):
            Name of the dihedral template from which this instance was created.

        particle_id1 ('int'):
            ID of the first particle in the quadruplet.

        particle_id2 ('int'):
            ID of the second particle in the quadruplet.

        particle_id3 ('int'):
            ID of the third particle in the quadruplet.

        particle_id4 ('int'):
            ID of the fourth particle in the quadruplet.

        added_to_engine ('bool'):
            Whether the dihedral has already been transferred to the simulation engine.
    """
    pmb_type: Literal["dihedral"] = "dihedral"
    dihedral_id: int
    name: str
    particle_id1: int
    particle_id2: int
    particle_id3: int
    particle_id4: int
    added_to_engine: bool = False

    @validator("dihedral_id", "particle_id1", "particle_id2", "particle_id3", "particle_id4")
    def validate_non_negative_int(cls, value, field):
        if value < 0:
            raise ValueError(f"{field.name} must be a non-negative integer.")
        return value
