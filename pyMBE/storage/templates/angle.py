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

from typing import Dict, Literal, Optional
from ..base_type import PMBBaseModel
from ..pint_quantity import PintQuantity
from pydantic import Field

class AngleTemplate(PMBBaseModel):
    """
    Template defining an angular potential in the pyMBE database.

    Attributes:
        pmb_type ('Literal["angle"]'):
            Fixed type identifier for this template. Always "angle".

        name ('str'):
            Unique name of the angle template, e.g., "A-B-C" where B is the central particle.

        angle_type ('str'):
            Type of angle potential. Examples: "harmonic", "cosine", "harmonic_cosine".

        side_particle1 ('Optional[str]'):
            Name of the first side particle in the angle triplet.

        central_particle ('Optional[str]'):
            Name of the central particle in the angle triplet.

        side_particle2 ('Optional[str]'):
            Name of the second side particle in the angle triplet.

        parameters ('Dict[str, PintQuantity]'):
            Dictionary of angle parameters (e.g., k, phi_0).

    Notes:
        - Values of the parameters are stored as PintQuantity objects for unit-aware calculations.
        - The canonical name sorts the two side particles alphabetically while keeping
          the central particle in the middle: ``"side_min-central-side_max"``.
    """
    pmb_type: Literal["angle"] = "angle"
    name: str = Field(default="default")
    angle_type: str
    side_particle1: Optional[str] = None
    central_particle: Optional[str] = None
    side_particle2: Optional[str] = None
    parameters: Dict[str, PintQuantity]

    @classmethod
    def make_angle_key(cls, side1, central, side2):
        """Return a canonical name for an angle between three particle names.

        The two side particles are sorted alphabetically, with the central
        particle kept in the middle position.

        Args:
            side1 ('str'):
                Name of the first side particle.

            central ('str'):
                Name of the central particle.

            side2 ('str'):
                Name of the second side particle.

        Returns:
            ('str'):
                Canonical angle name, e.g. "A-B-C".
        """
        sides = sorted([side1, side2])
        return f"{sides[0]}-{central}-{sides[1]}"

    def _make_name(self):
        """Create canonical name using particle names."""
        if not self.side_particle1 or not self.central_particle or not self.side_particle2:
            raise RuntimeError("Cannot generate angle name: side_particle1, central_particle, or side_particle2 missing.")
        self.name = self.make_angle_key(self.side_particle1, self.central_particle, self.side_particle2)

    def get_parameters(self, ureg):
        """
        Retrieve the angle parameters as Pint `Quantity` objects.

        Args:
            ureg ('pint.UnitRegistry'):
                Pint unit registry used to reconstruct physical quantities from storage.

        Returns:
            'Dict[str, pint.Quantity]':
                A dictionary mapping parameter names to their corresponding unit-aware Pint quantities.
        """
        pint_parameters = {}
        for parameter in self.parameters.keys():
            pint_parameters[parameter] = self.parameters[parameter].to_quantity(ureg)
        return pint_parameters
