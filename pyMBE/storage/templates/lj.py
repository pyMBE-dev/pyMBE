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

from pydantic import BaseModel, Field, root_validator
from ..pint_quantity import PintQuantity


class LJInteractionTemplate(BaseModel):
    """
    Template representing the Lennard–Jones (LJ) interaction parameters
    between two particle *states* stored in the pyMBE database.

    Attributes:
        pmb_type ('str'):
            Fixed identifier for the template type. Always ``"lj"``.

        state1 ('str'):
            Name of the first particle state in the pair.

        state2 ('str'):
            Name of the second particle state in the pair.

        sigma ('PintQuantity'):
            Lennard–Jones σ parameter (distance scale) after applying the combining rule.

        epsilon ('PintQuantity'):
            Lennard–Jones ε parameter (energy scale) after combining.

        cutoff ('PintQuantity'):
            Cutoff radius for the interaction.

        offset ('PintQuantity'):
            Offset applied to the potential (ESPResSo parameter).

        shift ('str | PintQuantity'):
            Shift applied at the cutoff. May be ``"auto"`` or a PintQuantity value.

        name ('str'):
            Auto-generated unique identifier for the interaction, built from ``state1`` and ``state2`` in alphabetical order. Cannot be set  manually by the user.

    Notes:
        - The order of ``state1`` and ``state2`` does **not** matter. The name is always generated as ``"min(state1, state2)-max(state1, state2)"``.
    """
    pmb_type: str = "lj"
    name: str = Field(default="", description="Automatically generated name")
    state1: str
    state2: str
    sigma: PintQuantity
    epsilon: PintQuantity
    cutoff: PintQuantity
    offset: PintQuantity
    shift: str | PintQuantity

    @classmethod
    def _make_name(cls, state1: str, state2: str) -> str:
        """
        Creates a canonical interaction name from two particle states.

        Args:
            state1 ('str'):
                Name of the first particle state.

            state2 ('str'):
                Name of the second particle state.

        Returns:
            ('str'):
                Canonical interaction name in the form ``"A-B"``,
                where ``A`` and ``B`` are sorted alphabetically.
        """
        s1, s2 = sorted([state1, state2])
        return f"{s1}-{s2}"

    # ------------------------------------------------------------------
    # Validators
    # ------------------------------------------------------------------

    @root_validator
    def _auto_generate_name(cls, values):
        """
        Automatically generates and enforces a standardized interaction name.

        The name is derived from ``state1`` and ``state2`` and overrides
        any manually provided value.

        Returns:
            ('dict'):
                Updated model values with the generated ``name`` field.
        """
        state1 = values.get("state1")
        state2 = values.get("state2")

        if state1 is not None and state2 is not None:
            values["name"] = cls._make_name(state1, state2)

        return values
