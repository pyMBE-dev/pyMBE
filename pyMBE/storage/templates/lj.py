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

from pydantic import BaseModel, Field, model_validator
from ..pint_quantity import PintQuantity


class LJInteractionTemplate(BaseModel):
    """
    Template representing the Lennard–Jones (LJ) interaction parameters
    between two particle *states*.

    The template **always generates the interaction name automatically**
    from the two provided state names, ensuring standardized naming and
    preventing inconsistencies between different LJ entries.

    The LJ parameters stored here correspond to the *final effective*
    values after applying the combining rule (e.g., Lorentz–Berthelot).
    This allows users to inspect, validate, or export the exact values
    that will be passed to the simulation engine.

    Attributes:
        pmb_type (str):
            Fixed identifier for the template type. Always ``"lj"``.
        state1 (str):
            Name of the first particle state in the pair.
        state2 (str):
            Name of the second particle state in the pair.
        sigma (PintQuantity):
            Lennard–Jones σ parameter (distance scale) after applying
            the combining rule.
        epsilon (PintQuantity):
            Lennard–Jones ε parameter (energy scale) after combining.
        cutoff (PintQuantity):
            Cutoff radius for the interaction.
        offset (PintQuantity):
            Offset applied to the potential (ESPResSo parameter).
        shift (str | PintQuantity):
            Shift applied at the cutoff. May be ``"auto"`` or a PintQuantity value.
        name (str):
            Auto-generated unique identifier for the interaction, built from
            ``state1`` and ``state2`` in alphabetical order. Cannot be set
            manually by the user.

    Notes:
        - The order of ``state1`` and ``state2`` does **not** matter.
          The name is always generated as ``"min(state1, state2)-max(state1, state2)"``.

    Examples:
        Creating an LJ interaction:

        >>> LJInteractionTemplate(
        ...     state1="HA",
        ...     state2="A-",
        ...     sigma=sigma,
        ...     epsilon=epsilon,
        ...     cutoff=cutoff,
        ...     offset=offset,
        ...     shift="auto",
        ... )
        
        Interaction between ``"L"`` and ``"W"`` results in:

        >>> LJInteractionTemplate(
        ...     state1="W",
        ...     state2="L",
        ...     ...
        ... ).name
        'L-W'
    """

    pmb_type: str = "lj"
    name: str = Field(default="", description="Automatically generated name")

    state1: str
    state2: str

    sigma: PintQuantity
    epsilon: PintQuantity
    cutoff: PintQuantity
    offset: PintQuantity
    shift: str | float


    @classmethod
    def _make_name(cls, state1: str, state2: str) -> str:
        """Create a canonical name from two states."""
        s1, s2 = sorted([state1, state2])
        return f"{s1}-{s2}"

    @model_validator(mode="after")
    def _auto_generate_name(self):
        """Enforce standardized automatic name."""
        object.__setattr__(self, "name", self._make_name(self.state1, self.state2))
        return self
