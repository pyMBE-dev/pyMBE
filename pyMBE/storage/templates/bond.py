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

from typing import Dict, Literal
from ..base_type import PMBBaseModel
from ..pint_quantity import PintQuantity
from pydantic import Field, model_validator


class BondTemplate(PMBBaseModel):
    """
    Template defining a bond in a pyMBE simulation.

    Attributes:
        pmb_type (Literal["bond"]): Fixed type identifier for this template. Always "bond".
        name (str): Unique name of the bond template, e.g., "HARMONIC_default".
        bond_type (str): Type of bond potential. Examples: "HARMONIC", "FENE".
        parameters (Dict[str, PintQuantity]): Dictionary of bond parameters.
            Common keys:
                - "k": Force constant (energy / distance^2)
                - "r0": Equilibrium bond length
                - "d_r_max": Maximum bond extension (for FENE)
        l0 (PintQuantity): Initial bond length when the bond is instantiated.
            
    Notes:
        Values are stored as PintQuantity objects for unit-aware calculations.
    """
    pmb_type: Literal["bond"] = "bond"
    name: str = Field(default="default")
    bond_type: str                      # "HARMONIC", "FENE"
    particle_name1: str | None = None
    particle_name2: str | None = None
    parameters: Dict[str, PintQuantity] # k, r0, d_r_max...

    def _make_name(self):
        """Create a canonical name for the bond."""
        if self.particle_name1 is None or self.particle_name2 is None:
            raise RuntimeError("The BondTemplate has no defined particle_name1 or particle_name2 and therefore the name could not be automatically generated")
        pn1, pn2 = sorted([self.particle_name1, self.particle_name2])
        self.name = f"{pn1}-{pn2}"

    