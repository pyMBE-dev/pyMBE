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
            
    Notes:
        Values are stored as PintQuantity objects for unit-aware calculations.
    """
    pmb_type: Literal["bond"] = "bond"
    name: str = Field(default="default")
    bond_type: str                      # "HARMONIC", "FENE"
    particle_name1: str | None = None
    particle_name2: str | None = None
    parameters: Dict[str, PintQuantity] # k, r0, d_r_max...

    @classmethod
    def make_bond_key(cls, pn1, pn2):
        """Return a canonical name for a bond between two particle names.

        Args:
            pn1 (str): Name of the first particle.
            pn2 (str): Name of the second particle.

        Returns:
            str: Canonical bond name, e.g. "A-B".
        """
        return "-".join(sorted([pn1, pn2]))

    def _make_name(self):
        """Create canonical name using particle names."""
        if not self.particle_name1 or not self.particle_name2:
            raise RuntimeError(
                "Cannot generate bond name: particle_name1 or particle_name2 missing."
            )

        self.name = self.make_bond_key(self.particle_name1, self.particle_name2)

    def get_parameters(self, ureg):
        """
        Retrieve the bond parameters as Pint `Quantity` objects.

        Args:
            ureg (pint.UnitRegistry) : Pint unit registry used to reconstruct physical quantities from storage.

        Returns:
            Dict[str, pint.Quantity]:
                A dictionary mapping parameter names to their corresponding
                unit-aware Pint quantities.

        Example:
            >>> bt = BondTemplate(
            ...     bond_type="harmonic",
            ...     particle_name1="A",
            ...     particle_name2="B",
            ...     parameters={"k": PintQuantity("100 kJ/mol/nm^2"),
            ...                 "r0": PintQuantity("0.3 nm")}
            ... )
            >>> bt.get_parameters()
            {'k': <Quantity(100, 'kilojoule / mole / nanometer ** 2')>,
             'r0': <Quantity(0.3, 'nanometer')>}
        """
        pint_parameters={}
        for parameter in self.parameters.keys():
            pint_parameters[parameter] = self.parameters[parameter].to_quantity(ureg)
        return pint_parameters