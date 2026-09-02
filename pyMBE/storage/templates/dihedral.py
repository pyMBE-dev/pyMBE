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

from typing import Dict, List, Literal, Optional
from ..base_type import PMBBaseModel
from ..pint_quantity import PintQuantity
from pydantic import Field

class DihedralTemplate(PMBBaseModel):
    """
    Template defining a dihedral (torsion) potential in the pyMBE database.

    Attributes:
        pmb_type ('Literal["dihedral"]'):
            Fixed type identifier for this template. Always "dihedral".

        name ('str'):
            Unique name of the dihedral template, e.g., "A-B-C-D" where the
            quadruplet is canonicalised under full reversal (see `make_dihedral_key`).

        dihedral_type ('str'):
            Type of dihedral potential. Supported: "periodic" (the ESPResSo
            ``Dihedral`` bond ``V(phi) = k*[1 - cos(n*phi - phi_0)]``). "tabulated"
            is reserved for force fields that need Ryckaert-Bellemans or improper
            dihedrals and is not implemented yet.

        particle_name1 ('Optional[str]'):
            Name of the first particle in the dihedral quadruplet.

        particle_name2 ('Optional[str]'):
            Name of the second particle in the dihedral quadruplet.

        particle_name3 ('Optional[str]'):
            Name of the third particle in the dihedral quadruplet.

        particle_name4 ('Optional[str]'):
            Name of the fourth particle in the dihedral quadruplet.

        terms ('List[Dict[str, PintQuantity]]'):
            List of dihedral terms. GROMACS-style force fields stack several
            ``[dihedrals]`` lines (different multiplicities) on the same
            quadruplet, so a dihedral template holds a list of terms. For
            ``dihedral_type == "periodic"`` each term is a dictionary with keys
            ``k`` (energy), ``phi_0`` (dimensionless, radians) and ``mult``
            (dimensionless, the integer multiplicity).

    Notes:
        - Values of the parameters are stored as PintQuantity objects for
          unit-aware calculations. The multiplicity ``mult`` is stored as a
          dimensionless PintQuantity so that all terms round-trip through the
          same serialisation path.
        - The canonical name is invariant under full reversal of the
          quadruplet: ``make_dihedral_key("A", "B", "C", "D")`` equals
          ``make_dihedral_key("D", "C", "B", "A")``.
        - GROMACS proper dihedral type 1, ``V(phi) = k*[1 + cos(n*phi - phi_s)]``,
          maps onto the ``"periodic"`` type with ``k`` unchanged, ``mult = n``
          and ``phi_0 = phi_s + pi`` (``phi_s`` converted from degrees to radians).
    """
    pmb_type: Literal["dihedral"] = "dihedral"
    name: str = Field(default="default")
    dihedral_type: str
    particle_name1: Optional[str] = None
    particle_name2: Optional[str] = None
    particle_name3: Optional[str] = None
    particle_name4: Optional[str] = None
    terms: List[Dict[str, PintQuantity]]

    @classmethod
    def make_dihedral_key(cls, name1, name2, name3, name4):
        """Return a canonical name for a dihedral between four particle names.

        A dihedral ``(a, b, c, d)`` describes the same torsion as its full
        reversal ``(d, c, b, a)``. The canonical name is the lexicographically
        smaller of the two orderings, joined by dashes.

        Args:
            name1 ('str'):
                Name of the first particle.

            name2 ('str'):
                Name of the second particle.

            name3 ('str'):
                Name of the third particle.

            name4 ('str'):
                Name of the fourth particle.

        Returns:
            ('str'):
                Canonical dihedral name, e.g. "A-B-C-D".
        """
        forward = (name1, name2, name3, name4)
        reverse = (name4, name3, name2, name1)
        n1, n2, n3, n4 = min(forward, reverse)
        return f"{n1}-{n2}-{n3}-{n4}"

    def _make_name(self):
        """Create canonical name using particle names."""
        if (not self.particle_name1 or not self.particle_name2
                or not self.particle_name3 or not self.particle_name4):
            raise RuntimeError("Cannot generate dihedral name: one or more of "
                               "particle_name1..4 is missing.")
        self.name = self.make_dihedral_key(self.particle_name1,
                                           self.particle_name2,
                                           self.particle_name3,
                                           self.particle_name4)

    def get_terms(self, ureg):
        """
        Retrieve the dihedral terms as dictionaries of Pint `Quantity` objects.

        Args:
            ureg ('pint.UnitRegistry'):
                Pint unit registry used to reconstruct physical quantities from storage.

        Returns:
            'List[Dict[str, pint.Quantity]]':
                A list with one dictionary per term, mapping parameter names to
                their corresponding unit-aware Pint quantities.
        """
        pint_terms = []
        for term in self.terms:
            pint_terms.append({key: value.to_quantity(ureg)
                               for key, value in term.items()})
        return pint_terms
