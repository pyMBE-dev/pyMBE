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


class ResidueInstance(PMBBaseModel):
    """
    Instance of a residue placed within a molecule during a simulation.

    ``ResidueInstance`` represents a concrete occurrence of a residue
    derived from a residue template. Each instance is uniquely indexed
    by ``residue_id`` and may optionally belong to a parent molecule,
    such as a peptide, protein, or generic molecule.

    Attributes:
        pmb_type (str):
            Fixed string identifying this object as a residue instance.
            Always ``"residue"``.
        name (str):
            Name of the residue template from which this instance is derived.
        residue_id (int):
            Unique non-negative integer identifying this residue instance
            within the database.
        molecule_id (int | None):
            Identifier of the parent molecule to which this residue belongs.
            ``None`` indicates that the residue is not assigned to any molecule.

    Notes:
        - ``ResidueInstance`` does not itself store particle-level
          information; instead, particles reference the residue via
          ``residue_id``.
        - Residues may be standalone (e.g., in coarse systems) or part of
          polymers, proteins, peptides, or hydrogels.
        - The sequence ordering and topology of residues are encoded at the
          molecule instance/template level, not here.
    """
    
    pmb_type: str = "residue"
    name: str            # residue template name
    residue_id: int
    molecule_id: int | None = None
    
    @field_validator("residue_id")
    def validate_residue_id(cls, rid):
        if rid < 0:
            raise ValueError("residue_id must be a non-negative integer.")
        return rid
