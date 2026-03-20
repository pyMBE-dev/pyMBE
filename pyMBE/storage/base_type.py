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

from pydantic import BaseModel, Field

class PMBBaseModel(BaseModel):
    """
    Base class for all pyMBE models.

    Provides common fields and validation behavior for pyMBE templates and instances.

    Attributes:
        pmb_type (str): Fixed type identifier. Subclasses must set this to a specific type.
        name (str): Unique name of the model instance or template.

    Config:
        validate_assignment (bool): Ensures that attribute assignments are validated.
        extra (str): Forbids extra attributes not defined in the model.
    """

    pmb_type: str = Field(frozen=True)
    name: str
    
    class Config:
        validate_assignment = True
        extra = "forbid"
