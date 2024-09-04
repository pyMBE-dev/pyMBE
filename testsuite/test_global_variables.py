#
# Copyright (C) 2024 pyMBE-dev team
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

import pyMBE
import numpy as np

# Test that two different instances of pyMBE do not share the same memory address for their global attributes
pmb1 = pyMBE.pymbe_library(seed=42)
pmb2 = pyMBE.pymbe_library(seed=43)
np.testing.assert_raises(AssertionError,
                        np.testing.assert_equal,
                        id(pmb1.units),
                        id(pmb2.units))

# Test that redefining the system of reduced units does not create a new pint.UnitRegistry

## Define a variables in the old unit registry
nm_length_1=pmb1.units.Quantity(1,"nm")

## Change the system of reduced units
pmb1.set_reduced_units(unit_length=0.4*pmb1.units.nm)

## Define variable in the new unit registry
nm_length_2=pmb1.units.Quantity(2,"nm")

## operations between old and new quantities should work normally 
np.testing.assert_equal((nm_length_1+nm_length_2).m_as("nm"), 
                        (pmb1.units.Quantity(3,"nm").m_as("nm")))
                        

# Test that set_reduced_units raises a ValueError if the wrong unit is provided
input_parameters={"unit_length": pmb1.units.Quantity(1, "J")}                      

np.testing.assert_raises(ValueError,
                        pmb1.set_reduced_units,
                        **input_parameters)
