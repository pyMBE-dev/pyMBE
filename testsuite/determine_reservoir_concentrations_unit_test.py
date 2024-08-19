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

# Import pyMBE and other libraries
import pyMBE
import numpy as np

def determine_reservoir_concentrations_test(pH_res, c_salt_res):

    # Create an instance of the pyMBE library
    pmb = pyMBE.pymbe_library(seed=42)

    # Determine the reservoir composition using pyMBE
    cH_res_pyMBE, cOH_res_pyMBE, cNa_res_pyMBE, cCl_res_pyMBE = pmb.determine_reservoir_concentrations(
            pH_res,
            c_salt_res * pmb.units.mol/ pmb.units.L,
            lambda x: 1.0)

    # Determine the reservoir concentrations without pyMBE
    cH_res = 10**(-pH_res) * pmb.units.mol/ pmb.units.L
    cOH_res = 10**(pH_res-14) * pmb.units.mol/ pmb.units.L
    c_salt_res = c_salt_res * pmb.units.mol/ pmb.units.L
    cNa_res = max(c_salt_res, c_salt_res + cOH_res - cH_res) 
    cCl_res = max(c_salt_res, c_salt_res + cH_res - cOH_res) 

    # Check the determine equilibrium constants
    np.testing.assert_allclose(cH_res, cH_res_pyMBE)
    np.testing.assert_allclose(cOH_res, cOH_res_pyMBE)
    np.testing.assert_allclose(cNa_res, cNa_res_pyMBE)
    np.testing.assert_allclose(cCl_res, cCl_res_pyMBE)

print("*** Unit test: check that pyMBE correctly calculates the reservoir composition when setting up the grand-reaction method. ***")
for pH_res in np.linspace(1.0, 13.0, 5):
    for c_salt_res in np.logspace(-5, 0, 5):
        determine_reservoir_concentrations_test(pH_res, c_salt_res)
print("*** Unit test passed ***")
