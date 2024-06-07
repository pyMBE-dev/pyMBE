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
import sys
import subprocess
import numpy as np
import pandas as pd

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)
script_path = pmb.get_resource("samples/peptide_mixture_grxmc_ideal.py")
data_path = pmb.get_resource("samples/data_peptide_grxmc.csv")

print("*** Grand reaction (G-RxMC) implementation tests ***\n")
print("*** Test that our implementation of the original G-RxMC method reproduces the Henderson-Hasselbalch equation corrected with the Donnan potential (HH+Don) for an ideal mixture of peptides ***")

run_command = [sys.executable, script_path, "--mode", "standard", "--test"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH_Donnan"], rtol=0.01, atol=0.05)
# Check if partition coefficients agree
np.testing.assert_allclose(data["xi_sim"], data["xi_HH_Donnan"], rtol=0.1, atol=0.1)

print("*** Test passed ***\n")


print("*** Test that our implementation of the G-RxMC method with unified ion types reproduces HH+Don for an ideal mixture of peptides ***")

run_command = [sys.executable, script_path, "--mode", "unified", "--test"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH_Donnan"], rtol=0.01, atol=0.05)
# Check if partition coefficients agree
np.testing.assert_allclose(data["xi_sim"], data["xi_HH_Donnan"], rtol=0.1, atol=0.1)

print("*** Test passed ***\n")
