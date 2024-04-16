# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import os
import tempfile
import subprocess
import numpy as np
import pandas as pd

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()
script_path = pmb.get_resource(f"samples/peptide_mixture_grxmc_ideal.py")
data_path = pmb.get_resource(f"samples/data_peptide_grxmc.csv")

print(f"*** Grand reaction (G-RxMC) implementation tests ***\n")
print(f"*** Test that our implementation of the original G-RxMC method reproduces the Henderson Hasselbalch equation corrected with the Donnan potential (HH+Don) for an ideal mixture of peptides ***\n")


run_command = ["python3", script_path, "--mode", "standard", "--test"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH_Donnan"], rtol=0.01, atol=0.05)
# Check if partition coefficients agree
np.testing.assert_allclose(data["xi_sim"], data["xi_HH_Donnan"], rtol=0.1, atol=0.1)

print(f"*** Test passed ***\n")

print(f"*** Test that our implementation of the G-RxMC method with unified ion types reproduces HH+Don for an ideal mixture of peptides ***\n")

run_command = ["python3", script_path, "--mode", "unified", "--test"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH_Donnan"], rtol=0.01, atol=0.05)
# Check if partition coefficients agree
np.testing.assert_allclose(data["xi_sim"], data["xi_HH_Donnan"], rtol=0.1, atol=0.1)

print(f"*** Test passed ***\n")
