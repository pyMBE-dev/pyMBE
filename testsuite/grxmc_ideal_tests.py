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

print(f"*** Running G-RxMC tests (ideal) ***\n")


print(f"*** Running test for standard G-RxMC ***\n")

run_command = ["python3", script_path, "--mode", "standard", "--no_plot"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
print(data)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH_Donnan"], rtol=0.01, atol=0.05)
# Check if partition coefficients agree
np.testing.assert_allclose(data["xi_sim"], data["xi_HH_Donnan"], rtol=0.1, atol=0.1)

print(f"*** Test passed ***\n")


print(f"*** Running test for unified G-RxMC ***\n")

run_command = ["python3", script_path, "--mode", "unified", "--no_plot"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
print(data)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH_Donnan"], rtol=0.01, atol=0.05)
# Check if partition coefficients agree
np.testing.assert_allclose(data["xi_sim"], data["xi_HH_Donnan"], rtol=0.1, atol=0.1)

print(f"*** Test passed ***\n")
