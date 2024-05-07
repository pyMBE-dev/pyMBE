# Import pyMBE and other libraries
import pyMBE
import subprocess
import numpy as np
import pandas as pd

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)
script_path = pmb.get_resource("samples/branched_polyampholyte.py")
data_path = pmb.get_resource("samples/data_polyampholyte_cph.csv")


print("*** Constant pH (cpH) implementation tests ***\n")
print("*** Test that our implementation of the cpH method reproduces the Henderson Hasselbalch equation for an ideal polyampholyte ***\n")

run_command = ["python3", script_path, "--test"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH"], rtol=0.15, atol=0.2)

print("*** Test passed ***\n")
