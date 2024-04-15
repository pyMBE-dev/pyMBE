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
script_path = pmb.get_resource(f"samples/branched_polyampholyte.py")
data_path = pmb.get_resource(f"samples/data_polyampholyte_cph.csv")


print(f"*** Running test for cpH (ideal) ***\n")

run_command = ["python3", script_path, "--no_plot"]
subprocess.check_output(run_command)

data = pd.read_csv(data_path)
# Check if charges agree
np.testing.assert_allclose(data["Z_sim"], data["Z_HH"], rtol=0.05, atol=0.1)

print(f"*** Test passed ***\n")
