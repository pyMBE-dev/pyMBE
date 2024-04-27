# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import os
import tempfile
import subprocess
import numpy as np
import pandas as pd
import argparse

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

# Command line arguments
parser = argparse.ArgumentParser(description='Test of the GCMC implementation in pyMBE.')
parser.add_argument('--mode',
                    type=str,
                    default= "ideal",
                    help='Set if an ideal or interacting system is simulated.')
args = parser.parse_args()

script_path=pmb.get_resource(f"samples/salt_solution_gcmc.py")
salt_concentrations=[0.0001, 0.001, 0.01, 0.1]

rtol=0.05 # relative tolerance
atol=0.0 # absolute tolerance

if args.mode == "ideal":
    print(f"*** Running test for GCMC of salt solution (ideal). ***")
elif args.mode == "interacting":
    print(f"*** Running test for GCMC of salt solution (interacting). ***")
with tempfile.TemporaryDirectory() as time_series_path:
    for c_salt_res in salt_concentrations:
        print(f"c_salt_res = {c_salt_res}")
        run_command=["python3", script_path, "--c_salt_res", str(c_salt_res), "--mode", "ideal", "--output", time_series_path, "--mode", args.mode, "--no_verbose"]
        print(subprocess.list2cmdline(run_command))
        subprocess.check_output(run_command)
    # Analyze all time series
    data=analysis.analyze_time_series(path_to_datafolder=time_series_path)

# Check concentration
test_concentration=np.sort(data["csalt","value"].to_numpy(dtype=float))
ref_concentration=np.sort(data["mean","c_salt"].to_numpy())
np.testing.assert_allclose(test_concentration, ref_concentration, rtol=rtol, atol=atol)
print(f"*** Test was successful ***")
