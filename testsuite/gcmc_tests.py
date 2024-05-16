# Import pyMBE and other libraries
import sys
import tempfile
import subprocess
import pyMBE
from lib import analysis
import numpy as np

# Template of the test

def gcmc_test(script_path, mode):
    if mode == "ideal":
        print("*** Running test for GCMC of salt solution (ideal). ***")
    elif mode == "interacting":
        print("*** Running test for GCMC of salt solution (interacting). ***")
    with tempfile.TemporaryDirectory() as time_series_path:
        for c_salt_res in salt_concentrations:
            print(f"c_salt_res = {c_salt_res}")
            run_command=[sys.executable, script_path, "--c_salt_res", str(c_salt_res), "--output", time_series_path, "--mode", mode, "--no_verbose"]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path)

    # Check concentration
    test_concentration=np.sort(data["csalt","value"].to_numpy(dtype=float))
    ref_concentration=np.sort(data["mean","c_salt"].to_numpy())
    np.testing.assert_allclose(test_concentration, ref_concentration, rtol=rtol, atol=atol)
    print("*** Test was successful ***")

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

script_path=pmb.get_resource("samples/salt_solution_gcmc.py")
salt_concentrations=[0.0001, 0.001, 0.01, 0.1]

rtol=0.05 # relative tolerance
atol=0.0 # absolute tolerance

# Ideal test
gcmc_test(script_path, "ideal")

# Interacting test
gcmc_test(script_path, "interacting")
