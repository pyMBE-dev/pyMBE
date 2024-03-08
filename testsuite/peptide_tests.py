# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import os
import tempfile
import subprocess
import numpy as np
import pandas as pd

# Template of the test

def run_peptide_test(script_path,test_pH_values,sequence,rtol,atol,mode="test"):
    """
    Runs a set of tests for a given peptide sequence.

    Args:
        script_path(`str`): Path to the script to run the test.
        test_pH_values(`lst`): List of pH values to be tested.
        sequence(`str`): Amino acid sequence of the peptide.
    """
    valid_modes=["test","save"]
    assert mode in valid_modes, f"Mode {mode} not supported, valid modes: {valid_modes}"

    print(f"Running tests for {sequence}")
    with tempfile.TemporaryDirectory() as time_series_path:
        for pH in test_pH_values:
            print(f"pH = {pH}")
            run_command=["python3", script_path, "--sequence", sequence, "--pH", str(pH), "--mode", "test", "--no_verbose", "--output", time_series_path]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path)
        data_path=pmb.get_resource(path="testsuite/peptide_tests_data")
    if mode == "test":
        # Get reference test data
        ref_data=pd.read_csv(data_path+f"/{sequence}.csv", header=[0, 1])
        # Check charge
        test_charge=np.sort(data["mean","charge"].to_numpy())
        ref_charge=np.sort(ref_data["mean","charge"].to_numpy())
        np.testing.assert_allclose(test_charge, ref_charge, rtol=rtol, atol=atol)
        # Check rg
        test_rg=np.sort(data["mean","rg"].to_numpy())
        ref_rg=np.sort(ref_data["mean","rg"].to_numpy())
        np.testing.assert_allclose(test_rg, ref_rg, rtol=rtol, atol=atol)
        print(f"Test for {sequence} was successful")
    elif mode == "save":
        # Save data for future testing
        data.to_csv(f"{data_path}/{sequence}.csv", index=False)
    else:
        raise RuntimeError

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

script_path=pmb.get_resource(f"samples/Beyer2024/peptide.py")
test_pH_values=[3,7,11]
rtol=0.1 # relative tolerance
atol=0.5 # absolute tolerance

# Run test for K_5-D_5 case
sequence="K"*5+"D"*5
run_peptide_test(script_path=script_path,
                    test_pH_values=test_pH_values,
                    sequence=sequence,
                    rtol=rtol,
                    atol=atol)

# Run test for E_5-H_5 case
sequence="E"*5+"H"*5
run_peptide_test(script_path=script_path,
                    test_pH_values=test_pH_values,
                    sequence=sequence,
                    rtol=rtol,
                    atol=atol)

# Run test for histatin-5 case
sequence="nDSHAKRHHGYKRKFHEKHHSHRGYc"
run_peptide_test(script_path=script_path,
                    test_pH_values=test_pH_values,
                    sequence=sequence,
                    rtol=rtol,
                    atol=atol)   
