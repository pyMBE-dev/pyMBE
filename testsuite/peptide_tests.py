# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import os 
import numpy as np
import pandas as pd
from tqdm import tqdm

# Template of the test

def run_peptide_test(script_path,test_pH_values,sequence,rtol,atol):
    """
    Runs a set of tests for a given peptide sequence.

    Args:
        script_path(`str`): Path to the script to run the test.
        test_pH_values(`lst`): List of pH values to be tested.
        sequence(`str`): Amino acid sequence of the peptide.
    """
    # Get data folder
    time_series_folder_path=pmb.get_resource(f"samples/Beyer2024/time_series")
    # clean up data folder
    os.system(f"rm {time_series_folder_path}/*")
    print(f"Running tests for {sequence}")
    for pH in (pbar := tqdm(test_pH_values)):
        pbar.set_description(f"pH = {pH}")
        run_command=f"python3 {script_path} --sequence {sequence} --pH {pH} --mode test"
        print(run_command)
        os.system(run_command)
    # Analyze all time series
    data=analysis.analyze_time_series(path_to_datafolder=time_series_folder_path)
    # Get reference test data
    data_path=pmb.get_resource(path="testsuite/peptide_tests_data")
    ref_data=pd.read_csv(data_path+f"/{sequence}.csv", header=[0, 1], index_col=0)
    # Check that the charge and radius of gyration are consistent
    np.testing.assert_allclose(data["mean","charge"].to_numpy(), 
                               ref_data["mean","charge"].to_numpy(), 
                               rtol=rtol, 
                               atol=atol)
    np.testing.assert_allclose(data["mean","rg"].to_numpy(), 
                               ref_data["mean","rg"].to_numpy(), 
                               rtol=rtol, 
                               atol=atol)
    print(f"Test for {sequence} succesful")

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