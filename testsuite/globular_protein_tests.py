# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import tempfile
import subprocess
import numpy as np
import pandas as pd

# Template of the test

def run_protein_test(script_path, test_pH_values, protein_pdb, rtol, atol,mode="test"):
    """
    Runs a set of tests for a given protein pdb.

    Args:
        script_path(`str`): Path to the script to run the test.
        test_pH_values(`lst`): List of pH values to be tested.
        protein_pdb(`str`): PDB code of the protein.
    """
    valid_modes=["test","save"]
    assert mode in valid_modes, f"Mode {mode} not supported, valid modes: {valid_modes}"

    print(f"Running tests for {protein_pdb}")
    with tempfile.TemporaryDirectory() as time_series_path:
       
        for pH in test_pH_values:
            print(f"pH = {pH}")
            
            if protein_pdb == '1f6s':
                  run_command=["python3", script_path, "--pdb", protein_pdb, "--pH", str(pH), "--path_to_cg", f"parameters/globular_proteins/{protein_pdb}.vtf", "--metal_ion_name","Ca", "--metal_ion_charge", "2","--mode", "test", "--no_verbose", "--output", time_series_path]
                  
            else:
                  run_command=["python3", script_path, "--pdb", protein_pdb, "--pH", str(pH), "--path_to_cg", f"parameters/globular_proteins/{protein_pdb}.vtf",  "--mode", "test", "--no_verbose", "--output", time_series_path]
                  
                  
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path)

        data_path=pmb.get_resource(path="testsuite/globular_protein_tests_data")
        
    if mode == "test":
        # Get reference test data
        ref_data=pd.read_csv(data_path+f"/{protein_pdb}.csv", header=[0, 1])
        # Check charge
        test_charge=np.sort(data["mean","charge"].to_numpy())
        ref_charge=np.sort(ref_data["mean","charge"].to_numpy())
        np.testing.assert_allclose(test_charge, ref_charge, rtol=rtol, atol=atol)
        print(f"Test for {protein_pdb} was successful")
        
    elif mode == "save":
        # Save data for future testing
        data.to_csv(f"{data_path}/{protein_pdb}.csv", index=False)
    else:
        raise RuntimeError

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

script_path=pmb.get_resource(f"samples/Beyer2024/globular_protein.py")
test_pH_values=[2,5,7]
rtol=0.1 # relative tolerance
atol=0.5 # absolute tolerance

# Run test for 1BEB case
protein_pdb = "1beb"
run_protein_test(script_path=script_path,
                    test_pH_values=test_pH_values,
                    protein_pdb=protein_pdb,
                    rtol=rtol,
                    atol=atol)

# Run test for 1F6S case
protein_pdb = "1f6s"
run_protein_test(script_path=script_path,
                    test_pH_values=test_pH_values,
                    protein_pdb=protein_pdb,
                    rtol=rtol,
                    atol=atol)   
