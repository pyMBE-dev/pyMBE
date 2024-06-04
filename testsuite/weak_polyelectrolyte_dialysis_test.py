# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import sys
import tempfile
import subprocess
import numpy as np
import pandas as pd

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

script_path=pmb.get_resource("samples/Beyer2024/weak_polyelectrolyte_dialysis.py")
test_pH_values=[3,5,7,9]
c_salt_res=0.01
c_mon_sys=0.435
pKa_value=4.0
rtol=0.1 # relative tolerance
atol=0.05 # absolute tolerance

print("*** Running test for weak polyelectrolyte dialysis with G-RxMC (interacting). ***")
with tempfile.TemporaryDirectory() as time_series_path:
    for pH in test_pH_values:
        print(f"pH = {pH}")
        run_command=[sys.executable, script_path, "--c_salt_res", str(c_salt_res), "--c_mon_sys", str(c_mon_sys), "--pKa_value", str(pKa_value), "--pH_res", str(pH), "--mode", "test", "--output", time_series_path, "--no_verbose"]
        print(subprocess.list2cmdline(run_command))
        subprocess.check_output(run_command)
    # Analyze all time series
    data=analysis.analyze_time_series(path_to_datafolder=time_series_path,
                                      filename_extension="_time_series.csv")
    data_path=pmb.get_resource(path="testsuite/weak_polyelectrolyte_dialysis_test_data")

# Get reference test data
ref_data=pd.read_csv(f"{data_path}/data.csv", header=[0, 1])
# Check charge
test_charge=np.sort(data["mean","alpha"].to_numpy())
ref_charge=np.sort(ref_data["mean","alpha"].to_numpy())
np.testing.assert_allclose(test_charge, ref_charge, rtol=rtol, atol=atol)
print("*** Test was successful ***")
