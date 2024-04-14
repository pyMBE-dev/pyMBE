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

print(f"*** Running G-RxMC tests (ideal) ***\n")


print(f"*** Running test for standard G-RxMC ***\n")

run_command = ["python3", script_path, "--mode", "standard-test"]
subprocess.check_output(run_command)

print(f"*** Test passed ***\n")


print(f"*** Running test for unified G-RxMC ***\n")

run_command = ["python3", script_path, "--mode", "unified-test"]
subprocess.check_output(run_command)

print(f"*** Test passed ***\n")
