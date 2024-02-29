# Load espresso, sugar and other necessary libraries
import sys
import os 
import inspect
import espressomd
import numpy as np
import pandas as pd
import argparse 
from tqdm import tqdm
from espressomd.io.writer import vtf
from espressomd import interactions
from espressomd import electrostatics

# Import pyMBE
import pyMBE
from lib import analysis
from lib import handy_functions as hf

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

parser = argparse.ArgumentParser(description='Script to create the data from Beyer2024')
parser.add_argument('--fig_label', 
                    type=str, 
                    required= True,  
                    help='Label of the corresponding figure in Beyer2024, currently supported: 6a, 6b, 6c.')
parser.add_argument('--mode', 
                    type=str, 
                    default= "long-run",  
                    help='sets for how long the simulation runs, valid modes are "short-run" and "long-run"')
args = parser.parse_args()

# Inputs
fig_label=args.fig_label
mode=args.mode

# Sanity checks
valid_fig_labels=["6a", "6b", "6c"]

if fig_label not in valid_fig_labels:
    raise ValueError(f"The figure label {fig_label} is not supported. Supported figure labels are {valid_fig_labels}")

valid_modes=["short-run","long-run"]
if mode not in valid_modes:
    raise ValueError(f"Mode {mode} is not currently supported, valid modes are {valid_modes}")

## Peptide plots (Fig. 6)
if fig_label in ["6a", "6b", "6c"]:
    script_path=pmb.get_resource(f"samples/Beyer2024/peptide.py")
    if fig_label == "6a":
        sequence="K"*5+"D"*5
    elif fig_label == "6b":
        sequence="E"*5+"H"*5
    elif fig_label == "6c":
        sequence="nDSHAKRHHGYKRKFHHSHRGYc"
    else:
        raise RuntimeError()
    pH_range = np.linspace(2, 12, num=21)

    for pH in pH_range:
        run_command=f"python3 {script_path} --sequence {sequence} --pH {pH} --mode {mode}"
        print(run_command)
        #os.system(run_command)

# Read all files in the subdir
data_files=[]
time_series_folder_path=pmb.get_resource(f"samples/Beyer2024/time_series")

data=pd.DataFrame()

with os.scandir(time_series_folder_path) as subdirectory:
    # Gather all data
    for subitem in subdirectory:
        if subitem.is_file():
            if 'time_series' in subitem.name:
                # Get parameters from the file name
                data_dict=analysis.get_params_from_dir_name(subitem.name.replace('_time_series.csv', ''))
                file_data=pd.DataFrame(data_dict, index=[0])
                # Get the observables for binning analysis
                time_series_data=analysis.read_csv_file(path=f"{time_series_folder_path}/{subitem.name}")
                analyzed_data=analysis.block_analyze(full_data=time_series_data)
                file_data=pd.concat([file_data,analyzed_data])
                data = pd.concat([data,file_data])

data_path=pmb.get_resource("samples/Beyer2024/")+"data"
if not os.path.exists(data_path):
    os.makedirs(data_path)

data.to_csv(f"{data_path}/fig{fig_label}.csv")







