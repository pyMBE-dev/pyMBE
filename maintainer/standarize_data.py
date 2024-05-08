# Import pyMBE and other libraries
import pyMBE
import numpy as np
import pandas as pd
import argparse 

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

# Expected inputs
supported_filenames=["data_landsgesell.csv",
                     "Glu-HisMSDE.csv",
                     "Lys-AspMSDE.csv",
                     "histatin5_SoftMatter.txt",
                     "1beb-10mM-torres.dat",
                     "1f6s-10mM-torres.dat"]


parser = argparse.ArgumentParser(description='Script to standarize the data from various authors')
parser.add_argument('--src_filename', 
                    type=str, 
                    required= True,  
                    help='File name from testsuite/data/src. Currently supported {supported_filenames}')
args = parser.parse_args()

# Inputs
filename=args.src_filename

# Outputs
output_filenames={"data_landsgesell.csv": "Landsgesell2020a.csv",
                  "Lys-AspMSDE.csv": "Lunkad2021a.csv",
                  "Glu-HisMSDE.csv": "Lunkad2021b.csv",
                  "histatin5_SoftMatter.txt": "Blanco2020a.csv",
                  "1beb-10mM-torres.dat": "Torres2017.csv",
                  "1f6s-10mM-torres.dat": "Torres2022.csv"}

# Sanity checks
if filename not in supported_filenames:
    ValueError(f"Filename {filename} not supported, supported files are {supported_filenames}")

# Extract the data from Ref.
ref_path=pmb.get_resource(f"testsuite/data/src/{filename}")
Refs_lunkad=["Glu-HisMSDE.csv","Lys-AspMSDE.csv"]
Ref_blanco=["histatin5_SoftMatter.txt"]

Ref_torres = ["1f6s-10mM-torres.dat","1beb-10mM-torres.dat" ]
Ref_landsgesell=["data_landsgesell.csv"]

if filename in Refs_lunkad:
    data=pd.read_csv(ref_path)
    Z_ref = 5*-1*data['aaa']+5*data['aab']
    # Error propagation calculation
    # 1/4 factor added to correct for bug in the original calculation of the error reported by the authors       
    Z_ref_err = 5/4*np.sqrt((data['eaa'])**2+(data['eab'])**2)   
    pH_range = np.linspace(2, 12, num=21)
    
elif filename in Ref_blanco:
    data=np.loadtxt(ref_path, delimiter=",")
    Z_ref=data[:,1]         
    Z_ref_err=data[:,2]
    pH_range = np.linspace(2, 12, num=21)
    
elif filename in Ref_torres:

    Z_ref = []
    Z_ref_err = []
    pH_range = []
    
    with open (ref_path,'r') as file:
        for line in file: 
            line_split = line.split ()
            pH = float (line_split[0])
            pH_range.append(pH)
            Z_ref.append (float(line_split[1]))
            Z_ref_err.append(float(line_split[2]))

elif filename in Ref_landsgesell:
    data = pd.read_csv(ref_path, sep="\t", index_col=False)
else:
    raise RuntimeError()


if filename in Refs_lunkad+Ref_blanco:
    pH_range = np.linspace(2, 12, num=21)

    # Store the data
    data=pd.DataFrame({"pH": pH_range,
                      "charge": Z_ref,
                      "charge_error": Z_ref_err})
if filename in Ref_torres:
    # Store the data
    data=pd.DataFrame({"pH": pH_range,
                    "charge": Z_ref,
                    "charge_error": Z_ref_err})

data_path=pmb.get_resource("testsuite/data")
data.to_csv(f"{data_path}/{output_filenames[filename]}", 
            index=False)
