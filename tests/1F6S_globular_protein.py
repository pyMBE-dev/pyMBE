import os
import sys
import inspect


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

#Directory of the protein model 
protein_filename = 'tests/coarse_grain_model_of_1f6s.vtf'
#Reads the VTF file of the protein model
pmb.load_protein_vtf_file(filename=protein_filename)

print(pmb.df.head(10))