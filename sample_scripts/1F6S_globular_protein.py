import os
import sys
import inspect
import espressomd

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')


#System Parameters 

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 

espresso_system = espressomd.System(box_l=[Box_L.to('reduced_length').magnitude] * 3)


#Directory of the protein model 
protein_filename = 'sample_scripts/coarse_grain_model_of_1f6s.vtf'
#Reads the VTF file of the protein model
pdb_code = '1f6s'

protein_positions = pmb.load_protein_vtf_in_df (name=pdb_code,filename=protein_filename)

pmb.create_protein_in_espresso(name=pdb_code,
                               number_of_proteins=1,
                               espresso_system=espresso_system,
                               positions=protein_positions)
print(pmb.df.head(100))