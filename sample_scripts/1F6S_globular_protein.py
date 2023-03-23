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

# Here you can adjust the width of the panda columns displayed when running the code 
pmb.pd.options.display.max_colwidth = 10

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

protein_name = '1f6s'
protein_filename = 'sample_scripts/coarse_grain_model_of_1f6s.vtf'

#Reads the VTF file of the protein model
protein_positions = pmb.load_protein_vtf_in_df (name=protein_name,filename=protein_filename)

#We define each aminoacid in the pyMBE data frame

pmb.load_pka_set (filename='reference_parameters/pka_sets/CRC1991.txt')

acidic_aminoacids = ['c','E','D','Y','C']
basic_aminoacids  = ['R','n','K','H']

bead_size = 0.4*pmb.units.nm
epsilon = 1*pmb.units('reduced_energy')

already_defined_AA=[]

protein_sequence = pmb.df.loc[pmb.df['name']== protein_name].sequence.values[0]

for aminoacid_key in pmb.protein_sequence_parser(sequence=protein_sequence):
    if aminoacid_key in already_defined_AA:
        continue
    if aminoacid_key in acidic_aminoacids:
        pmb.define_particle (name=aminoacid_key, acidity='acidic', diameter=bead_size, epsilon=epsilon)
    elif aminoacid_key in basic_aminoacids:
        pmb.define_particle (name=aminoacid_key, acidity='basic', diameter=bead_size, epsilon=epsilon)
    else:
        pmb.define_particle (name=aminoacid_key, q=0, diameter=bead_size, epsilon=epsilon)
    already_defined_AA.append(aminoacid_key)

pmb.define_particle(name='CA',q=0,diameter=bead_size,epsilon=epsilon)
pmb.define_particle(name='Ca',q=0,diameter=bead_size,epsilon=epsilon)



pmb.create_protein_in_espresso(name=protein_name,
                               number_of_proteins=1,
                               espresso_system=espresso_system,
                               positions=protein_positions)

print(pmb.df)



# from espressomd import visualization

# filename = 'protein.png'

# visualizer = visualization.openGLLive(
#     espresso_system, bond_type_radius=[0.3], particle_coloring='type', draw_axis=False, background_color=[1, 1, 1],
# particle_type_colors=[[1.02,0.51,0], # Brown
#                 [1,1,1],  # Grey
#                 [2.55,0,0], # Red
#                 [0,0,2.05],  # Blue
#                 [0,0,2.05],  # Blue
#                 [2.55,0,0], # Red
#                 [2.05,1.02,0]]) # Orange
# visualizer.screenshot(filename)

# from PIL import Image
# img = Image.open(filename)
# img.show()


