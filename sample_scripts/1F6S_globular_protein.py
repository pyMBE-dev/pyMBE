import os
import sys
import inspect
from tqdm import tqdm
import espressomd

from espressomd import interactions
from espressomd.io.writer import vtf
from espressomd import electrostatics 

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

from handy_scripts.handy_functions import  calculate_net_charge


# Here you can adjust the width of the panda columns displayed when running the code 
pmb.pd.options.display.max_colwidth = 10

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

#System Parameters 
pH_value = 2.0
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 

SEED = 77 
solvent_permitivity = 78.3
dt = 0.001

Samples_per_pH = 1000
MD_steps_per_sample = 1000
steps_eq = int(Samples_per_pH/3)
N_samples_print = 10  # Write the trajectory every 100 samples
probability_reaction = 0.5 

bead_size = 0.4*pmb.units.nm
epsilon = 1*pmb.units('reduced_energy')

espresso_system = espressomd.System(box_l=[Box_L.to('reduced_length').magnitude] * 3)

#Directory of the protein model 

protein_name = '1f6s'
protein_filename = 'sample_scripts/coarse_grain_model_of_1f6s.vtf'

#Reads the VTF file of the protein model
protein_positions = pmb.load_protein_vtf_in_df (name=protein_name,filename=protein_filename)


# Solution 

cation_name = 'Na'
anion_name = 'Cl'
pmb.define_particle(name=cation_name,  q=1, diameter=0.2*pmb.units.nm, epsilon=epsilon)
pmb.define_particle(name=anion_name,  q=-1, diameter=0.36*pmb.units.nm,  epsilon=epsilon)

#We define each aminoacid in the pyMBE data frame

pmb.load_pka_set (filename='reference_parameters/pka_sets/CRC1991.txt')

acidic_aminoacids = ['c','E','D','Y','C']
basic_aminoacids  = ['R','n','K','H']

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

pmb.center_pmb_object_in_the_simulation_box (name=protein_name,espresso_system=espresso_system)

pmb.create_counterions_in_espresso (pmb_object='particle',cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system)

c_salt_calculated = pmb.create_added_salt_in_espresso (espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)


basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.drop_duplicates().to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.drop_duplicates().to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print('The box length of the system is', Box_L.to('reduced_length'), Box_L.to('nm'))
print('The ionisable groups in the protein are ', list_ionisible_groups)

RE, sucessfull_reactions_labels = pmb.setup_constantpH_reactions_in_espresso (counter_ion=cation_name, constant_pH=pH_value, SEED = SEED )
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy
pmb.setup_lj_interactions_in_espresso (espresso_system=espresso_system)

print('\nMinimazing system energy\n')
espresso_system.cell_system.skin = 0.4
espresso_system.time_step = dt 
print('steepest descent')
espresso_system.integrator.set_steepest_descent(f_max=0, gamma=0.1, max_displacement=0.1)
espresso_system.integrator.run(1000)
print('velocity verlet')
espresso_system.integrator.set_vv()  # to switch back to velocity Verlet
espresso_system.integrator.run(1000)
espresso_system.thermostat.turn_off()
print('\nMinimization finished \n')

#Electrostatic energy setup 
print ('Electrostatics setup')
bjerrum_length = pmb.e.to('reduced_charge')**2 / (4 *pmb.units.pi*pmb.units.eps0* solvent_permitivity * pmb.kT.to('reduced_energy'))
coulomb_prefactor = bjerrum_length.to('reduced_length') * pmb.kT.to('reduced_energy')
coulomb = espressomd.electrostatics.P3M ( prefactor = coulomb_prefactor.magnitude, accuracy=1e-3)

print('\nBjerrum length ', bjerrum_length.to('nm'), '=', bjerrum_length.to('reduced_length'))

espresso_system.time_step = dt
espresso_system.actors.add(coulomb)

# save the optimal parameters and add them by hand
p3m_params = coulomb.get_params()
espresso_system.actors.remove(coulomb)

coulomb = espressomd.electrostatics.P3M (
                            prefactor = coulomb_prefactor.magnitude,
                            accuracy = 1e-3,
                            mesh = p3m_params['mesh'],
                            alpha = p3m_params['alpha'] ,
                            cao = p3m_params['cao'],
                            r_cut = p3m_params['r_cut'],
                            tune = False,
                                )

espresso_system.actors.add(coulomb)

print("\nElectrostatics successfully added to the system \n")

#Save the initial state 

n_frame = 0
pmb.write_output_vtf_file(espresso_system=espresso_system,n_frame=n_frame)

print (f'Optimizing skin\n')
espresso_system.time_step = dt 
espresso_system.integrator.set_vv()
espresso_system.thermostat.set_langevin(kT=pmb.kT.to('reduced_energy').magnitude, gamma=0.1, seed=SEED)

espresso_system.cell_system.tune_skin ( min_skin = 10, 
                                        max_skin = espresso_system.box_l[0]/2., tol=1e-3, 
                                        int_steps=1000, adjust_max_skin=True)

print('Optimized skin value: ', espresso_system.cell_system.skin, '\n')


observables_df = pmb.pd.DataFrame()
time_step = []
net_charge_list = []

net_charge_amino_save = {}

for step in tqdm(range(Samples_per_pH+steps_eq)):
        
        if pmb.np.random.random() > probability_reaction:
            espresso_system.integrator.run(steps=MD_steps_per_sample)
        else:
            RE.reaction(steps=total_ionisible_groups)

        calculated_net_charge = calculate_net_charge (espresso_system=espresso_system,pmb_df = pmb.df, name =protein_name)

        net_charge = calculated_net_charge['net_charge']
        net_charge_aminoacids = calculated_net_charge ['net_charge_aminoacids']

        time_step.append (str(espresso_system.time))
        net_charge_list.append (net_charge)

        if len(net_charge_amino_save.keys()) == 0:
            for amino in net_charge_aminoacids.keys():
                net_charge_amino_save [amino] = []
        for amino in net_charge_aminoacids.keys():            
            net_charge_amino_save [amino].append (net_charge_aminoacids[amino])

        if (step % N_samples_print == 0 ):
            n_frame +=1
            pmb.write_output_vtf_file(espresso_system=espresso_system,n_frame=n_frame)

observables_df['time'] = time_step 
observables_df['Znet'] = net_charge_list

for amino in net_charge_amino_save.keys():
    observables_df[amino] = net_charge_amino_save[amino]

observables_df.to_csv(f'pH-{pH_value}_observables_.csv',index=False)


print(pmb.df)
print (observables_df)


from espressomd import visualization

espresso_system.time_step=1e-3
espresso_system.thermostat.set_langevin(kT=1, gamma=1.0, seed=24)
espresso_system.cell_system.tune_skin(min_skin=0.01, max_skin =4.0, tol=0.1, int_steps=1000)
visualizer = espressomd.visualization.openGLLive(espresso_system, bond_type_radius=[0.3], background_color=[1, 1, 1],particle_type_colors=[[1.02,0.51,0], # Brown
                [1,1,1],  # Grey
                [2.55,0,0], # Red
                [0,0,2.05],  # Blue
                [0,0,2.05],  # Blue
                [2.55,0,0], # Red
                [2.05,1.02,0]])     
visualizer.run(1)

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
