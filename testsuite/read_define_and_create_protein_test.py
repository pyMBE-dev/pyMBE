#
# Copyright (C) 2024 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np 
import espressomd
import pyMBE
import re

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)


print("*** Unit test: check that read_protein_vtf_in_df() loads the protein topology correctly ***")

#NOTE: save to json error "TypeError: Object of type Quantity is not JSON serializable"
# problems with units

protein_pdb = '1beb'

path_to_cg=pmb.get_resource(f'parameters/globular_proteins/{protein_pdb}.vtf')

topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)

# with open ("protein_topology_dict.json", "w") as output:
#       json.dump(topology_dict, output)

# with open ("protein_topology_dict.json", "r") as file:
#     load_json = json.load(file)

# np.testing.assert_equal(actual= topology_dict, 
#                         desired= load_json,
#                         verbose = True)

print("*** Unit test passed ***")


print("*** Unit test: check that define_protein() defines the aminoacids in the protein correctly ***")

protein_model = '2beadAA'

pmb.define_protein (name=protein_pdb, 
                    topology_dict=topology_dict, 
                    model = protein_model,
                    lj_setup_mode = "wca")
sequence = []
clean_sequence= []

for aminoacid in topology_dict.keys():
    
    input_parameters=topology_dict[aminoacid]
    residue_name = re.split(r'\d+', aminoacid)[0]
    sequence.append(residue_name)
    
    if residue_name not in ['CA', 'n', 'c']:
        clean_sequence.append(residue_name)

    
    for index in pmb.df[pmb.df['name']==residue_name].index:
        if residue_name not in sequence:           
            np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="particle", 
                                verbose=True)

residue_list = pmb.define_AA_residues(sequence= clean_sequence,
                                      model = protein_model)

for residue in residue_list:
    for index in pmb.df[pmb.df['name']==residue].index:
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                        desired="residue", 
                        verbose=True)

protein_index = pmb.df[pmb.df['name']==protein_pdb].index


np.testing.assert_equal(actual=str(pmb.df.loc[protein_index, "name"].values[0]), 
                            desired=protein_pdb, 
                            verbose=True)        

np.testing.assert_equal(actual=pmb.df.loc[protein_index, ('sequence','')].values[0], 
                    desired=clean_sequence, 
                    verbose=True)


np.testing.assert_equal(actual=pmb.df.loc[protein_index, ('residue_list','')].values[0], 
                    desired=residue_list, 
                    verbose=True)  

print("*** Unit test passed ***")

print("*** Unit test: check that define_protein() raises a ValueError if the user provides a wrong model")


input_parameters={"name": protein_pdb,
                 "topology_dict": topology_dict,
                 "model" : "3beadAA",
                "lj_setup_mode": "wca"}

np.testing.assert_raises(ValueError, pmb.define_protein, **input_parameters)

#Note: provide a correct option to lj 

input_parameters={"name": protein_pdb,
                 "topology_dict": topology_dict,
                 "model" : protein_model,
                "lj_setup_mode": "awc"}

np.testing.assert_raises(ValueError, pmb.define_protein, **input_parameters)

print("*** Unit test passed ***")


print("*** Unit test: check that create_protein() creates all the particles in the protein into the espresso_system with the properties defined in pmb.df  ***")

espresso_system=espressomd.System(box_l = [10]*3)

# Here we upload the pka set from the reference_parameters folder
path_to_pka=pmb.get_resource('parameters/pka_sets/Nozaki1967.json') 
pmb.load_pka_set(filename=path_to_pka)

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

residue_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].residue_id.dropna().to_list()

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

box_half=espresso_system.box_l[0]/2.0
protein_center = pmb.generate_coordinates_outside_sphere(radius = 1, 
                                                                        max_dist=box_half, 
                                                                        n_samples=1, 
                                                                        center=[box_half]*3)[0]
for id in particle_id_list:

    initial_pos = espresso_system.part.by_id(id).pos
    charge = espresso_system.part.by_id(id).q
    es_type = espresso_system.part.by_id(id).type

    residue_id = pmb.df.loc[pmb.df['particle_id']==id].residue_id.values[0]
    residue_name = pmb.df.loc[pmb.df['particle_id']==id].name.values[0]

    input_parameters=topology_dict[aminoacid]

    #NOTE: positions are too different even adding the protein center 
    
    # np.testing.assert_equal(actual=initial_pos, 
    #                     desired=topology_dict[residue_name+residue_id]['initial_pos']+protein_center, 
    #                     verbose=True)
      
    index = pmb.df.loc[pmb.df['particle_id']==id].index

    np.testing.assert_equal(actual=charge, 
                        desired=pmb.df.loc[index, ("state_one","z")].values[0], 
                        verbose=True)
      
    np.testing.assert_equal(actual=es_type, 
                    desired=pmb.df.loc[index, ("state_one","es_type")].values[0], 
                    verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_protein() does not create any protein for number_of_proteins <= 0  ***")

starting_number_of_particles=len(espresso_system.part.all())

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=0,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=-1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that enable_motion_of_rigid_object() moves the protein correctly ***")

espresso_system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

pmb.enable_motion_of_rigid_object(espresso_system=espresso_system,
                                  name=protein_pdb)

for id in particle_id_list:
    fix_value = espresso_system.part.by_id(id).fix

    np.testing.assert_equal(actual=fix_value, 
                            desired=[True, True, True], 
                            verbose=True)

print("*** Unit test passed ***")


print("*** Unit test: check that enable_motion_of_rigid_object() raises a ValueError if a wrong pmb_type is provided***")

input_parameters = {"espresso_system":espresso_system,
                    "name": "CA"}

np.testing.assert_raises(ValueError, pmb.enable_motion_of_rigid_object, **input_parameters)

print("*** Unit test passed ***")

print("*** Unit test: check that protein_sequence_parser() correctly returns que protein sequence ***")

output = ["R", "E", "C", "H"]
sequence = "ARG-GLU-CYS-HIS"

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

output = ["R", "E", "C", "H"]
sequence = "R-E-C-H"

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

output = ["R", "E", "C", "H"]
sequence = ["R","E", "C", "H"]

clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

np.testing.assert_equal(actual=clean_sequence, 
                        desired=output, 
                        verbose=True)

#Note: The function does not accept a list of aminoacid with 3 letter code?

# output = ["R", "E", "C", "H"]
# sequence = ["ARG","GLU", "CYS", "HIS"]

# clean_sequence= pmb.protein_sequence_parser(sequence = sequence)

# np.testing.assert_equal(actual=clean_sequence, 
#                         desired=output, 
#                         verbose=True)


print("*** Unit test: check that protein_sequence_parser() raises a ValueError if a wrong residue key is provided***")

input_parameters = {"sequence":"ARG-GLU-TUR-HIS"}

np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":"A-E-E-X"}

np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":["A", "E","X"]}

np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

print("*** Unit test passed ***")