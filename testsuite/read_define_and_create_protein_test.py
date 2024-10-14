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

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)


print("*** Unit test: check that read_protein_vtf_in_df() loads the protein topology correctly ***")

protein_pdb = '1beb'

path_to_cg=pmb.get_resource(f'parameters/globular_proteins/{protein_pdb}.vtf')

topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)


# np.testing.assert_equal(actual= topology_dict, 
#                         desired= {},
#                         verbose = True)


print("*** Unit test passed ***")


print("*** Unit test: check that define_protein() defines the aminoacids in the protein correctly ***")

pmb.define_protein (name=protein_pdb, 
                    topology_dict=topology_dict, 
                    model = '2beadAA',
                    lj_setup_mode = "wca")


for aminoacid in topology_dict.keys():
    input_parameters=topology_dict[aminoacid]
    for index in pmb.df[pmb.df['name']==aminoacid].index:
          
            np.testing.assert_equal(actual=str(protein_pdb), 
                                desired="protein", 
                                verbose=True)
          
            np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="residue", 
                                verbose=True)
          
            np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                                desired="particle", 
                                verbose=True)

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

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

np.testing.assert_equal(actual = len(particle_id_list), 
                        desired = len(topology_dict.keys()), 
                        verbose = True)


#Check that creates proteins correctly creates all the particles in the topology_dict()

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



