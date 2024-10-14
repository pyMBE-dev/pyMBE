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


np.testing.assert_equal(actual= str(pmb.df.loc[index, "pmb_type"].values[0]), 
                        desired= 'protein',
                        verbose = True)


# print("*** Unit test: check that define_protein() does not setup any particle if no parameters are provided ***")


# output = pmb.define_protein(name="", model="2beadAA", topology_dict={})
# np.testing.assert_equal(actual=output, 
#                         desired=0, 
#                         verbose=True)



# Simulation parameters



# pmb.read_protein_vtf_in_df

# pmb.create_protein(name= ,
#                    number_of_proteins=,
#                    espresso_system=,
#                    topology_dict=, )