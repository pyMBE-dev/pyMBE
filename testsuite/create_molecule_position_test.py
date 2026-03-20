#
# Copyright (C) 2024-2026 pyMBE-dev team
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

import espressomd
import pyMBE
import unittest as ut

pmb = pyMBE.pymbe_library(seed=42)

pmb.define_particle(name='central_mon',
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))

pmb.define_particle(name='side_mon',
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))

pmb.define_residue(name = 'res1',
                   central_bead = 'central_mon',
                   side_chains = ['side_mon', 'side_mon'])

bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant,
                 }

pmb.define_default_bond(bond_type = bond_type, bond_parameters = harmonic_bond)

molecule_name = 'generic_molecule'
pmb.define_molecule(name=molecule_name, residue_list = ['res1']*5)

# Create an instance of an espresso system
L = 52
espresso_system=espressomd.System(box_l = [L]*3)
pos_list = [[10,10,10], [20,20,20], [30,30,30]]

class Test(ut.TestCase):
    def test_create_molecule_at_position(self):
        """
        Check that the positions of the central bead of the first residue in the generated molecules are equal to the input positions
        """
        molecule_ids = pmb.create_molecule(name=molecule_name,
                                           number_of_molecules= 3,
                                           espresso_system=espresso_system,
                                           use_default_bond=True,
                                           list_of_first_residue_positions = pos_list)
        particle_id_map = pmb.get_particle_id_map(object_name=molecule_name)
        central_bead_pos = []
        for molecule_id in molecule_ids:
            pids = particle_id_map["molecule_map"][molecule_id]
            central_bead_id = min(pids)
            central_bead_pos.append(espresso_system.part.by_id(central_bead_id).pos.tolist())
        self.assertListEqual(pos_list,
                             central_bead_pos)
        for molid in molecule_ids:
            pmb.delete_instances_in_system(instance_id=molid,
                                           pmb_type="molecule",
                                           espresso_system=espresso_system)
        
    def test_sanity_create_molecule(self):
        """
        Sanity tests for input positions in create_molecule
        """

        # Check that create_molecule raises a ValueError if the user does not provide a nested list for list_of_first_residue_positions 
        input_parameters={"name": "generic_molecule",
                        "number_of_molecules": 1,
                        "espresso_system": espresso_system,
                        "list_of_first_residue_positions": [1,2,3]}
        self.assertRaises(ValueError, 
                          pmb.create_molecule, 
                          **input_parameters)
        # Check that create_molecule raises a ValueError if the user does not provide a nested list with three coordinates
        input_parameters={"name": "generic_molecule",
                        "number_of_molecules": 1,
                        "espresso_system": espresso_system,
                        "list_of_first_residue_positions": [[1,2]]}
        self.assertRaises(ValueError, 
                          pmb.create_molecule, 
                          **input_parameters)
        # Check that create_molecule raises a ValueError if the user does not provide a the same number of first_residue_positions as number_of_molecules
        input_parameters={"name": "generic_molecule",
                        "number_of_molecules": 2,
                        "espresso_system": espresso_system,
                        "list_of_first_residue_positions": [[1,2,3]]}
        self.assertRaises(ValueError, 
                          pmb.create_molecule, 
                          **input_parameters)
        
    def test_center_molecule_in_simulation_box(self):
        """
        Unit tests for center_molecule_in_simulation_box
        """
        molecule_ids = pmb.create_molecule(name=molecule_name,
                                           number_of_molecules= 3,
                                           espresso_system=espresso_system,
                                           use_default_bond=True,
                                           list_of_first_residue_positions = pos_list)

        # Check that center_molecule_in_simulation_box works correctly for cubic boxes
        pmb.center_object_in_simulation_box(instance_id=molecule_ids[0], 
                                            espresso_system=espresso_system,
                                            pmb_type="molecule")
        center_of_mass = pmb.calculate_center_of_mass(instance_id=molecule_ids[0],
                                                      pmb_type="molecule", 
                                                      espresso_system=espresso_system)
        center_of_mass_ref = [L/2]*3
        for ind in range(len(center_of_mass)):
            self.assertAlmostEqual(center_of_mass[ind], 
                                center_of_mass_ref[ind])
        #Check that center_molecule_in_simulation_box works correctly for non-cubic boxes
        espresso_system.change_volume_and_rescale_particles(d_new=3*L, dir="z")
        
        pmb.center_object_in_simulation_box(instance_id=molecule_ids[2],
                                            pmb_type="molecule", 
                                            espresso_system=espresso_system)
        center_of_mass = pmb.calculate_center_of_mass(instance_id=molecule_ids[2],
                                                      pmb_type="molecule", 
                                                      espresso_system=espresso_system)
        center_of_mass_ref = [L/2, L/2, 1.5*L]
        for ind in range(len(center_of_mass)):
            self.assertAlmostEqual(center_of_mass[ind], 
                                center_of_mass_ref[ind])

    def test_sanity_center_object_in_simulation_box(self):
        """
        Sanity tests for center_molecule_in_simulation_box
        """
        # Check that center_molecule_in_simulation_box raises a Value Error if a wrong molecule_id is provided

        input_parameters = {"instance_id": 20 ,
                            "pmb_type": "molecule",
                            "espresso_system":espresso_system}

        self.assertRaises(ValueError, 
                          pmb.center_object_in_simulation_box, 
                           **input_parameters)
        


if __name__ == "__main__":
    ut.main()

