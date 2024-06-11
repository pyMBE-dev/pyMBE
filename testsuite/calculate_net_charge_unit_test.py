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
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

print("*** Unit test: check that calculate_net_charge calculates the charge in a molecule properly ***")

pmb.define_particle(name='0P',
                        q=0)

pmb.define_particle(name='+1p',
                    q=+1)

pmb.define_particle(name='-1p',
                    q=-1)

pmb.define_residue(
    name = 'R1',
    central_bead = '+1p',
    side_chains = ['0P']
    )

pmb.define_residue(
    name = 'R2',
    central_bead = '-1p',
    side_chains = ['R1']
    )

bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant,
                }

pmb.define_default_bond(bond_type = bond_type, 
                        bond_parameters = harmonic_bond)

molecule_name = 'generic_molecule'
pmb.define_molecule(name=molecule_name, 
                    residue_list = ['R1']*2+['R2']*3)

# Create an instance of an espresso system
espresso_system=espressomd.System(box_l = [10]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
molecules = pmb.create_molecule(name=molecule_name,
                        number_of_molecules= 2,
                        espresso_system=espresso_system,
                        use_default_bond=True,)

charge_map=pmb.calculate_net_charge(molecule_name=molecule_name,
                                    espresso_system=espresso_system)

# Check mean charge
np.testing.assert_equal(charge_map["mean"],2.0)
# Check molecule charge map
np.testing.assert_equal(charge_map["molecules"],{0: 2.0, 1: 2.0})
# Check residue charge map
np.testing.assert_equal(charge_map["residues"],{0: 1.0, 1: 1.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 1.0, 6: 1.0, 7: 0.0, 8: 0.0, 9: 0.0})

print("*** Unit test passed ***")
print("*** Unit test: check that calculate_net_charge raises a ValueError if one provides the name of an object that is not a molecule ***")
input_parameters={"molecule_name":"R1", 
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.calculate_net_charge, **input_parameters)
print("*** Unit test passed ***")

