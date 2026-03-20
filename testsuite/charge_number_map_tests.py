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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import unittest as ut

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)
pmb.define_particle(name="I", 
                    pka = np.nan,
                    sigma = 1* pmb.units.nm,
                    epsilon= 1 *pmb.units.reduced_energy,
                    z =5)

pmb.define_particle(name = "A", 
                    acidity = "acidic",
                    sigma = 1* pmb.units.nm,
                    epsilon = 1 *pmb.units.reduced_energy,
                    pka =4)

pmb.define_particle(name ="B", 
                    acidity = "basic",
                    sigma = 1* pmb.units.nm,
                    epsilon = 1 *pmb.units.reduced_energy,
                    pka =   4)

charge_map = pmb.get_charge_number_map()
type_map = pmb.get_type_map()

class Test(ut.TestCase):
    def test_inert_particle(self):
        """
        Check that get_charge_number_map works correctly for inert particles
        """
        self.assertEqual(charge_map[type_map["I"]],
                        5)

    def test_acidic_particle(self):
        """
        Check that get_charge_number_map works correctly for acidic particles
        """
        self.assertEqual(charge_map[type_map["AH"]],
                         0)
        self.assertEqual(charge_map[type_map["A"]],
                         -1)

    def test_basic_particle(self):
        """
        Check that get_charge_number_map works correctly for basic particles
        """
        self.assertEqual(charge_map[type_map["BH"]],
                         1)
        self.assertEqual(charge_map[type_map["B"]],
                         0)
        


if __name__ == '__main__':
    ut.main()