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
import pyMBE
import pandas as pd
import numpy as np
import unittest as ut

class Test(ut.TestCase):
    def test_pka_set_format(self):
        """
        Check that the different pKa sets are correctly formatted
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pka_root=pmb.root / "parameters" / "pka_sets" 
        for path in pka_root.glob("*.json"):
            pmb.load_pka_set(path)
            pmb.db.delete_reactions()

    def test_sanity_load_datasets(self):
        """
        Check that the order to execute load_pka_set() and load_databaasedoes not change the resulting parameters in pyMBE database
        """

        # First order of loading parameters
        pmb1 = pyMBE.pymbe_library(seed=42)
        path_to_interactions=pmb1.root / "parameters" / "peptides" / "Lunkad2021"
        path_to_pka=pmb1.root / "parameters" / "pka_sets" / "Hass2015.json"
        pmb1.load_database (folder=path_to_interactions)
        pmb1.load_pka_set(filename=path_to_pka)
        

        # Second order of loading parameters
        pmb2 = pyMBE.pymbe_library(seed=23)
        path_to_interactions=pmb2.root / "parameters" / "peptides" / "Lunkad2021"
        path_to_pka=pmb2.root / "parameters" / "pka_sets" / "Hass2015.json"
        pmb2.load_pka_set(filename=path_to_pka)
        pmb2.load_database(folder=path_to_interactions)

        pmb_types_to_test = ["particle_state",
                             "particle",
                             "bond"]
        for pmb_type in pmb_types_to_test:
            pd.testing.assert_frame_equal(pmb1.get_templates_df(pmb_type=pmb_type),
                                          pmb2.get_templates_df(pmb_type=pmb_type))
        pd.testing.assert_frame_equal(pmb1.get_reactions_df(),
                                      pmb2.get_reactions_df())
        
    def test_sanity_check_pka_set(self):
        """
        Check that  check_pka_set raises a ValueError if data is missing important fields
        """
        pmb = pyMBE.pymbe_library(seed=42)
        np.testing.assert_raises(ValueError, pmb.check_pka_set, {"name" : {}})
        np.testing.assert_raises(ValueError, pmb.check_pka_set, {"name" : {"pka_value": 1.}})
        np.testing.assert_raises(ValueError, pmb.check_pka_set, {"name" : {"acidity": 1.}})

if __name__ == "__main__":
    ut.main()