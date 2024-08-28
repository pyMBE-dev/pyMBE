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

import unittest as ut
import numpy as np
import pathlib
import pyMBE


class Test(ut.TestCase):
    data_root = pathlib.Path(__file__).parent.resolve() / "henderson_hasselbalch_tests_data"

    def test(self):
        pmb = pyMBE.pymbe_library(seed=42)
        print("*** Running HH tests ***\n")

        # Peptide parameters
        sequence1 = 5 * "D" + 8 * "H"
        sequence2 = 3 * "E" + 7 * "R"
        pep1_concentration = 1e-2 * pmb.units.mol/pmb.units.L
        pep2_concentration = 1e-2 * pmb.units.mol/pmb.units.L
        c_salt=5e-3 * pmb.units.mol/ pmb.units.L
        model = '1beadAA'

        # Load pKa-values
        path_to_pka=pmb.get_resource("parameters/pka_sets/Nozaki1967.json")
        pmb.load_pka_set(path_to_pka)

        # Define the peptides in the pyMBE data frame
        pmb.define_peptide(name = "peptide_1",
                sequence = sequence1,
                model = model)

        pmb.define_peptide(name = "peptide_2",
                sequence = sequence2,
                model = model)


        print("*** Check that Henderson-Hasselbalch equation works correctly ***")

        # Calculate charge according to Henderson-Hasselbalch equation
        pH_range = np.linspace(2, 12, num=200)
        Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1",
                                pH_list = pH_range)
        Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2",
                                pH_list = pH_range)

        """
        with open(self.data_root / "HH.csv", "wb") as f:
            np.savetxt(f, np.asarray(Z_HH_1).reshape(1,-1), delimiter=",")
            np.savetxt(f, np.asarray(Z_HH_2).reshape(1,-1), delimiter=",")
        """

        data_path = pmb.get_resource(path=self.data_root)
        ref_data_HH = np.loadtxt(f"{data_path}/HH.csv", delimiter=",")
        np.testing.assert_allclose(Z_HH_1, ref_data_HH[0,:])
        np.testing.assert_allclose(Z_HH_2, ref_data_HH[1,:])

        print("*** Test passed ***\n")


        print("*** Check that Henderson-Hasselbalch equation + Donnan works correctly ***")

        HH_Donnan_dict = pmb.calculate_HH_Donnan(
                c_macro = {"peptide_1": pep1_concentration,
                           "peptide_2": pep2_concentration},
                c_salt = c_salt,
                pH_list = pH_range)

        """
        with open(self.data_root / "HH_Donnan.csv", "wb") as f:
            np.savetxt(f, np.asarray(HH_Donnan_dict["charges_dict"]["peptide_1"]).reshape(1,-1), delimiter=",")
            np.savetxt(f, np.asarray(HH_Donnan_dict["charges_dict"]["peptide_2"]).reshape(1,-1), delimiter=",")
        """

        ref_data_HH_Donnan = np.loadtxt(f"{data_path}/HH_Donnan.csv", delimiter=",")
        np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_1"], ref_data_HH_Donnan[0,:])
        np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_2"], ref_data_HH_Donnan[1,:])

        print("*** Test passed ***\n")


        print("*** Check that HH and HH_Don are consistent ***")

        Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1",
                                pH_list = HH_Donnan_dict["pH_system_list"])
        Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2",
                                pH_list = HH_Donnan_dict["pH_system_list"])

        np.testing.assert_allclose(Z_HH_1, HH_Donnan_dict["charges_dict"]["peptide_1"])
        np.testing.assert_allclose(Z_HH_2, HH_Donnan_dict["charges_dict"]["peptide_2"])


        print("*** Test passed***")

if __name__ == "__main__":
    ut.main()
