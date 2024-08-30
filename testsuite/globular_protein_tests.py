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

from lib import analysis
import sys
import pathlib
import tempfile
import subprocess
import multiprocessing
import numpy as np
import pandas as pd
import unittest as ut


root = pathlib.Path(__file__).parent.parent.resolve()
data_root = root / "testsuite" / "globular_protein_tests_data"
script_path = root / "samples" / "Beyer2024" / "globular_protein.py"
test_pH_values = [2, 5, 7]
tasks = ["1beb", "1f6s"]
mode = "test"


def kernel(protein_pdb):
    """
    Runs a set of tests for a given protein pdb.

    Args:
        protein_pdb(`str`): PDB code of the protein.
    """
    with tempfile.TemporaryDirectory() as time_series_path:
        for pH in test_pH_values:
            print(f"pH = {pH}")
            run_command=[sys.executable, script_path, "--pdb", protein_pdb, "--pH", str(pH),
                         "--path_to_cg", f"parameters/globular_proteins/{protein_pdb}.vtf",
                         "--mode", "test", "--no_verbose", "--output", time_series_path]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path,
                                          filename_extension="_time_series.csv")
    return (protein_pdb, data)


class Test(ut.TestCase):

    def test_globular_protein(self):
        with multiprocessing.Pool(processes=2) as pool:
            results = dict(pool.map(kernel, tasks, chunksize=1))

        rtol=0.1 # relative tolerance
        atol=0.5 # absolute tolerance
        for protein_pdb, data in results.items():
            # Save data for future testing
            if mode == "save":
                data.to_csv(data_root / f"{protein_pdb}.csv", index=False)
                continue
            assert mode == "test", f"Mode {mode} not supported, valid modes: ['save', 'test']"
            with self.subTest(msg=f"Protein {protein_pdb}"):
                # Get reference test data
                ref_data=pd.read_csv(data_root / f"{protein_pdb}.csv", header=[0, 1])
                # Check charge
                test_charge=np.sort(data["mean","charge"].to_numpy())
                ref_charge=np.sort(ref_data["mean","charge"].to_numpy())
                np.testing.assert_allclose(
                    test_charge, ref_charge, rtol=rtol, atol=atol)

if __name__ == "__main__":
    ut.main()
