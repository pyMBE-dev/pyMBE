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


# Tests that samples/peptide.py and analysis.py work properly

import sys
import pathlib
import tempfile
import subprocess
import multiprocessing
from lib import analysis
import numpy as np
import pandas as pd
import unittest as ut

root = pathlib.Path(__file__).parent.parent.resolve()
data_root = root / "samples" / "time_series" / "peptide"
sample_path = root / "samples" / "peptide.py"
analysis_path = root / "samples" / "analysis.py"
reference_path = root / "testsuite" / "peptide_tests_data"
sequences=["EEEEDDDD"]
test_pH_values=[3,7]
mode = "save"

def kernel(sequence):
    """
    Runs a set of tests for a given peptide sequence.

    Args:
        sequence(`str`): Amino acid sequence of the peptide.
    """
    with tempfile.TemporaryDirectory() as time_series_path:
        for pH in test_pH_values:
            print(f"pH = {pH}")
            run_command=[sys.executable, sample_path, "--sequence", sequence,
                         "--pH", str(pH), "--output", time_series_path]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        run_command=[sys.executable, analysis_path, "--data_folder", time_series_path]
        subprocess.check_output(run_command)
        analyzed_data=pd.read_csv(time_series_path + "/analyzed_data.csv", header=[0, 1])
        # Save data for future testing
        if mode == "save":
            analyzed_data.to_csv(reference_path / f"peptide_sample.csv", index=False)
        else:
            assert mode == "test", f"Mode {mode} not supported, valid modes: ['save', 'test']"

        
    return analyzed_data


class Test(ut.TestCase):

    def test_peptide(self):
        with multiprocessing.Pool(processes=2) as pool:
            analyzed_data=pool.map(kernel, sequences, chunksize=2)[0]
        rtol=0.1 # relative tolerance
        atol=0.5 # absolute tolerance
        reference_data=pd.read_csv(reference_path / f"peptide_sample.csv", header=[0, 1])
        pd.testing.assert_series_equal (analyzed_data[("mean","charge")], 
                                       reference_data[("mean","charge")], 
                                       rtol=rtol, 
                                       atol=atol)
        print(analyzed_data[("mean","charge")])
if __name__ == "__main__":
    ut.main()