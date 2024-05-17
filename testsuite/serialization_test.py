import io
import json
import contextlib
import unittest as ut
import numpy as np
import pandas as pd
import pyMBE
import lib.analysis


class Serialization(ut.TestCase):

    def test_json_encoder(self):
        encoder = pyMBE.pymbe_library.NumpyEncoder
        # Python types
        self.assertEqual(json.dumps(1, cls=encoder), "1")
        self.assertEqual(json.dumps([1, 2], cls=encoder), "[1, 2]")
        self.assertEqual(json.dumps((1, 2), cls=encoder), "[1, 2]")
        self.assertEqual(json.dumps({1: 2}, cls=encoder), """{"1": 2}""")
        # NumPy types
        self.assertEqual(json.dumps(np.array([1, 2]), cls=encoder), "[1, 2]")
        self.assertEqual(json.dumps(np.array(1), cls=encoder), "1")
        self.assertEqual(json.dumps(np.int32(1), cls=encoder), "1")
        # Pandas types
        with self.assertRaisesRegex(TypeError, "Object of type Series is not JSON serializable"):
            json.dumps(pd.Series([1, 2]), cls=encoder)

    def test_parameters_to_path(self):
        params = {"kT": 2., "phi": -np.pi, "n": 3, "fene": True, "name": "pep"}
        name = lib.analysis.built_output_name(params)
        self.assertEqual(name, "kT-2_phi--3.14_n-3_fene-True_name-pep")
        params_out = lib.analysis.get_params_from_dir_name(name)
        params_ref = {"kT": "2", "phi": "-3.14", "n": "3",
                      "fene": "True", "name": "pep"}
        self.assertEqual(params_out, params_ref)

    def test_pint_units(self):
        ref_output = [
            "Current set of reduced units:",
            "0.355 nanometer = 1 reduced_length",
            "4.1164e-21 joule = 1 reduced_energy",
            "1.6022e-19 coulomb = 1 reduced_charge",
            "Temperature: 298.15 kelvin",
        ]
        pmb = pyMBE.pymbe_library(SEED=42)
        with contextlib.redirect_stdout(io.StringIO()) as f:
            pmb.print_reduced_units()
        self.assertEqual(f.getvalue().strip("\n").split("\n"), ref_output)


if __name__ == "__main__":
    ut.main()
