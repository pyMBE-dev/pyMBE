import io
import json
import contextlib
import unittest as ut
import numpy as np
import pandas as pd
import pyMBE
import lib.analysis as ana


class Serialization(ut.TestCase):

    def test_get_dt(self):
        data = pd.DataFrame.from_dict( {'time': [0, 1, 2,], 'obs': ['1.0', '2.0', '4.0']} )
        dt = ana.get_dt(data)
        self.assertAlmostEqual(dt, 1.0, delta = 1e-7)

if __name__ == "__main__":
    ut.main()
