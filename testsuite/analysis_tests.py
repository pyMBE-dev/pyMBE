import unittest as ut
import pandas as pd
import lib.analysis as ana


class Serialization(ut.TestCase):

    def test_get_dt(self):
        print("*** Unit test: test that analysis.get_dt returns the right time step ***")
        data = pd.DataFrame.from_dict( {'time': [0, 1, 2,], 'obs': ['1.0', '2.0', '4.0']} )
        dt, n_warnings = ana.get_dt(data)
        self.assertAlmostEqual(dt, 1.0, delta = 1e-7)
        self.assertEqual(n_warnings, 0)
        print("*** Unit passed ***")

        print("*** Unit test: test that analysis.get_dt prints a warning if there are values with repeated time steps ***")
        data = pd.DataFrame.from_dict( {'time': [0, 1, 1,], 'obs': ['1.0', '2.0', '4.0']} )
        dt, n_warnings = ana.get_dt(data,verbose=True)
        self.assertAlmostEqual(dt, 1.0, delta = 1e-7)
        self.assertEqual(n_warnings, 1)
        print("*** Unit passed ***")

        print("*** Unit test: test that analysis.get_dt raises a ValueError if the column with the time is not found ***")
        data = pd.DataFrame.from_dict( {'ns': [0, 1, 2,], 'obs': ['1.0', '2.0', '4.0']} )
        inputs = {"data": data}
        self.assertRaises(ValueError, ana.get_dt, **inputs)
        
        print("*** Unit passed ***")

        print("*** Unit test: test that analysis.get_dt raises a ValueError if the time is stored at uneven intervals ***")
        data = pd.DataFrame.from_dict( {'time': [0, 1, 4,], 'obs': ['1.0', '2.0', '4.0']} )
        inputs = {"data": data}
        self.assertRaises(ValueError, ana.get_dt, **inputs)
        
        print("*** Unit passed ***")

if __name__ == "__main__":
    print("*** lib.analysis unit tests ***")
    ut.main()
