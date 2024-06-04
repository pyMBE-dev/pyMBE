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

    def test_get_params_from_file_name(self):
        print("*** Unit test: test that get_params_from_file_name parses a filename without minus separtor ***")
        filename = 'density_0.001_N_1000_T_2.00.csv'
        correct_params = {'density': '0.001', 'N': '1000', 'T': '2.00'}
        params = ana.get_params_from_file_name(filename, 
                                                minus_separator = False)
        self.assertEqual(correct_params,params)
        print("*** Unit passed ***")

        print("*** Unit test: test that get_params_from_file_name parses a filename with minus separtor ***")
        filename = 'N-064_Solvent-good_Init-coil_observables.csv'
        correct_params = {'N': 64, 'Solvent': 'good', 'Init': 'coil'}
        params = ana.get_params_from_file_name(filename, 
                                                minus_separator = True)
        self.assertEqual(correct_params,params)
        print("*** Unit passed ***")

        print("*** Unit test: test that get_params_from_file_name parses a filename with a different extension ***")
        filename = 'density_0.001_N_1000_T_2.00_observables.txt'
        correct_params = {'density': '0.001', 'N': '1000', 'T': '2.00'}
        params = ana.get_params_from_file_name(filename, 
                                                minus_separator = False,
                                                filename_extension="_observables.txt")
        self.assertEqual(correct_params,params)
        print("*** Unit passed ***")

        print("*** Unit test: test that get_params_from_file_name raises a ValueError if a filename with a wrong formating is provided ***")
        inputs = {"file_name": 'density_0.001_N_1000_T_f_2.00_observables.txt',
                  "filename_extension": "_observables.txt"}
        self.assertRaises(ValueError, ana.get_params_from_file_name, **inputs)
        print("*** Unit passed ***")

if __name__ == "__main__":
    print("*** lib.analysis unit tests ***")
    ut.main()
