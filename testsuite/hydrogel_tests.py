import sys
import pathlib
import tempfile
import subprocess
import multiprocessing
import numpy as np
import pandas as pd
from scipy import interpolate
import unittest as ut
from lib import analysis
import pyMBE

pmb = pyMBE.pymbe_library(seed=42)

# Values from Landsgesell2022
unit_length=0.355*pmb.units.nm
temperature=300*pmb.units.K
pmb.set_reduced_units(unit_length=unit_length,
                      temperature=temperature)

root = pathlib.Path(__file__).parent.parent.resolve()
data_root = root / "testsuite" / "hydrogel_tests_data"
script_path = root / "samples" / "weak_polyacid_hydrogel_grxmc.py"
mode = "test"
mpc = 40
ref_max_L = 89.235257606948
ref_cs = pmb.units.Quantity(0.00134770889, "1/reduced_length**3")
boltzmann_constant = 1.38065e-23
temperature_in_kelvin = 300
kT_SI = boltzmann_constant * temperature_in_kelvin
sigma_in_meter = 3.55e-10
conv_p = kT_SI / (sigma_in_meter ** 3) / 1e5

# Test cases separated into pressure and titration tests
pressure_test_cases = {
    "c_salt_res": 0.05,
    "pH": [5, 5],
    "L_fractions": [0.30259595, 0.50963529],
}

titration_test_cases = {
    "c_salt_res": 0.01,
    "swelling_eq": {
        4: 0.36630036630074175,
        6: 0.5574136008913612,
    },
}

def run_simulation(single_case, test_type):
    """
    Runs a simulation for a given single set of parameters and returns analyzed data.
    """
    with tempfile.TemporaryDirectory() as time_series_path:
        if test_type == "pressure":
            pH = single_case["pH"]
            L_fraction = single_case["L_fraction"]
            c_salt_res = single_case["c_salt_res"]
        elif test_type == "titration":
            pH = single_case["pH"]
            L_fraction = single_case["L_fraction"]
            c_salt_res = single_case["c_salt_res"]
        else:
            raise ValueError(f"Unknown test_type: {test_type}")

        print(f"Running simulation for pH={pH}, L_fraction={L_fraction}, c_salt_res={c_salt_res}")
        run_command = [
            sys.executable, script_path,
            "--mpc", str(mpc),
            "--csalt_res", str(c_salt_res),
            "--L_fraction", str(L_fraction),
            "--pH_res", str(pH),
            "--mode", "test",
            "--output", time_series_path
        ]
        subprocess.check_output(run_command)

        result_key = frozenset({
            "pH": pH,
            "c_salt_res": c_salt_res,
            "L_fraction": L_fraction
        }.items())

        data = analysis.analyze_time_series(path_to_datafolder=time_series_path)
        return (result_key, data)

class HydrogelTest(ut.TestCase):
    
    def test_pressure(self):
        test_cases = [
            {"c_salt_res": pressure_test_cases["c_salt_res"],
            "pH": pH,
            "L_fraction": L_frac}
            for pH, L_frac in zip(pressure_test_cases["pH"], pressure_test_cases["L_fractions"])
        ]
        with multiprocessing.Pool(processes=2) as pool:
            results = dict(pool.starmap(run_simulation, [(tc, "pressure") for tc in test_cases]))
        rtol = 0.4  # Relative tolerance
        atol = 0.4  # Absolute tolerance

        data_path = pmb.get_resource("testsuite/data")
        data_ref = pd.read_csv(f"{data_path}/Landsgesell2022a.csv")

        # Compare pressure values
        with self.subTest(msg=f"Testing pressure for c_salt_res={pressure_test_cases['c_salt_res']}"):
            pressures_ref = pmb.units.Quantity((data_ref[(data_ref["pH"] == 5) & (data_ref["cs_bulk"] == 0.00134770889)])["total_isotropic_pressures"].values,"reduced_energy/reduced_length**3")
            box_l_ref = data_ref[(data_ref["pH"] == 5) & (data_ref["cs_bulk"] == 0.00134770889)]["#final_box_l"].values
            data_path = pmb.get_resource("parameters/salt")
            monovalent_salt_ref_data=pd.read_csv(f"{data_path}/excess_chemical_potential_excess_pressure.csv")
            cs_bulk = pmb.units.Quantity(monovalent_salt_ref_data['cs_bulk_[1/sigma^3]'].values,"1/reduced_length**3")
            excess_press = pmb.units.Quantity(monovalent_salt_ref_data['excess_pressure_[kT/sigma^3]'].values, "reduced_energy/reduced_length**3")
            excess_press = interpolate.interp1d(cs_bulk.m_as("1/L"),
                                    excess_press.m_as("bar"))
            p_id = 2*ref_cs*pmb.kT
            p_res = p_id + excess_press(ref_cs.m_as("1/L"))*pmb.units.bar
            for case in test_cases:
                L_fraction = case["L_fraction"]
                pH = case["pH"]
                key = frozenset({
                   "pH": pH,
                   "c_salt_res": case["c_salt_res"],
                   "L_fraction": L_fraction
                }.items())
                index = np.where(np.isclose(box_l_ref/ref_max_L, L_fraction))[0][0]
                pressure_ref = (pressures_ref[index]-p_res).m_as("bar")
                results[key][("Lfraction", "value")] = pd.to_numeric(results[key][("Lfraction", "value")], errors='coerce')
                mask = np.isclose(results[key][("Lfraction","value")].astype(float), float(L_fraction),atol=0.01)
                test_pressure = results[key][mask][("mean","pressure")]
                test_pressure_value = test_pressure.iloc[0]  # or test_pressure.item()
                test_pressure = pmb.units.Quantity(test_pressure_value, "reduced_energy/reduced_length**3")
                p_sys_minus_p_res = test_pressure.m_as("bar") - p_res.m_as("bar")
                np.testing.assert_allclose(p_sys_minus_p_res, pressure_ref, rtol=rtol, atol=atol)
            
    def test_titration(self):
        test_cases = [
            {
                "c_salt_res": titration_test_cases["c_salt_res"],
                "pH": pH,
                "L_fraction": L_frac
            }
            for pH, L_frac in titration_test_cases["swelling_eq"].items()
        ]
        with multiprocessing.Pool(processes=2) as pool:
            results = dict(pool.starmap(run_simulation, [(tc, "titration") for tc in test_cases]))
        rtol = 0.05
        atol = 0.05
        data_path = pmb.get_resource("testsuite/data")
        data_ref = pd.read_csv(f"{data_path}/Landsgesell2022a.csv")

        with self.subTest(msg=f"Testing titration curve for c_salt_res={titration_test_cases['c_salt_res']}"):
            for case in test_cases:
                pH = case["pH"]
                key = frozenset({
                    "pH": pH,
                    "c_salt_res": case["c_salt_res"],
                    "L_fraction": case["L_fraction"]
                }.items())
                data_ref_filtered = data_ref[np.isclose(data_ref['cs'], case["c_salt_res"])].sort_values(by="pH")
                ref_alpha = data_ref_filtered[data_ref_filtered["pH"] == pH]["alpha"].values[0]
                results[key][("pH", "value")] = pd.to_numeric(results[key][("pH", "value")], errors='coerce')
                mask = np.isclose(results[key][("pH", "value")].astype(float), float(pH), atol=1e-5)
                test_alpha = results[key][mask][("mean", "alpha")]
                test_alpha_val = test_alpha.iloc[0]
                np.testing.assert_allclose(test_alpha_val, ref_alpha, rtol=rtol, atol=atol)

if __name__ == "__main__":
    ut.main()

