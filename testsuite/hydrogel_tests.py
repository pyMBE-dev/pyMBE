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
                      temperature=temperature,
                      verbose=False)

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
print(conv_p)
# Test cases separated into pressure and titration tests
#pressure_test_cases = {
#    "c_salt_res": 0.05,
#    "pH": 5,
#    "L_fractions": [0.30259595, 0.50963529],
#}

pressure_test_cases = {
    "c_salt_res": 0.05,
    "pH": 5,
    "L_fractions": [0.30259595],
}

#titration_test_cases = {
#    "c_salt_res": 0.01,
#    "swelling_eq": {
#        4: 0.36630036630074175,
#        6: 0.5574136008913612,
#    },
#}

titration_test_cases = {
    "c_salt_res": 0.01,
    "swelling_eq": {
        4: 0.36630036630074175}
}

def run_simulation(params, test_type):
    """
    Runs a simulation for a given set of parameters and returns analyzed data.
    """
    with tempfile.TemporaryDirectory() as time_series_path:
        if test_type == "pressure":
            pH_values = [params["pH"]]
            L_fractions = params["L_fractions"]
        else:
            pH_values = params["swelling_eq"].keys()
            L_fractions = params["swelling_eq"].values()

        for pH, L_fraction in zip(pH_values, L_fractions):
            print(f"Running simulation for pH={pH}, L_fraction={L_fraction}, c_salt_res={params['c_salt_res']}")
            run_command = [
                sys.executable, script_path,
                "--mpc", str(mpc),
                "--csalt_res", str(params["c_salt_res"]),
                "--L_fraction", str(L_fraction),
                "--pH_res", str(pH),
                "--mode", "test",
                "--no_verbose",
                "--output", time_series_path
            ]
            subprocess.check_output(run_command)

        # Analyze the data
        print(time_series_path)
        data = analysis.analyze_time_series(path_to_datafolder=time_series_path)
    params_hashable = {
        k: tuple(v) if isinstance(v, list) else
           tuple(v.items()) if isinstance(v, dict) else v
        for k, v in params.items()
    }

    return (frozenset(params_hashable.items()), data)


class HydrogelTest(ut.TestCase):
    def test_hydrogel(self):
        with multiprocessing.Pool(processes=2) as pool:
            results = dict(pool.starmap(run_simulation, [(pressure_test_cases, "pressure"), (titration_test_cases, "titration")]))
        rtol = 0.1  # Relative tolerance
        atol = 0.05  # Absolute tolerance

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
            pressure_key = frozenset(
            (k, tuple(v)) if isinstance(v, list) else
            (k, tuple(v.items())) if isinstance(v, dict) else
            (k, v)
            for k, v in pressure_test_cases.items())

            for L_fraction in pressure_test_cases["L_fractions"]:
                index = np.where(np.isclose(box_l_ref/ref_max_L, L_fraction))[0][0]
                pressure_ref = (pressures_ref[index]-p_res).m_as("bar")
                print("pressure_key.columns",results[pressure_key].columns)
                results[pressure_key][("Lfraction", "value")] = pd.to_numeric(results[pressure_key][("Lfraction", "value")], errors='coerce')
                mask = np.isclose(results[pressure_key][("Lfraction","value")].astype(float), float(L_fraction),atol=0.01)
                print(results[pressure_key][("Lfraction","value")].astype(float), float(L_fraction))
                print("p mask", mask)
                print("results[pressure_key][mask]",results[pressure_key][mask])
                test_pressure = results[pressure_key][mask][("mean","pressure")]
                print("test_pressure",test_pressure*conv_p)
                p_sys_minus_p_res = test_pressure*conv_p - p_res
                np.testing.assert_allclose(p_sys_minus_p_res, pressure_ref, rtol=rtol, atol=atol)

        # Compare titration curve values
        with self.subTest(msg=f"Testing titration curve for c_salt_res={titration_test_cases['c_salt_res']}"):
            for pH, ref_L in titration_test_cases["swelling_eq"].items():
                data_ref_filtered = data_ref[np.isclose(data_ref['cs'], titration_test_cases["c_salt_res"])].sort_values(by="pH")
                ref_alpha = data_ref_filtered[data_ref_filtered["pH"] == pH]["alpha"].values[0]
                titration_key = frozenset({('c_salt_res', titration_test_cases['c_salt_res']),
                           ('swelling_eq', tuple(titration_test_cases['swelling_eq'].items()))})
                print("ref_alpha",ref_alpha)
                results[titration_key][("pH", "value")] = pd.to_numeric(results[titration_key][("pH", "value")], errors='coerce')
                mask = np.isclose(results[titration_key][("pH", "value")].astype(float), float(pH), atol=1e-5)
                print("mask",mask)
                print("results[titration_key][mask]",results[titration_key][mask])
                test_alpha = results[titration_key][mask][("mean", "alpha")]
                print("test_alpha",test_alpha)
                np.testing.assert_allclose(test_alpha, ref_alpha, rtol=rtol, atol=atol)


if __name__ == "__main__":
    ut.main()

