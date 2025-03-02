import numpy as np
import pandas as pd
from lib import analysis
import pyMBE
from scipy import interpolate

pmb = pyMBE.pymbe_library(seed=42)

boltzmann_constant = 1.38065e-23
temperature_in_kelvin = 300
kT_SI = boltzmann_constant * temperature_in_kelvin
sigma_in_meter = 3.55e-10
m_to_nm = 1e9
sigma_in_nm = sigma_in_meter * m_to_nm
conv_p = kT_SI / (sigma_in_meter ** 3) / 1e5

def p_id_bar(c_molar, T=300):
    c_m3_to_c_liter = 1e3
    pa_to_bar = 1e-5
    R = 8.314  # J K^-1 mol^{-1}
    return pa_to_bar * c_m3_to_c_liter * c_molar * R * T

L_target = np.array([0.35, 0.55])
csalt = 0.05
max_L = 89.235257606948  # Maximum chain length for MPC=39
pH = 5
pKa = 4

p_id_bar_value = p_id_bar(csalt) * 2

path_to_ex_pot = pmb.get_resource("testsuite/data/src/")
data = pd.read_csv(f'{path_to_ex_pot}/excess_chemical_potential_excess_pressure.csv')
excess_press = interpolate.interp1d(data["#cs_bulk_[1/sigma^3]"] * 37.115, 
                                    data["excess_pressure_[kT/sigma^3]"] * conv_p)
p_res = p_id_bar_value + excess_press(csalt)
p_res_err = data["std_error_excess_pressure_[kT/sigma^3]"][11] * conv_p

pressures = []
err_pressures = []

for l_des in L_target:
    data_path = pmb.get_resource("testsuite/data/src")
    data = pd.read_csv(f"{data_path}/csalt_{csalt}_L_target_{l_des}_pH_{pH}_pKa_{pKa}_time_series.csv")
    analyzed_data = analysis.block_analyze(data)
    d_mean = analyzed_data["mean"]
    d_err = analyzed_data["err_mean"]
    pressures.append(d_mean["pressure"])
    err_pressures.append(d_err["pressure"])

plotting_data = pd.DataFrame({
    "L/L_max": L_target,
    "P_sys - P_res": [p * conv_p - p_res for p in pressures],
    "P_err_total": [p_res_err + p for p in err_pressures]
})
plotting_data.to_csv("plot_pv_curve.csv", index=False)

path_to_wk_gel_data = pmb.get_resource("testsuite/data/src/")
df = pd.read_csv(f"{path_to_wk_gel_data}/weak-gel_total_data.csv")
pressures_ref = df[(df["pH"] == 5) & (df["cs_bulk"] == 0.00134770889)]["total_isotropic_pressures"]
final_box_l = df[(df["pH"] == 5) & (df["cs_bulk"] == 0.00134770889)]["#final_box_l"]

p_ref = np.array([p * conv_p - p_res for p in pressures_ref])
box_l_ref = np.array(final_box_l) / max_L

plotting_ref = pd.DataFrame({
    "box_l": box_l_ref,
    "p_ref": p_ref
})
plotting_ref.to_csv("plot_ref.csv", index=False)


