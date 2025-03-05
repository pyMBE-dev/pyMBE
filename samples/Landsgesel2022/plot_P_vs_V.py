import matplotlib.pyplot as plt
import pandas as pd

# Load the processed data
plot_data_1 = pd.read_csv("plot_pv_curve.csv")
plot_data_2 = pd.read_csv("plot_ref.csv")

L_target = plot_data_1["L/L_max"]
p_sys_minus_p_res = plot_data_1["P_sys - P_res"]
p_err_total = plot_data_1["P_err_total"]

box_l_ref = plot_data_2["box_l"]
p_ref = plot_data_2["p_ref"]

# Plot
plt.errorbar(L_target, p_sys_minus_p_res, yerr=p_err_total, fmt="o", capsize=3, elinewidth=2,
             ecolor='tab:blue', label="pH=5, cs=0.05M")
plt.xlabel(r"$L/L_{max}$", fontsize=17)
plt.ylabel(r"$P_{sys}-P_{res}$ [bar]", fontsize=17)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.plot(box_l_ref, p_ref, label="Landsgesell et al")
plt.legend(fontsize=17)
plt.show()
