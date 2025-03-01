import numpy as np
import matplotlib.pyplot as plt
import pyMBE
import pandas as pd

pmb = pyMBE.pymbe_library(seed=42)
pH = np.array([4,6])
alphas = []
err_alphas = []
L_target=[0.366,0.557]
for l_des,ph in zip(L_target,pH):
    data_path = pmb.get_resource("testsuite/data/src")
    data = pd.read_csv(f"{data_path}/csalt_0.01_L_target_{l_des}_pH_{ph}_pKa_4_time_series.csv")
    alphas.append(np.mean(data["alpha"]))
    err_alphas.append(np.std(data["alpha"]))

x = [1,2,3,4,5,6,7,8,9,10,11,12,13]
path_to_eq_gel_MD = pmb.get_resource("testsuite/data/src/")
df = pd.read_csv(f"{path_to_eq_gel_MD}/equilibrium_values_gel_MD.csv")
df['cs'] = pd.to_numeric(df['cs'], errors='coerce')
alphas_ref = df[np.isclose(df['cs'], 0.01)]["alpha"]

plt.errorbar(pH, alphas, yerr = err_alphas, fmt='o', capsize=3, elinewidth=2, ecolor='tab:blue', label=r"$c_{\mathrm{salt}}$=0.01M")
plt.xlabel("pH in the reservoir",fontsize=17)
plt.ylabel(r"degree of ionization $\alpha$", fontsize=17)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.plot(x,alphas_ref,marker="*",label="Landsgesell et al.")
plt.legend(fontsize=17)
plt.show()
