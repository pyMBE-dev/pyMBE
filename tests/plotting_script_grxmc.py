import numpy as np
import os
import pickle
import itertools
import math
import matplotlib.pyplot as plt
import pandas as pd

width = 3.033
height = 0.9*3.033

plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r"\usepackage{mathptmx}")
plt.rcParams["font.family"] = "serif"
plt.tight_layout()

import matplotlib as mpl
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.linewidth'] = 1.0

def ideal_alpha_acid(pH, pKa):
    return 1. / (1 + 10**(pKa - pH))

def calculate_mean(x):
    return np.mean(x)

def autocorrelation(data, normalized=True):
    """
    Compute autocorrelation using FFT
    """
    nobs = len(data)
    corr_data = data - data.mean()
    n = 2**int(math.log(nobs, 2))
    corr_data = corr_data[:n]
    Frf = np.fft.fft(corr_data)
    acf = np.fft.ifft(Frf * np.conjugate(Frf))/corr_data.shape[0]
    if normalized:
        acf /= acf[0]
    acf = np.real(acf)
    # only return half of the ACF 
    # (see 4.3.1 "Kreuzkorrelationsfunktion" 
    # of https://github.com/arnolda/padc)
    return acf[:int(corr_data.shape[0]/2)]

def calculate_uncorrelated_samples(data):
    """
    Error estimation for time series of simulation observables and take into
    account that these series are correlated (which
    enhances the estimated statistical error).
    """
    # calculate the normalized autocorrelation function of data
    acf = autocorrelation(data)
    # calculate the integrated correlation time tau_int
    # (Janke, Wolfhard. "Statistical analysis of simulations: Data correlations
    # and error estimation." Quantum Simulations of Complex Many-Body Systems:
    # From Theory to Algorithms 10 (2002): 423-445.)
    tau_int = 0.5
    for i in range(len(acf)):
        tau_int += acf[i]
        if i >= 6 * tau_int:
            break
    # mean value of the time series
    data_mean = np.mean(data)
    # calculate the so called effective length of the time series N_eff
    if tau_int > 0.5:
        N_eff = int(len(data) / (2.0 * tau_int))
        # finally the error is sqrt(var(data)/N_eff)
        #stat_err = np.sqrt(np.var(data) / N_eff)
    else:
        #stat_err = np.sqrt(np.var(data) / len(data))
        N_eff = len(data)

    delta = int(len(data)//N_eff)
    uncorrelated_mean_values = [np.mean(data[i*delta:(i+1)*delta]) for i in range(N_eff)]
    return uncorrelated_mean_values 

def calculate_uncorrelated_error(array):
    N = len(array)
    if N == 1:
        return np.nan
    else:
        return np.sqrt((np.dot(array, array) - np.mean(array) ** 2 * N ) / (N * (N - 1)))

if __name__ == "__main__":

    # Conversion between SI and simulation units
    T_SI = 300
    K_BOLTZMANN_SI = 1.3806503 * 10**-23
    SIGMA_SI = 3.55 * 10**-10
    CONVERSION_FACTOR_PRESSURE = T_SI * K_BOLTZMANN_SI / (10**5 * SIGMA_SI**3)
    sigma = 3.55e-10 # Sigma in SI units
    avo = 6.022e+23 # Avogadro's number in SI units
    pref = 1/(10**3 * avo * sigma**3) # Prefactor to mol/L
    N_MONOMERS = 40
    LJ_SIGMA = 1.0
    BOND_LENGTH = 0.96 * LJ_SIGMA
    INITIAL_BOX_LENGTH = (N_MONOMERS + 1) * BOND_LENGTH / (0.25 * np.sqrt(3))

    print("pref:", pref)

    # Parameters
    list_c_salt = [0.01] 
    list_c_monomer = [0.435] 
    list_pKa_values = [4.0]
    list_pH_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0]

    index = pd.MultiIndex.from_product([list_c_salt, list_c_monomer, list_pH_values, list_pKa_values], names = ["c_salt", "c_monomer", "pH", "pKa"])
    data = pd.DataFrame(index=index).reset_index()

    data["alpha"] = np.nan
    data["error_alpha"] = np.nan

    for c_salt, c_monomer, pH, pKa in itertools.product(list_c_salt, list_c_monomer, list_pH_values, list_pKa_values):

        filename = os.path.abspath('.') + "/data/data_grxmc/pKa" + str(pKa) + "/salt_concentration" + str(c_salt) + "/monomer_concentration" + str(c_monomer) + "/pH" + str(pH) +  "/data.pkl"
        print(filename)
        if os.path.isfile(filename):
            simulation_data = pickle.load(open(filename, "rb"))
            print('Reading', filename)

            # Calculate mean values of the various observables
            alpha_samples = calculate_uncorrelated_samples(np.asarray(simulation_data["alphas"]))
            partition_coefficient_samples = calculate_uncorrelated_samples(np.asarray(simulation_data["partition_coefficient"]))

            data.loc[(data["c_salt"] == c_salt) & (data["c_monomer"] == c_monomer) & (data["pH"] == pH) & (data["pKa"] == pKa), ["alpha"]] = np.mean(alpha_samples)
            data.loc[(data["c_salt"] == c_salt) & (data["c_monomer"] == c_monomer) & (data["pH"] == pH) & (data["pKa"] == pKa), ["error_alpha"]] = calculate_uncorrelated_error(alpha_samples)

            data.loc[(data["c_salt"] == c_salt) & (data["c_monomer"] == c_monomer) & (data["pH"] == pH) & (data["pKa"] == pKa), ["partition_coefficient"]] = np.mean(partition_coefficient_samples)
            data.loc[(data["c_salt"] == c_salt) & (data["c_monomer"] == c_monomer) & (data["pH"] == pH) & (data["pKa"] == pKa), ["error_partition_coefficient"]] = calculate_uncorrelated_error(partition_coefficient_samples)

            pH_range = simulation_data["pH_range"]
            alpha_HH = simulation_data["alpha_HH"]
            alpha_HH_Donnan = simulation_data["alpha_HH_Donnan"]


    reference_data = pd.read_csv(open(os.path.abspath('..') + '/reference_data/data_landsgesell.csv'), sep="\t", index_col=False)

    for c_salt, c_monomer, pKa in itertools.product(list_c_salt, list_c_monomer, list_pKa_values):
        filtered_simulation_data = data.loc[(data["c_salt"] == c_salt) & (data["c_monomer"] == c_monomer) & (data["pKa"] == pKa)]

        filtered_reference_data = reference_data.loc[np.isclose(reference_data["cs_bulk"], c_salt/pref, rtol=1e-03) & np.isclose(reference_data["c_poly_[mol/l]"], c_monomer/50, rtol=1e-03) & np.isclose(reference_data["Kcideal_in_mol_per_l"], 10**(-pKa), rtol=1e-03)]
        
        plt.plot(pH_range, alpha_HH, label=r"HH", color="black")
        plt.plot(pH_range, alpha_HH_Donnan, label=r"HH+Don", color="#666666", linestyle="dashed")

        plt.errorbar(filtered_reference_data["pH"], filtered_reference_data["degree_of_dissociation"], yerr=filtered_reference_data["err_degree_of_dissociation"], linestyle="none", marker="s", label="Landsgesell et al.", color="C0")
        plt.errorbar(filtered_simulation_data["pH"], filtered_simulation_data["alpha"], yerr=filtered_simulation_data["error_alpha"], linestyle="none", marker="o", label="pyMBE", color="C1", fillstyle="none", markeredgewidth=1.5)
        plt.ylabel(r"Degree of ionization $\alpha$")
        plt.xlabel(r"pH in the reservoir")
        #plt.legend()
        plt.legend(frameon=False, loc="lower left", fontsize=9, bbox_to_anchor=(0,1.02,1,0.2), mode="expand", borderaxespad=0, ncol=2)
        plt.xticks([2,4,6,8,10,12])
        plt.grid(which='major', color='#CCCCCC', linestyle='--', linewidth=0.6)
        plt.tight_layout()
        #plt.show()
        plt.savefig('./comparison_landsgesell.pdf', bbox_inches="tight")
        plt.close()
