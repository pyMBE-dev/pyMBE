import numpy as np
import matplotlib.pyplot as plt
import pyMBE
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Plots alpha vs pH  from weak_gel.py and the corresponding reference data from Landsgesell2022.')
parser.add_argument('--path_to_data',
                    type=str,
                    required= False,
                    default="samples/Landsgesel2022/time_series/analyzed_data.csv",
                    help='path to the analyzed data')
args = parser.parse_args()

pmb = pyMBE.pymbe_library(seed=42)

# Read the reference data
data_path = pmb.get_resource("testsuite/data")
data_ref = pd.read_csv(f"{data_path}/Landsgesell2022a.csv")
data_ref['cs'] = pd.to_numeric(data_ref['cs'], errors='coerce')
data_ref = data_ref[np.isclose(data_ref['cs'], 0.01)].sort_values(by="pH")
# Read the analyzed time series
time_series_folder_path=pmb.get_resource(args.path_to_data)
analyzed_data = pd.read_csv(time_series_folder_path, header=[0,1])
analyzed_data_to_plot={"pH" : [],
                       "alpha": [],
                       "alpha_err": []    
                       }
for pH_value in data_ref["pH"].drop_duplicates():
    sorted_ref_df = data_ref[data_ref["pH"] == pH_value]
    reference_L_fraction = sorted_ref_df["best_L"].values[0]/sorted_ref_df["max_L"].values[0]
    sorted_analyzed_df = analyzed_data[
        (analyzed_data["pH"]["value"] == pH_value) & (abs(analyzed_data["Ltarget"]["value"] - reference_L_fraction) < 1)
    ]
    if not sorted_analyzed_df.empty:
        analyzed_data_to_plot["pH"].append(sorted_analyzed_df["pH"]["value"].values[0])
        analyzed_data_to_plot["alpha"].append(sorted_analyzed_df["mean"]["alpha"].values[0])
        analyzed_data_to_plot["alpha_err"].append(sorted_analyzed_df["err_mean"]["alpha"].values[0])


plt.xlabel("pH in the reservoir",fontsize=12)
plt.ylabel(r"degree of ionization $\alpha$", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.errorbar(analyzed_data_to_plot["pH"], 
             analyzed_data_to_plot["alpha"], 
             yerr = analyzed_data_to_plot["alpha_err"], 
             fmt='o', 
             capsize=3, 
             elinewidth=2, 
             ecolor='tab:blue', 
             label=r"pyMBE")
plt.plot(data_ref["pH"].drop_duplicates(),
         data_ref["alpha"].drop_duplicates(),
         marker="*",
         label="Landsgesell et al. 2022")

plt.legend(fontsize=12, loc="lower right")
plt.show()
