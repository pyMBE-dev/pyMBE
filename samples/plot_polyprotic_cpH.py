#
# Copyright (C) 2024-2026 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Plots the charge of a polyprotic acid or base from cpH simulation data and the corresponding analytical HH solution

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pathlib
import pandas as pd
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

parser = argparse.ArgumentParser(description='Plots the titration data from polyprotic_cpH.py and the corresponding analytical solution.')
parser.add_argument('--acidity',
                    type=str,
                    default='acidic',
                    choices=['acidic', 'basic'],
                    help='whether the particle is acidic or basic (default: acidic)')
parser.add_argument('--pka',
                    type=float,
                    nargs='+',
                    default=[2.16, 7.21, 12.32],
                    help='pKa values for each deprotonation step (default: 2.16 7.21 12.32, phosphoric acid)')
parser.add_argument('--path_to_data',
                    type=pathlib.Path,
                    required=False,
                    default=None,
                    help='path to the analyzed data (default: time_series/polyprotic_<acidity>_cpH/analyzed_data.csv)')
parser.add_argument('--output',
                    type=pathlib.Path,
                    required=False,
                    default=None,
                    help='output directory (default: time_series/polyprotic_<acidity>_cpH)')
parser.add_argument('--mode',
                    type=str,
                    required=False,
                    default="plot",
                    choices=["plot", "store_HH"],
                    help='mode to execute the script; available options are: store_HH (stores the analytical HH solution, used for testing) and plot (produces the plots)')
args = parser.parse_args()

acidity = args.acidity
pka_list = args.pka
n = len(pka_list)

if args.output is None:
    args.output = pathlib.Path(__file__).parent / "time_series" / f"polyprotic_{acidity}_cpH"
if args.path_to_data is None:
    args.path_to_data = args.output / "analyzed_data.csv"

# Define the molecule (necessary to calculate the HH analytical solution with pyMBE)
pmb.define_polyprotic_particle(name="particle",
                               sigma=1 * pmb.units('reduced_length'),
                               epsilon=1 * pmb.units('reduced_energy'),
                               n=n,
                               acidity=acidity,
                               pka_list=pka_list)

pmb.define_residue(name="Res_particle",
                   central_bead="particle",
                   side_chains=[])

pmb.define_molecule(name="polyprotic_molecule",
                    residue_list=["Res_particle"])

# Calculate the ideal titration curve with Henderson-Hasselbalch equation
if args.mode == "plot":
    pH_range_HH = np.linspace(0, 14, num=200)
elif args.mode == "store_HH":
    pH_range_HH = [1.0, 5.0, 10.0, 13.0]
Z_HH = pmb.calculate_HH(template_name="polyprotic_molecule",
                         pH_list=pH_range_HH)

if args.mode == "plot":
    # Read the analyzed data produced with polyprotic_cpH.py + analyze_time_series.py
    analyzed_data = pd.read_csv(args.path_to_data, header=[0, 1])

    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.errorbar(analyzed_data["pH"]["value"],
                analyzed_data["mean"]["charge"],
                yerr=analyzed_data["err_mean"]["charge"],
                fmt='o',
                capsize=3,
                color="red",
                label='cpH simulation')

    ax.plot(pH_range_HH,
            Z_HH,
            "--r",
            label='HH analytical solution')

    plt.legend()
    plt.xlabel('pH')
    plt.ylabel(f'Net charge of the {n}-protic {acidity} particle / e')
    plt.show()
    plt.close()

elif args.mode == "store_HH":
    HH_data = pd.DataFrame({"pH": pH_range_HH,
                             "Z_HH": Z_HH})
    HH_data.to_csv(args.output / "HH_data.csv",
                    index=False)
