"""
Compare pyMBE's polyprotic Henderson-Hasselbalch implementation against:
  1. The analytical calculate_Z formula from pbs_simulations
  2. Ideal constant-pH simulation data from pbs_simulations sample_data

Produces comparison plots for monoprotic, diprotic, and triprotic acids.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pyMBE

# Add pbs_simulations analysis to path
PBS_ROOT = os.path.expanduser("~/ntnu/pbs_simulations")
sys.path.insert(0, os.path.join(PBS_ROOT, "scripts"))
from analysis import analyze_time_series

# ---------- Reference analytical formula from pbs_simulations ----------

def calculate_Z_reference(pK, pH_array):
    """Exact polyprotic HH formula from pbs_simulations/scripts/post_processing.py"""
    pK = np.asarray(pK)
    pH_array = np.atleast_1d(pH_array).astype(float)
    y = pH_array[None, :]  # (1, N_pH)
    i_range = np.arange(1, len(pK) + 1)[:, None]  # (n, 1)
    cumsum_pK = np.cumsum(pK)[:, None]  # (n, 1)
    fractions = 10.0 ** (-cumsum_pK + i_range * y)
    numerator = np.sum(i_range * fractions, axis=0)
    denominator = 1.0 + np.sum(fractions, axis=0)
    return -numerator / denominator


def compute_sim_Z(sample_dir, n):
    """Compute average charge Z from simulation sample data using block binning."""
    analyzed = analyze_time_series(path_to_datafolder=sample_dir,
                                   ignore_files=["analyzed_data.csv"])
    pH_vals = []
    Z_vals = []
    Z_err = []
    for idx in range(len(analyzed)):
        row = analyzed.iloc[idx]
        pH_val = float(row[("pH", "value")])
        n_ha = float(row[("mean", "n_ha")])
        Z_num = -sum(float(row[("mean", f"n_a{i+1}")]) * (i + 1) for i in range(n))
        Z_den = n_ha + sum(float(row[("mean", f"n_a{i+1}")]) for i in range(n))
        pH_vals.append(pH_val)
        Z_vals.append(Z_num / Z_den if Z_den > 0 else 0.0)
    order = np.argsort(pH_vals)
    return np.array(pH_vals)[order], np.array(Z_vals)[order]


def compute_pyMBE_Z(pka_list, pH_array):
    """Compute charge using pyMBE's calculate_HH."""
    pmb = pyMBE.pymbe_library(seed=42)
    n = len(pka_list)
    if n == 1:
        pmb.define_particle(name="acid",
                            sigma=0.35 * pmb.units.nm,
                            epsilon=1 * pmb.units("reduced_energy"),
                            acidity="acidic",
                            pka=pka_list[0])
    else:
        pmb.define_polyprotic_particle(name="acid",
                                       sigma=0.35 * pmb.units.nm,
                                       epsilon=1 * pmb.units("reduced_energy"),
                                       n=n,
                                       acidity="acidic",
                                       pka_list=pka_list)
    pmb.define_residue(name="res", central_bead="acid", side_chains=[])
    pmb.define_molecule(name="mol", residue_list=["res"])
    return pmb.calculate_HH(template_name="mol", pH_list=list(pH_array))


# ---------- Configuration ----------

pK_values = {
    "monoprotic": [2.16],
    "diprotic": [2.16, 7.21],
    "triprotic": [2.16, 7.21, 12.32],
}

sample_dirs = {
    "monoprotic": os.path.join(PBS_ROOT, "sample_data", "ideal-monoprotic"),
    "diprotic": os.path.join(PBS_ROOT, "sample_data", "ideal-diprotic"),
    "triprotic": os.path.join(PBS_ROOT, "sample_data", "ideal-triprotic"),
}

# ---------- Plot ----------

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

for ax, label in zip(axes, ["monoprotic", "diprotic", "triprotic"]):
    pK = pK_values[label]
    n = len(pK)
    pH_fine = np.linspace(0.5, 14.5, 200)

    # 1. Reference analytical formula
    Z_ref = calculate_Z_reference(pK, pH_fine)

    # 2. pyMBE calculate_HH
    Z_pyMBE = compute_pyMBE_Z(pK, pH_fine)

    # 3. Simulation data
    pH_sim, Z_sim = compute_sim_Z(sample_dirs[label], n)

    # Plot
    ax.plot(pH_fine, Z_ref, "k-", lw=2, label="Analytical (pbs_simulations)")
    ax.plot(pH_fine, Z_pyMBE, "r--", lw=2, label="pyMBE calculate_HH")
    ax.scatter(pH_sim, Z_sim, c="blue", s=30, zorder=5, label="Simulation data (ideal)")

    # pKa markers
    for i, pk in enumerate(pK):
        ax.axvline(x=pk, color="gray", ls=":", alpha=0.6)
        ax.text(pk + 0.2, -n + 0.3, f"pK{i+1}={pk}", fontsize=8, color="gray")

    ax.set_xlabel("pH")
    ax.set_title(f"{n}-protic acid")
    ax.grid(alpha=0.3)

    # MSE between pyMBE and reference
    mse = np.mean((np.array(Z_pyMBE) - Z_ref) ** 2)
    ax.text(0.05, 0.05, f"MSE(pyMBE vs analytical) = {mse:.2e}",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

axes[0].set_ylabel("Average charge Z")
axes[0].legend(loc="lower left", fontsize=8)

plt.suptitle("Polyprotic HH: pyMBE vs Analytical vs Simulation", fontsize=14)
plt.tight_layout()
plt.savefig("tests/polyprotic_HH_comparison.png", dpi=150, bbox_inches="tight")
print("Plot saved to tests/polyprotic_HH_comparison.png")
plt.show()
