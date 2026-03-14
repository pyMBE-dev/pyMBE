"""
Compare pyMBE's polyprotic HH for basic particles against the analytical formula.
No simulation data available for bases — analytical comparison only.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pyMBE


def calculate_Z_base_reference(pK, pH_array):
    """
    Analytical polyprotic HH for a basic particle.
    Fully protonated has charge +n, each deprotonation removes +1.
    Z = n - sum(i * 10^(-cumsum(pK)[i] + i*pH)) / (1 + sum(10^(-cumsum(pK)[i] + i*pH)))
    """
    pK = np.asarray(pK)
    n = len(pK)
    pH_array = np.atleast_1d(pH_array).astype(float)
    y = pH_array[None, :]
    i_range = np.arange(1, n + 1)[:, None]
    cumsum_pK = np.cumsum(pK)[:, None]
    fractions = 10.0 ** (-cumsum_pK + i_range * y)
    numerator = np.sum(i_range * fractions, axis=0)
    denominator = 1.0 + np.sum(fractions, axis=0)
    return n - numerator / denominator


def compute_pyMBE_Z_base(pka_list, pH_array):
    """Compute charge using pyMBE's calculate_HH for a basic particle."""
    pmb = pyMBE.pymbe_library(seed=42)
    n = len(pka_list)
    if n == 1:
        pmb.define_particle(name="base",
                            sigma=0.35 * pmb.units.nm,
                            epsilon=1 * pmb.units("reduced_energy"),
                            acidity="basic",
                            pka=pka_list[0])
    else:
        pmb.define_polyprotic_particle(name="base",
                                       sigma=0.35 * pmb.units.nm,
                                       epsilon=1 * pmb.units("reduced_energy"),
                                       n=n,
                                       acidity="basic",
                                       pka_list=pka_list)
    pmb.define_residue(name="res", central_bead="base", side_chains=[])
    pmb.define_molecule(name="mol", residue_list=["res"])
    return pmb.calculate_HH(template_name="mol", pH_list=list(pH_array))


pK_values = {
    "monoprotic": [9.25],
    "diprotic": [6.0, 10.0],
    "triprotic": [4.0, 8.0, 12.0],
}

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
pH_fine = np.linspace(0.5, 14.5, 200)

for ax, label in zip(axes, ["monoprotic", "diprotic", "triprotic"]):
    pK = pK_values[label]
    n = len(pK)

    Z_ref = calculate_Z_base_reference(pK, pH_fine)
    Z_pyMBE = compute_pyMBE_Z_base(pK, pH_fine)

    ax.plot(pH_fine, Z_ref, "k-", lw=2, label="Analytical")
    ax.plot(pH_fine, Z_pyMBE, "r--", lw=2, label="pyMBE calculate_HH")

    for i, pk in enumerate(pK):
        ax.axvline(x=pk, color="gray", ls=":", alpha=0.6)
        ax.text(pk + 0.2, 0.3, f"pK{i+1}={pk}", fontsize=8, color="gray")

    ax.set_xlabel("pH")
    ax.set_title(f"{n}-protic base")
    ax.grid(alpha=0.3)

    mse = np.mean((np.array(Z_pyMBE) - Z_ref) ** 2)
    ax.text(0.05, 0.05, f"MSE = {mse:.2e}",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

axes[0].set_ylabel("Average charge Z")
axes[0].legend(loc="upper right", fontsize=8)

plt.suptitle("Polyprotic HH (basic): pyMBE vs Analytical", fontsize=14)
plt.tight_layout()
plt.savefig("tests/polyprotic_HH_comparison_base.png", dpi=150, bbox_inches="tight")
print("Plot saved to tests/polyprotic_HH_comparison_base.png")
plt.show()
