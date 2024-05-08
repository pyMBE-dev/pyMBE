# Plots the charge of a peptide with a given sequence using either a custom or reference set of pKa values

import os
import sys
import inspect
import numpy as np
import matplotlib.pyplot as plt

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(SEED=42)

# Input parameters

sequence="nDEHKc"
pH_values = np.linspace(2, 12, num=21)
load_pka_set_from_file=True   # If set to false, uses custom_pka_set
path_to_pka=pmb.get_resource("parameters/pka_sets/Nozaki1967.json")

custom_pka_set={"D" : {"pka_value": 4.0, "acidity": "acidic"},
                "E" : {"pka_value": 4.4, "acidity": "acidic"},
                "H" : {"pka_value": 6.8, "acidity": "basic"},
                "K" : {"pka_value": 10.4, "acidity": "basic"},
                "n" : {"pka_value": 8.0, "acidity": "basic"},
                "c" : {"pka_value": 3.6, "acidity": "acidic"}}

pmb.define_peptide(name="example_pep",
                   sequence=sequence,
                   model="1beadAA")

if load_pka_set_from_file:
    pka_set=pmb.load_pka_set(filename=path_to_pka)
    print('pka_set stored in pyMBE: ', pmb.get_pka_set())
else:
    pka_set=custom_pka_set

# pka_set is an optional argument, if it is not provided sugar will use the one stored in pmb.pka_set

Z_HH = pmb.calculate_HH(molecule_name="example_pep", 
                        pH_list=pH_values, 
                        pka_set=pka_set)

# Plot

plt.figure(figsize=[11, 9])
plt.suptitle('Peptide sequence: '+sequence)
plt.plot(pH_values, Z_HH, "-k", label='Henderson-Hasselbach')
plt.axhline(y=0.0, color="gray", linestyle="--")
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.legend()

plt.show()
