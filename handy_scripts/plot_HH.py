# Plots the charge of a molecule with a given sequence using either a custom or reference set of pKa values
# Developed by:
# Dr. Pablo M. Blanco (Charles University) 

import sys
import numpy as np
import matplotlib.pyplot as plt

# Load the sugar library
sugar_path='.'
sys.path.insert(0, sugar_path) 

import sugar

# Create an instance of sugar

sg=sugar.sugar_library()

# Input parameters

sequence="nDEHKc"
pH_values = np.linspace(2, 12, num=21)
load_pka_set_from_file=False   # If set to false, uses custom_pka_set
path_to_pka_set_file='reference_parameters/pka_sets/Nozaki1967.txt' 

custom_pka_set={"D" : {"pka_value": 4.0, "acidity": "acidic"},
                "E" : {"pka_value": 4.4, "acidity": "acidic"},
                "H" : {"pka_value": 6.8, "acidity": "basic"},
                "K" : {"pka_value": 10.4, "acidity": "basic"},
                "n" : {"pka_value": 8.0, "acidity": "basic"},
                "c" : {"pka_value": 3.6, "acidity": "acidic"}}

if load_pka_set_from_file:
    pka_set=sg.load_pka_set(filename=path_to_pka_set_file)
    print('pka_set stored in sugar: ', sg.pka_set)
else:
    pka_set=custom_pka_set

# pka_set is an optional argument, if it is not provided sugar will use the one stored in sg.pka_set

Z_HH = sg.calculate_HH(sequence=sequence, pH=pH_values, pka_set=pka_set)

# Plot

plt.figure(figsize=[11, 9])
plt.suptitle('Peptide sequence: '+sequence)
plt.plot(pH_values, Z_HH, "-k", label='Henderson-Hasselbach')
plt.axhline(y=0.0, color="gray", linestyle="--")
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.legend()

plt.show()