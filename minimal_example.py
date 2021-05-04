import espressomd
import protein_library as pl
import numpy as np
import matplotlib.pyplot as plt

from espressomd import reaction_ensemble, interactions
from protein_library import protein, molecule

# Minimal example on how to use the protein library
# The program simulates the titration of a peptide in a given pH range
# Simulation done in ideal conditions using a 1 bead per aminoacid model

# Simulation parameters

pH = np.linspace(2, 12, num=20)
Steps_per_sim= 1000
steps_eq=int(Steps_per_sim/3)
SEED=12345
sigma=1

# Write the protein sequence in a list

seq = ["COOH",
        "G",
        "G",
        "G",
        "P",
        "Q",
        "V",
        "T",
        "R",
        "G",
        "D",
        "V",
        "F",
        "T",
        "M",
        "P",
        "NH2"]

# Create a protein object with the input sequence 

peptide = protein(sequence=seq)

# The protein object stores information of the protein and contains dictionaries with specific information of each aminoacid type (pKa-values, charge, radii and default types). By default, types 10-39 are reserved for aminoacid groups

peptide.N = 2 # Number of protein chains in the system
peptide.Nm = len(peptide.sequence) # Number of groups in the protein chain

print("Input protein sequence: ") 
print(peptide.sequence)

# Internally, the protein library works with the 3-letter notation for aminoacid groups in capital letters.
# However, the one-letter notation can be used as input if the sequence is checked.

pl.check_protein_sequence(peptide) # It also corrects capital letter missprints


print()
print("Corrected protein sequence: ") 
print(peptide.sequence)

# Count the number of titrable groups in the peptide

N_ti=pl.count_titrable_groups(peptide)
print()
print("Number of titrable groups in the sequence: ", N_ti)

# Create counter-ion objects 
# Remember that types 10-39 are reserved for protein groups and should be avoided 

cation=molecule(type=0, q=1, radi=sigma/2)
anion=molecule(type=1, q=-1, radi=sigma/2)

# Set-up the reactions of the acid/base groups

RE = reaction_ensemble.ConstantpHEnsemble(temperature=1, exclusion_radius=1, seed=SEED)
pl.setup_protein_acidbase_reactions(RE, peptide, cation)

# Set-up the harmonic bond for the peptide

harmonic_bond = interactions.HarmonicBond(k=100, r_0=sigma)

# Create the system

system = espressomd.System(box_l=[len(peptide.sequence)*3] * 3)

# Add potentials to the system

system.bonded_inter.add(harmonic_bond)

# Create the peptide

peptide.ids=pl.create_protein(system, peptide, harmonic_bond)

# Create the counter-ions

cation.ids, anion.ids = pl.create_counterions(system,cation,anion)

# Set-up espresso to track the ionizable groups

pl.track_ionization(system, peptide)

# Average charge list

Z_pH=[]

for pH_value in pH:

    Z_sim=[]
    RE.constant_pH = pH_value

    for step in range(Steps_per_sim+steps_eq):
    
        RE.reaction(N_ti)

        if ( step > steps_eq):

            Z, Z2=pl.calculate_protein_charge(system,peptide)
            Z_sim.append(Z)
    Z_sim=np.array(Z_sim)
    Z_pH.append(Z_sim.mean())
    print("pH = {:6.4g} done".format(pH_value))

# Calculate the Henderson-Hasselbach prediction in the given pH-range

Z_HH = pl.calculate_HH(pH, peptide)


fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(pH, Z_pH, "ro", label='Simulation-ideal')
ax.plot(pH, Z_HH, "-k", label='Henderson-Hasselbach')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Charge of the peptide / e')
plt.show()

