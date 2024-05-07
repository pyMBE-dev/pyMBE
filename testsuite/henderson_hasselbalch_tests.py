import numpy as np
import pyMBE
pmb = pyMBE.pymbe_library(SEED=42)

print(f"*** Running HH tests ***\n")

# Peptide parameters
sequence1 = 5 * "D" + 8 * "H"
sequence2 = 3 * "E" + 7 * "R"
pep1_concentration = 1e-2 * pmb.units.mol/pmb.units.L
pep2_concentration = 1e-2 * pmb.units.mol/pmb.units.L
c_salt=5e-3 * pmb.units.mol/ pmb.units.L
model = '1beadAA' 

# Load pKa-values
path_to_pka=pmb.get_resource("parameters/pka_sets/Nozaki1967.json")
pmb.load_pka_set(path_to_pka)

# Define the peptides in the pyMBE data frame
pmb.define_peptide(name = "peptide_1",
        sequence = sequence1,
        model = model)

pmb.define_peptide(name = "peptide_2",
        sequence = sequence2,
        model = model)


print(f"*** Check that Henderson-Hasselbalch equation works correctly ***")

# Calculate charge according to Henderson-Hasselbalch equation
pH_range = np.linspace(2, 12, num=200)
Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1", 
                        pH_list = pH_range)
Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2", 
                        pH_list = pH_range)

"""
with open("henderson_hasselbalch_tests_data/HH.csv", "wb") as f:
    np.savetxt(f, np.asarray(Z_HH_1).reshape(1,-1), delimiter=",")
    np.savetxt(f, np.asarray(Z_HH_2).reshape(1,-1), delimiter=",")
"""

data_path = pmb.get_resource(path="testsuite/henderson_hasselbalch_tests_data")
ref_data_HH = np.loadtxt(data_path+f"/HH.csv", delimiter=",")
np.testing.assert_allclose(Z_HH_1, ref_data_HH[0,:])
np.testing.assert_allclose(Z_HH_2, ref_data_HH[1,:])

print(f"*** Test passed ***\n")


print(f"*** Check that Henderson-Hasselbalch equation + Donnan works correctly ***")

HH_Donnan_dict = pmb.calculate_HH_Donnan(
        c_macro = {"peptide_1": pep1_concentration,
                   "peptide_2": pep2_concentration}, 
        c_salt = c_salt, 
        pH_list = pH_range)

"""
with open("henderson_hasselbalch_tests_data/HH_Donnan.csv", "wb") as f:
    np.savetxt(f, np.asarray(HH_Donnan_dict["charges_dict"]["peptide_1"]).reshape(1,-1), delimiter=",")
    np.savetxt(f, np.asarray(HH_Donnan_dict["charges_dict"]["peptide_2"]).reshape(1,-1), delimiter=",")
"""

ref_data_HH_Donnan = np.loadtxt(data_path+f"/HH_Donnan.csv", delimiter=",")
np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_1"], ref_data_HH_Donnan[0,:])
np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_2"], ref_data_HH_Donnan[1,:])

print(f"*** Test passed ***\n")


print(f"*** Check that HH and HH_Don are consistent ***")

Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1", 
                        pH_list = HH_Donnan_dict["pH_system_list"])
Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2", 
                        pH_list = HH_Donnan_dict["pH_system_list"])

np.testing.assert_allclose(Z_HH_1, HH_Donnan_dict["charges_dict"]["peptide_1"])
np.testing.assert_allclose(Z_HH_2, HH_Donnan_dict["charges_dict"]["peptide_2"])


print(f"*** Test passed***")
