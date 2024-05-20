import numpy as np 
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(SEED=42)

print("*** Unit test: Check  that check_aminoacid_key returns True for any latter valid in the one letter amino acid code***")
valid_AA_keys=['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
for key in valid_AA_keys:
    np.testing.assert_equal(actual=pmb.check_aminoacid_key(key=key), 
                        desired=True, 
                        verbose=True)
print("*** Unit test passed ***\n")
print("*** Unit test: Check  that check_aminoacid_key returns False for a key not valid in the one letter amino acid code ***")
np.testing.assert_equal(actual=pmb.check_aminoacid_key(key="B"), 
                        desired=False, 
                        verbose=True)
print("*** Unit test passed ***\n")

print("*** Unit test: Check  that check_if_metal_ion returns True for any key corresponding to a supported metal ion ***")
for key in pmb.get_metal_ions_charge_map().keys():
    np.testing.assert_equal(actual=pmb.check_if_metal_ion(key=key), 
                        desired=True, 
                        verbose=True)
print("*** Unit test passed ***\n")
print("*** Unit test: Check  that check_if_metal_ion returns False for a key not corresponding to a supported metal ion ***")
np.testing.assert_equal(actual=pmb.check_if_metal_ion(key="B"), 
                        desired=False, 
                        verbose=True)
print("*** Unit test passed ***\n")

print("*** Unit test: Check  that get_metal_ions_charge_map returns the correct charge map for metals ***")
metal_charge_map = {"Ca": 2}
pmb_metal_charge_map = pmb.get_metal_ions_charge_map()

np.testing.assert_equal(actual=pmb_metal_charge_map, 
                    desired=metal_charge_map, 
                    verbose=True)
print("*** Unit test passed ***\n")

metal_charge_map = {"Ca": 2}
