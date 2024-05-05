import numpy as np 
import espressomd
import warnings
from espressomd import interactions
import pyMBE

pmb = pyMBE.pymbe_library(SEED=42)


print(f"*** Check that the different pKa sets are correctly formatted ***")

list_of_pka_sets = ["Bienkiewicz1999", 
                    "CRC1991", 
                    "Dobrev2020", 
                    "Hass2015", 
                    "Nozaki1967", 
                    "Platzer2014", 
                    "Thurlkill2006"]
for pka_set in list_of_pka_sets:
    print("Checking", pka_set)
    path_to_pka=pmb.get_resource("parameters/pka_sets/" + pka_set + ".json")
    pmb.load_pka_set(path_to_pka)

print(f"*** Test passed ***")
