# Import pyMBE and other libraries
import numpy as np
import pyMBE

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

print("*** Particle acidity unit tests ***")
print(f"*** Unit test: check that set_particle_acidity raises a ValueError if pKa is not provided and pKa is acidic or basic  ***")
input_parametersA={"name":"A", 
                   "acidity": "acidic" }

input_parametersB= {"name": "B",
                   "acidity": "basic"}
np.testing.assert_raises(ValueError, pmb.set_particle_acidity,**input_parametersA)
np.testing.assert_raises(ValueError, pmb.set_particle_acidity, **input_parametersB)
print(f"*** Unit test passed ***")
print(f"*** All unit tests passed ***")