# Import pyMBE and other libraries
import pyMBE
import numpy as np

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

print(f"*** Unit test: check that all LJ input parameters in define_particle are correctly stored in pmb.df***")

input_parameters={"name":"A", 
                    "sigma":1*pmb.units.nm, 
                    "epsilon":pmb.units.Quantity(1,"reduced_energy"), 
                    "cutoff":2*pmb.units.nm, 
                    "offset":3*pmb.units.nm}

pmb.define_particle(**input_parameters)

for parameter_key in input_parameters.keys():
    np.testing.assert_equal(actual=pmb.df[parameter_key].values[0], 
                            desired=input_parameters[parameter_key], 
                            verbose=True)
print("*** Unit test passed ***")
print(f"*** Unit test: check that `cutoff` and `offset` default to 0***")
# Clean pmb.df
pmb.setup_df()
# Define dummy particle
pmb.define_particle(name="A")
for parameter_key in ["offset","cutoff"]:
    np.testing.assert_equal(actual=pmb.df[parameter_key].values[0], 
                            desired=pmb.units.Quantity(0,"reduced_length"), 
                            verbose=True)
print("*** Unit test passed ***")

print(f"*** Unit test: check that define_particle raises a ValueError if sigma is provided with the wrong dimensionality ***")
input_parameters={"name":"B", 
                   "sigma":1*pmb.units.ns }
np.testing.assert_raises(ValueError, pmb.define_particle, **input_parameters)
print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_particle raises a ValueError if offset is provided with the wrong dimensionality ***")
input_parameters={"name":"B", 
                   "offset":1*pmb.units.ns }
np.testing.assert_raises(ValueError, pmb.define_particle, **input_parameters)
print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_particle raises a ValueError if cutoff is provided with the wrong dimensionality ***")
input_parameters={"name":"B", 
                   "cutoff":1*pmb.units.ns }
np.testing.assert_raises(ValueError, pmb.define_particle, **input_parameters)
print(f"*** Unit test passed ***")

print(f"*** Unit test: check that define_particle raises a ValueError if epsilon is provided with the wrong dimensionality ***")
input_parameters={"name":"B", 
                   "epsilon":1*pmb.units.ns }
np.testing.assert_raises(ValueError, pmb.define_particle, **input_parameters)
print(f"*** Unit test passed ***")

print(f"*** Unit test: test that setup_lj_interactions uses the parameters provided in define_particle to setup the LJ interactions ***")
# Clean pmb.df
pmb.setup_df()
# Define particles
input_parameters={"name":"A", 
                    "sigma":1*pmb.units.nm, 
                    "epsilon":pmb.units.Quantity(1,"reduced_energy"), 
                    "cutoff":2**(1./6.)*pmb.units.nm, 
                    "offset":1*pmb.units.nm}
pmb.define_particle(**input_parameters)
input_parameters={"name":"B", 
                    "sigma":2*pmb.units.nm, 
                    "epsilon":pmb.units.Quantity(2,"reduced_energy"), 
                    "cutoff":2*2**(1./6.)*pmb.units.nm, 
                    "offset":2*pmb.units.nm}
pmb.define_particle(**input_parameters)
input_parameters={"name":"C", 
                    "sigma":0*pmb.units.nm, 
                    "epsilon":pmb.units.Quantity(2,"reduced_energy"), 
                    "cutoff":2*2**(1./6.)*pmb.units.nm, 
                    "offset":2*pmb.units.nm}
pmb.define_particle(**input_parameters)

# Create a dummy instance of an espresso system
import espressomd
espresso_system=espressomd.System(box_l = [10]*3)

pmb.setup_lj_interactions(espresso_system=espresso_system)

print(pmb.filter_df("LennardJones"))
