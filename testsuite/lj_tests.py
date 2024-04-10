# Import pyMBE and other libraries
import pyMBE
import numpy as np

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library()

print("*** LJ unit tests ***")
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
print(f"*** Unit test: check that `offset` defaults to 0***")
# Clean pmb.df
pmb.setup_df()
# Define dummy particle
pmb.define_particle(name="A")
np.testing.assert_equal(actual=pmb.df["offset"].values[0], 
                        desired=pmb.units.Quantity(0,"reduced_length"), 
                        verbose=True)
print("*** Unit test passed ***")

print(f"*** Unit test: check that `cutoff` defaults to `2**(1./6.) reduced_length` ***")
# Clean pmb.df
pmb.setup_df()
# Define dummy particle
pmb.define_particle(name="A")
np.testing.assert_equal(actual=pmb.df["cutoff"].values[0], 
                        desired=pmb.units.Quantity(2**(1./6.),"reduced_length"), 
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

print(f"*** Unit test: test that setup_lj_interactions sets up inert particles correctly ***")

# Clean pmb.df
pmb.setup_df()
# Define particles
A_input_parameters={"name":"A", 
                    "sigma":1*pmb.units.nm, 
                    "epsilon":pmb.units.Quantity(1,"reduced_energy"), 
                    "cutoff":2**(1./6.)*pmb.units.nm, 
                    "offset":1*pmb.units.nm}

B_input_parameters={"name":"B", 
                    "sigma":2*pmb.units.nm, 
                    "epsilon":pmb.units.Quantity(2,"reduced_energy"), 
                    "cutoff":2*2**(1./6.)*pmb.units.nm, 
                    "offset":2*pmb.units.nm,
                    "acidity": "acidic",
                    "pka": 3}
C_input_parameters={"name":"C", 
                "sigma":0*pmb.units.nm, 
                "epsilon":pmb.units.Quantity(2,"reduced_energy"), 
                "cutoff":2*2**(1./6.)*pmb.units.nm, 
                "offset":2*pmb.units.nm}
pmb.define_particle(**A_input_parameters)
pmb.define_particle(**B_input_parameters)
pmb.define_particle(**C_input_parameters)

# Create a dummy instance of an espresso system
import espressomd
espresso_system=espressomd.System(box_l = [50]*3)

# Setup LJ interactions
pmb.setup_lj_interactions(espresso_system=espresso_system)

# Check A-A LJ setup
setup_AA_lj_parameters=pmb.df[pmb.df['name']=="LJ: A-A"].parameters_of_the_potential.values[0]

for parameter_key in ["sigma","offset","cutoff"]:
    np.testing.assert_equal(actual=setup_AA_lj_parameters[parameter_key], 
                            desired=A_input_parameters[parameter_key].to("reduced_length").magnitude, 
                            verbose=True)
np.testing.assert_equal(actual=setup_AA_lj_parameters["epsilon"], 
                            desired=A_input_parameters["epsilon"].to("reduced_energy").magnitude, 
                            verbose=True)

print(f"*** Unit test passed ***")
print(f"*** Unit test: test that setup_lj_interactions sets up acid/base particles correctly ***")


# Check B-B, B-BH, BH-BH setup
labels=["B-B", "BH-B", "BH-BH"]

for label in labels:
    setup_lj_parameters=pmb.df[pmb.df['name']==f"LJ: {label}"].parameters_of_the_potential.values[0]
    for parameter_key in ["sigma","offset","cutoff"]:
        np.testing.assert_equal(actual=setup_lj_parameters[parameter_key], 
                                desired=B_input_parameters[parameter_key].to("reduced_length").magnitude, 
                                verbose=True)
    np.testing.assert_equal(actual=setup_lj_parameters["epsilon"], 
                                desired=B_input_parameters["epsilon"].to("reduced_energy").magnitude, 
                                verbose=True)

print(f"*** Unit test passed ***")
print(f"*** Unit test: test that setup_lj_interactions sets up LJ interaction between different particles correctly ***")


# Calculate the reference parameters
# Assuming Lorentz-Berthelot combining rule
# Check A-BH, A-B, setup
labels=["A-BH", "A-B"]

ref_lj_parameters={}
for parameter_key in ["sigma","offset","cutoff"]:
    ref_lj_parameters[parameter_key]=(A_input_parameters[parameter_key]+B_input_parameters[parameter_key])/2
ref_lj_parameters["epsilon"]=np.sqrt(A_input_parameters["epsilon"]*B_input_parameters["epsilon"])

# Check the parameters set up by pyMBE against the reference parameters
for label in labels:
    setup_lj_parameters=pmb.df[pmb.df['name']==f"LJ: {label}"].parameters_of_the_potential.values[0]
    for parameter_key in ["sigma","offset","cutoff"]:
        np.testing.assert_equal(actual=setup_lj_parameters[parameter_key], 
                                desired=ref_lj_parameters[parameter_key].to("reduced_length").magnitude, 
                                verbose=True)
    np.testing.assert_equal(actual=setup_lj_parameters["epsilon"], 
                                desired=ref_lj_parameters["epsilon"].to("reduced_energy").magnitude, 
                                verbose=True)
print(f"*** Unit test passed ***")

print(f"*** Unit test: test that setup_lj_interactions does not set up any LJ interactions for particles with sigma = 0 ***")

lj_labels=pmb.filter_df("LennardJones")["name"].values
# Check that no interaction between particle C and any other particle has been set up
# Particle C has sigma = 0 (ideally behaving particle)

for label in lj_labels:
    if "C" in label:
        raise Exception("*** Unit Test failed ***")

print(f"*** Unit test passed ***")
print(f"*** All unit tests passed ***")
