# Import pyMBE and other libraries
import numpy as np
import pyMBE

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

def check_acid_base_setup(input_parameters,acidity_setup):
    """
    Checks if pyMBE stores in the pmb.df the input parameters for acid/base particles correctly.

    Args:
        input_parameters(`dict`): dictionary with the input parameters for define_particle.
        acidity_setup(`dict`): dictionary with the expected setup that pyMBE should do in the pmb.df foor acid/base particles.

    """
    pmb.define_particle(**input_parameters)
    if input_parameters["acidity"] == "inert":
        input_parameters.pop("q")
    # Checks that the input parameters are stored properly
    for parameter_key in input_parameters.keys():
        np.testing.assert_equal(actual=pmb.df[parameter_key].values[0], 
                                    desired=input_parameters[parameter_key], 
                                    verbose=True)
    # Checks that the setup of the acid base properties is done correctly
    for state in ["state_one","state_two"]:
        for state_atribute in ["label","charge"]:
            np.testing.assert_equal(actual=pmb.df[state][state_atribute].values[0], 
                                    desired=acidity_setup[state][state_atribute], 
                                    verbose=True)
    # checks that pyMBE assigns different espresso type to each state
    np.testing.assert_raises(AssertionError, np.testing.assert_equal, pmb.df["state_one"]["es_type"].values[0], pmb.df["state_two"]["es_type"].values[0])


print("*** Particle acidity unit tests ***")
print("*** Unit test: check that all acid/base input parameters in define_particle for an inert particle are correctly stored in pmb.df***")
# Clean pmb.df
pmb.setup_df()
input_parameters={"name":"I", 
                  "acidity": "inert",
                  "pka": np.nan,
                  "q":2}
acidity_setup={"state_one":{"label":f"{input_parameters['name']}",
                         "charge":2},
            "state_two":{"label": np.nan,
                         "charge":np.nan},}

check_acid_base_setup(input_parameters=input_parameters,
                      acidity_setup=acidity_setup)

print("*** Unit test passed ***")
print("*** Unit test: check that all acid/base input parameters in define_particle for an acid are correctly stored in pmb.df***")
# Clean pmb.df
pmb.setup_df()
input_parameters={"name":"A", 
                  "acidity": "acidic",
                  "pka":4}
acidity_setup={"state_one":{"label":f"{input_parameters['name']}H",
                         "charge":0},
            "state_two":{"label":f"{input_parameters['name']}",
                         "charge":-1},}

check_acid_base_setup(input_parameters=input_parameters,
                      acidity_setup=acidity_setup)
print("*** Unit test passed ***")
print("*** Unit test: check that all acid/base input parameters in define_particle for a base are correctly stored in pmb.df***")
# Clean pmb.df
pmb.setup_df()
input_parameters={"name":"B", 
                  "acidity": "basic",
                  "pka":9}
acidity_setup={"state_one":{"label":f"{input_parameters['name']}H",
                         "charge":1},
            "state_two":{"label":f"{input_parameters['name']}",
                         "charge":0},}

check_acid_base_setup(input_parameters=input_parameters,
                      acidity_setup=acidity_setup)
print("*** Unit test passed ***")
print("*** Unit test: check that set_particle_acidity raises a ValueError if pKa is not provided and pKa is acidic or basic  ***")
input_parametersA={"name":"A", 
                   "acidity": "acidic" }

input_parametersB= {"name": "B",
                   "acidity": "basic"}
np.testing.assert_raises(ValueError, pmb.set_particle_acidity,**input_parametersA)
np.testing.assert_raises(ValueError, pmb.set_particle_acidity, **input_parametersB)
print("*** Unit test passed ***")
print("*** Unit test: check that set_particle_acidity raises a ValueError if a non-supported acidity is provided  ***")
input_parametersA={"name":"A", 
                   "acidity": "random" }
np.testing.assert_raises(ValueError, pmb.set_particle_acidity,**input_parametersA)
print("*** Unit test passed ***")
print("*** All unit tests passed ***")
