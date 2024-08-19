#
# Copyright (C) 2024 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import espressomd

def reaction_method_test_template(parameters):

    # Create an instance of the pyMBE library
    pmb = pyMBE.pymbe_library(seed=42)

    if parameters["method"] in ["cpH", "grxmc", "grxmc_unified"]:
        # Define the acidic particle
        pmb.define_particle(
            name = "A",
            acidity = "acidic",
            pka = parameters["pK_acid"],
            sigma = 1*pmb.units('reduced_length'),
            epsilon = 1*pmb.units('reduced_energy'))

        # Define the basic particle
        pmb.define_particle(
            name = "B",
            acidity = "basic",
            pka = parameters["pK_base"],
            sigma = 1*pmb.units('reduced_length'),
            epsilon = 1*pmb.units('reduced_energy'))


    # Define the ions
    pmb.define_particle(
            name="Na", 
            z=1,
            sigma = 1*pmb.units('reduced_length'),
            epsilon = 1*pmb.units('reduced_energy'))

    pmb.define_particle(
            name="Cl", 
            z=-1,
            sigma = 1*pmb.units('reduced_length'),
            epsilon = 1*pmb.units('reduced_energy'))

    pmb.define_particle(
            name="H", 
            z=1,
            sigma = 1*pmb.units('reduced_length'),
            epsilon = 1*pmb.units('reduced_energy'))

    pmb.define_particle(
            name="OH", 
            z=-1,
            sigma = 1*pmb.units('reduced_length'),
            epsilon = 1*pmb.units('reduced_energy'))
    
    if parameters["method"] == "cpH":
        # Add the reactions using pyMBE
        cpH, labels = pmb.setup_cpH(counter_ion="H", constant_pH=parameters["pH"])

        # Check the number of reactions
        np.testing.assert_equal(len(cpH.reactions), 4)

        # Check the equilibrium constants
        np.testing.assert_allclose(cpH.reactions[0].gamma, 10**(-parameters["pK_acid"]))
        np.testing.assert_allclose(cpH.reactions[1].gamma, 10**parameters["pK_acid"])

        np.testing.assert_allclose(cpH.reactions[2].gamma, 10**(-parameters["pK_base"]))
        np.testing.assert_allclose(cpH.reactions[3].gamma, 10**parameters["pK_base"])


    elif parameters["method"] == "gcmc": 
        # Add the reactions using pyMBE
        gcmc = pmb.setup_gcmc(
                c_salt_res=parameters["c_salt_res"] * pmb.units.mol/ pmb.units.L, 
                salt_cation_name="Na", 
                salt_anion_name="Cl", 
                activity_coefficient=lambda x: 1.0)

        # Check the number of reactions
        np.testing.assert_equal(len(gcmc.reactions), 2)

        # Check the equilibrium constants
        K_NaCl = (parameters["c_salt_res"] * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude**2
        np.testing.assert_allclose(gcmc.reactions[0].gamma, K_NaCl)
        np.testing.assert_allclose(gcmc.reactions[1].gamma, 1/K_NaCl)

    elif parameters["method"] == "grxmc": 
        grxmc, labels, ionic_strength_res = pmb.setup_grxmc_reactions(
                pH_res=parameters["pH_res"], 
                c_salt_res=parameters["c_salt_res"] * pmb.units.mol/ pmb.units.L, 
                proton_name="H", 
                hydroxide_name="OH", 
                salt_cation_name="Na", 
                salt_anion_name="Cl", 
                activity_coefficient=lambda x: 1.0)

        # Check the number of reactions
        np.testing.assert_equal(len(grxmc.reactions), 28)

        # Determine the reservoir concentrations independent from pyMBE
        cH_res = (10**(-parameters["pH_res"]) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        cOH_res = (10**(parameters["pH_res"]-14) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        c_salt_res = (parameters["c_salt_res"] * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        cNa_res = max(c_salt_res, c_salt_res + cOH_res - cH_res) 
        cCl_res = max(c_salt_res, c_salt_res + cH_res - cOH_res) 
        Ka_acid = (10**(-parameters["pK_acid"]) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        Ka_base = (10**(-parameters["pK_base"]) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude

        # Check the equilibrium constants
        np.testing.assert_allclose(grxmc.reactions[0].gamma, cH_res*cOH_res)
        np.testing.assert_allclose(grxmc.reactions[1].gamma, 1/(cH_res*cOH_res))
        np.testing.assert_allclose(grxmc.reactions[2].gamma, cNa_res*cCl_res)
        np.testing.assert_allclose(grxmc.reactions[3].gamma, 1/(cNa_res*cCl_res))
        np.testing.assert_allclose(grxmc.reactions[4].gamma, cNa_res*cOH_res)
        np.testing.assert_allclose(grxmc.reactions[5].gamma, 1/(cNa_res*cOH_res))
        np.testing.assert_allclose(grxmc.reactions[6].gamma, cH_res*cCl_res)
        np.testing.assert_allclose(grxmc.reactions[7].gamma, 1/(cH_res*cCl_res))

        np.testing.assert_allclose(grxmc.reactions[8].gamma, cNa_res/cH_res)
        np.testing.assert_allclose(grxmc.reactions[9].gamma, cH_res/cNa_res)
        np.testing.assert_allclose(grxmc.reactions[10].gamma, cCl_res/cOH_res)
        np.testing.assert_allclose(grxmc.reactions[11].gamma, cOH_res/cCl_res)

        np.testing.assert_allclose(grxmc.reactions[12].gamma, Ka_acid)
        np.testing.assert_allclose(grxmc.reactions[13].gamma, 1/Ka_acid)
        np.testing.assert_allclose(grxmc.reactions[14].gamma, Ka_acid*cNa_res/cH_res)
        np.testing.assert_allclose(grxmc.reactions[15].gamma, cH_res/(cNa_res*Ka_acid))
        np.testing.assert_allclose(grxmc.reactions[16].gamma, Ka_acid/(cH_res*cOH_res))
        np.testing.assert_allclose(grxmc.reactions[17].gamma, cH_res*cOH_res/Ka_acid)
        np.testing.assert_allclose(grxmc.reactions[18].gamma, Ka_acid/(cH_res*cCl_res))
        np.testing.assert_allclose(grxmc.reactions[19].gamma, cH_res*cCl_res/Ka_acid)

        np.testing.assert_allclose(grxmc.reactions[20].gamma, Ka_base)
        np.testing.assert_allclose(grxmc.reactions[21].gamma, 1/Ka_base)
        np.testing.assert_allclose(grxmc.reactions[22].gamma, Ka_base*cNa_res/cH_res)
        np.testing.assert_allclose(grxmc.reactions[23].gamma, cH_res/(cNa_res*Ka_base))
        np.testing.assert_allclose(grxmc.reactions[24].gamma, Ka_base/(cH_res*cOH_res))
        np.testing.assert_allclose(grxmc.reactions[25].gamma, cH_res*cOH_res/Ka_base)
        np.testing.assert_allclose(grxmc.reactions[26].gamma, Ka_base/(cH_res*cCl_res))
        np.testing.assert_allclose(grxmc.reactions[27].gamma, cH_res*cCl_res/Ka_base)

    elif parameters["method"] == "grxmc_unified": 
        grxmc, labels, ionic_strength_res = pmb.setup_grxmc_unified(
                pH_res=parameters["pH_res"], 
                c_salt_res=parameters["c_salt_res"] * pmb.units.mol/ pmb.units.L, 
                cation_name="H", 
                anion_name="OH", 
                activity_coefficient=lambda x: 1.0)

        # Check the number of reactions
        np.testing.assert_equal(len(grxmc.reactions), 10)

        # Determine the reservoir concentrations independent from pyMBE
        cH_res = (10**(-parameters["pH_res"]) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        cOH_res = (10**(parameters["pH_res"]-14) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        c_salt_res = (parameters["c_salt_res"] * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        cNa_res = max(c_salt_res, c_salt_res + cOH_res - cH_res) 
        cCl_res = max(c_salt_res, c_salt_res + cH_res - cOH_res) 
        c_cation_res = cH_res + cNa_res
        c_anion_res = cOH_res + cCl_res
        Ka_acid = (10**(-parameters["pK_acid"]) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude
        Ka_base = (10**(-parameters["pK_base"]) * pmb.units.mol/ pmb.units.L).to('1/(N_A * reduced_length**3)').magnitude

        # Check the equilibrium constants
        np.testing.assert_allclose(grxmc.reactions[0].gamma, c_cation_res*c_anion_res)
        np.testing.assert_allclose(grxmc.reactions[1].gamma, 1/(c_cation_res*c_anion_res))

        np.testing.assert_allclose(grxmc.reactions[2].gamma, Ka_acid*c_cation_res/cH_res)
        np.testing.assert_allclose(grxmc.reactions[3].gamma, cH_res/(Ka_acid*c_cation_res))
        np.testing.assert_allclose(grxmc.reactions[4].gamma, Ka_acid/(cH_res*c_anion_res))
        np.testing.assert_allclose(grxmc.reactions[5].gamma, cH_res*c_anion_res/Ka_acid)

        np.testing.assert_allclose(grxmc.reactions[6].gamma, Ka_base*c_cation_res/cH_res)
        np.testing.assert_allclose(grxmc.reactions[7].gamma, cH_res/(Ka_base*c_cation_res))
        np.testing.assert_allclose(grxmc.reactions[8].gamma, Ka_base/(cH_res*c_anion_res))
        np.testing.assert_allclose(grxmc.reactions[9].gamma, cH_res*c_anion_res/Ka_base)
        

# Set up the espresso system
espresso_system=espressomd.System(box_l = [10.0]*3)

# cpH test
print("*** Unit test: check that reactions are correctly set up in the cpH method. ***")
parameters = {
        "method": "cpH",
        "pK_acid": 4.0,
        "pK_base": 8.0,
        "pH": 7.0
        }
reaction_method_test_template(parameters)
print("*** Unit test passed ***")

# gcmc test
print("*** Unit test: check that reactions are correctly set up in the GCMC method. ***")
parameters = {
        "method": "gcmc",
        "c_salt_res": 1,
        }
reaction_method_test_template(parameters)
print("*** Unit test passed ***")

# grxmc test
print("*** Unit test: check that reactions are correctly set up in the G-RxMC method. ***")
parameters = {
        "method": "grxmc",
        "pK_acid": 4.0,
        "pK_base": 9.0,
        "c_salt_res": 1,
        "pH_res": 5.0
        }
reaction_method_test_template(parameters)
print("*** Unit test passed ***")

# grxmc unified test
print("*** Unit test: check that reactions are correctly set up in the unified G-RxMC method. ***")
parameters = {
        "method": "grxmc_unified",
        "pK_acid": 4.0,
        "pK_base": 9.0,
        "c_salt_res": 1,
        "pH_res": 5.0
        }
reaction_method_test_template(parameters)
print("*** Unit test passed ***")
