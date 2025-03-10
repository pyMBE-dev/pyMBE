#
# Copyright (C) 2025 pyMBE-dev team
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

# Tests that the utility functions in lib/handy_functions work properly
# Import pyMBE and other libraries

import unittest as ut
import espressomd
import pyMBE
import lib.handy_functions as hf
import logging
import io

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)] )
# Create instances of espresso and pyMBE
espresso_system = espressomd.System(box_l=[10,10,10])
seed = 23
pmb = pyMBE.pymbe_library(seed=seed)
kT = pmb.kT

langevin_inputs={"espresso_system":espresso_system, 
                "kT" : kT, 
                "seed" : seed,
                "time_step" : 1e-2, 
                "gamma": 1, 
                "tune_skin" : False, 
                "min_skin": 1,
                "max_skin": None,
                "tolerance": 1e-3,
                "int_steps": 200,
                "adjust_max_skin": True}

class Test(ut.TestCase):
    def test_exceptions_langevin_setup(self):
        print("\n*** Testing exceptions in lib.handy_functions.setup_langevin_dynamics ***")
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["seed"] = "pyMBE"
        self.assertRaises(TypeError, hf.setup_langevin_dynamics, **broken_inputs)
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["time_step"] = -1
        self.assertRaises(ValueError, hf.setup_langevin_dynamics, **broken_inputs)
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["gamma"] = -1
        self.assertRaises(ValueError, hf.setup_langevin_dynamics, **broken_inputs)
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["min_skin"] = 10
        broken_inputs["max_skin"] = 1
        self.assertRaises(ValueError, hf.setup_langevin_dynamics, **broken_inputs)
        print("*** Unit test passed ***")
    def test_langevin_setup(self):
        import re
        print("\n*** Testing setup in lib.handy_functions.setup_langevin_dynamics ***")
        hf.setup_langevin_dynamics(**langevin_inputs)
        ## Test setup of the integrator
        self.assertEqual(first=langevin_inputs["time_step"],
                         second=espresso_system.time_step,
                         msg="The input time step in `lib.handy_functions.setup_langevin_dynamics` is not consisted with the one in the espresso simulation System")
        ## Test setup of the thermostat
        thermostat_setup=espresso_system.thermostat.get_state()[0]
        self.assertEqual(first="LANGEVIN",
                         second=thermostat_setup["type"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different thermostat than Langevin")
        self.assertEqual(first=langevin_inputs["kT"].m_as("reduced_energy"),
                         second=thermostat_setup["kT"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different kT than the input one")
        self.assertEqual(first=langevin_inputs["seed"],
                         second=thermostat_setup["seed"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different seed than the input one")
        self.assertEqual(first=[langevin_inputs["gamma"]]*3,
                         second=thermostat_setup["gamma"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different seed than the input one")
        print("*** Unit test passed ***")
        print("\n*** Testing optimization of skin in lib.handy_functions.setup_langevin_dynamics ***")
        espresso_system.thermostat.turn_off()
        espresso_system.part.add(pos=[1,1,1])
        espresso_system.part.add(pos=[2,2,2])
        espresso_system.non_bonded_inter[0,0].lennard_jones.set_params(epsilon = 1, 
                                                                       sigma = 1, 
                                                                       cutoff = 3,
                                                                       shift = "auto")
        langevin_inputs_opt_skin  = langevin_inputs.copy() 
        langevin_inputs_opt_skin["tune_skin"] = True
        hf.setup_langevin_dynamics(**langevin_inputs_opt_skin)
        log_contents = log_stream.getvalue()
        optimized_skin = float(re.search(r"Optimized skin value: ([\d.]+)", log_contents).group(1))
        self.assertEqual(first=optimized_skin,
                         second=espresso_system.cell_system.skin,
                         msg="The optimized skin has not been set in espresso")
        print("*** Unit test passed ***")
if __name__ == "__main__":
    ut.main()