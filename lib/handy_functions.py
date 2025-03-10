#
# Copyright (C) 2024-2025 pyMBE-dev team
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

import logging
# Create a logger for this module
logger = logging.getLogger(__name__)

def setup_electrostatic_interactions (units, espresso_system, kT, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3,verbose=True):
    """
    Sets up electrostatic interactions in espressomd. 

    Args:
        units(`pint.UnitRegistry`): Unit registry.
        espresso_system: instance of espressmd system class
        kT(`float`): Thermal energy.
        c_salt: Added salt concentration. If provided, the program outputs the debye screening length. It is a mandatory parameter for the Debye-Huckel method.
        solvent_permittivity: Solvent relative permitivity, by default chosen per water at 298.15 K
        method: method prefered for computing the electrostatic interactions. Currently only P3M (label = p3m) and Debye-Huckel (label = DH) are implemented
        tune_p3m: If true (default), tunes the p3m parameters to improve efficiency
        accuracy: desired accuracy for electrostatics, by default 1e-3
        verbose (`bool`): switch to activate/deactivate verbose. Defaults to True.
    """
    import espressomd.electrostatics
    import espressomd.version
    import numpy as np
    import scipy.constants

    # Initial checks
    valid_methods_list=['p3m', 'DH']
    if method not in valid_methods_list:
        raise ValueError('provided an unknown label for method, valid values are', valid_methods_list)
    if c_salt is None and method == 'DH':
        raise ValueError('Please provide the added salt concentration c_salt to setup the Debye-Huckel potential')

    e = scipy.constants.e * units.C
    N_A = scipy.constants.N_A / units.mol

    BJERRUM_LENGTH = e.to('reduced_charge')**2 / (4 * units.pi * units.eps0 * solvent_permittivity * kT.to('reduced_energy'))
    if verbose:
        print(f"\n Bjerrum length {BJERRUM_LENGTH.to('nm')} = {BJERRUM_LENGTH.to('reduced_length')}")

    COULOMB_PREFACTOR=BJERRUM_LENGTH.to('reduced_length') * kT.to('reduced_energy') 
    
    if c_salt is not None:
        if c_salt.check('[substance] [length]**-3'):
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*N_A*c_salt)
        elif c_salt.check('[length]**-3'):
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*c_salt)
        else:
            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)
        if verbose:
            print(f"Debye kappa {KAPPA.to('nm')} = {KAPPA.to('reduced_length')}")

    if method == 'p3m':

        coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.magnitude, 
                                                accuracy=accuracy,
                                                verbose=verbose,
                                                tune=tune_p3m)

        if tune_p3m:
            espresso_system.time_step=0.01
            if espressomd.version.friendly() == "4.2":
                espresso_system.actors.add(coulomb)
            else:
                espresso_system.electrostatics.solver = coulomb

            # save the optimal parameters and add them by hand

            p3m_params = coulomb.get_params()
            if espressomd.version.friendly() == "4.2":
                espresso_system.actors.remove(coulomb)
            else:
                espresso_system.electrostatics.solver = None
            coulomb = espressomd.electrostatics.P3M(
                                        prefactor = COULOMB_PREFACTOR.magnitude,
                                        accuracy = accuracy,
                                        mesh = p3m_params['mesh'],
                                        alpha = p3m_params['alpha'] ,
                                        cao = p3m_params['cao'],
                                        r_cut = p3m_params['r_cut'],
                                        tune = False
                                        )

    elif method == 'DH':

        coulomb = espressomd.electrostatics.DH(prefactor = COULOMB_PREFACTOR.magnitude, 
                                            kappa = (1./KAPPA).to('1/ reduced_length').magnitude, 
                                            r_cut = KAPPA.to('reduced_length').magnitude)

    
    if espressomd.version.friendly() == "4.2":
        espresso_system.actors.add(coulomb)
    else:
        espresso_system.electrostatics.solver = coulomb
    if verbose:
        print("\n Electrostatics successfully added to the system \n")

    return

def minimize_espresso_system_energy(espresso_system, skin=1, gamma=1, Nsteps=10000, time_step=1e-5, max_displacement=0.1, verbose=True, reset=True):
    """
    Does a steppest descent minimization to relax the system energy

    Inputs:
    espresso_system: instance of espressmd system class
    skin: skin parameter for verlet list (default 2 reduced length)
    gamma: dammping constant (Default=1 reduced length)
    Nsteps: total number of steps of the minimization (Default=10000)
    time_step: Time step used for the energy minimization (Default=1e-2)
    max_displacement: maximum particle displacement allowed (Default=0.1 reduced length)
    verbose (`bool`): switch to activate/deactivate verbose. Defaults to True.
    """
    if verbose:
        print("\n*** Minimizing system energy... ***\n")
    espresso_system.cell_system.skin = skin
    espresso_system.time_step=time_step
    if verbose:
        print("steepest descent")
    espresso_system.integrator.set_steepest_descent(f_max=1e-3, gamma=gamma, max_displacement=max_displacement)
    espresso_system.integrator.run(int(Nsteps/2))
    if verbose:
        print("velocity verlet")
    espresso_system.integrator.set_vv()  # to switch back to velocity Verlet
    espresso_system.integrator.run(int(Nsteps/2))
    espresso_system.thermostat.turn_off()
    # Reset the time of the system to 0
    if reset:
        espresso_system.time = 0.
    if verbose:
        print("\n Minimization finished \n")
    return

def setup_langevin_dynamics(espresso_system, kT, seed,time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
    """
    Sets up Langevin Dynamics for an ESPResSo simulation system.

    Parameters:
    espresso_system(`espressomd.system.System`): system object of espressomd library.
    kT (`pint.Quantity`): Target temperature in reduced energy units.
    seed (int): Seed for the random number generator for the thermostat.
    time_step (float, optional): Integration time step. Defaults to 1e-2.
    gamma (float, optional): Damping coefficient for the Langevin thermostat. Defaults to 1.
    tune_skin (bool, optional): Whether to optimize the skin parameter. Defaults to True.
    min_skin (float, optional): Minimum skin value for optimization. Defaults to 1.
    max_skin (float, optional): Maximum skin value for optimization. Defaults to None, which is handled by setting its value to box length / 2.
    tolerance (float, optional): Tolerance for skin optimization. Defaults to 1e-3.
    int_steps (int, optional): Number of integration steps for tuning. Defaults to 200.
    adjust_max_skin (bool, optional): Whether to adjust the maximum skin value during tuning. Defaults to True.
    """        
    if not isinstance(seed, int):
        raise TypeError("SEED must be an integer.")
    if not isinstance(time_step, (float, int)) or time_step <= 0:
        raise ValueError("time_step must be a positive number.")
    if not isinstance(gamma, (float, int)) or gamma <= 0:
        raise ValueError("gamma must be a positive number.")
    if max_skin is None:
            max_skin=espresso_system.box_l[0]/2
    if min_skin >= max_skin:
        raise ValueError("min_skin must be smaller than max_skin.")
    espresso_system.time_step=time_step
    espresso_system.integrator.set_vv()
    espresso_system.thermostat.set_langevin(kT= kT.to('reduced_energy').magnitude, 
                                            gamma= gamma, 
                                            seed= seed)
    # Optimize the value of skin
    if tune_skin:
        logging.info("*** Optimizing skin ... ***")
        espresso_system.cell_system.tune_skin(min_skin=min_skin, 
                                              max_skin=max_skin, 
                                              tol=tolerance, 
                                              int_steps=int_steps, 
                                              adjust_max_skin=adjust_max_skin)
        logging.info(f"Optimized skin value: {espresso_system.cell_system.skin}")
    return

def get_number_of_particles(espresso_system, ptype):
    import espressomd.version
    if espressomd.version.friendly() == "4.2":
        args = (ptype,)
        kwargs = {}
    else:
        args = ()
        kwargs = {"type": ptype}
    return espresso_system.number_of_particles(*args, **kwargs)

def do_reaction(algorithm, steps):
    import espressomd.version
    if espressomd.version.friendly() == '4.2':
        algorithm.reaction(reaction_steps=steps)
    else:
        algorithm.reaction(steps=steps)
