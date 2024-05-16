
def setup_electrostatic_interactions (units, espresso_system, kT, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3,verbose=True):
    """
    Sets up electrostatic interactions in espressomd. 
    Inputs:
    system: instance of espressmd system class
    c_salt: Added salt concentration. If provided, the program outputs the debye screening length. It is a mandatory parameter for the Debye-Huckel method. 
    solvent_permittivity: Solvent relative permitivity, by default chosen per water at 298.15 K
    method: method prefered for computing the electrostatic interactions. Currently only P3M (label = p3m) and Debye-Huckel (label = DH) are implemented
    tune_p3m: If true (default), tunes the p3m parameters to improve efficiency
    accuracy: desired accuracy for electrostatics, by default 1e-3
    verbose (`bool`): switch to activate/deactivate verbose. Defaults to True.
    """
    import espressomd.electrostatics
    import numpy as np
    #Initial checks

    valid_methods_list=['p3m', 'DH']

    if method not in valid_methods_list:

        raise ValueError('provided an unknown label for method, valid values are', valid_methods_list)

    if c_salt is None and method == 'DH':

        raise ValueError('Please provide the added salt concentration c_salt to settup the Debye-Huckel potential')
        

    e = 1.60217662e-19 *units.C
    N_A=6.02214076e23    / units.mol

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
    print()

    if method == 'p3m':

        coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.magnitude, 
                                                accuracy=accuracy,
                                                verbose=verbose,
                                                tune=tune_p3m)

        if tune_p3m:
            espresso_system.time_step=0.01
            espresso_system.actors.add(coulomb)

            # save the optimal parameters and add them by hand

            p3m_params = coulomb.get_params()
            espresso_system.actors.remove(coulomb)
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

    
    espresso_system.actors.add(coulomb)
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

def setup_langevin_dynamics(espresso_system, kT, SEED,time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
    """
    Sets up Langevin Dynamics in espressomd.
    espresso_system: instance of espressmd system class
    time_step: time s

    """        
        
    espresso_system.time_step=time_step
    espresso_system.integrator.set_vv()
    espresso_system.thermostat.set_langevin(kT=kT.to('reduced_energy').magnitude, gamma=gamma, seed= SEED)

    # Optimize the value of skin

    if tune_skin:

        print("\n*** Optimizing skin ... ***")

        if max_skin is None:

            max_skin=espresso_system.box_l[0]/2

        espresso_system.cell_system.tune_skin(min_skin=min_skin, max_skin=max_skin, tol=tolerance, int_steps=int_steps, adjust_max_skin=adjust_max_skin)

        print("Optimized skin value: ", espresso_system.cell_system.skin, "\n")

    return

def create_random_seed():
    """
    Generates a seed for the random number generator using the system time in seconds.
    """
    import time 
    SEED=int(time.time())
    print('\n The chosen seed for the random number generator is ', SEED)
    return SEED

def visualize_espresso_system(espresso_system):
    """
    Uses espresso visualizator for displaying the current state of the espresso_system
    """ 
    
    import threading
    from espressomd import visualization
        
    visualizer = visualization.openGLLive(espresso_system)
    
    def main_thread():
        while True:
            espresso_system.integrator.run(1)
            visualizer.update()

    t = threading.Thread(target=main_thread)
    t.daemon = True
    t.start()
    visualizer.start()
    return

def do_snapshot_espresso_system(espresso_system, filename):
    """
    Uses espresso visualizator for creating a snapshot of the current state of the espresso_system
    """ 
    
    from espressomd import visualization
    

    visualizer = visualization.openGLLive(
            espresso_system, bond_type_radius=[0.3], particle_coloring='type', draw_axis=False, background_color=[1, 1, 1],
    particle_type_colors=[[1.02,0.51,0], # Brown
                        [1,1,1],  # Grey
                        [2.55,0,0], # Red
                        [0,0,2.05],  # Blue
                        [0,0,2.05],  # Blue
                        [2.55,0,0], # Red
                        [2.05,1.02,0]]) # Orange
    visualizer.screenshot(filename)

    return
