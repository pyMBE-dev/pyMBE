"""
NOTE: Many of this functions rely on a depracted version of the sugar library and need to be fixed
"""

def setup_electrostatic_interactions_in_espresso(units, espresso_system, kT, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3):
    """
    Setups electrostatic interactions in espressomd. 
    Inputs:
    system: instance of espressmd system class
    c_salt: Added salt concentration. If provided, the program outputs the debye screening length. It is a mandatory parameter for the Debye-Huckel method. 
    solvent_permittivity: Solvent relative permitivity, by default chosen per water at 298.15 K
    method: method prefered for computing the electrostatic interactions. Currently only P3M (label = p3m) and Debye-Huckel (label = DH) are implemented
    tune_p3m: If true (default), tunes the p3m parameters to improve efficiency
    accuracy: desired accuracy for electrostatics, by default 1e-3
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

    print('\n Bjerrum length ', BJERRUM_LENGTH.to('nm'), '=', BJERRUM_LENGTH.to('reduced_length'))

    COULOMB_PREFACTOR=BJERRUM_LENGTH.to('reduced_length') * kT.to('reduced_energy') 
    
    if c_salt is not None:

        if c_salt.check('[substance] [length]**-3'):

            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*N_A*c_salt)

        elif c_salt.check('[length]**-3'):
            
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*c_salt)

        else:

            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)


        print('Debye kappa ', KAPPA.to('nm'), '=', KAPPA.to('reduced_length'), )
    print()

    if method == 'p3m':

        coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.magnitude, accuracy=accuracy)

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
    """

    if verbose:

        print("\n*** Minimazing system energy... ***\n")
    
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

def setup_langevin_dynamics_in_espresso(espresso_system, kT, SEED,time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
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
    Creates the seed for the random number generator from the system hour
    """
    import time 
    SEED=int(time.time())
    print('\n The chosen seed for the random number generator is ', SEED)
    return SEED

def block_analyze(input_data, n_blocks=16):
    '''         
    Performs a binning analysis of input_data. 
    Divides the samples in ``n_blocks`` equispaced blocks
    and returns the mean, its uncertainty, the correlation time 
    and the block size        
    '''
    # NOTE: Depracted function, check soft_matter_wiki
    import numpy as np
    data = np.asarray(input_data)
    block = 0
    # this number of blocks is recommended by Janke as a reasonable compromise
    # between the conflicting requirements on block size and number of blocks
    block_size = int(data.shape[1] / n_blocks)
    print(f"block_size: {block_size}")
    # initialize the array of per-block averages
    block_average = np.zeros((n_blocks, data.shape[0]))
    # calculate averages per each block
    for block in range(n_blocks):
        block_average[block] = np.average(data[:, block * block_size: (block + 1) * block_size], axis=1)
    # calculate the average and average of the square
    av_data = np.average(data, axis=1)
    av2_data = np.average(data * data, axis=1)
    # calculate the variance of the block averages
    block_var = np.var(block_average, axis=0)
    # calculate standard error of the mean
    err_data = np.sqrt(block_var / (n_blocks - 1))
    # estimate autocorrelation time using the formula given by Janke
    # this assumes that the errors have been correctly estimated
    tau_data = np.zeros(av_data.shape)
    for val in range(av_data.shape[0]):
        if av_data[val] == 0:
            # unphysical value marks a failure to compute tau
            tau_data[val] = -1.0
        else:
            tau_data[val] = 0.5 * block_size * n_blocks / (n_blocks - 1) * block_var[val] \
                / (av2_data[val] - av_data[val] * av_data[val])

    # check if the blocks contain enough data for reliable error estimates
    print("uncorrelated samples per block:\nblock_size/tau = ",
        block_size/tau_data)
    threshold = 10.  # block size should be much greater than the correlation time
    if np.any(block_size / tau_data < threshold):
        print("\nWarning: some blocks may contain less than ", threshold, "uncorrelated samples."
        "\nYour error estimated may be unreliable."
        "\nPlease, check them using a more sophisticated method or run a longer simulation.")
        print("? block_size/tau > threshold ? :", block_size/tau_data > threshold)
    else:
        print("\nAll blocks seem to contain more than ", threshold, "uncorrelated samples.\
        Error estimates should be OK.")

    return av_data, err_data, tau_data, block_size

def write_progress(step, total_steps, initial_simulation_time, units):
    """
        Writes the progress of the loop and estimates the time for its completion
        
        Args:
            step(int): current step in the production loop
            total_steps(int): total number of steps in the production loop
            initial_simulation_time(float): computer time when the simulation started, in seconds
            units(obj): UnitRegistry object from the pint library

        Note:
            It assumes that the simulation starts with step = 0
    """
    import time 
    
    time_act=time.time()*units.s
    perc_sim=100 *(step+1) / (total_steps)
    time_per_step= (time_act - initial_simulation_time)/(step+1)
    remaining_time=(total_steps - step +1)*time_per_step
    elapsed_time=time_act-initial_simulation_time

    def find_right_time_units(time):
        """
        Given a pint variable with time units, it returns in which time scale it is
        """

        if (time.to('s').magnitude/60 < 1):

            time_unit='s'

        elif (time.to('s').magnitude/3600 < 1):

            time_unit='min'

        elif (time.to('s').magnitude/(3600*24) < 1):

            time_unit='hour'

        else:

            time_unit='day'

        return time_unit

    time_unit_elapsed_time=find_right_time_units(elapsed_time)
    time_unit_remaining_time=find_right_time_units(remaining_time)

    print("{0:.2g}% done, elapsed time {1:.2g}s; estimated completion in {2:.2g}s".format(perc_sim,elapsed_time.to(time_unit_elapsed_time),remaining_time.to(time_unit_remaining_time)))

    return

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