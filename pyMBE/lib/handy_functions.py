#
# Copyright (C) 2024-2026 pyMBE-dev team
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
import re
import numpy as np
import scipy

def calculate_initial_bond_length(bond_parameters, bond_type, lj_parameters):
    """
    Calculates the initial bond length that is used when setting up molecules,
    based on the minimum of the sum of bonded and short-range (LJ) interactions.
    
    Args:
        bond_object(`espressomd.interactions.BondedInteractions`): instance of a bond object from espressomd library
        bond_type(`str`): label identifying the used bonded potential
        epsilon(`pint.Quantity`): LJ epsilon of the interaction between the particles
        sigma(`pint.Quantity`): LJ sigma of the interaction between the particles
        cutoff(`pint.Quantity`): cutoff-radius of the LJ interaction 
        offset(`pint.Quantity`): offset of the LJ interaction
    """    
    def truncated_lj_potential(x, epsilon, sigma, cutoff,offset):
        if x>cutoff:
            return 0.0
        else:
            return 4*epsilon*((sigma/(x-offset))**12-(sigma/(x-offset))**6) - 4*epsilon*((sigma/cutoff)**12-(sigma/cutoff)**6)
    epsilon=lj_parameters["epsilon"].m_as("reduced_energy")
    sigma=lj_parameters["sigma"].m_as("reduced_length")
    cutoff=lj_parameters["cutoff"].m_as("reduced_length")
    offset=lj_parameters["offset"].m_as("reduced_length")
    if bond_type == "harmonic":
        r_0 = bond_parameters['r_0'].m_as("reduced_length")
        k = bond_parameters['k'].m_as("reduced_energy/reduced_length**2")
        l0 = scipy.optimize.minimize(lambda x: 0.5*k*(x-r_0)**2 + truncated_lj_potential(x, epsilon, sigma, cutoff, offset), x0=r_0).x
    elif bond_type == "FENE":
        r_0 = bond_parameters['r_0'].m_as("reduced_length")
        k = bond_parameters['k'].m_as("reduced_energy/reduced_length**2")
        d_r_max = bond_parameters['d_r_max'].m_as("reduced_length")
        l0 = scipy.optimize.minimize(lambda x: -0.5*k*(d_r_max**2)*np.log(1-((x-r_0)/d_r_max)**2) + truncated_lj_potential(x, epsilon, sigma, cutoff,offset), x0=1.0).x
    return l0

def check_aminoacid_key(key):
    """
    Checks if `key` corresponds to a valid aminoacid letter code.

    Args:
        key(`str`): key to be checked.

    Returns:
        `bool`: True if `key` is a valid aminoacid letter code, False otherwise.
    """
    valid_AA_keys=['V', #'VAL'
                    'I', #'ILE'
                    'L', #'LEU'
                    'E', #'GLU'
                    'Q', #'GLN'
                    'D', #'ASP'
                    'N', #'ASN'
                    'H', #'HIS'
                    'W', #'TRP'
                    'F', #'PHE'
                    'Y', #'TYR'
                    'R', #'ARG' 
                    'K', #'LYS'
                    'S', #'SER'
                    'T', #'THR'
                    'M', #'MET'
                    'A', #'ALA'
                    'G', #'GLY'
                    'P', #'PRO'
                    'C', #'CYS'
                    "n", # n terminus
                    "c", # c terminus
                    ] 
    if key in valid_AA_keys:
        return True
    else:
        return False

def check_if_metal_ion(key):
    """
    Checks if `key` corresponds to a label of a supported metal ion.

    Args:
        key(`str`): key to be checked

    Returns:
        (`bool`): True if `key`  is a supported metal ion, False otherwise.
    """
    if key in get_metal_ions_charge_number_map().keys():
        return True
    else:
        return False

def define_protein_AA_particles(topology_dict, pmb, pka_set,  lj_setup_mode="wca"):
    """
    Defines particle templates in pyMBE for all unique residue/atom types appearing
    in a protein topology dictionary.

    The Lennard-Jones parameters (σ, ε, offset) are generated according to the
    selected setup mode (currently only the WCA scheme is supported).  
    
    Metal ions are automatically assigned their correct valence charge.

    Args:
        topology_dict (dict):
            Dictionary defining the structure of a protein.
            Keys must be residue/particle identifiers such as `"ALA1"`, `"LYS2"`,
            `"ZN3"`, etc., where the alphabetical prefix encodes the residue/
            particle type.

            Each entry must contain:
                - `"radius"` (float): Effective radius of the bead, used to
                  compute the Lennard-Jones offset.

            Example:
                {
                    "ALA1": {"radius": 0.5, ...},
                    "GLY2": {"radius": 0.4, ...},
                    "ZN3":  {"radius": 0.2, ...},
                }

        pmb (pyMBE.pymbe_library):
            Instance of the pyMBE library.

        Dictionary of the form:
                {"particle_name": {"pka_value": float,
                                   "acidity": "acidic" | "basic"}}

        lj_setup_mode (str, optional):
            Determines how Lennard-Jones parameters are assigned. Defaults to `"wca"`.           

    Raises:
        ValueError:
            If `lj_setup_mode` is not supported.

    Notes:
        - Particle names are extracted by stripping trailing digits
          (e.g., `"ALA1"` → `"ALA"`).
        - For metal ions (identified via `check_if_metal_ion()`), the correct
          ionic charge is retrieved from the metal-ion charge map.
        - The Lennard-Jones offset is computed as:
                offset = 2 * radius - sigma
    """
    valid_lj_setups = ["wca"]
    if lj_setup_mode not in valid_lj_setups:
            raise ValueError('Invalid key for the lj setup, supported setup modes are {valid_lj_setups}')
    if lj_setup_mode == "wca":
        sigma = 1*pmb.units.Quantity("reduced_length")
        epsilon = 1*pmb.units.Quantity("reduced_energy")
    part_dict={}
    metal_ions_charge_number_map=get_metal_ions_charge_number_map()
    defined_particles=[]
    for particle in topology_dict.keys():
        particle_name = re.split(r'\d+', particle)[0] 
        if particle_name not in defined_particles:
            part_dict = {"name" : particle_name}
            if lj_setup_mode == "wca":
                part_dict["sigma"] = sigma
                part_dict["offset"]= topology_dict[particle]['radius']*2-sigma
                part_dict["epsilon"] = epsilon
            if particle_name in pka_set.keys():
                part_dict["acidity"] = pka_set[particle_name]["acidity"]
            else:
                if check_if_metal_ion(key=particle_name):
                    z=metal_ions_charge_number_map[particle_name]
                else:
                    z=0
                part_dict["z"]=z
        if particle_name not in defined_particles:
            pmb.define_particle(**part_dict)
            defined_particles.append(particle_name)

def define_protein_AA_residues(sequence, model, pmb):
    """
    Define residue templates in the pyMBE database for a protein topology dict.

    Args:
        topology_dict (dict):
                Dictionary defining the internal structure of the protein.
                Expected format:
                    {
                        "ResidueName1": {
                            "initial_pos": np.ndarray,
                            "chain_id": int,
                            "radius": float
                        },
                        "ResidueName2": { ... },
                        ...
                    }
                The `"initial_pos"` entry is required and represents the residue’s
                reference coordinates before shifting to the protein's center-of-mass.

        model (str):
            Coarse-grained representation to use. Supported options:
                - `"1beadAA"`
                - `"2beadAA"`

        pmb (pyMBE.pymbe_library):
            Instance of the pyMBE library.
    Return:
        (list of str): List of the defined residue names

    Notes:
        - Supported models:
            - `"1beadAA"`: Each amino acid is represented by a single bead.  
                The central bead is the amino-acid name itself, and no side chains are used.
            - `"2beadAA"`: Each amino acid is represented by two beads, except for terminal or special residues:
                * `"c"`, `"n"`, and `"G"` (glycine) are treated as single-bead residues.
                * All other residues use `"CA"` (central bead) plus one side-chain bead named after the amino acid.

        - Residue names are constructed as `"AA-<residue>"`, e.g., `"AA-A"`, `"AA-L"`.
    """
    residue_list = []
    for item in sequence:
        if model == '1beadAA':
            central_bead = item
            side_chains = []
        elif model == '2beadAA':
            if item in ['c','n', 'G']: 
                central_bead = item
                side_chains = []
            else:
                central_bead = 'CA'              
                side_chains = [item]
        residue_name='AA-'+item
        if residue_name not in residue_list:   
            pmb.define_residue(name = residue_name, 
                                central_bead = central_bead,
                                side_chains = side_chains)              
        residue_list.append(residue_name)
    return residue_list

def define_peptide_AA_residues(sequence,model, pmb):
    """
    Define residue templates in the pyMBE database for a given model.

    Args:
        sequence (list of str):
            Ordered amino-acid sequence of the peptide or protein. Each element must
            be a residue identifier compatible with the selected model.

        model (str):
            Coarse-grained representation to use. Supported options:
                - `"1beadAA"`
                - `"2beadAA"`

        pmb (pyMBE.pymbe_library):
            Instance of the pyMBE library.

    Notes:
        - Supported models:
            - `"1beadAA"`: Each amino acid is represented by a single bead.  
                The central bead is the amino-acid name itself, and no side chains are used.
            - `"2beadAA"`: Each amino acid is represented by two beads, except for terminal or special residues:
                * `"c"`, `"n"`, and `"G"` (glycine) are treated as single-bead residues.
                * All other residues use `"CA"` (central bead) plus one side-chain bead named after the amino acid.

        - Residue names are constructed as `"AA-<residue>"`, e.g., `"AA-A"`, `"AA-L"`.
    """
    defined_residues = []
    for residue_name in sequence:
        if model == '1beadAA':
            central_bead = residue_name
            side_chains = []
        elif model == '2beadAA':
            if residue_name in ['c','n', 'G']: 
                central_bead = residue_name
                side_chains = []
            else:
                central_bead = 'CA'              
                side_chains = [residue_name]
        residue_name='AA-'+residue_name
        if residue_name not in defined_residues:   
            pmb.define_residue(name = residue_name, 
                                central_bead = central_bead,
                                side_chains = side_chains)              
            defined_residues.append(residue_name)

def do_reaction(algorithm, steps):
    """
    Executes reaction steps using an ESPResSo reaction algorithm with
    version-compatible calling semantics.

    This function wraps the `reaction` method of an ESPResSo reaction
    algorithm to account for differences in the method signature between
    ESPResSo versions.

    Args:
        algorithm:
            ESPResSo reaction algorithm object (e.g. constant pH,
            reaction ensemble, or similar).
        steps (int):
            Number of reaction steps to perform.

    Notes:
        - In ESPResSo 4.2, the `reaction` method expects the number of steps
          to be passed as the keyword argument `reaction_steps`.
        - In newer ESPResSo versions, the keyword argument is `steps`.
        - This helper function provides a stable interface across ESPResSo
          versions by dispatching to the appropriate keyword internally.
    """
    import espressomd.version
    if espressomd.version.friendly() == '4.2':
        algorithm.reaction(reaction_steps=steps)
    else:
        algorithm.reaction(steps=steps)

def get_number_of_particles(espresso_system, ptype):
    """
    Returns the number of particles of a given ESPResSo particle type.

    This function provides a compatibility wrapper around
    `espresso_system.number_of_particles`, which has a different calling
    signature depending on the ESPResSo version.

    Args:
        espresso_system (espressomd.system.System):
            ESPResSo system object from which the particle count is queried.
        ptype (int):
            ESPResSo particle type identifier.

    Returns:
        int:
            Number of particles in `espresso_system` with particle type `ptype`.

    Notes:
        - In ESPResSo 4.2, `number_of_particles` expects the particle type
          as a positional argument.
        - In later ESPResSo versions, the particle type must be passed as a
          keyword argument (`type=ptype`).
        - This helper function hides these API differences and provides
          a uniform interface across ESPResSo versions.
    """
    import espressomd.version
    if espressomd.version.friendly() == "4.2":
        args = (ptype,)
        kwargs = {}
    else:
        args = ()
        kwargs = {"type": ptype}
    return espresso_system.number_of_particles(*args, **kwargs)

def get_residues_from_topology_dict(topology_dict, model):
    """
    Groups beads from a topology dictionary into residues and assigns residue names.

    Args:
        topology_dict (dict):
            Dictionary describing the molecular topology, where keys are bead
            identifiers (e.g. "CA12", "SC12") that encode both residue type and
            residue index.
        model (str):
            Protein model identifier. Supported values are:
            - `"1beadAA"`: single-bead-per-amino-acid model.
            - `"2beadAA"`: two-bead-per-amino-acid model, where CA beads are excluded
              from residue name assignment.

    Returns:
        dict:
            Dictionary mapping residue indices (as strings) to residue data:
            {
                resid: {
                    "beads": [bead_id1, bead_id2, ...],
                    "resname": residue_name
                },
                ...
            }

    Notes:
        - Bead identifiers are parsed by separating alphabetic prefixes
          (residue or bead type) from numeric residue indices.
        - For the `"2beadAA"` model, beads named `"CA"` are excluded when
          determining the residue name.
        - Residues that only contain CA beads (i.e., no side-chain beads)
          are assigned the residue name `"G"` (glycine).
        - Residue indices are returned as strings, consistent with the parsed
          bead identifiers.
    """
    if model not in {"1beadAA", "2beadAA"}:
        raise ValueError(f"Unknown protein model '{model}'")    
    if model == "1beadAA":
        excluded_residue_names = []
    elif model == "2beadAA":
        excluded_residue_names = ["CA"]
    # GROUP BEADS BY RESIDUE
    residues = {}
    for bead_id in topology_dict.keys():
        # extract prefix and index number
        prefix = re.split(r'\d+', bead_id)[0]         
        index_match = re.findall(r'\d+', bead_id)
        if not index_match:
            raise ValueError(f"Topology key '{bead_id}' does not contain a residue index.")
        resid = index_match[0]
        if resid not in residues:
            residues[resid] = {"beads": []}
        residues[resid]["beads"].append(bead_id)
        if prefix not in excluded_residue_names:
            residues[resid]["resname"] = prefix
    
    # Assign name to glycine residues (only with CA beads)
    for bead_id in residues:
        if "resname" not in residues[bead_id]:
            residues[bead_id]["resname"] = "G"
    return residues

def get_metal_ions_charge_number_map():
    """
    Gets a map with the charge numbers of all the metal ions supported.

    Returns:
        metal_charge_number_map(dict): Has the structure {"metal_name": metal_charge_number}

    """
    metal_charge_number_map = {"Ca": 2}
    return metal_charge_number_map

def protein_sequence_parser(sequence):
        '''
        Parses `sequence` to the one letter code for amino acids.
        
        Args:
            sequence(`str` or `lst`): Sequence of the amino acid. 

        Returns:
            clean_sequence(`lst`): `sequence` using the one letter code.
        
        Note:
            - Accepted formats for `sequence` are:
                - `lst` with one letter or three letter code of each aminoacid in each element
                - `str` with the sequence using the one letter code
                - `str` with the squence using the three letter code, each aminoacid must be separated by a hyphen "-"
        
        '''
        # Aminoacid key
        keys={"ALA": "A",
                "ARG": "R",
                "ASN": "N",
                "ASP": "D",
                "CYS": "C",
                "GLU": "E",
                "GLN": "Q",
                "GLY": "G",
                "HIS": "H",
                "ILE": "I",
                "LEU": "L",
                "LYS": "K",
                "MET": "M",
                "PHE": "F",
                "PRO": "P",
                "SER": "S",
                "THR": "T",
                "TRP": "W",
                "TYR": "Y",
                "VAL": "V",
                "PSER": "J",
                "PTHR": "U",
                "PTyr": "Z",
                "NH2": "n",
                "COOH": "c"}
        clean_sequence=[]
        if isinstance(sequence, str):
            if sequence.find("-") != -1:
                splited_sequence=sequence.split("-")
                for residue in splited_sequence:
                    if len(residue) == 1:
                        if residue in keys.values():
                            residue_ok=residue
                        else:
                            if residue.upper() in keys.values():
                                residue_ok=residue.upper()
                            else:
                                raise ValueError("Unknown one letter code for a residue given: ", residue, " please review the input sequence")
                        clean_sequence.append(residue_ok)
                    else:
                        if residue in keys.keys():
                            clean_sequence.append(keys[residue])
                        else:
                            if residue.upper() in keys.keys():
                                clean_sequence.append(keys[residue.upper()])
                            else:
                                raise ValueError("Unknown  code for a residue: ", residue, " please review the input sequence")
            else:
                for residue in sequence:
                    if residue in keys.values():
                        residue_ok=residue
                    else:
                        if residue.upper() in keys.values():
                            residue_ok=residue.upper()
                        else:
                            raise ValueError("Unknown one letter code for a residue: ", residue, " please review the input sequence")
                    clean_sequence.append(residue_ok)
        if isinstance(sequence, list):
            for residue in sequence:
                if residue in keys.values():
                    residue_ok=residue
                else:
                    if residue.upper() in keys.values():
                        residue_ok=residue.upper()
                    elif (residue.upper() in keys.keys()):
                        residue_ok= keys[residue.upper()]
                    else:
                        raise ValueError("Unknown code for a residue: ", residue, " please review the input sequence")
                clean_sequence.append(residue_ok)
        return clean_sequence


def relax_espresso_system(espresso_system, seed, gamma=1e-3, Nsteps_steepest_descent=5000, max_displacement=0.01, Nsteps_iter_relax=500):
    """
    Relaxes the energy of the given ESPResSo system by performing the following steps:
    (1) Steepest descent energy minimization, to remove large forces and relax the system to a local minimum.
    (2) A Langevin Dynamics run, to further relax the system and ensure that it is in thermal equilibrium.

    This function is useful to avoid code repetition in the sample scripts of pyMBE, but it is by no means general-purpose.
    Similarly, the default parameters are not universal and should be adapted to the specific system at hand.
    In general, system relaxation is a complex procedure and should be adapted for each particular application.
    If you experience crashes or unexpected behavior, please consider using your own relaxation procedure.

    Args:
        espresso_system (`espressomd.system.System`): system object of espressomd library.
        seed (`int`): Seed for the random number generator for the thermostat.
        gamma (`float`, optional): Starting damping constant for Langevin dynamics. Defaults to  1e-3 reduced time**-1.
        Nsteps_steepest_descent (`int`, optional): Total number of steps for steepest descent minimization. Defaults to 5000.
        max_displacement (`float`, optional): Maximum particle displacement allowed during minimization. Defaults to 0.01 reduced length.
        Nsteps_iter_relax (`int`, optional): Number of steps per iteration for Langevin dynamics relaxation. Defaults to 500.

    Return:
        (`float`): minimum distance between particles in the system after the relaxation

    Note:
        The thermostat is turned off by the end of the procedure. 
        Make sure the system is initialized properly before calling this function.
    """
    # Sanity checks
    if gamma <= 0:
        raise ValueError("The damping constant 'gamma' must be positive.")
    if Nsteps_steepest_descent <= 0 or Nsteps_iter_relax <= 0:
        raise ValueError("Step counts must be positive integers.")
    if max_displacement <= 0:
        raise ValueError("'max_displacement' must be positive.")
    logging.debug("*** Relaxing the energy of the system... ***")
    logging.debug("*** Starting steepest descent minimization ***")
    espresso_system.thermostat.turn_off()
    espresso_system.integrator.set_steepest_descent(f_max=0,
                                                    gamma=gamma, 
                                                    max_displacement=max_displacement)
    espresso_system.integrator.run(Nsteps_steepest_descent)
    logging.debug("*** Finished steepest descent minimization ***")
    logging.debug("*** Starting Langevin Dynamics relaxation ***")
    
    espresso_system.integrator.set_vv()
    espresso_system.thermostat.set_langevin(kT=1., gamma=gamma, seed=seed)
    espresso_system.integrator.run(Nsteps_iter_relax)
    espresso_system.thermostat.turn_off()

    logging.debug("*** Finished Langevin Dynamics relaxation ***")
    logging.info(f"*** Minimum particle distance after relaxation: {espresso_system.analysis.min_dist()} ***")
    logging.debug("*** Relaxation finished ***")
    return espresso_system.analysis.min_dist()

def setup_langevin_dynamics(espresso_system, kT, seed,time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
    """
    Sets up Langevin Dynamics for an ESPResSo simulation system.

    Args:
        espresso_system (`espressomd.system.System`): system object of espressomd library.
        kT (`pint.Quantity`): Target temperature in reduced energy units.
        seed (`int`): Seed for the random number generator for the thermostat.
        time_step (`float`, optional): Integration time step. Defaults to 1e-2.
        gamma (`float`, optional): Damping coefficient for the Langevin thermostat. Defaults to 1.
        tune_skin (`bool`, optional): Whether to optimize the skin parameter. Defaults to True.
        min_skin (`float`, optional): Minimum skin value for optimization. Defaults to 1.
        max_skin (`float`, optional): Maximum skin value for optimization. Defaults to None, which is handled by setting its value to box length / 2.
        tolerance (`float`, optional): Tolerance for skin optimization. Defaults to 1e-3.
        int_steps (`int`, optional): Number of integration steps for tuning. Defaults to 200.
        adjust_max_skin (`bool`, optional): Whether to adjust the maximum skin value during tuning. Defaults to True.
    """        
    if not isinstance(seed, int):
        raise TypeError("seed must be an integer.")
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
        logging.debug("*** Optimizing skin ... ***")
        espresso_system.cell_system.tune_skin(min_skin=min_skin, 
                                              max_skin=max_skin, 
                                              tol=tolerance, 
                                              int_steps=int_steps, 
                                              adjust_max_skin=adjust_max_skin)
        logging.info(f"Optimized skin value: {espresso_system.cell_system.skin}")

def setup_electrostatic_interactions(units, espresso_system, kT, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3, params=None, verbose=False):
    """
    Sets up electrostatic interactions in an ESPResSo system.

    Args:
        units(`pint.UnitRegistry`): Unit registry for handling physical units.
        espresso_system(`espressomd.system.System`): system object of espressomd library.
        kT(`pint.Quantity`): Thermal energy.
        c_salt(`pint.Quantity`): Added salt concentration. If provided, the program outputs the debye screening length. It is a mandatory parameter for the Debye-Hückel method.
        solvent_permittivity (`float`): Solvent relative permittivity. Defaults to 78.5, correspoding to its value in water at 298.15 K.
        method(`str`): Method for computing electrostatic interactions. Defaults to "p3m". 
        tune_p3m(`bool`): If True, tunes P3M parameters for efficiency. Defaults to True. 
        accuracy(`float`): Desired accuracy for electrostatics. Defaults to 1e-3.
        params(`dict`): Additional parameters for the electrostatic method. For P3M, it can include 'mesh', 'alpha', 'cao' and `r_cut`. For Debye-Hückel, it can include 'r_cut'.
        verbose(`bool`): If True, enables verbose output for P3M tuning. Defaults to False.

    Note:
        `c_salt` is a mandatory argument for setting up the Debye-Hückel electrostatic potential.
        The calculated Bjerrum length is ouput to the log. If `c_salt` is provided, the calculated Debye screening length is also output to the log.
        Currently, the only supported electrostatic methods are P3M ("p3m") and Debye-Hückel ("dh").
    """
    import espressomd.electrostatics
    import espressomd.version
    import numpy as np
    import scipy.constants
    logging.debug("*** Starting electrostatic interactions setup... ***")
    # Initial sanity checks
    if not hasattr(units, 'Quantity'):
        raise TypeError("Invalid 'units' argument: Expected a pint.UnitRegistry object")
    valid_methods_list=['p3m', 'dh']
    if method not in valid_methods_list:
        raise ValueError('Method not supported, supported methods are', valid_methods_list)
    if c_salt is None and method == 'dh':
        raise ValueError('Please provide the added salt concentration c_salt to setup the Debye-Huckel potential')
    e = scipy.constants.e * units.C
    N_A = scipy.constants.N_A / units.mol
    BJERRUM_LENGTH = e**2 / (4 * units.pi * units.eps0 * solvent_permittivity * kT)
    logging.info(f" Bjerrum length {BJERRUM_LENGTH.to('nm')} = {BJERRUM_LENGTH.to('reduced_length')}")
    COULOMB_PREFACTOR=BJERRUM_LENGTH * kT 
    if c_salt is not None:
        if c_salt.check('[substance] [length]**-3'):
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*N_A*c_salt)
        elif c_salt.check('[length]**-3'):
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*c_salt)
        else:
            raise ValueError('Unknown units for c_salt, supported units for salt concentration are [mol / volume] or [particle / volume]', c_salt)
        
        logging.info(f"Debye kappa {KAPPA.to('nm')} = {KAPPA.to('reduced_length')}")

    if params is None:
        params = {}

    if method == 'p3m':
        logging.debug("*** Setting up Coulomb electrostatics using the P3M method ***")
        coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.m_as("reduced_length * reduced_energy"), 
                                                accuracy=accuracy,
                                                verbose=verbose,
                                                tune=tune_p3m,
                                                **params)

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
            coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.m_as("reduced_length * reduced_energy"),
                                                    accuracy = accuracy,
                                                    mesh = p3m_params['mesh'],
                                                    alpha = p3m_params['alpha'] ,
                                                    cao = p3m_params['cao'],
                                                    r_cut = p3m_params['r_cut'],
                                                    tune = False)

    elif method == 'dh':
        logging.debug("*** Setting up Debye-Hückel electrostatics ***")
        if params:
            r_cut = params['r_cut']
        else:
            r_cut = 3*KAPPA.to('reduced_length').magnitude
            
        coulomb = espressomd.electrostatics.DH(prefactor = COULOMB_PREFACTOR.m_as("reduced_length * reduced_energy"), 
                                               kappa = (1./KAPPA).to('1/ reduced_length').magnitude, 
                                               r_cut = r_cut)
    if espressomd.version.friendly() == "4.2":
        espresso_system.actors.add(coulomb)
    else:
        espresso_system.electrostatics.solver = coulomb
    logging.debug("*** Electrostatics successfully added to the system ***")