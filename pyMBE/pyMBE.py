#
# Copyright (C) 2023-2025 pyMBE-dev team
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

import re
import json
import pint
import numpy as np
import pandas as pd
import scipy.constants
import scipy.optimize
import logging
import importlib.resources

# Database
from pyMBE.storage.manager import Manager
from pyMBE.storage.pint_quantity import PintQuantity
## Templates
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleState
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.templates.hydrogel import HydrogelTemplate, HydrogelNode, HydrogelChain
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.templates.lj import LJInteractionTemplate
## Instances
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.instances.hydrogel import HydrogelInstance
## Reactions
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
# Utilities
import pyMBE.lib.handy_functions as hf
import pyMBE.storage.io as io

class pymbe_library():
    """
    Core library for the Molecular Builder for ESPResSo (pyMBE).

    Provides access to functions to define templates to build coarse-grained models of macromolecules.

    Attributes:
        N_A (pint.Quantity): Avogadro number.
        kB (pint.Quantity): Boltzmann constant.
        e (pint.Quantity): Elementary charge.
        kT (pint.Quantity): Thermal energy at the set temperature.
        Kw (pint.Quantity): Ionic product of water. Used in G-RxMC method setup.
        db (Manager): Database manager instance for pyMBE objects.
        rng (np.random.Generator): Random number generator initialized with the provided seed.
        units (pint.UnitRegistry): Pint UnitRegistry for unit-aware calculations.
        lattice_builder (optional): Placeholder for lattice builder object, initially None.
        root (importlib.resources.Path): Root path to the pyMBE package resources.
    """

    def __init__(self, seed, temperature=None, unit_length=None, unit_charge=None, Kw=None):
        """
        Initialize the pyMBE library.

        Sets up the reduced unit system using temperature, unit length, and unit charge,
        initializes the pyMBE database, and sets default physical constants.

        Args:
            seed (int): Seed for the random number generator.
            temperature (pint.Quantity, optional): Simulation temperature. Defaults to 298.15 K if None.
            unit_length (pint.Quantity, optional): Reference length for reduced units. Defaults to 0.355 nm if None.
            unit_charge (pint.Quantity, optional): Reference charge for reduced units. Defaults to 1 elementary charge if None.
            Kw (pint.Quantity, optional): Ionic product of water in mol^2/l^2. Defaults to 1e-14 mol^2/l^2 if None.

        Notes:
            - Initializes `self.rng` for random number generation.
            - Sets fundamental constants: Avogadro number (`N_A`), Boltzmann constant (`kB`), elementary charge (`e`).
            - Initializes the reduced units via `set_reduced_units`.
            - Prepares an empty database (`self.db`) for pyMBE objects.
            - Initializes placeholders for `lattice_builder` and package resource path (`root`).
        """
        # Seed and RNG
        self.seed=seed
        self.rng = np.random.default_rng(seed)
        self.units=pint.UnitRegistry()
        self.N_A=scipy.constants.N_A / self.units.mol
        self.kB=scipy.constants.k * self.units.J / self.units.K
        self.e=scipy.constants.e * self.units.C
        self.set_reduced_units(unit_length=unit_length, 
                               unit_charge=unit_charge,
                               temperature=temperature, 
                               Kw=Kw)
        
        self.db = Manager(units=self.units)
        self.lattice_builder = None
        self.root = importlib.resources.files(__package__)
        self._bond_instances={}   

    def _create_espresso_bond_instance(self, bond_type, bond_parameters):
        """
        Creates an ESPResSo bond instance.

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.

        Note:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain:
                - k (`Pint.Quantity`)      : Magnitude of the bond. It should have units of energy/length**2 
                using the `pmb.units` UnitRegistry.
                - r_0 (`Pint.Quantity`)    : Equilibrium bond length. It should have units of length using 
                the `pmb.units` UnitRegistry.
           
            For a FENE bond the dictionary must additionally contain:
                - d_r_max (`Pint.Quantity`): Maximal stretching length for FENE. It should have 
                units of length using the `pmb.units` UnitRegistry. Default 'None'.

        Returns:
            (`espressomd.interactions`): instance of an ESPResSo bond object
        """
        from espressomd import interactions
        valid_bond_types   = ["harmonic", "FENE"] 
        if bond_type not in valid_bond_types:
            raise NotImplementedError(f"Bond type '{bond_type}' currently not implemented in pyMBE, accepted types are {valid_bond_types}")
        required_parameters = {"harmonic": ["r_0","k"],
                                "FENE": ["r_0","k","d_r_max"]}
        for required_parameter in required_parameters[bond_type]:
            if required_parameter not in bond_parameters.keys():
                raise ValueError(f"Missing required parameter {required_parameter} for {bond_type} bond")
        if bond_type == 'harmonic':
            bond_instance = interactions.HarmonicBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                      r_0 = bond_parameters["r_0"].m_as("reduced_length"))
        elif bond_type == 'FENE':
            bond_instance    = interactions.FeneBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                      r_0 = bond_parameters["r_0"].m_as("reduced_length"),
                                                      d_r_max = bond_parameters["d_r_max"].m_as("reduced_length"))    
        return bond_instance

    def _create_hydrogel_chain(self, hydrogel_chain, nodes, espresso_system):
        """
        Creates a chain between two nodes of a hydrogel.

        Args:
            hydrogel_chain(HydrogelChain): template of a hydrogel chain
            nodes(dict): {node_index: {"name": node_particle_name, "pos": node_position, "id": node_particle_instance_id}}
            espresso_system (espressomd.system.System): ESPResSo system object where the hydrogel chain will be created.

        Return:
            (int): molecule_id of the created hydrogel chian.

        Note:
            - For example, if the chain is defined between node_start = ``[0 0 0]`` and node_end = ``[1 1 1]``, the chain will be placed between these two nodes.
            - The chain will be placed in the direction of the vector between `node_start` and `node_end`. 
            - This function does not support default bonds.
        """
        if self.lattice_builder is None:
            raise ValueError("LatticeBuilder is not initialized. Use `initialize_lattice_builder` first.")
        molecule_tpl = self.db.get_template(pmb_type="molecule",
                                            name=hydrogel_chain.molecule_name)
        residue_list = molecule_tpl.residue_list
        molecule_name = molecule_tpl.name
        node_start = hydrogel_chain.node_start
        node_end = hydrogel_chain.node_end
        node_start_label = self.lattice_builder._create_node_label(node_start)
        node_end_label = self.lattice_builder._create_node_label(node_end)
        _, reverse = self.lattice_builder._get_node_vector_pair(node_start, node_end)
        if node_start != node_end or residue_list == residue_list[::-1]:
            RuntimeError(f"Aborted creation because hydrogel chain between '{node_start}' and '{node_end}' because pyMBE could not resolve a unique topology for that chain")
        if reverse:
            residue_list = residue_list[::-1]
        start_node_id = nodes[node_start_label]["id"]
        end_node_id = nodes[node_end_label]["id"]
        # Finding a backbone vector between node_start and node_end
        vec_between_nodes = np.array(nodes[node_start_label]["pos"]) - np.array(nodes[node_end_label]["pos"])
        vec_between_nodes = vec_between_nodes - self.lattice_builder.box_l * np.round(vec_between_nodes/self.lattice_builder.box_l)
        backbone_vector = np.array((vec_between_nodes/(self.lattice_builder.mpc + 1)))
        backbone_vector = backbone_vector / np.linalg.norm(backbone_vector)
        # Calculate the start position of the chain
        chain_residues = self.db.get_template(pmb_type="molecule",
                                              name=molecule_name).residue_list
        part_start_chain_name = self.db.get_template(pmb_type="residue",
                                                     name=chain_residues[0]).central_bead
        part_end_chain_name = self.db.get_template(pmb_type="residue",
                                                   name=chain_residues[-1]).central_bead
        lj_parameters = self.get_lj_parameters(particle_name1=nodes[node_start_label]["name"],
                                               particle_name2=part_start_chain_name)
        bond_tpl = self.get_bond_template(particle_name1=nodes[node_start_label]["name"],
                                        particle_name2=part_start_chain_name)
        l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                              bond_type=bond_tpl.bond_type,
                                              bond_parameters=bond_tpl.get_parameters(ureg=self.units))
        first_bead_pos = np.array((nodes[node_start_label]["pos"])) + np.array(backbone_vector)*l0
        mol_id = self.create_molecule(name=molecule_name,  # Use the name defined earlier
                                      number_of_molecules=1,  # Creating one chain
                                      espresso_system=espresso_system,
                                      list_of_first_residue_positions=[first_bead_pos.tolist()],#Start at the first node
                                      backbone_vector=np.array(backbone_vector)/l0,
                                      use_default_bond=False)
        # Bond chain to the hydrogel nodes
        chain_pids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                             attribute="molecule_id",
                                                             value=mol_id)
        start_bond_instance = self.get_espresso_bond_instance(particle_name1=nodes[node_start_label]["name"],
                                                              particle_name2=part_start_chain_name,
                                                              espresso_system=espresso_system)    
        end_bond_instance = self.get_espresso_bond_instance(particle_name1=nodes[node_end_label]["name"],
                                                              particle_name2=part_end_chain_name,
                                                              espresso_system=espresso_system)
        espresso_system.part.by_id(start_node_id).add_bond((start_bond_instance, chain_pids[0]))
        espresso_system.part.by_id(chain_pids[-1]).add_bond((end_bond_instance, end_node_id))
        return mol_id

    def _create_hydrogel_node(self, node_index, node_name, espresso_system):
        """
        Set a node residue type.
        
        Args:
            node_index(`str`): Lattice node index in the form of a string, e.g. "[0 0 0]".
            node_name(`str`): name of the node particle defined in pyMBE.
            espresso_system (espressomd.system.System): ESPResSo system object where the hydrogel node will be created.

        Returns:
            node_position(`list`): Position of the node in the lattice.
            p_id(`int`): Particle ID of the node.
        """
        if self.lattice_builder is None:
            raise ValueError("LatticeBuilder is not initialized. Use `initialize_lattice_builder` first.")

        node_position = np.array(node_index)*0.25*self.lattice_builder.box_l
        p_id = self.create_particle(name = node_name,
                                    espresso_system=espresso_system,
                                    number_of_particles=1,
                                    position = [node_position])
        key = self.lattice_builder._get_node_by_label(f"[{node_index[0]} {node_index[1]} {node_index[2]}]")
        self.lattice_builder.nodes[key] = node_name

        return node_position.tolist(), p_id[0]

    def _get_residue_list_from_sequence(self, sequence):
        """
        Convinience function to get a `residue_list` from a protein or peptide `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein.

        Returns:
            residue_list(`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.             
        """
        residue_list = []
        for item in sequence:
            residue_name='AA-'+item
            residue_list.append(residue_name)
        return residue_list
    
    def _delete_particles_from_espresso(self, particle_ids, espresso_system):
        """
        Remove a list of particles from an ESPResSo simulation system.

        Args:
        particle_ids : Iterable[int]
            A list (or other iterable) of ESPResSo particle IDs to remove.

        espresso_system : espressomd.system.System
            The ESPResSo simulation system from which the particles
            will be removed.

        Note:
        - This method removes particles only from the ESPResSo simulation,
        **not** from the pyMBE database. Database cleanup must be handled
        separately by the caller.
        - Attempting to remove a non-existent particle ID will raise
        an ESPResSo error.
        """
        for pid in particle_ids:
            espresso_system.part.by_id(pid).remove()


    def calculate_center_of_mass_of_molecule(self, molecule_id, espresso_system):
        """
        Calculates the center of the molecule with a given molecule_id.

        Args:
            molecule_id(`int`): Id of the molecule whose center of mass is to be calculated.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
        
        Returns:
            center_of_mass(`lst`): Coordinates of the center of mass.
        """
        center_of_mass = np.zeros(3)
        axis_list = [0,1,2]
        mol_inst = self.db.get_instance(pmb_type="molecule",
                                        instance_id=molecule_id)
        molecule_name = mol_inst.name
        particle_id_list = self.get_particle_id_map(object_name=molecule_name)["all"]
        for pid in particle_id_list:
            for axis in axis_list:
                center_of_mass [axis] += espresso_system.part.by_id(pid).pos[axis]
        center_of_mass = center_of_mass / len(particle_id_list)
        return center_of_mass

    def calculate_HH(self, molecule_name, pmb_type, pH_list=None, pka_set=None):
        """
        Calculates the charge per molecule according to the ideal Henderson-Hasselbalch titration curve 
        for molecules with the name `molecule_name`.

        Args:
            molecule_name(`str`): name of the molecule to calculate the ideal charge for
            pH_list(`lst`): pH-values to calculate. 
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}

        Returns:
            Z_HH(`lst`): Henderson-Hasselbalch prediction of the charge of `sequence` in `pH_list`

        Note:
            - This function supports objects with pmb types: "molecule", "peptide" and "protein".
            - If no `pH_list` is given, 50 equispaced pH-values ranging from 2 to 12 are calculated
            - If no `pka_set` is given, the pKa values are taken from `pmb.df`
            - This function should only be used for single-phase systems. For two-phase systems `pmb.calculate_HH_Donnan`  should be used.
        """
        _DFm._check_if_name_is_defined_in_df(name = molecule_name,
                                             df = self.df)
        self._check_supported_molecule(molecule_name = molecule_name,
                                       valid_pmb_types = ["molecule","peptide","protein"])
        if pH_list is None:
            pH_list=np.linspace(2,12,50)
        if pka_set is None:
            pka_set=self.get_pka_set()
        index = self.df.loc[self.df['name'] == molecule_name].index[0].item() 
        residue_list = self.df.at [index,('residue_list','')].copy()
        particles_in_molecule = []
        for residue in residue_list:
            list_of_particles_in_residue = self.search_particles_in_residue(residue)
            if len(list_of_particles_in_residue) == 0:
                logging.warning(f"The residue {residue} has no particles defined in the pyMBE DataFrame, it will be ignored.")
                continue
            particles_in_molecule += list_of_particles_in_residue
        if len(particles_in_molecule) == 0:
            return [None]*len(pH_list)
        self.check_pka_set(pka_set=pka_set)
        charge_number_map = self.get_charge_number_map()
        Z_HH=[]
        for pH_value in pH_list:    
            Z=0
            for particle in particles_in_molecule:
                if particle in pka_set.keys():
                    if pka_set[particle]['acidity'] == 'acidic':
                        psi=-1
                    elif pka_set[particle]['acidity']== 'basic':
                        psi=+1
                    Z+=psi/(1+10**(psi*(pH_value-pka_set[particle]['pka_value'])))                      
                else:
                    state_one_type = self.df.loc[self.df['name']==particle].state_one.es_type.values[0]
                    Z+=charge_number_map[state_one_type]
            Z_HH.append(Z)
        return Z_HH

    def calculate_HH_Donnan(self, c_macro, c_salt, pH_list=None, pka_set=None):
        """
        Calculates the charge on the different molecules according to the Henderson-Hasselbalch equation coupled to the Donnan partitioning.

        Args:
            c_macro('dict'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
            c_salt('float'): Salt concentration in the reservoir.
            pH_list('lst'): List of pH-values in the reservoir. 
            pka_set('dict'): {"name": {"pka_value": pka, "acidity": acidity}}.

        Returns:
            {"charges_dict": {"name": charges}, "pH_system_list": pH_system_list, "partition_coefficients": partition_coefficients_list}
            pH_system_list ('lst'): List of pH_values in the system.
            partition_coefficients_list ('lst'): List of partition coefficients of cations.

        Note:
            - If no `pH_list` is given, 50 equispaced pH-values ranging from 2 to 12 are calculated
            - If no `pka_set` is given, the pKa values are taken from `pmb.df`
            - If `c_macro` does not contain all charged molecules in the system, this function is likely to provide the wrong result.
        """
        if pH_list is None:
            pH_list=np.linspace(2,12,50)
        if pka_set is None:
            pka_set=self.get_pka_set() 
        self.check_pka_set(pka_set=pka_set)

        partition_coefficients_list = []
        pH_system_list = []
        Z_HH_Donnan={}
        for key in c_macro:
            Z_HH_Donnan[key] = []

        def calc_charges(c_macro, pH):
            """
            Calculates the charges of the different kinds of molecules according to the Henderson-Hasselbalch equation.

            Args:
                c_macro('dic'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
                pH('float'): pH-value that is used in the HH equation.

            Returns:
                charge('dict'): {"molecule_name": charge}
            """
            charge = {}
            for name in c_macro:
                charge[name] = self.calculate_HH(name, [pH], pka_set)[0]
            return charge

        def calc_partition_coefficient(charge, c_macro):
            """
            Calculates the partition coefficients of positive ions according to the ideal Donnan theory.

            Args:
                charge('dict'): {"molecule_name": charge}
                c_macro('dict'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
            """
            nonlocal ionic_strength_res
            charge_density = 0.0
            for key in charge:
                charge_density += charge[key] * c_macro[key]
            return (-charge_density / (2 * ionic_strength_res) + np.sqrt((charge_density / (2 * ionic_strength_res))**2 + 1)).magnitude

        for pH_value in pH_list:    
            # calculate the ionic strength of the reservoir
            if pH_value <= 7.0:
                ionic_strength_res = 10 ** (-pH_value) * self.units.mol/self.units.l + c_salt 
            elif pH_value > 7.0:
                ionic_strength_res = 10 ** (-(14-pH_value)) * self.units.mol/self.units.l + c_salt

            #Determine the partition coefficient of positive ions by solving the system of nonlinear, coupled equations
            #consisting of the partition coefficient given by the ideal Donnan theory and the Henderson-Hasselbalch equation.
            #The nonlinear equation is formulated for log(xi) since log-operations are not supported for RootResult objects.
            equation = lambda logxi: logxi - np.log10(calc_partition_coefficient(calc_charges(c_macro, pH_value - logxi), c_macro))
            logxi = scipy.optimize.root_scalar(equation, bracket=[-1e2, 1e2], method="brentq")
            partition_coefficient = 10**logxi.root

            charges_temp = calc_charges(c_macro, pH_value-np.log10(partition_coefficient))
            for key in c_macro:
                Z_HH_Donnan[key].append(charges_temp[key])

            pH_system_list.append(pH_value - np.log10(partition_coefficient))
            partition_coefficients_list.append(partition_coefficient)

        return {"charges_dict": Z_HH_Donnan, "pH_system_list": pH_system_list, "partition_coefficients": partition_coefficients_list}

    def calculate_net_charge(self, espresso_system, molecule_name, dimensionless=False):
        '''
        Calculates the net charge per molecule of molecules with `name` = molecule_name. 
        Returns the net charge per molecule and a maps with the net charge per residue and molecule.

        Args:
            espresso_system(`espressomd.system.System`): system information 
            molecule_name(`str`): name of the molecule to calculate the net charge
            dimensionless(`bool'): sets if the charge is returned with a dimension or not

        Returns:
            (`dict`): {"mean": mean_net_charge, "molecules": {mol_id: net_charge_of_mol_id, }, "residues": {res_id: net_charge_of_res_id, }}

        Note:
            - The net charge of the molecule is averaged over all molecules of type `name` 
            - The net charge of each particle type is averaged over all particle of the same type in all molecules of type `name`
        '''        
        self._check_supported_molecule(molecule_name=molecule_name,
                                        valid_pmb_types=["molecule","protein","peptide"])

        id_map = self.get_particle_id_map(object_name=molecule_name)
        def create_charge_map(espresso_system,id_map,label):
            charge_number_map={}
            for super_id in id_map[label].keys():
                if dimensionless:
                    net_charge=0
                else:
                    net_charge=0 * self.units.Quantity(1,'reduced_charge')
                for pid in id_map[label][super_id]:
                    if dimensionless:
                        net_charge+=espresso_system.part.by_id(pid).q
                    else:
                        net_charge+=espresso_system.part.by_id(pid).q * self.units.Quantity(1,'reduced_charge')
                charge_number_map[super_id]=net_charge
            return charge_number_map
        net_charge_molecules=create_charge_map(label="molecule_map",
                                                espresso_system=espresso_system,
                                                id_map=id_map)
        net_charge_residues=create_charge_map(label="residue_map",
                                                espresso_system=espresso_system,
                                                id_map=id_map)       
        if dimensionless:
            mean_charge=np.mean(np.array(list(net_charge_molecules.values())))
        else:
            mean_charge=np.mean(np.array([value.magnitude for value in net_charge_molecules.values()]))*self.units.Quantity(1,'reduced_charge')
        return {"mean": mean_charge, "molecules": net_charge_molecules, "residues": net_charge_residues}

    def center_molecule_in_simulation_box(self, molecule_id, espresso_system):
        """
        Centers the pmb object matching `molecule_id` in the center of the simulation box in `espresso_md`.
        
        Args:
            molecule_id(`int`): Id of the molecule to be centered.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
        """
        if len(self.df.loc[self.df['molecule_id']==molecule_id].pmb_type) == 0:
            raise ValueError("The provided molecule_id is not present in the pyMBE dataframe.")      
        center_of_mass = self.calculate_center_of_mass_of_molecule(molecule_id=molecule_id,espresso_system=espresso_system)
        box_center = [espresso_system.box_l[0]/2.0,
                      espresso_system.box_l[1]/2.0,
                      espresso_system.box_l[2]/2.0]
        molecule_name = self.df.loc[(self.df['molecule_id']==molecule_id) & (self.df['pmb_type'].isin(["molecule","protein"]))].name.values[0]
        particle_id_list = self.get_particle_id_map(object_name=molecule_name)["all"]
        for pid in particle_id_list:
            es_pos = espresso_system.part.by_id(pid).pos
            espresso_system.part.by_id(pid).pos = es_pos - center_of_mass + box_center
        return 

    def check_dimensionality(self, variable, expected_dimensionality):
        """
        Checks if the dimensionality of `variable` matches `expected_dimensionality`.

        Args:
            variable(`pint.Quantity`): Quantity to be checked.
            expected_dimensionality(`str`): Expected dimension of the variable.

        Returns:
            (`bool`): `True` if the variable if of the expected dimensionality, `False` otherwise.

        Note:
            - `expected_dimensionality` takes dimensionality following the Pint standards [docs](https://pint.readthedocs.io/en/0.10.1/wrapping.html?highlight=dimensionality#checking-dimensionality).
            - For example, to check for a variable corresponding to a velocity `expected_dimensionality = "[length]/[time]"`    
        """
        correct_dimensionality=variable.check(f"{expected_dimensionality}")      
        if not correct_dimensionality:
            raise ValueError(f"The variable {variable} should have a dimensionality of {expected_dimensionality}, instead the variable has a dimensionality of {variable.dimensionality}")
        return correct_dimensionality   

    def check_pka_set(self, pka_set):
        """
        Checks that `pka_set` has the formatting expected by the pyMBE library.
       
        Args:
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_name, pka_entry in pka_set.items():
                if required_key not in pka_entry:
                    raise ValueError(f'missing a required key "{required_key}" in entry "{pka_name}" of pka_set ("{pka_entry}")')
        return

    def create_added_salt(self, espresso_system, cation_name, anion_name, c_salt):    
        """
        Creates a `c_salt` concentration of `cation_name` and `anion_name` ions into the `espresso_system`.

        Args:
            espresso_system(`espressomd.system.System`): instance of an espresso system object.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            c_salt(`float`): Salt concentration.
            
        Returns:
            c_salt_calculated(`float`): Calculated salt concentration added to `espresso_system`.
        """ 
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=cation_name)
        cation_charge = cation_tpl.states[cation_tpl.initial_state].z
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=anion_name)
        anion_charge = anion_tpl.states[anion_tpl.initial_state].z
        if cation_charge <= 0:
            raise ValueError(f'ERROR cation charge must be positive, charge {cation_charge}')
        if anion_charge >= 0:
            raise ValueError(f'ERROR anion charge must be negative, charge {anion_charge}')
        # Calculate the number of ions in the simulation box
        volume=self.units.Quantity(espresso_system.volume(), 'reduced_length**3')
        if c_salt.check('[substance] [length]**-3'):
            N_ions= int((volume*c_salt.to('mol/reduced_length**3')*self.N_A).magnitude)
            c_salt_calculated=N_ions/(volume*self.N_A)
        elif c_salt.check('[length]**-3'):
            N_ions= int((volume*c_salt.to('reduced_length**-3')).magnitude)
            c_salt_calculated=N_ions/volume
        else:
            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)
        N_cation = N_ions*abs(anion_charge)
        N_anion = N_ions*abs(cation_charge)
        self.create_particle(espresso_system=espresso_system, 
                             name=cation_name, 
                             number_of_particles=N_cation)
        self.create_particle(espresso_system=espresso_system, 
                             name=anion_name, 
                             number_of_particles=N_anion)
        if c_salt_calculated.check('[substance] [length]**-3'):
            logging.info(f"added salt concentration of {c_salt_calculated.to('mol/L')} given by {N_cation} cations and {N_anion} anions")
        elif c_salt_calculated.check('[length]**-3'):
            logging.info(f"added salt concentration of {c_salt_calculated.to('reduced_length**-3')} given by {N_cation} cations and {N_anion} anions")
        return c_salt_calculated

    def create_bond(self, particle_id1, particle_id2, espresso_system, use_default_bond=False):
        """
        Creates a bond between two particle instances in an ESPResSo system and registers it in the pyMBE database.

        This method performs the following steps:
            1. Retrieves the particle instances corresponding to `particle_id1` and `particle_id2` from the database.
            2. Retrieves or creates the corresponding ESPResSo bond instance using the bond template.
            3. Adds the ESPResSo bond instance to the ESPResSo system if it was newly created.
            4. Adds the bond to the first particle's bond list in ESPResSo.
            5. Creates a `BondInstance` in the database and registers it.

        Args:
            particle_id1 (int): pyMBE and ESPResSo ID of the first particle.
            particle_id2 (int): pyMBE and ESPResSo ID of the second particle.
            espresso_system (espressomd.system.System): ESPResSo system object where the bond will be created.
            use_default_bond (bool, optional): If True, use a default bond template if no specific template exists. Defaults to False.
        """
        particle_inst_1 = self.db.get_instance(pmb_type="particle",
                                               instance_id=particle_id1)
        particle_inst_2 = self.db.get_instance(pmb_type="particle",
                                               instance_id=particle_id2)

        bond_inst = self.get_espresso_bond_instance(particle_name1=particle_inst_1.name,
                                                particle_name2=particle_inst_2.name,
                                                espresso_system=espresso_system,
                                                use_default_bond=use_default_bond)
        espresso_system.part.by_id(particle_id1).add_bond((bond_inst, particle_id2))
        if use_default_bond:
            bond_name = "default"
        else:
            bond_name = BondTemplate.make_bond_key(pn1=particle_inst_1.name, 
                                                   pn2=particle_inst_2.name)
        pmb_bond_instance = BondInstance(bond_id=self.db._propose_instance_id(pmb_type="bond"),
                                         name=bond_name,
                                         particle_id1=particle_id1,
                                         particle_id2=particle_id2)
        self.db._register_instance(instance=pmb_bond_instance)

    def create_counterions(self, object_name, cation_name, anion_name, espresso_system):
        """
        Creates particles of `cation_name` and `anion_name` in `espresso_system` to counter the net charge of `object_name`.
        
        Args:
            object_name(`str`): `name` of a pyMBE object.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.

        Returns: 
            counterion_number(`dict`): {"name": number}

        Note:
            This function currently does not support the creation of counterions for hydrogels.
        """ 
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=cation_name)
        cation_charge = cation_tpl.states[cation_tpl.initial_state].z
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=anion_name)
        anion_charge = anion_tpl.states[anion_tpl.initial_state].z
        object_ids = self.get_particle_id_map(object_name=object_name)["all"]
        counterion_number={}
        object_charge={}
        for name in ['positive', 'negative']:
            object_charge[name]=0
        for id in object_ids:
            if espresso_system.part.by_id(id).q > 0:
                object_charge['positive']+=1*(np.abs(espresso_system.part.by_id(id).q ))
            elif espresso_system.part.by_id(id).q < 0:
                object_charge['negative']+=1*(np.abs(espresso_system.part.by_id(id).q ))
        if object_charge['positive'] % abs(anion_charge) == 0:
            counterion_number[anion_name]=int(object_charge['positive']/abs(anion_charge))
        else:
            raise ValueError('The number of positive charges in the pmb_object must be divisible by the  charge of the anion')
        if object_charge['negative'] % abs(cation_charge) == 0:
            counterion_number[cation_name]=int(object_charge['negative']/cation_charge)
        else:
            raise ValueError('The number of negative charges in the pmb_object must be divisible by the  charge of the cation')
        if counterion_number[cation_name] > 0: 
            self.create_particle(espresso_system=espresso_system, 
                                 name=cation_name, 
                                 number_of_particles=counterion_number[cation_name])
        else:
            counterion_number[cation_name]=0
        if counterion_number[anion_name] > 0:
            self.create_particle(espresso_system=espresso_system, 
                                 name=anion_name, 
                                 number_of_particles=counterion_number[anion_name])
        else:
            counterion_number[anion_name] = 0
        logging.info('the following counter-ions have been created: ')
        for name in counterion_number.keys():
            logging.info(f'Ion type: {name} created number: {counterion_number[name]}')
        return counterion_number

    def create_hydrogel(self, name, espresso_system):
        """ 
        Creates a hydrogel in espresso_system using a pyMBE hydrogel template given by `name`

        Args:
            name(str): name of the hydrogel template in the pyMBE database.
            espresso_system (espressomd.system.System): ESPResSo system object where the hydrogel will be created.

        Returns:
            (int): id of the hydrogel instance created.
        """

        hydrogel_tpl = self.db.get_template(pmb_type="hydrogel",
                                            name=name)
        assembly_id = self.db._propose_instance_id(pmb_type="hydrogel")
        # Create the nodes
        nodes = {}
        node_topology = hydrogel_tpl.node_map
        for node in node_topology:
            node_index = node.lattice_index
            node_name = node.particle_name
            node_pos, node_id = self._create_hydrogel_node(node_index=node_index,
                                                          node_name=node_name,
                                                          espresso_system=espresso_system)
            node_label = self.lattice_builder._create_node_label(node_index=node_index)
            nodes[node_label] = {"name": node_name, "id": node_id, "pos": node_pos} 
            self.db._update_instance(instance_id=node_id,
                                     pmb_type="particle",
                                     attribute="assembly_id",
                                     value=assembly_id)
        for hydrogel_chain in hydrogel_tpl.chain_map:
            molecule_id = self._create_hydrogel_chain(hydrogel_chain=hydrogel_chain,
                                                      nodes=nodes, 
                                                      espresso_system=espresso_system)
            self.db._update_instance(instance_id=molecule_id,
                                     pmb_type="molecule",
                                     attribute="assembly_id",
                                     value=assembly_id)
        self.db._propagate_id(root_type="hydrogel", 
                                root_id=assembly_id, 
                                attribute="assembly_id", 
                                value=assembly_id)
        # Register an hydrogel instance in the pyMBE databasegit 
        self.db._register_instance(HydrogelInstance(name=name,
                                                    assembly_id=assembly_id))
        return assembly_id

    def create_molecule(self, name, number_of_molecules, espresso_system, list_of_first_residue_positions=None, backbone_vector=None, use_default_bond=False):
        """
        Creates `number_of_molecules` molecule of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the molecule type to be created. `name` must be defined in `pmb.df`
            espresso_system(`espressomd.system.System`): Instance of a system object from espressomd library.
            number_of_molecules(`int`): Number of molecules or peptides of type `name` to be created.
            list_of_first_residue_positions(`list`, optional): List of coordinates where the central bead of the first_residue_position will be created, random by default.
            backbone_vector(`list` of `float`): Backbone vector of the molecule, random by default. Central beads of the residues in the `residue_list` are placed along this vector. 
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle with undefined bonds in `pymbe.df`

        Returns:
            created_molecule_id_list(`list` of `int`): List with the `molecule_id` of the pyMBE molecule instances created into `espresso_system`.

        Note:
            Despite its name, this function can be used to create both molecules and peptides.    
        """
        supported_pmb_types = ["molecule", "peptide"]
        if number_of_molecules <= 0:
            return {}
        if list_of_first_residue_positions is not None:
            for item in list_of_first_residue_positions:
                if not isinstance(item, list):
                    raise ValueError("The provided input position is not a nested list. Should be a nested list with elements of 3D lists, corresponding to xyz coord.")
                elif len(item) != 3:
                    raise ValueError("The provided input position is formatted wrong. The elements in the provided list does not have 3 coordinates, corresponding to xyz coord.")

            if len(list_of_first_residue_positions) != number_of_molecules:
                raise ValueError(f"Number of positions provided in {list_of_first_residue_positions} does not match number of molecules desired, {number_of_molecules}")
        # Sanity tests, this function should work for both molecules and peptides
        registered_pmb_types_with_name = self.db._find_template_types(name=name)
        if len(registered_pmb_types_with_name) > 1: 
            raise KeyError(f"Detected multiple templates with the same name '{name}' in the pyMBE database, pmb_types: {registered_pmb_types_with_name}. Molecule creation aborted to avoid ambiguity.")  
        elif len(registered_pmb_types_with_name) == 0:
            logging.warning(f"No template with name '{name}'  defined in the pyMBE database, nothing will be created.")   
            return
        pmb_type = registered_pmb_types_with_name[0]
        if pmb_type not in supported_pmb_types:
            raise KeyError(f"Unsupported template type {pmb_type} for template {name}. Supported template types are {supported_pmb_types}")
        # Generate an arbitrary random unit vector
        if backbone_vector is None:
            backbone_vector = self.generate_random_points_in_a_sphere(center=[0,0,0],
                                                    radius=1, 
                                                    n_samples=1,
                                                    on_surface=True)[0]
        else:
            backbone_vector = np.array(backbone_vector)
        first_residue = True
        molecule_tpl = self.db.get_template(pmb_type=pmb_type,
                                            name=name)
        residue_list = molecule_tpl.residue_list
        pos_index = 0 
        molecule_ids = []
        for _ in range(number_of_molecules):        
            molecule_id = self.db._propose_instance_id(pmb_type=pmb_type)
            for residue in residue_list:
                if first_residue:
                    if list_of_first_residue_positions is None:
                        central_bead_pos = None
                    else:
                        for item in list_of_first_residue_positions:
                            central_bead_pos = [np.array(list_of_first_residue_positions[pos_index])]
                            
                    residue_id = self.create_residue(name=residue,
                                                     espresso_system=espresso_system, 
                                                     central_bead_position=central_bead_pos,  
                                                     use_default_bond= use_default_bond, 
                                                     backbone_vector=backbone_vector)
                    
                    # Add molecule_id to the residue instance and all particles associated
                    self.db._propagate_id(root_type="residue", 
                                           root_id=residue_id,
                                           attribute="molecule_id", 
                                           value=molecule_id)
                    particle_ids_in_residue = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                                      attribute="residue_id",
                                                                                      value=residue_id)
                    prev_central_bead_id = particle_ids_in_residue[0]
                    prev_central_bead_name = self.db.get_instance(pmb_type="particle", 
                                                                  instance_id=prev_central_bead_id).name
                    prev_central_bead_pos = espresso_system.part.by_id(prev_central_bead_id).pos
                    first_residue = False          
                else:
                    
                    # Calculate the starting position of the new residue
                    residue_tpl = self.db.get_template(pmb_type="residue",
                                                       name=residue)
                    lj_parameters = self.get_lj_parameters(particle_name1=prev_central_bead_name,
                                                           particle_name2=residue_tpl.central_bead)
                    bond_tpl = self.get_bond_template(particle_name1=prev_central_bead_name,
                                                      particle_name2=residue_tpl.central_bead,
                                                      use_default_bond=use_default_bond)
                    l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                                          bond_type=bond_tpl.bond_type,
                                                          bond_parameters=bond_tpl.get_parameters(ureg=self.units))
                    central_bead_pos = prev_central_bead_pos+backbone_vector*l0
                    # Create the residue
                    residue_id = self.create_residue(name=residue, 
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=[central_bead_pos],
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    # Add molecule_id to the residue instance and all particles associated
                    self.db._propagate_id(root_type="residue", 
                                          root_id=residue_id, 
                                          attribute="molecule_id", 
                                          value=molecule_id)
                    particle_ids_in_residue = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                                      attribute="residue_id",
                                                                                      value=residue_id)
                    central_bead_id = particle_ids_in_residue[0]

                    # Bond the central beads of the new and previous residues
                    self.create_bond(particle_id1=prev_central_bead_id,
                                     particle_id2=central_bead_id,
                                     espresso_system=espresso_system,
                                     use_default_bond=use_default_bond)
                    
                    prev_central_bead_id = central_bead_id                    
                    prev_central_bead_name = self.db.get_instance(pmb_type="particle", instance_id=central_bead_id).name
                    prev_central_bead_pos =central_bead_pos
            # Create a Peptide or Molecule instance and register it on the pyMBE database
            if pmb_type == "molecule":
                inst = MoleculeInstance(molecule_id=molecule_id,
                                        name=name)
            elif pmb_type == "peptide":
                inst = PeptideInstance(name=name,
                                       molecule_id=molecule_id)
            self.db._register_instance(inst)
            first_residue = True
            pos_index+=1
            molecule_ids.append(molecule_id)
            
        return molecule_id
    
    def create_particle(self, name, espresso_system, number_of_particles, position=None, fix=False):
        """
        Creates one or more particles in an ESPResSo system based on the particle template in the pyMBE database.
        
        Args:
            name(`str`): Label of the particle template in the pyMBE database. 
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            number_of_particles(`int`): Number of particles to be created.
            position(list of [`float`,`float`,`float`], optional): Initial positions of the particles. If not given, particles are created in random positions. Defaults to None.
            fix(`bool`, optional): Controls if the particle motion is frozen in the integrator, it is used to create rigid objects. Defaults to False.
        Returns:
            created_pid_list(`list` of `int`): List with the ids of the particles created into `espresso_system`.
        """       
        if number_of_particles <=0:
            return []
        if not self.db._has_template(name=name, pmb_type="particle"):
            logging.warning(f"Particle template with name '{name}' is not defined in the pyMBE database, no particle will be created.")
            return []
        
        part_tpl = self.db.get_template(pmb_type="particle",
                                        name=name)
        initial_state = part_tpl.states[part_tpl.initial_state]
        z = initial_state.z
        es_type = initial_state.es_type
        
        # Create the new particles into  ESPResSo 
        created_pid_list=[]
        for index in range(number_of_particles):
            if position is None:
                particle_position = self.rng.random((1, 3))[0] *np.copy(espresso_system.box_l)
            else:
                particle_position = position[index]
            
            particle_id = self.db._propose_instance_id(pmb_type="particle")
            created_pid_list.append(particle_id)
            kwargs = dict(id=particle_id, pos=particle_position, type=es_type, q=z)
            if fix:
                kwargs["fix"] = 3 * [fix]
            espresso_system.part.add(**kwargs)
            part_inst = ParticleInstance(name=name,
                                         particle_id=particle_id,
                                         initial_state=initial_state.name)
            self.db._register_instance(part_inst)
                              
        return created_pid_list

    def create_protein(self, name, number_of_proteins, espresso_system, topology_dict):
        """
        Creates one or more protein molecules in an ESPResSo system based on the 
        protein template in the pyMBE database and a provided topology.

        Args:
            name (str):
                Name of the protein template stored in the pyMBE database.
            
            number_of_proteins (int):
                Number of protein molecules to generate.  
            
            espresso_system (espressomd.system.System):
                The ESPResSo simulation system where the protein molecules will be created.
            
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

        Notes:
            - Particles are created using `create_particle()` with `fix=True`,
            meaning they are initially immobilized.
            - The function assumes all residues in `topology_dict` correspond to
            particle templates already defined in the pyMBE database.
            - Bonds between residues are not created here; it assumes a rigid body representation of the protein.
        """
        if number_of_proteins <= 0:
            return

        protein_tpl = self.db.get_template(pmb_type="protein", name=name)
        box_half = espresso_system.box_l[0] / 2.0

        residues = hf.get_residues_from_topology_dict(topology_dict=topology_dict,
                                                      model=protein_tpl.model)
        # Create protein
        for _ in range(number_of_proteins):
            # create a molecule identifier in pyMBE
            molecule_id = self.db._propose_instance_id(pmb_type="protein")
            # place protein COM randomly
            protein_center = self.generate_coordinates_outside_sphere(radius=1,
                                                                      max_dist=box_half,
                                                                      n_samples=1,
                                                                      center=[box_half]*3)[0]
            # CREATE RESIDUES + PARTICLES
            for _, rdata in residues.items():
                base_resname = rdata["resname"]  
                residue_name = f"AA-{base_resname}"
                # residue instance ID
                residue_id = self.db._propose_instance_id("residue")
                # register ResidueInstance
                self.db._register_instance(ResidueInstance(name=residue_name,
                                                           residue_id=residue_id,
                                                           molecule_id=molecule_id))

                # PARTICLE CREATION
                for bead_id in rdata["beads"]:
                    bead_type = re.split(r'\d+', bead_id)[0]
                    relative_pos = topology_dict[bead_id]["initial_pos"]
                    absolute_pos = relative_pos + protein_center
                    particle_id = self.create_particle(name=bead_type,
                                                       espresso_system=espresso_system,
                                                       number_of_particles=1,
                                                       position=[absolute_pos],
                                                       fix=True)[0]

                    # update metadata
                    self.db._update_instance(instance_id=particle_id,
                                             pmb_type="particle",
                                             attribute="molecule_id",
                                             value=molecule_id)
                    self.db._update_instance(instance_id=particle_id,
                                             pmb_type="particle",
                                             attribute="residue_id",
                                             value=residue_id)
            protein_inst = ProteinInstance(name=name,
                                           molecule_id=molecule_id)
            self.db._register_instance(protein_inst)

    def create_residue(self, name, espresso_system, central_bead_position=None,use_default_bond=False, backbone_vector=None):
        """
        Creates a residue of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the residue type to be created. `name` must be defined in `pmb.df`
            espresso_system(`espressomd.system.System`): Instance of a system object from espressomd library.
            central_bead_position(`list` of `float`): Position of the central bead.
            use_default_bond(`bool`): Switch to control if a bond of type `default` is used to bond a particle whose bonds types are not defined in `pmb.df`
            backbone_vector(`list` of `float`): Backbone vector of the molecule. All side chains are created perpendicularly to `backbone_vector`.

        Returns:
            (int) : residue_id of the residue created.
        """
        if not self.db._has_template(name=name, pmb_type="residue"):
            logging.warning(f"Residue template with name '{name}' is not defined in the pyMBE database, no residue will be created.")
            return
        res_tpl = self.db.get_template(pmb_type="residue",
                                       name=name)
        # Assign a residue_id
        residue_id = self.db._propose_instance_id(pmb_type="residue")
        res_inst = ResidueInstance(name=name,
                                   residue_id=residue_id)
        self.db._register_instance(res_inst)
        # create the principal bead   
        central_bead_name = res_tpl.central_bead 
        central_bead_id = self.create_particle(name=central_bead_name,
                                               espresso_system=espresso_system,
                                               position=central_bead_position,
                                               number_of_particles = 1)[0]
        
        central_bead_position=espresso_system.part.by_id(central_bead_id).pos
        # Assigns residue_id to the central_bead particle created.
        self.db._update_instance(pmb_type="particle",
                                 instance_id=central_bead_id,
                                 attribute="residue_id",
                                 value=residue_id)
        
        # create the lateral beads  
        side_chain_list = res_tpl.side_chains
        side_chain_beads_ids = []
        for side_chain_name in side_chain_list:
            pmb_type_list = self.db._find_template_types(name=side_chain_name)
            if len(pmb_type_list) > 1:
                raise KeyError(f"Detected multiple templates with the same name '{side_chain_name}' in the pyMBE database, pmb_types: {pmb_type_list}. Residue creation aborted to avoid ambiguity.")  
            elif not pmb_type_list:
                logging.warning(f"Element in side chain with name '{name}' is not defined in the pyMBE database, nothing will be created.")
                continue
            pmb_type = pmb_type_list[0]
            if pmb_type == 'particle':
                lj_parameters = self.get_lj_parameters(particle_name1=central_bead_name,
                                                       particle_name2=side_chain_name)
                bond_tpl = self.get_bond_template(particle_name1=central_bead_name,
                                                  particle_name2=side_chain_name,
                                                  use_default_bond=use_default_bond)
                l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                                      bond_type=bond_tpl.bond_type,
                                                      bond_parameters=bond_tpl.get_parameters(ureg=self.units))               
                if backbone_vector is None:
                    bead_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                radius=l0, 
                                                                n_samples=1,
                                                                on_surface=True)[0]
                else:
                    bead_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=np.array(backbone_vector),
                                                                                                magnitude=l0)
                    
                side_bead_id = self.create_particle(name=side_chain_name, 
                                                    espresso_system=espresso_system,
                                                    position=[bead_position], 
                                                    number_of_particles=1)[0]
                side_chain_beads_ids.append(side_bead_id)
                self.db._update_instance(pmb_type="particle",
                                         instance_id=side_bead_id,
                                         attribute="residue_id",
                                         value=residue_id)
                self.create_bond(particle_id1=central_bead_id,
                                 particle_id2=side_bead_id,
                                 espresso_system=espresso_system,
                                 use_default_bond=use_default_bond)
            elif pmb_type == 'residue':
                side_residue_tpl = self.db.get_template(name=side_chain_name,
                                                        pmb_type=pmb_type)
                central_bead_side_chain = side_residue_tpl.central_bead
                lj_parameters = self.get_lj_parameters(particle_name1=central_bead_name,
                                                       particle_name2=central_bead_side_chain)
                bond_tpl = self.get_bond_template(particle_name1=central_bead_name,
                                                  particle_name2=central_bead_side_chain,
                                                  use_default_bond=use_default_bond)
                l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                                      bond_type=bond_tpl.bond_type,
                                                      bond_parameters=bond_tpl.get_parameters(ureg=self.units))
                if backbone_vector is None:
                    residue_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                radius=l0, 
                                                                n_samples=1,
                                                                on_surface=True)[0]
                else:
                    residue_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                                    magnitude=l0)
                side_residue_id = self.create_residue(name=side_chain_name, 
                                                      espresso_system=espresso_system,
                                                      central_bead_position=[residue_position],
                                                      use_default_bond=use_default_bond)
                # Find particle ids of the inner residue
                side_chain_beads_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                               attribute="residue_id",
                                                                               value=side_residue_id)
                # Change the residue_id of the residue in the side chain to the one of the outer residue
                for particle_id in side_chain_beads_ids:
                    self.db._update_instance(instance_id=particle_id,
                                             pmb_type="particle",
                                             attribute="residue_id",
                                             value=residue_id)
                # Remove the instance of the inner residue
                self.db.delete_instance(pmb_type="residue",
                                        instance_id=side_residue_id)
                self.create_bond(particle_id1=central_bead_id,
                                 particle_id2=side_chain_beads_ids[0],
                                 espresso_system=espresso_system,
                                 use_default_bond=use_default_bond)        
        return  residue_id  

    def define_bond(self, bond_type, bond_parameters, particle_pairs):
        """
        Defines bond templates for each particle pair in `particle_pairs` in the pyMBE database.

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.
            particle_pairs(`lst`): list of the `names` of the `particles` to be bonded.

        Note:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain the following parameters:
                - k (`pint.Quantity`)      : Magnitude of the bond. It should have units of energy/length**2 
                using the `pmb.units` UnitRegistry.
                - r_0 (`pint.Quantity`)    : Equilibrium bond length. It should have units of length using 
                the `pmb.units` UnitRegistry.
           
            For a FENE bond the dictionary must contain the same parameters as for a HARMONIC bond and:              
                - d_r_max (`pint.Quantity`): Maximal stretching length for FENE. It should have 
                units of length using the `pmb.units` UnitRegistry. Default 'None'.
        """

        parameters_expected_dimensions={"r_0": "length",
                                        "k": "energy/length**2",
                                        "d_r_max": "length"}

        parameters_tpl = {}
        for key in bond_parameters.keys():
            parameters_tpl[key]= PintQuantity.from_quantity(q=bond_parameters[key],
                                                            expected_dimension=parameters_expected_dimensions[key],
                                                            ureg=self.units)

        bond_names=[]
        for particle_name1, particle_name2 in particle_pairs:
            
            tpl = BondTemplate(particle_name1=particle_name1,
                               particle_name2=particle_name2,
                               parameters=parameters_tpl,
                               bond_type=bond_type)
            tpl._make_name()
            if tpl.name in bond_names:
                raise RuntimeError(f"Bond {tpl.name} has already been defined, please check the list of particle pairs")
            else:
                self.db._register_template(tpl)
                bond_names.append(tpl.name)

    
    def define_default_bond(self, bond_type, bond_parameters):
        """
        Defines a bond template as a "default" template in the pyMBE database.
        
        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.
            
        Note:
            - Currently, only harmonic and FENE bonds are supported. 
        """
        parameters_expected_dimensions={"r_0": "length",
                                        "k": "energy/length**2",
                                        "d_r_max": "length"}
        parameters_tpl = {}
        for key in bond_parameters.keys():
            parameters_tpl[key]= PintQuantity.from_quantity(q=bond_parameters[key],
                                                            expected_dimension=parameters_expected_dimensions[key],
                                                            ureg=self.units)
        tpl = BondTemplate(parameters=parameters_tpl,
                               bond_type=bond_type)
        tpl.name = "default"
        self.db._register_template(tpl)
    
    def define_hydrogel(self, name, node_map, chain_map):
        """
        Defines a hydrogel template in the pyMBE database.

        Args:
            name(`str`): Unique label that identifies the `hydrogel`.
            node_map(`list of dict`): [{"particle_name": , "lattice_index": }, ... ]
            chain_map(`list of dict`): [{"node_start": , "node_end": , "residue_list": , ... ]
        """
        
        # Sanity tests
        node_indices = {tuple(entry['lattice_index']) for entry in node_map}
                
        chain_map_connectivity = set()
        for entry in chain_map:
            start = self.lattice_builder.node_labels[entry['node_start']]
            end = self.lattice_builder.node_labels[entry['node_end']]
            chain_map_connectivity.add((start,end))

        if self.lattice_builder.lattice.connectivity != chain_map_connectivity:
            raise ValueError("Incomplete hydrogel: A diamond lattice must contain correct 16 lattice index pairs")

        
        diamond_indices = {tuple(row) for row in self.lattice_builder.lattice.indices}
        if node_indices != diamond_indices:
            raise ValueError(f"Incomplete hydrogel: A diamond lattice must contain exactly 8 lattice indices, {diamond_indices} ")
        
        # Register information in the pyMBE database
        nodes=[]
        for entry in node_map:
            nodes.append(HydrogelNode(particle_name=entry["particle_name"],
                                      lattice_index=entry["lattice_index"]))
        chains=[]
        for chain in chain_map:
            chains.append(HydrogelChain(node_start=chain["node_start"],
                                        node_end=chain["node_end"],
                                        molecule_name=chain["molecule_name"]))
        tpl = HydrogelTemplate(name=name,
                               node_map=nodes,
                               chain_map=chains)
        self.db._register_template(tpl)

    def define_molecule(self, name, residue_list):
        """
        Defines a molecule template in the pyMBE database.

        Args:
            name(`str`): Unique label that identifies the `molecule`.
            residue_list(`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.  
        """
        tpl = MoleculeTemplate(name=name,
                               residue_list=residue_list)
        self.db._register_template(tpl)

    def define_particle(self, name,  sigma, epsilon, z=0, acidity=pd.NA, pka=pd.NA, cutoff=pd.NA, offset=pd.NA):
        """
        Defines a particle template in the pyMBE database.

        Args:
            name(`str`): Unique label that identifies this particle type.  
            sigma(`pint.Quantity`): Sigma parameter used to set up Lennard-Jones interactions for this particle type. 
            epsilon(`pint.Quantity`): Epsilon parameter used to setup Lennard-Jones interactions for this particle tipe.
            z(`int`, optional): Permanent charge number of this particle type. Defaults to 0.
            acidity(`str`, optional): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to pd.NA.
            pka(`float`, optional): If `particle` is an acid or a base, it defines its  pka-value. Defaults to pd.NA.
            cutoff(`pint.Quantity`, optional): Cutoff parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
            offset(`pint.Quantity`, optional): Offset parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
            
        Note:
            - `sigma`, `cutoff` and `offset` must have a dimensitonality of `[length]` and should be defined using pmb.units.
            - `epsilon` must have a dimensitonality of `[energy]` and should be defined using pmb.units.
            - `cutoff` defaults to `2**(1./6.) reduced_length`. 
            - `offset` defaults to 0.
            - For more information on `sigma`, `epsilon`, `cutoff` and `offset` check `pmb.setup_lj_interactions()`.
        """ 
        
        # If `cutoff` and `offset` are not defined, default them to the following values
        if pd.isna(cutoff):
            cutoff=self.units.Quantity(2**(1./6.), "reduced_length")
        if pd.isna(offset):
            offset=self.units.Quantity(0, "reduced_length")

        tpl = ParticleTemplate(name=name, 
                               sigma=PintQuantity.from_quantity(q=sigma, expected_dimension="length", ureg=self.units), 
                               epsilon=PintQuantity.from_quantity(q=epsilon, expected_dimension="energy", ureg=self.units),
                               cutoff=PintQuantity.from_quantity(q=cutoff, expected_dimension="length", ureg=self.units), 
                               offset=PintQuantity.from_quantity(q=offset, expected_dimension="length", ureg=self.units))
                
        # Define particle acid/base properties
        self.set_particle_acidity(particle_template=tpl, 
                                  acidity=acidity, 
                                  default_charge_number=z, 
                                  pka=pka)
        return 
    
    def define_peptide(self, name, sequence, model):
        """
        Defines a peptide template in the pyMBE database.

        Args:
            name (`str`): Unique label that identifies the peptide.
            sequence (`str`): Sequence of the peptide.
            model (`str`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
        """
        valid_keys = ['1beadAA','2beadAA']
        if model not in valid_keys:
            raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
        clean_sequence = self.protein_sequence_parser(sequence=sequence)    
        residue_list = self._get_residue_list_from_sequence(sequence=clean_sequence)
        tpl = PeptideTemplate(name=name,
                              residue_list=residue_list,
                              model=model,
                              sequence=sequence)
        self.db._register_template(tpl)        
    
    def define_protein(self, name, sequence, model):
        """
        Defines a protein template in the pyMBE database.

        Args:
            name (`str`): Unique label that identifies the protein.
            sequence (`str`): Sequence of the protein.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.

        Note:
            - Currently, only `lj_setup_mode="wca"` is supported. This corresponds to setting up the WCA potential.
        """
        valid_model_keys = ['1beadAA','2beadAA']
        if model not in valid_model_keys:
            raise ValueError('Invalid key for the protein model, supported models are {valid_model_keys}')
        
        residue_list = self._get_residue_list_from_sequence(sequence=sequence)
        tpl = ProteinTemplate(name=name,
                              model=model,
                              residue_list=residue_list,
                              sequence=sequence)
        self.db._register_template(tpl)
    
    def define_residue(self, name, central_bead, side_chains):
        """
        Defines a residue template in the pyMBE database.

        Args:
            name(`str`): Unique label that identifies the residue.
            central_bead(`str`): `name` of the `particle` to be placed as central_bead of the residue.
            side_chains(`list` of `str`): List of `name`s of the pmb_objects to be placed as side_chains of the residue. Currently, only pyMBE objects of type `particle` or `residue` are supported.
        """
        tpl = ResidueTemplate(name=name,
                              central_bead=central_bead,
                              side_chains=side_chains)
        self.db._register_template(tpl)
        return    

    def delete_instances_in_system(self, instance_id, pmb_type, espresso_system):
        """
        Deletes the instance with instance_id from the ESPResSo system. 
        Related assembly, molecule, residue, particles and bond instances will also be deleted from the pyMBE dataframe.

        Args:
            instance_id (int): id of the assembly to be deleted. 
            pmb_type (str): the instance type to be deleted. 
            espresso_system (espressomd.system.System): Instance of a system class from espressomd library.

        """
        if pmb_type == "particle":
            instance_identifier = "particle_id"
        elif pmb_type == "residue":
            instance_identifier = "residue_id"
        elif pmb_type in self.db._molecule_like_types:
            instance_identifier = "molecule_id"
        elif pmb_type in self.db._assembly_like_types:
            instance_identifier = "assembly_id"
        particle_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                               attribute=instance_identifier,
                                                               value=instance_id)
        self._delete_particles_from_espresso(particle_ids=particle_ids,
                                             espresso_system=espresso_system)
        
        self.db.delete_instance(pmb_type=pmb_type,
                                instance_id=instance_id,
                                cascade=True)

    def determine_reservoir_concentrations(self, pH_res, c_salt_res, activity_coefficient_monovalent_pair, max_number_sc_runs=200):
        """
        Determines the concentrations of the various species in the reservoir for given values of the pH and salt concentration.
        To do this, a system of nonlinear equations involving the pH, the ionic product of water, the activity coefficient of an
        ion pair and the concentrations of the various species is solved numerically using a self-consistent approach.
        More details can be found in chapter 5.3 of Landsgesell (doi.org/10.18419/opus-10831).
        This is a modified version of the code by Landsgesell et al. (doi.org/10.18419/darus-2237).

        Args:
            pH_res('float'): pH-value in the reservoir.
            c_salt_res('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            activity_coefficient_monovalent_pair('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.

        Returns:
            cH_res('pint.Quantity'): Concentration of H+ ions.
            cOH_res('pint.Quantity'): Concentration of OH- ions.
            cNa_res('pint.Quantity'): Concentration of Na+ ions.
            cCl_res('pint.Quantity'): Concentration of Cl- ions.
        """

        self_consistent_run = 0
        cH_res = 10**(-pH_res) * self.units.mol/self.units.l #initial guess for the concentration of H+

        def determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res):
            #Calculate and initial guess for the concentrations of various species based on ideal gas estimate
            cOH_res = self.Kw / cH_res 
            cNa_res = None
            cCl_res = None
            if cOH_res>=cH_res:
                #adjust the concentration of sodium if there is excess OH- in the reservoir:
                cNa_res = c_salt_res + (cOH_res-cH_res)
                cCl_res = c_salt_res
            else:
                # adjust the concentration of chloride if there is excess H+ in the reservoir
                cCl_res = c_salt_res + (cH_res-cOH_res)
                cNa_res = c_salt_res
                
            def calculate_concentrations_self_consistently(cH_res, cOH_res, cNa_res, cCl_res):
                nonlocal max_number_sc_runs, self_consistent_run
                if self_consistent_run<max_number_sc_runs:
                    self_consistent_run+=1
                    ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
                    cOH_res = self.Kw / (cH_res * activity_coefficient_monovalent_pair(ionic_strength_res))
                    if cOH_res>=cH_res:
                        #adjust the concentration of sodium if there is excess OH- in the reservoir:
                        cNa_res = c_salt_res + (cOH_res-cH_res)
                        cCl_res = c_salt_res
                    else:
                        # adjust the concentration of chloride if there is excess H+ in the reservoir
                        cCl_res = c_salt_res + (cH_res-cOH_res)
                        cNa_res = c_salt_res
                    return calculate_concentrations_self_consistently(cH_res, cOH_res, cNa_res, cCl_res)
                else:
                    return cH_res, cOH_res, cNa_res, cCl_res
            return calculate_concentrations_self_consistently(cH_res, cOH_res, cNa_res, cCl_res)

        cH_res, cOH_res, cNa_res, cCl_res = determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_pH = -np.log10(cH_res.to('mol/L').magnitude * np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))

        while abs(determined_pH-pH_res)>1e-6:
            if determined_pH > pH_res:
                cH_res *= 1.005
            else:
                cH_res /= 1.003
            cH_res, cOH_res, cNa_res, cCl_res = determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res)
            ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
            determined_pH = -np.log10(cH_res.to('mol/L').magnitude * np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))
            self_consistent_run=0

        return cH_res, cOH_res, cNa_res, cCl_res

    def enable_motion_of_rigid_object(self, name, espresso_system):
        '''
        Enables the motion of the rigid object `name` in the `espresso_system`.

        Args:
            name(`str`): Label of the object.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.

        Note:
            - It requires that espressomd has the following features activated: ["VIRTUAL_SITES_RELATIVE", "MASS"].
        '''
        logging.info('enable_motion_of_rigid_object requires that espressomd has the following features activated: ["VIRTUAL_SITES_RELATIVE", "MASS"]')
        self._check_supported_molecule(molecule_name=name,
                                        valid_pmb_types= ['protein'])
        molecule_ids_list = self.df.loc[self.df['name']==name].molecule_id.to_list()
        for molecule_id in molecule_ids_list:    
            particle_ids_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
            center_of_mass = self.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)
            rigid_object_center = espresso_system.part.add(pos=center_of_mass,
                                                           rotation=[True,True,True], 
                                                           type=self.propose_unused_type())
            
            rigid_object_center.mass = len(particle_ids_list)
            momI = 0
            for pid in particle_ids_list:
                momI += np.power(np.linalg.norm(center_of_mass - espresso_system.part.by_id(pid).pos), 2)

            rigid_object_center.rinertia = np.ones(3) * momI
            
            for particle_id in particle_ids_list:
                pid = espresso_system.part.by_id(particle_id)
                pid.vs_auto_relate_to(rigid_object_center.id)
        return

    def find_value_from_es_type(self, es_type, column_name):
        """
        Finds a value in `pmb.df` for a `column_name` and `es_type` pair.

        Args:
            es_type(`int`): value of the espresso type
            column_name(`str`): name of the column in `pymbe.df`

        Returns:
            Value in `pymbe.df` matching  `column_name` and `es_type`
        """
        idx = pd.IndexSlice
        for state in ['state_one', 'state_two']:            
            index = self.df.loc[self.df[(state, 'es_type')] == es_type].index
            if len(index) > 0:
                if column_name == 'label':
                    label = self.df.loc[idx[index[0]], idx[(state,column_name)]]
                    return label
                else: 
                    column_name_value = self.df.loc[idx[index[0]], idx[(column_name,'')]]
                    return column_name_value

    def format_node(self, node_list):
        return "[" + " ".join(map(str, node_list)) + "]"


    def generate_coordinates_outside_sphere(self, center, radius, max_dist, n_samples):
        """
        Generates coordinates outside a sphere centered at `center`.

        Args:
            center(`lst`): Coordinates of the center of the sphere.
            radius(`float`): Radius of the sphere.
            max_dist(`float`): Maximum distance from the center of the spahre to generate coordinates.
            n_samples(`int`): Number of sample points.

        Returns:
            coord_list(`lst`): Coordinates of the sample points.
        """
        if not radius > 0: 
            raise ValueError (f'The value of {radius} must be a positive value')
        if not radius < max_dist:
            raise ValueError(f'The min_dist ({radius} must be lower than the max_dist ({max_dist}))')
        coord_list = []
        counter = 0
        while counter<n_samples:
            coord = self.generate_random_points_in_a_sphere(center=center, 
                                            radius=max_dist,
                                            n_samples=1)[0]
            if np.linalg.norm(coord-np.asarray(center))>=radius:
                coord_list.append (coord)
                counter += 1
        return coord_list
    
    def generate_random_points_in_a_sphere(self, center, radius, n_samples, on_surface=False):
        """
        Uniformly samples points from a hypersphere. If on_surface is set to True, the points are
        uniformly sampled from the surface of the hypersphere.
        
        Args:
            center(`lst`): Array with the coordinates of the center of the spheres.
            radius(`float`): Radius of the sphere.
            n_samples(`int`): Number of sample points to generate inside the sphere.
            on_surface (`bool`, optional): If set to True, points will be uniformly sampled on the surface of the hypersphere.

        Returns:
            samples(`list`): Coordinates of the sample points inside the hypersphere.
        """
        # initial values
        center=np.array(center)
        d = center.shape[0]
        # sample n_samples points in d dimensions from a standard normal distribution
        samples = self.rng.normal(size=(n_samples, d))
        # make the samples lie on the surface of the unit hypersphere
        normalize_radii = np.linalg.norm(samples, axis=1)[:, np.newaxis]
        samples /= normalize_radii
        if not on_surface:
            # make the samples lie inside the hypersphere with the correct density
            uniform_points = self.rng.uniform(size=n_samples)[:, np.newaxis]
            new_radii = np.power(uniform_points, 1/d)
            samples *= new_radii
        # scale the points to have the correct radius and center
        samples = samples * radius + center
        return samples 

    def generate_trial_perpendicular_vector(self,vector,magnitude):
        """
        Generates an orthogonal vector to the input `vector`.

        Args:
            vector(`lst`): arbitrary vector.
            magnitude(`float`): magnitude of the orthogonal vector.
            
        Returns:
            (`lst`): Orthogonal vector with the same magnitude as the input vector.
        """ 
        np_vec = np.array(vector) 
        if np.all(np_vec == 0):
            raise ValueError('Zero vector')
        np_vec /= np.linalg.norm(np_vec) 
        # Generate a random vector 
        random_vector = self.generate_random_points_in_a_sphere(radius=1, 
                                                                center=[0,0,0],
                                                                n_samples=1, 
                                                                on_surface=True)[0]
        # Project the random vector onto the input vector and subtract the projection
        projection = np.dot(random_vector, np_vec) * np_vec
        perpendicular_vector = random_vector - projection
        # Normalize the perpendicular vector to have the same magnitude as the input vector
        perpendicular_vector /= np.linalg.norm(perpendicular_vector) 
        return perpendicular_vector*magnitude
            
    def get_bond_template(self, particle_name1, particle_name2, use_default_bond=False) :
        """
        Searches for bond template linking particle templates with `particle_name1` and `particle_name2` names in the pyMBE database and returns it.
        If `use_default_bond` is activated and a "default" bond is defined, returns the default bond template instead.

        Args:
            particle_name1(`str`): label of the type of the first particle type of the bonded particles.
            particle_name2(`str`): label of the type of the second particle type of the bonded particles.
            use_default_bond(`bool`, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            bond(`espressomd.interactions.BondedInteractions`): bond object from the espressomd library.
        
        Note:
            - If `use_default_bond`=True and no bond is defined between `particle_name1` and `particle_name2`, it returns the default bond defined in `pmb.df`.
        """
        if use_default_bond:
            bond_key = "default"
        else:
            bond_key = BondTemplate.make_bond_key(pn1=particle_name1,
                                              pn2=particle_name2)
        bond_tpl = self.db.get_template(name=bond_key,
                                        pmb_type="bond")
        return bond_tpl
    
    def get_charge_number_map(self):
        """
        Construct a mapping from ESPResSo particle types to their charge numbers.

        Returns:
            dict[int, float]:
                Dictionary mapping ESPResSo particle types to charge numbers,
                ``{es_type: z}``.

        Notes:
            - The mapping is built from particle *states*, not instances.
            - If multiple templates define states with the same ``es_type``,
            the last encountered definition will overwrite previous ones.
            This behavior is intentional and assumes database consistency.
            - Neutral particles (``z = 0``) are included in the map.
        """
        charge_number_map = {}
        particle_templates = self.db.get_templates("particle")
        for tpl in particle_templates.values():
            for state in tpl.states.values():
                if state.es_type is None:
                    continue
                charge_number_map[state.es_type] = state.z
        return charge_number_map


    def get_espresso_bond_instance(self, particle_name1, particle_name2, espresso_system, use_default_bond=False):
        """
        Retrieve or create a bond instance in an ESPResSo system for a given pair of particle names.

        This method checks whether a bond instance already exists in the database for the
        specified particle pair. If it exists, it retrieves the corresponding ESPResSo bond
        instance. Otherwise, it creates a new ESPResSo bond instance using the bond template.

        Args:
            particle_name1 (str): Name of the first particle involved in the bond.
            particle_name2 (str): Name of the second particle involved in the bond.
            espresso_system: An ESPResSo system object where the bond will be added or retrieved.
            use_default_bond (bool, optional): If True, use a default bond template when no 
                specific template exists for the particle pair. Defaults to False.

        Returns:
            (espressomd.interactions.BondedInteraction): The ESPResSo bond instance object.

        Raises:
            KeyError: If no bond template is found for the particle pair and `use_default_bond` is False.

        Note:
            When a new bond instance is created, it is not added to the ESPResSo system.
        """
        bond_tpl = self.get_bond_template(particle_name1=particle_name1,
                                          particle_name2=particle_name2,
                                          use_default_bond=use_default_bond)        
        if bond_tpl.name in self._bond_instances.keys():
            bond_inst = self._bond_instances[bond_tpl.name]
        else:   
            # Create an instance of the bond 
            bond_inst = self._create_espresso_bond_instance(bond_type=bond_tpl.bond_type,
                                                             bond_parameters=bond_tpl.get_parameters(self.units))
            self._bond_instances[bond_tpl.name]= bond_inst
            espresso_system.bonded_inter.add(bond_inst)
        return bond_inst

    def get_instances_df(self, pmb_type):
        """
        Returns a dataframe with all instances of type `pmb_type` in the pyMBE database.

        Args:
            pmb_type(`str`): pmb type to search instances in the pyMBE database.
        
        Returns:
            instances_df(`Pandas.Dataframe`): Dataframe with all instances of type `pmb_type`.
        """
        return self.db._get_instances_df(pmb_type=pmb_type)

    def get_lj_parameters(self, particle_name1, particle_name2, combining_rule='Lorentz-Berthelot'):
        """
        Returns the Lennard-Jones parameters for the interaction between the particle types given by 
        `particle_name1` and `particle_name2` in `pymbe.df`, calculated according to the provided combining rule.

        Args:
            particle_name1 (str): label of the type of the first particle type
            particle_name2 (str): label of the type of the second particle type
            combining_rule (`string`, optional): combining rule used to calculate `sigma` and `epsilon` for the potential betwen a pair of particles. Defaults to 'Lorentz-Berthelot'.

        Returns:
            {"epsilon": epsilon_value, "sigma": sigma_value, "offset": offset_value, "cutoff": cutoff_value}

        Note:
            - Currently, the only `combining_rule` supported is Lorentz-Berthelot.
            - If the sigma value of `particle_name1` or `particle_name2` is 0, the function will return an empty dictionary. No LJ interactions are set up for particles with sigma = 0.
        """
        supported_combining_rules=["Lorentz-Berthelot"]
        if combining_rule not in supported_combining_rules:
            raise ValueError(f"Combining_rule {combining_rule} currently not implemented in pyMBE, valid keys are {supported_combining_rules}")
        part_tpl1 = self.db.get_template(name=particle_name1,
                                         pmb_type="particle")
        part_tpl2 = self.db.get_template(name=particle_name2,
                                         pmb_type="particle")
        lj_parameters1 = part_tpl1.get_lj_parameters(ureg=self.units)
        lj_parameters2 = part_tpl2.get_lj_parameters(ureg=self.units)

        # If one of the particle has sigma=0, no LJ interations are set up between that particle type and the others    
        if part_tpl1.sigma.magnitude == 0 or part_tpl2.sigma.magnitude == 0:
            return {}
        # Apply combining rule
        if combining_rule == 'Lorentz-Berthelot':
            sigma=(lj_parameters1["sigma"]+lj_parameters2["sigma"])/2
            cutoff=(lj_parameters1["cutoff"]+lj_parameters2["cutoff"])/2
            offset=(lj_parameters1["offset"]+lj_parameters2["offset"])/2
            epsilon=np.sqrt(lj_parameters1["epsilon"]*lj_parameters2["epsilon"])
        return {"sigma": sigma, "cutoff": cutoff, "offset": offset, "epsilon": epsilon}    

    def get_particle_id_map(self, object_name):
        """
        Collect all particle IDs associated with an object of given name in the
        pyMBE database. Works for particles, residues, molecules, proteins,
        peptides, and assemblies.

        Relies in the internal method Manager.get_particle_id_map, see method for the detailed code.

        Args:
            object_name (str): Name of the object.

        Returns:
            dict: {"all": [particle_ids],
                   "residue_map": {residue_id: [particle_ids]},
                   "molecule_map": {molecule_id: [particle_ids]},
                   "assembly_map": {assembly_id: [particle_ids]},}
        """
        return self.db.get_particle_id_map(object_name=object_name)

    def get_pka_set(self):
        '''
        Gets the pka-values and acidities of the particles with acid/base properties in `pmb.df`
        
        Returns:
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}
        '''
        titratables_AA_df = self.df[[('name',''),('pka',''),('acidity','')]].drop_duplicates().dropna()
        pka_set = {}
        for index in titratables_AA_df.name.keys():
            name = titratables_AA_df.name[index]
            pka_value = titratables_AA_df.pka[index]
            acidity = titratables_AA_df.acidity[index]   
            pka_set[name] = {'pka_value':pka_value,'acidity':acidity}
        return pka_set 
    
    def get_radius_map(self, dimensionless=True):
        '''
        Gets the effective radius of each `espresso type` in `pmb.df`. 

        Args:
            dimensionless('bool'): controlls if the returned radii have a dimension. Defaults to False.
        
        Returns:
            radius_map(`dict`): {espresso_type: radius}.

        Note:
            The radius corresponds to (sigma+offset)/2
        '''
        if "particle" not in self.db._templates:
            return {}          
        result = {}
        for _, tpl in self.db._templates["particle"].items():
            radius = (tpl.sigma.to_quantity(self.units) + tpl.offset.to_quantity(self.units))/2.0
            if dimensionless:
                radius = radius.magnitude
            for _, state in tpl.states.items():
                result[state.es_type] = radius
        return result
        
    def get_reduced_units(self):
        """
        Returns the  current set of reduced units defined in pyMBE.

        Returns:
            reduced_units_text(`str`): text with information about the current set of reduced units.

        """
        unit_length=self.units.Quantity(1,'reduced_length')
        unit_energy=self.units.Quantity(1,'reduced_energy')
        unit_charge=self.units.Quantity(1,'reduced_charge')
        reduced_units_text = "\n".join(["Current set of reduced units:",
                                       f"{unit_length.to('nm'):.5g} = {unit_length}",
                                       f"{unit_energy.to('J'):.5g} = {unit_energy}",
                                       f"{unit_charge.to('C'):.5g} = {unit_charge}",
                                       f"Temperature: {(self.kT/self.kB).to('K'):.5g}"
                                        ])   
        return reduced_units_text

    def get_templates_df(self, pmb_type):
        """
        Returns a dataframe with all templates of type `pmb_type` in the pyMBE database.

        Args:
            pmb_type(`str`): pmb type to search templates in the pyMBE database.
        
        Returns:
            templates_df(`Pandas.Dataframe`): Dataframe with all templates of type `pmb_type`.
        """
        return self.db._get_templates_df(pmb_type=pmb_type)

    def get_type_map(self):
        """
        Return the mapping of ESPResSo types for all particle states defined in the pyMBE database.
        
        Returns:
            dict[str, int]:
                A dictionary mapping each particle state to its corresponding ESPResSo type:{
                            state_name: es_type,
                            ...
                        }
        """
        
        return self.db.get_es_types_map()

    def initialize_lattice_builder(self, diamond_lattice):
        """
        Initialize the lattice builder with the DiamondLattice object.

        Args:
            diamond_lattice(`DiamondLattice`): DiamondLattice object from the `lib/lattice` module to be used in the LatticeBuilder.
        """
        from .lib.lattice import LatticeBuilder, DiamondLattice
        if not isinstance(diamond_lattice, DiamondLattice):
            raise TypeError("Currently only DiamondLattice objects are supported.")
        self.lattice_builder = LatticeBuilder(lattice=diamond_lattice)
        logging.info(f"LatticeBuilder initialized with mpc={diamond_lattice.mpc} and box_l={diamond_lattice.box_l}")
        return self.lattice_builder

    def load_database(self, folder, format='csv'):
        """
        Loads a pyMBE database stored in `folder`.

        Args:
            folder (str or Path): Path to the folder where the pyMBE database was stored.
            format (str, optional): Format of the database to be loaded. Defaults to 'csv'.

        Note:
            - The folder must contain the files generated by `pmb.save_database()`.
            - Currently, only 'csv' format is supported.
        """
        supported_formats = ['csv']
        if format not in supported_formats:
            raise ValueError(f"Format {format} not supported. Supported formats are {supported_formats}")
        if format == 'csv':
            io._load_database_csv(self.db, 
                                folder=folder)
            

    def load_interaction_parameters(self, filename, overwrite=False):
        """
        Loads the interaction parameters stored in `filename` into `pmb.df`
        
        Args:
            filename(`str`): name of the file to be read
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
        """
        without_units = ['z','es_type']
        with_units = ['sigma','epsilon','offset','cutoff']

        with open(filename, 'r') as f:
            interaction_data = json.load(f)
            interaction_parameter_set = interaction_data["data"]

        for key in interaction_parameter_set:
            param_dict=interaction_parameter_set[key]
            object_type=param_dict.pop('object_type')
            if object_type == 'particle':
                not_required_attributes={}    
                for not_required_key in without_units+with_units:
                    if not_required_key in param_dict.keys():
                        if not_required_key in with_units:
                            not_required_attributes[not_required_key] = _DFm._create_variable_with_units(variable=param_dict.pop(not_required_key), 
                                                                                                         units_registry=self.units)
                        elif not_required_key in without_units:
                            not_required_attributes[not_required_key]=param_dict.pop(not_required_key)
                    else:
                        not_required_attributes[not_required_key]=pd.NA
                self.define_particle(name=param_dict.pop('name'),
                                z=not_required_attributes.pop('z'),
                                sigma=not_required_attributes.pop('sigma'),
                                offset=not_required_attributes.pop('offset'),
                                cutoff=not_required_attributes.pop('cutoff'),
                                epsilon=not_required_attributes.pop('epsilon'),
                                overwrite=overwrite)
            elif object_type == 'residue':
                self.define_residue(**param_dict)
            elif object_type == 'molecule':
                self.define_molecule(**param_dict)
            elif object_type == 'peptide':
                self.define_peptide(**param_dict)
            elif object_type == 'bond':
                particle_pairs = param_dict.pop('particle_pairs')
                bond_parameters = param_dict.pop('bond_parameters')
                bond_type = param_dict.pop('bond_type')
                if bond_type == 'harmonic':
                    k =  _DFm._create_variable_with_units(variable=bond_parameters.pop('k'), 
                                                          units_registry=self.units)
                    r_0 = _DFm._create_variable_with_units(variable=bond_parameters.pop('r_0'), 
                                                          units_registry=self.units)
                    bond = {'r_0'    : r_0,
                            'k'      : k,
                            }

                elif bond_type == 'FENE':
                    k = _DFm._create_variable_with_units(variable=bond_parameters.pop('k'), 
                                                         units_registry=self.units)
                    r_0 = _DFm._create_variable_with_units(variable=bond_parameters.pop('r_0'), 
                                                           units_registry=self.units)
                    d_r_max = _DFm._create_variable_with_units(variable=bond_parameters.pop('d_r_max'), 
                                                               units_registry=self.units)
                    bond =  {'r_0'    : r_0,
                             'k'      : k,
                             'd_r_max': d_r_max,
                             }
                else:
                    raise ValueError("Current implementation of pyMBE only supports harmonic and FENE bonds")
                self.define_bond(bond_type=bond_type, 
                                bond_parameters=bond, 
                                particle_pairs=particle_pairs)
            else:
                raise ValueError(object_type+' is not a known pmb object type')
            
        return
    
    def load_pka_set(self, filename, overwrite=True):
        """
        Loads the pka_set stored in `filename` into `pmb.df`.
        
        Args:
            filename(`str`): name of the file with the pka set to be loaded. Expected format is {name:{"acidity": acidity, "pka_value":pka_value}}.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to True. 
        """
        with open(filename, 'r') as f:
            pka_data = json.load(f)
            pka_set = pka_data["data"]

        self.check_pka_set(pka_set=pka_set)

        for key in pka_set:
            acidity = pka_set[key]['acidity']
            pka_value = pka_set[key]['pka_value']
            self.set_particle_acidity(name=key, 
                                      acidity=acidity, 
                                      pka=pka_value, 
                                      overwrite=overwrite)
        return


    def propose_unused_type(self):
        """
        Propose an unused ESPResSo particle type.

        This method scans the full `type_map` produced by `get_type_map()`,
        which contains all particle templates and their associated state `es_type`.
        It extracts all integer `es_type` values and returns the next available
        integer type, ensuring no collisions with existing ones.

        Returns:
            int: The next available integer ESPResSo type. Returns ``0`` if no
            integer types are currently defined.

    
        """
        type_map = self.get_type_map()
        # Flatten all es_type values across all particles and states
        all_types = []
        for es_type in type_map.values():
                all_types.append(es_type)
        # If no es_types exist, start at 0
        if not all_types:
            return 0
        return max(all_types) + 1

    def protein_sequence_parser(self, sequence):
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
    
    def read_pmb_df (self,filename):
        """
        Reads a pyMBE's Dataframe stored in `filename` and stores the information into pyMBE.

        Args:
            filename(`str`): path to file.

        Note:
            This function only accepts files with CSV format. 
        """
        if filename.rsplit(".", 1)[1] != "csv":
            raise ValueError("Only files with CSV format are supported")
        df = pd.read_csv (filename,header=[0, 1], index_col=0)
        self.df = _DFm._setup_df()
        columns_names = pd.MultiIndex.from_frame(self.df)
        columns_names = columns_names.names
        multi_index = pd.MultiIndex.from_tuples(columns_names)
        df.columns = multi_index
        _DFm._convert_columns_to_original_format(df=df, 
                                                 units_registry=self.units)
        self.df = df            
        self.df.fillna(pd.NA, 
                       inplace=True)
        return self.df
    
    def read_protein_vtf_in_df (self,filename,unit_length=None):
        """
        Loads a coarse-grained protein model in a vtf file `filename` into the `pmb.df` and it labels it with `name`.

        Args:
            filename(`str`): Path to the vtf file with the coarse-grained model.
            unit_length(`obj`): unit of length of the the coordinates in `filename` using the pyMBE UnitRegistry. Defaults to None.

        Returns:
            topology_dict(`dict`): {'initial_pos': coords_list, 'chain_id': id, 'sigma': sigma_value}

        Note:
            - If no `unit_length` is provided, it is assumed that the coordinates are in Angstrom.
        """

        logging.info(f'Loading protein coarse grain model file: {filename}')

        coord_list = []
        particles_dict = {}

        if unit_length is None:
            unit_length = 1 * self.units.angstrom 

        with open (filename,'r') as protein_model:
            for line in protein_model :
                line_split = line.split()
                if line_split:
                    line_header = line_split[0]
                    if line_header == 'atom':
                        atom_id  = line_split[1]
                        atom_name = line_split[3]
                        atom_resname = line_split[5]
                        chain_id = line_split[9]
                        radius = float(line_split [11])*unit_length 
                        particles_dict [int(atom_id)] = [atom_name , atom_resname, chain_id, radius]
                    elif line_header.isnumeric(): 
                        atom_coord = line_split[1:] 
                        atom_coord = [(float(i)*unit_length).to('reduced_length').magnitude for i in atom_coord]
                        coord_list.append (atom_coord)

        numbered_label = []
        i = 0   
        
        for atom_id in particles_dict.keys():
    
            if atom_id == 1:
                atom_name = particles_dict[atom_id][0]
                numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2],particles_dict[atom_id][3]]
                numbered_label.append(numbered_name)

            elif atom_id != 1: 
            
                if particles_dict[atom_id-1][1] != particles_dict[atom_id][1]:
                    i += 1                    
                    count = 1
                    atom_name = particles_dict[atom_id][0]
                    numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2],particles_dict[atom_id][3]]
                    numbered_label.append(numbered_name)
                    
                elif particles_dict[atom_id-1][1] == particles_dict[atom_id][1]:
                    if count == 2 or particles_dict[atom_id][1] == 'GLY':
                        i +=1  
                        count = 0
                    atom_name = particles_dict[atom_id][0]
                    numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2],particles_dict[atom_id][3]]
                    numbered_label.append(numbered_name)
                    count +=1

        topology_dict = {}

        for i in range (0, len(numbered_label)):   
            topology_dict [numbered_label[i][0]] = {'initial_pos': coord_list[i] ,
                                                    'chain_id':numbered_label[i][1],
                                                    'radius':numbered_label[i][2] }

        return topology_dict

    def save_database(self, folder, format='csv'):
        """
        Saves the current pyMBE database into a file `filename`.

        Args:
            folder (str or Path): Path to the folder where the database files will be saved.

        """
        supported_formats = ['csv']
        if format not in supported_formats:
            raise ValueError(f"Format {format} not supported. Supported formats are: {supported_formats}")
        if format == 'csv':
            io._save_database_csv(self.db, 
                                folder=folder)

    def search_particles_in_residue(self, residue_name):
        '''
        Searches for all particles in a given residue of name `residue_name`.

        Args:
            residue_name (`str`): name of the residue to be searched

        Returns:
            list_of_particles_in_residue (`lst`): list of the names of all particles in the residue

        Note:
            - The function returns a name per particle in residue, i.e. if there are multiple particles with the same type `list_of_particles_in_residue` will have repeated items.
            - The function will return an empty list if the residue is not defined in `pmb.df`.
            - The function will return an empty list if the particles are not defined in the pyMBE DataFrame.
        '''
        index_residue = self.df.loc[self.df['name'] == residue_name].index[0].item() 
        central_bead = self.df.at [index_residue, ('central_bead', '')]
        list_of_side_chains = self.df.at[index_residue, ('side_chains', '')]
        list_of_particles_in_residue = []
        if central_bead is not pd.NA:
                    list_of_particles_in_residue.append(central_bead)
        if list_of_side_chains is not pd.NA:
            for side_chain in list_of_side_chains:
                if _DFm._check_if_name_is_defined_in_df(name=side_chain, df=self.df): 
                    object_type = self.df[self.df['name']==side_chain].pmb_type.values[0]
                else:
                    continue
                if object_type == "residue":
                    list_of_particles_in_side_chain_residue = self.search_particles_in_residue(side_chain)
                    list_of_particles_in_residue += list_of_particles_in_side_chain_residue
                elif object_type == "particle":
                    if side_chain is not pd.NA:
                        list_of_particles_in_residue.append(side_chain)
        return list_of_particles_in_residue        

    def set_particle_acidity(self, particle_template, acidity=pd.NA, default_charge_number=0, pka=pd.NA):
        """
        Sets the particle acidity including the charges in each of its possible states. 

        Args:
            name(`str`): Unique label that identifies the `particle`. 
            acidity(`str`): Identifies whether the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to None.
            default_charge_number (`int`): Charge number of the particle. Defaults to 0.
            pka(`float`, optional):  If `particle` is an acid or a base, it defines its pka-value. Defaults to pandas.NA.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
     
        Note:
            - For particles with  `acidity = acidic` or `acidity = basic`, `state_one` and `state_two` correspond to the protonated and 
        deprotonated states, respectively. 
            - For particles without an acidity `acidity = pandas.NA`, only `state_one` is defined.
            - Each state has the following properties as sub-indexes: `label`,`charge` and `es_type`. 
        """
        acidity_valid_keys = ['inert','acidic', 'basic']
        if not pd.isna(acidity):
            if acidity not in acidity_valid_keys:
                raise ValueError(f"Acidity {acidity} provided for particle name  {particle_template.name} is not supported. Valid keys are: {acidity_valid_keys}")
            if acidity in ['acidic', 'basic'] and pd.isna(pka):
                raise ValueError(f"pKa not provided for particle with name {particle_template.name} with acidity {acidity}. pKa must be provided for acidic or basic particles.")   
            if acidity == "inert":
                acidity = pd.NA
                logging.warning("the keyword 'inert' for acidity has been replaced by setting acidity = pd.NA. For backwards compatibility, acidity has been set to pd.NA. Support for `acidity = 'inert'` may be deprecated in future releases of pyMBE")

        
        # Define the first state
        if pka is pd.NA:
            # Inert particle with a single state
            z_state_one = default_charge_number
            name_state_one = particle_template.name
        else:
            if acidity == "acidic":
                z_state_one = 0
            elif acidity == "basic":
                z_state_one = 1
            name_state_one = particle_template.name + "H"

        particle_template.add_state(ParticleState(name=name_state_one,
                                    z=z_state_one,
                                    es_type=self.propose_unused_type()))
        self.db._register_template(particle_template)

        # For monoprotic acid/base particles, define the second state
        if pka is not pd.NA:
            if acidity == "acidic":
                z_state_two = -1
            elif acidity == "basic":
                z_state_two = 0
            name_state_two = particle_template.name
            particle_template.add_state(ParticleState(name=name_state_two,
                                        z=z_state_two,
                                        es_type=self.propose_unused_type()))
        
            reaction = Reaction(participants=[ReactionParticipant(particle_name=particle_template.name,
                                                                  state_name=name_state_one, 
                                                                  coefficient=-1),
                                              ReactionParticipant(particle_name=particle_template.name,
                                                                  state_name=name_state_two,
                                                                  coefficient=1)],
                                reaction_type="acid/base",
                                pK=pka)
            self.db._register_reaction(reaction)

    
    def set_reduced_units(self, unit_length=None, unit_charge=None, temperature=None, Kw=None):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args:
            unit_length(`pint.Quantity`,optional): Reduced unit of length defined using the `pmb.units` UnitRegistry. Defaults to None. 
            unit_charge(`pint.Quantity`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            temperature(`pint.Quantity`,optional): Temperature of the system, defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw(`pint.Quantity`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 

        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
        """
        if unit_length is None:
            unit_length= 0.355*self.units.nm
        if temperature is None:
            temperature = 298.15 * self.units.K
        if unit_charge is None:
            unit_charge = scipy.constants.e * self.units.C
        if Kw is None:
            Kw = 1e-14
        # Sanity check
        variables=[unit_length,temperature,unit_charge]
        dimensionalities=["[length]","[temperature]","[charge]"]
        for variable,dimensionality in zip(variables,dimensionalities):
            self.check_dimensionality(variable,dimensionality)
        self.Kw=Kw*self.units.mol**2 / (self.units.l**2)
        self.kT=temperature*self.kB
        self.units._build_cache()
        self.units.define(f'reduced_energy = {self.kT} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')
        logging.info(self.get_reduced_units())
        return

    def setup_cpH (self, counter_ion, constant_pH, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up the Acid/Base reactions for acidic/basic `particles` defined in `pmb.df` to be sampled in the constant pH ensemble. 

        Args:
            counter_ion(`str`): `name` of the counter_ion `particle`.
            constant_pH(`float`): pH-value.
            exclusion_range(`pint.Quantity`, optional): Below this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.

        Returns:
            RE(`reaction_methods.ConstantpHEnsemble`): Instance of a reaction_methods.ConstantpHEnsemble object from the espressomd library.
            sucessfull_reactions_labels(`lst`): Labels of the reactions set up by pyMBE.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if pka_set is None:
            pka_set=self.get_pka_set()    
        self.check_pka_set(pka_set=pka_set)
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ConstantpHEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    constant_pH=constant_pH,
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )
        sucessfull_reactions_labels=[]
        charge_number_map = self.get_charge_number_map()
        for name in pka_set.keys():
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f'The acid-base reaction of {name} has not been set up because its particle type is not defined in the pyMBE DataFrame.')
                continue
            gamma=10**-pka_set[name]['pka_value']
            state_one_type   = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type   = self.df.loc[self.df['name']==name].state_two.es_type.values[0]
            counter_ion_type = self.df.loc[self.df['name']==counter_ion].state_one.es_type.values[0]
            RE.add_reaction(gamma=gamma,
                            reactant_types=[state_one_type],
                            product_types=[state_two_type, counter_ion_type],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            counter_ion_type: charge_number_map[counter_ion_type]})
            sucessfull_reactions_labels.append(name)
        return RE, sucessfull_reactions_labels

    def setup_gcmc(self, c_salt_res, salt_cation_name, salt_anion_name, activity_coefficient, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up grand-canonical coupling to a reservoir of salt.
        For reactive systems coupled to a reservoir, the grand-reaction method has to be used instead.

        Args:
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable'): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`pint.Quantity`, optional): For distances shorter than this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.

        Returns:
            RE (`reaction_methods.ReactionEnsemble`): Instance of a reaction_methods.ReactionEnsemble object from the espressomd library.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        determined_activity_coefficient = activity_coefficient(c_salt_res)
        K_salt = (c_salt_res.to('1/(N_A * reduced_length**3)')**2) * determined_activity_coefficient

        salt_cation_es_type = self.df.loc[self.df['name']==salt_cation_name].state_one.es_type.values[0]
        salt_anion_es_type = self.df.loc[self.df['name']==salt_anion_name].state_one.es_type.values[0]     

        salt_cation_charge = self.df.loc[self.df['name']==salt_cation_name].state_one.z.values[0]
        salt_anion_charge = self.df.loc[self.df['name']==salt_anion_name].state_one.z.values[0]     

        if salt_cation_charge <= 0:
            raise ValueError('ERROR salt cation charge must be positive, charge ', salt_cation_charge)
        if salt_anion_charge >= 0:
            raise ValueError('ERROR salt anion charge must be negative, charge ', salt_anion_charge)

        # Grand-canonical coupling to the reservoir
        RE.add_reaction(
            gamma = K_salt.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ salt_cation_es_type, salt_anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                salt_cation_es_type: salt_cation_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        return RE

    def setup_grxmc_reactions(self, pH_res, c_salt_res, proton_name, hydroxide_name, salt_cation_name, salt_anion_name, activity_coefficient, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up Acid/Base reactions for acidic/basic 'particles' defined in 'pmb.df', as well as a grand-canonical coupling to a 
        reservoir of small ions. 
        This implementation uses the original formulation of the grand-reaction method by Landsgesell et al. [1].

        [1] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.

        Args:
            pH_res ('float): pH-value in the reservoir.
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            proton_name ('str'): Name of the proton (H+) particle.
            hydroxide_name ('str'): Name of the hydroxide (OH-) particle.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable'): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`pint.Quantity`, optional): For distances shorter than this value, no particles will be inserted.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.

        Returns:
            RE (`obj`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
            sucessful_reactions_labels(`lst`): Labels of the reactions set up by pyMBE.
            ionic_strength_res ('pint.Quantity'): Ionic strength of the reservoir (useful for calculating partition coefficients).
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if pka_set is None:
            pka_set=self.get_pka_set()    
        self.check_pka_set(pka_set=pka_set)
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        cH_res, cOH_res, cNa_res, cCl_res = self.determine_reservoir_concentrations(pH_res, c_salt_res, activity_coefficient)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_activity_coefficient = activity_coefficient(ionic_strength_res)
        K_W = cH_res.to('1/(N_A * reduced_length**3)') * cOH_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        K_NACL = cNa_res.to('1/(N_A * reduced_length**3)') * cCl_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        K_HCL = cH_res.to('1/(N_A * reduced_length**3)') * cCl_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient

        proton_es_type = self.df.loc[self.df['name']==proton_name].state_one.es_type.values[0]
        hydroxide_es_type = self.df.loc[self.df['name']==hydroxide_name].state_one.es_type.values[0]     
        salt_cation_es_type = self.df.loc[self.df['name']==salt_cation_name].state_one.es_type.values[0]
        salt_anion_es_type = self.df.loc[self.df['name']==salt_anion_name].state_one.es_type.values[0]     

        proton_charge = self.df.loc[self.df['name']==proton_name].state_one.z.values[0]
        hydroxide_charge = self.df.loc[self.df['name']==hydroxide_name].state_one.z.values[0]     
        salt_cation_charge = self.df.loc[self.df['name']==salt_cation_name].state_one.z.values[0]
        salt_anion_charge = self.df.loc[self.df['name']==salt_anion_name].state_one.z.values[0]     

        if proton_charge <= 0:
            raise ValueError('ERROR proton charge must be positive, charge ', proton_charge)
        if salt_cation_charge <= 0:
            raise ValueError('ERROR salt cation charge must be positive, charge ', salt_cation_charge)
        if hydroxide_charge >= 0:
            raise ValueError('ERROR hydroxide charge must be negative, charge ', hydroxide_charge)
        if salt_anion_charge >= 0:
            raise ValueError('ERROR salt anion charge must be negative, charge ', salt_anion_charge)

        # Grand-canonical coupling to the reservoir
        # 0 = H+ + OH-
        RE.add_reaction(
            gamma = K_W.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ proton_es_type, hydroxide_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                proton_es_type: proton_charge, 
                hydroxide_es_type: hydroxide_charge, 
            }
        )

        # 0 = Na+ + Cl-
        RE.add_reaction(
            gamma = K_NACL.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ salt_cation_es_type, salt_anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                salt_cation_es_type: salt_cation_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        # 0 = Na+ + OH-
        RE.add_reaction(
            gamma = (K_NACL * K_W / K_HCL).magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ salt_cation_es_type, hydroxide_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                salt_cation_es_type: salt_cation_charge, 
                hydroxide_es_type: hydroxide_charge, 
            }
        )

        # 0 = H+ + Cl-
        RE.add_reaction(
            gamma = K_HCL.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ proton_es_type, salt_anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                proton_es_type: proton_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        # Annealing moves to ensure sufficient sampling
        # Cation annealing H+ = Na+
        RE.add_reaction(
            gamma = (K_NACL / K_HCL).magnitude,
            reactant_types = [proton_es_type],
            reactant_coefficients = [ 1 ],
            product_types = [ salt_cation_es_type ],
            product_coefficients = [ 1 ],
            default_charges = {
                proton_es_type: proton_charge, 
                salt_cation_es_type: salt_cation_charge, 
            }
        )

        # Anion annealing OH- = Cl- 
        RE.add_reaction(
            gamma = (K_HCL / K_W).magnitude,
            reactant_types = [hydroxide_es_type],
            reactant_coefficients = [ 1 ],
            product_types = [ salt_anion_es_type ],
            product_coefficients = [ 1 ],
            default_charges = {
                hydroxide_es_type: hydroxide_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        sucessful_reactions_labels=[]
        charge_number_map = self.get_charge_number_map()
        for name in pka_set.keys():
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f'The acid-base reaction of {name} has not been set up because its particle type is not defined in the dataframe.')
                continue

            Ka = (10**-pka_set[name]['pka_value'] * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')

            state_one_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type = self.df.loc[self.df['name']==name].state_two.es_type.values[0]

            # Reaction in terms of proton: HA = A + H+
            RE.add_reaction(gamma=Ka.magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, proton_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            proton_es_type: proton_charge})

            # Reaction in terms of salt cation: HA = A + Na+
            RE.add_reaction(gamma=(Ka * K_NACL / K_HCL).magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, salt_cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            salt_cation_es_type: salt_cation_charge})

            # Reaction in terms of hydroxide: OH- + HA = A
            RE.add_reaction(gamma=(Ka / K_W).magnitude,
                            reactant_types=[state_one_type, hydroxide_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            hydroxide_es_type: hydroxide_charge})

            # Reaction in terms of salt anion: Cl- + HA = A
            RE.add_reaction(gamma=(Ka / K_HCL).magnitude,
                            reactant_types=[state_one_type, salt_anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            salt_anion_es_type: salt_anion_charge})

            sucessful_reactions_labels.append(name)
        return RE, sucessful_reactions_labels, ionic_strength_res

    def setup_grxmc_unified(self, pH_res, c_salt_res, cation_name, anion_name, activity_coefficient, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up Acid/Base reactions for acidic/basic 'particles' defined in 'pmb.df', as well as a grand-canonical coupling to a 
        reservoir of small ions. 
        This implementation uses the formulation of the grand-reaction method by Curk et al. [1], which relies on "unified" ion types X+ = {H+, Na+} and X- = {OH-, Cl-}. 
        A function that implements the original version of the grand-reaction method by Landsgesell et al. [2] is also available under the name 'setup_grxmc_reactions'.

        [1] Curk, T., Yuan, J., & Luijten, E. (2022). Accelerated simulation method for charge regulation effects. The Journal of Chemical Physics, 156(4).
        [2] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.

        Args:
            pH_res ('float'): pH-value in the reservoir.
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            cation_name ('str'): Name of the cationic particle.
            anion_name ('str'): Name of the anionic particle.
            activity_coefficient ('callable'): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`pint.Quantity`, optional): Below this value, no particles will be inserted.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius per each espresso_type. Defaults to `False`.

        Returns:
            RE (`reaction_ensemble.ReactionEnsemble`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
            sucessful_reactions_labels(`lst`): Labels of the reactions set up by pyMBE.
            ionic_strength_res ('float'): Ionic strength of the reservoir (useful for calculating partition coefficients).
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if pka_set is None:
            pka_set=self.get_pka_set()    
        self.check_pka_set(pka_set=pka_set)
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        cH_res, cOH_res, cNa_res, cCl_res = self.determine_reservoir_concentrations(pH_res, c_salt_res, activity_coefficient)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_activity_coefficient = activity_coefficient(ionic_strength_res)
        a_hydrogen = (10 ** (-pH_res) * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')
        a_cation = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * np.sqrt(determined_activity_coefficient)
        a_anion = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * np.sqrt(determined_activity_coefficient)
        K_XX = a_cation * a_anion

        cation_es_type = self.df.loc[self.df['name']==cation_name].state_one.es_type.values[0]
        anion_es_type = self.df.loc[self.df['name']==anion_name].state_one.es_type.values[0]     
        cation_charge = self.df.loc[self.df['name']==cation_name].state_one.z.values[0]
        anion_charge = self.df.loc[self.df['name']==anion_name].state_one.z.values[0]     
        if cation_charge <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ', cation_charge)
        if anion_charge >= 0:
            raise ValueError('ERROR anion charge must be negative, charge ', anion_charge)

        # Coupling to the reservoir: 0 = X+ + X-
        RE.add_reaction(
            gamma = K_XX.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ cation_es_type, anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                cation_es_type: cation_charge, 
                anion_es_type: anion_charge, 
            }
        )

        sucessful_reactions_labels=[]
        charge_number_map = self.get_charge_number_map()
        for name in pka_set.keys():
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f'The acid-base reaction of {name} has not been set up because its particle type is not defined in the dataframe.')
                continue

            Ka = 10**-pka_set[name]['pka_value'] * self.units.mol/self.units.l
            gamma_K_AX = Ka.to('1/(N_A * reduced_length**3)').magnitude * a_cation / a_hydrogen

            state_one_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type = self.df.loc[self.df['name']==name].state_two.es_type.values[0]

            # Reaction in terms of small cation: HA = A + X+
            RE.add_reaction(gamma=gamma_K_AX.magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            cation_es_type: cation_charge})

            # Reaction in terms of small anion: X- + HA = A
            RE.add_reaction(gamma=gamma_K_AX.magnitude / K_XX.magnitude,
                            reactant_types=[state_one_type, anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            anion_es_type: anion_charge})

            sucessful_reactions_labels.append(name)
        return RE, sucessful_reactions_labels, ionic_strength_res

    def setup_lj_interactions(self, espresso_system, shift_potential=True, combining_rule='Lorentz-Berthelot'):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle states defined in the pyMBE database.

        Args:
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            shift_potential(`bool`, optional): If True, a shift will be automatically computed such that the potential is continuous at the cutoff radius. Otherwise, no shift will be applied. Defaults to True.
            combining_rule(`string`, optional): combining rule used to calculate `sigma` and `epsilon` for the potential between a pair of particles. Defaults to 'Lorentz-Berthelot'.
            warning(`bool`, optional): switch to activate/deactivate warning messages. Defaults to True.

        Note:
            - Currently, the only `combining_rule` supported is Lorentz-Berthelot.
            - Check the documentation of ESPResSo for more info about the potential https://espressomd.github.io/doc4.2.0/inter_non-bonded.html

        """
        from itertools import combinations_with_replacement

        particle_templates = self.db.get_templates("particle")

        shift = "auto" if shift_potential else 0

        # Flatten states with template context
        state_entries = []
        for tpl in particle_templates.values():
            for state in tpl.states.values():
                state_entries.append((tpl, state))

        # Iterate over all unique state pairs
        for (tpl1, state1), (tpl2, state2) in combinations_with_replacement(state_entries, 2):

            lj_parameters = self.get_lj_parameters(particle_name1=tpl1.name,
                                                   particle_name2=tpl2.name,
                                                   combining_rule=combining_rule)
            if not lj_parameters:
                continue

            espresso_system.non_bonded_inter[state1.es_type, state2.es_type].lennard_jones.set_params(
                epsilon=lj_parameters["epsilon"].to("reduced_energy").magnitude,
                sigma=lj_parameters["sigma"].to("reduced_length").magnitude,
                cutoff=lj_parameters["cutoff"].to("reduced_length").magnitude,
                offset=lj_parameters["offset"].to("reduced_length").magnitude,
                shift=shift)

            lj_template = LJInteractionTemplate(state1=state1.name,
                                                state2=state2.name,
                                                sigma=PintQuantity.from_quantity(q=lj_parameters["sigma"],
                                                                                 expected_dimension="length",
                                                                                 ureg=self.units),
                                                epsilon=PintQuantity.from_quantity(q=lj_parameters["epsilon"],
                                                                                   expected_dimension="energy",
                                                                                   ureg=self.units),
                                                cutoff=PintQuantity.from_quantity(q=lj_parameters["cutoff"],
                                                                                  expected_dimension="length",
                                                                                  ureg=self.units),
                                                offset=PintQuantity.from_quantity(q=lj_parameters["offset"],
                                                                                  expected_dimension="length",
                                                                                  ureg=self.units),
                                                shift=shift)
            self.db._register_template(lj_template)

    def write_pmb_df (self, filename):
        '''
        Writes the pyMBE dataframe into a csv file
        
        Args:
            filename(`str`): Path to the csv file 
        '''

        columns_with_list_or_dict = ['residue_list','side_chains', 'parameters_of_the_potential','sequence']
        df = self.df.copy(deep=True)
        for column_name in columns_with_list_or_dict:
            df[column_name] = df[column_name].apply(lambda x: json.dumps(x) if isinstance(x, (np.ndarray, tuple, list, dict)) or pd.notnull(x) else x)
        df['bond_object'] = df['bond_object'].apply(lambda x: f'{x.__class__.__name__}({json.dumps({**x.get_params(), "bond_id": x._bond_id})})' if pd.notnull(x) else x)
        df.fillna(pd.NA, inplace=True)
        df.to_csv(filename)
        return
