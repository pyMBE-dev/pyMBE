#
# Copyright (C) 2023-2024 pyMBE-dev team
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
import sys
import ast
import json
import pint
import numpy as np
import pandas as pd
import scipy.optimize


class pymbe_library():

    """
    The library for the Molecular Builder for ESPResSo (pyMBE)

    Attributes:
        N_A(`obj`): Avogadro number using the `pmb.units` UnitRegistry.
        Kb(`obj`): Boltzmann constant using the `pmb.units` UnitRegistry.
        e(`obj`): Elemental charge using the `pmb.units` UnitRegistry.
        df(`obj`): PandasDataframe used to bookkeep all the information stored in pyMBE. Typically refered as `pmb.df`. 
        kT(`obj`): Thermal energy using the `pmb.units` UnitRegistry. It is used as the unit of reduced energy.
        Kw(`obj`): Ionic product of water using the `pmb.units` UnitRegistry. Used in the setup of the G-RxMC method.
    """
    units = pint.UnitRegistry()
    N_A=6.02214076e23    / units.mol
    Kb=1.38064852e-23    * units.J / units.K
    e=1.60217662e-19 *units.C
    df=None
    kT=None
    Kw=None
    SEED=None
    rng=None


    class NumpyEncoder(json.JSONEncoder):
        """
        Custom JSON encoder that converts NumPy arrays to Python lists
        and NumPy scalars to Python scalars.
        """
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.generic):
                return obj.item()
            return super().default(obj)

    def __init__(self, SEED, temperature=None, unit_length=None, unit_charge=None, Kw=None):
        """
        Initializes the pymbe_library by setting up the reduced unit system with `temperature` and `reduced_length` 
        and sets up  the `pmb.df` for bookkeeping.

        Args:
            temperature(`pint.Quantity`,optional): Value of the temperature in the pyMBE UnitRegistry. Defaults to None.
            unit_length(`pint.Quantity`, optional): Value of the unit of length in the pyMBE UnitRegistry. Defaults to None.
            unit_charge (`pint.Quantity`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw (`pint.Quantity`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 
        
        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
        """
        # Seed and RNG
        self.SEED=SEED
        self.rng = np.random.default_rng(SEED)
        self.set_reduced_units(unit_length=unit_length, unit_charge=unit_charge,
                               temperature=temperature, Kw=Kw, verbose=False)
        self.setup_df()
        return

    def add_bond_in_df(self, particle_id1, particle_id2, use_default_bond=False):
        """
        Adds a bond entry on the `pymbe.df` storing the particle_ids of the two bonded particles.

        Args:
            particle_id1(`int`): particle_id of the type of the first particle type of the bonded particles
            particle_id2(`int`): particle_id of the type of the second particle type of the bonded particles
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle whose bond types are not defined in `pmb.df`. Defaults to False.

        Returns:
            index(`int`): Row index where the bond information has been added in pmb.df.
        """
        particle_name1 = self.df.loc[(self.df['particle_id']==particle_id1) & (self.df['pmb_type']=="particle")].name.values[0]
        particle_name2 = self.df.loc[(self.df['particle_id']==particle_id2) & (self.df['pmb_type']=="particle")].name.values[0]
        
        bond_key = self.find_bond_key(particle_name1=particle_name1,
                                    particle_name2=particle_name2, 
                                    use_default_bond=use_default_bond)
        if not bond_key:
            return
        self.copy_df_entry(name=bond_key,column_name='particle_id2',number_of_copies=1)
        indexs = np.where(self.df['name']==bond_key)
        index_list = list (indexs[0])
        used_bond_df = self.df.loc[self.df['particle_id2'].notnull()]
        #without this drop the program crashes when dropping duplicates because the 'bond' column is a dict
        used_bond_df = used_bond_df.drop([('bond_object','')],axis =1 )
        used_bond_index = used_bond_df.index.to_list()
        for index in index_list:
            if index not in used_bond_index:
                self.clean_df_row(index=int(index))
                self.df.at[index,'particle_id'] = particle_id1
                self.df.at[index,'particle_id2'] = particle_id2
                break
        return index

    def add_bonds_to_espresso(self, espresso_system) :
        """
        Adds all bonds defined in `pmb.df` to `espresso_system`.

        Args:
            espresso_system(`espressomd.system.System`): system object of espressomd library
        """

        if 'bond' in self.df.values:
            bond_df = self.df.loc[self.df ['pmb_type'] == 'bond']
            bond_list = bond_df.bond_object.values.tolist()
            for bond in bond_list:
                espresso_system.bonded_inter.add(bond)
        else:
            print ('WARNING: There are no bonds defined in pymbe.df')
        
        return

    def add_value_to_df(self,index,key,new_value, verbose=True, non_standard_value=False, overwrite=False):
        """
        Adds a value to a cell in the `pmb.df` DataFrame.

        Args:
            index(`int`): index of the row to add the value to.
            key(`str`): the column label to add the value to.
            verbose(`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.
            non_standard_value(`bool`, optional): Switch to enable insertion of non-standard values, such as `dict` objects. Defaults to False.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False.
        """

        token = "#protected:"

        def protect(obj):
            if non_standard_value:
                return token + json.dumps(obj, cls=self.NumpyEncoder)
            return obj

        def deprotect(obj):
            if non_standard_value and isinstance(obj, str) and obj.startswith(token):
                return json.loads(obj.removeprefix(token))
            return obj

        # Make sure index is a scalar integer value
        index = int (index)
        assert isinstance(index, int), '`index` should be a scalar integer value.'
        idx = pd.IndexSlice
        if self.check_if_df_cell_has_a_value(index=index,key=key):
            old_value= self.df.loc[index,idx[key]]
            if protect(old_value) != protect(new_value):
                name=self.df.loc[index,('name','')]
                pmb_type=self.df.loc[index,('pmb_type','')]
                if verbose:
                    print(f"WARNING: you are attempting to redefine the properties of {name} of pmb_type {pmb_type}")    
                if overwrite and verbose:
                    print(f'WARNING: overwritting the value of the entry `{key}`: old_value = {old_value} new_value = {new_value}')
                if not overwrite:
                    if verbose:
                        print(f"WARNING: pyMBE has preserved of the entry `{key}`: old_value = {old_value}. If you want to overwrite it with new_value = {new_value}, activate the switch overwrite = True ")
                    return
        self.df.loc[index,idx[key]] = protect(new_value)
        if non_standard_value:
            self.df[key] = self.df[key].apply(deprotect)
        return
    
    def assign_molecule_id(self, name, molecule_index, pmb_type, used_molecules_id):
        """
        Assigns the `molecule_id` of the pmb object given by `pmb_type`
        
        Args:
            name(`str`): Label of the molecule type to be created. `name` must be defined in `pmb.df`
            pmb_type(`str`): pmb_object_type to assign the `molecule_id` 
            molecule_index(`int`): index of the current `pmb_object_type` to assign the `molecule_id`
            used_molecules_id(`lst`): list with the `molecule_id` values already used.
        
        Returns:
            molecule_id(`int`): Id of the molecule
        """

        self.clean_df_row(index=int(molecule_index))
        
        if self.df['molecule_id'].isnull().values.all():
            molecule_id = 0        
        else:
            # check if a residue is part of another molecule
            check_residue_name = self.df[self.df['residue_list'].astype(str).str.contains(name)]
            mol_pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]
            if not check_residue_name.empty and mol_pmb_type == pmb_type:
                for value in check_residue_name.index.to_list():                  
                    if value not in used_molecules_id:                              
                        molecule_id = self.df.loc[value].molecule_id.values[0]                    
                        break
            else:
                molecule_id = self.df['molecule_id'].max() +1

        self.add_value_to_df (key=('molecule_id',''),
                                index=int(molecule_index),
                                new_value=molecule_id, 
                                verbose=False)

        return molecule_id
    
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
        molecule_name = self.df.loc[(self.df['molecule_id']==molecule_id) & (self.df['pmb_type'].isin(["molecule","protein"]))].name.values[0]
        particle_id_list = self.get_particle_id_map(object_name=molecule_name)["all"]
        for pid in particle_id_list:
            for axis in axis_list:
                center_of_mass [axis] += espresso_system.part.by_id(pid).pos[axis]
        center_of_mass = center_of_mass / len(particle_id_list)
        return center_of_mass

    def calculate_HH(self, molecule_name, pH_list=None, pka_set=None):
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
        if pH_list is None:
            pH_list=np.linspace(2,12,50)
        if pka_set is None:
            pka_set=self.get_pka_set() 
        self.check_pka_set(pka_set=pka_set)
        charge_map = self.get_charge_map()
        Z_HH=[]
        for pH_value in pH_list:    
            Z=0
            index = self.df.loc[self.df['name'] == molecule_name].index[0].item() 
            residue_list = self.df.at [index,('residue_list','')]
            sequence = self.df.at [index,('sequence','')]
            if np.any(pd.isnull(sequence)):
                # Molecule has no sequence
                for residue in residue_list:
                    list_of_particles_in_residue = self.search_particles_in_residue(residue)
                    for particle in list_of_particles_in_residue:
                        if particle in pka_set.keys():
                            if pka_set[particle]['acidity'] == 'acidic':
                                psi=-1
                            elif pka_set[particle]['acidity']== 'basic':
                                psi=+1
                            else:
                                psi=0
                            Z+=psi/(1+10**(psi*(pH_value-pka_set[particle]['pka_value'])))                      
                Z_HH.append(Z)
            else:
                # Molecule has a sequence
                if not isinstance(sequence, list):
                    # If the df has been read by file, the sequence needs to be parsed.
                    sequence = self.parse_sequence_from_file(sequence=sequence)
                for name in sequence:
                    if name in pka_set.keys():
                        if pka_set[name]['acidity'] == 'acidic':
                            psi=-1
                        elif pka_set[name]['acidity']== 'basic':
                            psi=+1
                        else:
                            psi=0
                        Z+=psi/(1+10**(psi*(pH_value-pka_set[name]['pka_value'])))
                    else:
                        state_one_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
                        Z+=charge_map[state_one_type]
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

    def calculate_initial_bond_length(self, bond_object, bond_type, epsilon, sigma, cutoff, offset):
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

        epsilon_red=epsilon.to('reduced_energy').magnitude
        sigma_red=sigma.to('reduced_length').magnitude
        cutoff_red=cutoff.to('reduced_length').magnitude
        offset_red=offset.to('reduced_length').magnitude

        if bond_type == "harmonic":
            r_0 = bond_object.params.get('r_0')
            k = bond_object.params.get('k')
            l0 = scipy.optimize.minimize(lambda x: 0.5*k*(x-r_0)**2 + truncated_lj_potential(x, epsilon_red, sigma_red, cutoff_red, offset_red), x0=r_0).x
        elif bond_type == "FENE":
            r_0 = bond_object.params.get('r_0')
            k = bond_object.params.get('k')
            d_r_max = bond_object.params.get('d_r_max')
            l0 = scipy.optimize.minimize(lambda x: -0.5*k*(d_r_max**2)*np.log(1-((x-r_0)/d_r_max)**2) + truncated_lj_potential(x, epsilon_red, sigma_red, cutoff_red,offset_red), x0=1.0).x
        return l0

    def calculate_net_charge (self, espresso_system, molecule_name):
        '''
        Calculates the net charge per molecule of molecules with `name` = molecule_name. 
        Returns the net charge per molecule and a maps with the net charge per residue and molecule.

        Args:
            espresso_system(`espressomd.system.System`): system information 
            molecule_name(`str`): name of the molecule to calculate the net charge

        Returns:
            (`dict`): {"mean": mean_net_charge, "molecules": {mol_id: net_charge_of_mol_id, }, "residues": {res_id: net_charge_of_res_id, }}

        Note:
            - The net charge of the molecule is averaged over all molecules of type `name` 
            - The net charge of each particle type is averaged over all particle of the same type in all molecules of type `name`
        '''        
        valid_pmb_types = ["molecule", "protein"]
        pmb_type=self.df.loc[self.df['name']==molecule_name].pmb_type.values[0]
        if pmb_type not in valid_pmb_types:
            raise ValueError("The pyMBE object with name {molecule_name} has a pmb_type {pmb_type}. This function only supports pyMBE types {valid_pmb_types}")      

        id_map = self.get_particle_id_map(object_name=molecule_name)
        def create_charge_map(espresso_system,id_map,label):
            charge_map={}
            for super_id in id_map[label].keys():
                net_charge=0
                for pid in id_map[label][super_id]:
                    net_charge+=espresso_system.part.by_id(pid).q
                charge_map[super_id]=net_charge
            return charge_map
        net_charge_molecules=create_charge_map(label="molecule_map",
                                                espresso_system=espresso_system,
                                                id_map=id_map)
        net_charge_residues=create_charge_map(label="residue_map",
                                                espresso_system=espresso_system,
                                                id_map=id_map)       
        mean_charge=np.mean(np.array(list(net_charge_molecules.values())))
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

    def check_aminoacid_key(self, key):
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
                       'C'] #'CYS'
        if key in valid_AA_keys:
            return True
        else:
            return False

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

    def check_if_df_cell_has_a_value(self, index,key):
        """
        Checks if a cell in the `pmb.df` at the specified index and column has a value.

        Args:
            index(`int`): Index of the row to check.
            key(`str`): Column label to check.

        Returns:
            `bool`: `True` if the cell has a value, `False` otherwise.
        """
        idx = pd.IndexSlice
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return not pd.isna(self.df.loc[index, idx[key]])

    def check_if_name_is_defined_in_df(self, name, pmb_type_to_be_defined):
        """
        Checks if `name` is defined in `pmb.df`.

        Args:
            name(`str`): label to check if defined in `pmb.df`.
            pmb_type_to_be_defined(`str`): pmb object type corresponding to `name`.

        Returns:
            `bool`: `True` for success, `False` otherwise.
        """
        if name in self.df['name'].unique():
            current_object_type = self.df[self.df['name']==name].pmb_type.values[0]
            if current_object_type != pmb_type_to_be_defined:
                raise ValueError (f"The name {name} is already defined in the df with a pmb_type = {current_object_type}, pymMBE does not support objects with the same name but different pmb_types")
            return True            
        else:
            return False

    def check_if_metal_ion(self,key):
        """
        Checks if `key` corresponds to a label of a supported metal ion.

        Args:
            key(`str`): key to be checked

        Returns:
            (`bool`): True if `key`  is a supported metal ion, False otherwise.
        """
        if key in self.get_metal_ions_charge_map().keys():
            return True
        else:
            return False

    def check_pka_set(self, pka_set):
        """
        Checks that `pka_set` has the formatting expected by the pyMBE library.
       
        Args:
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_entry in pka_set.values():
                if required_key not in pka_entry.keys():
                    raise ValueError(f'missing a required key "{required_key}" in the following entry of pka_set "{pka_entry}"')
        return

    def clean_df_row(self, index, columns_keys_to_clean=("particle_id", "particle_id2", "residue_id", "molecule_id")):
        """
        Cleans the columns of `pmb.df` in `columns_keys_to_clean` of the row with index `index` by assigning them a np.nan value.

        Args:
            index(`int`): Index of the row to clean.
            columns_keys_to_clean(`list` of `str`, optional): List with the column keys to be cleaned. Defaults to [`particle_id`, `particle_id2`, `residue_id`, `molecule_id`].
        """   
        for column_key in columns_keys_to_clean:
            self.add_value_to_df(key=(column_key,''),index=index,new_value=np.nan, verbose=False)
        return

    def convert_columns_to_original_format (self, df):
        """
        Converts the columns of the Dataframe to the original format in pyMBE.
        
        Args:
            df(`DataFrame`): dataframe with pyMBE information as a string  
        
        """

        columns_dtype_int = ['particle_id','particle_id2', 'residue_id','molecule_id', 'model',('state_one','es_type'),('state_two','es_type'),('state_one','charge'),('state_two','charge') ]  

        columns_with_units = ['sigma', 'epsilon', 'cutoff', 'offset']

        columns_with_list_or_dict = ['residue_list','side_chains', 'parameters_of_the_potential','sequence']

        for column_name in columns_dtype_int:
            df[column_name] = df[column_name].astype(object)
            
        for column_name in columns_with_list_or_dict:
            if df[column_name].isnull().all():
                df[column_name] = df[column_name].astype(object)
            else:
                df[column_name] = df[column_name].apply(lambda x: ast.literal_eval(str(x)) if pd.notnull(x) else x)

        for column_name in columns_with_units:
            df[column_name] = df[column_name].apply(lambda x: self.create_variable_with_units(x) if pd.notnull(x) else x)

        df['bond_object'] = df['bond_object'].apply(lambda x: self.convert_str_to_bond_object(x) if pd.notnull(x) else x)

        return df
    
    def convert_str_to_bond_object (self, bond_str):
        
        """
        Convert a row read as a `str` to the corresponding bond object. There are two supported bonds: HarmonicBond and FeneBond

        Args:
            bond_str(`str`): string with the information of a bond object

        Returns:
            bond_object(`obj`): EsPRESSo bond object 
        """
        
        from espressomd.interactions import HarmonicBond
        from espressomd.interactions import FeneBond

        supported_bonds = ['HarmonicBond', 'FeneBond']

        for bond in supported_bonds:
            variable = re.subn(f'{bond}', '', bond_str)
            if variable[1] == 1: 
                params = ast.literal_eval(variable[0])
                if bond == 'HarmonicBond':
                    bond_object = HarmonicBond(r_cut =params['r_cut'], k = params['k'], r_0=params['r_0'])
                elif bond == 'FeneBond':
                    bond_object = FeneBond(k = params['k'], d_r_max =params['d_r_max'], r_0=params['r_0'])

        return bond_object 

    def copy_df_entry(self, name, column_name, number_of_copies):
        '''
        Creates 'number_of_copies' of a given 'name' in `pymbe.df`.

        Args:
            name(`str`): Label of the particle/residue/molecule type to be created. `name` must be defined in `pmb.df`
            column_name(`str`): Column name to use as a filter. 
            number_of_copies(`int`): number of copies of `name` to be created.
        
        Note:
            - Currently, column_name only supports "particle_id", "particle_id2", "residue_id" and "molecule_id" 
        '''

        valid_column_names=["particle_id", "residue_id", "molecule_id", "particle_id2" ]
        if column_name not in valid_column_names:
            raise ValueError(f"{column_name} is not a valid column_name, currently only the following are supported: {valid_column_names}")
        df_by_name = self.df.loc[self.df.name == name]
        if number_of_copies != 1:           
            if df_by_name[column_name].isnull().values.any():       
                df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies-1), ignore_index=True)
            else:
                df_by_name = df_by_name[df_by_name.index == df_by_name.index.min()] 
                df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies), ignore_index=True)
                df_by_name_repeated[column_name] =np.NaN
            # Concatenate the new particle rows to  `df`
            self.df = pd.concat ([self.df,df_by_name_repeated], ignore_index=True)
        else:
            if not df_by_name[column_name].isnull().values.any():     
                df_by_name = df_by_name[df_by_name.index == df_by_name.index.min()] 
                df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies), ignore_index=True)
                df_by_name_repeated[column_name] =np.NaN
                self.df = pd.concat ([self.df,df_by_name_repeated], ignore_index=True)
        return

    def create_added_salt (self, espresso_system, cation_name, anion_name, c_salt, verbose=True):    
        """
        Creates a `c_salt` concentration of `cation_name` and `anion_name` ions into the `espresso_system`.

        Args:
            espresso_system(`espressomd.system.System`): instance of an espresso system object.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            c_salt(`float`): Salt concentration.
            verbose(`bool`): switch to activate/deactivate verbose. Defaults to True.
            
        Returns:
            c_salt_calculated(`float`): Calculated salt concentration added to `espresso_system`.
        """
        cation_name_charge = self.df.loc[self.df['name']==cation_name].state_one.charge.values[0]
        anion_name_charge = self.df.loc[self.df['name']==anion_name].state_one.charge.values[0]     
        if cation_name_charge <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ',cation_name_charge)
        if anion_name_charge >= 0:
            raise ValueError('ERROR anion charge must be negative, charge ', anion_name_charge)
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
        N_cation = N_ions*abs(anion_name_charge)
        N_anion = N_ions*abs(cation_name_charge)
        self.create_particle(espresso_system=espresso_system, name=cation_name, number_of_particles=N_cation)
        self.create_particle(espresso_system=espresso_system, name=anion_name, number_of_particles=N_anion)
        if verbose:
            if c_salt_calculated.check('[substance] [length]**-3'):
                print(f"\n Added salt concentration of {c_salt_calculated.to('mol/L')} given by {N_cation} cations and {N_anion} anions")
            elif c_salt_calculated.check('[length]**-3'):
                print(f"\n Added salt concentration of {c_salt_calculated.to('reduced_length**-3')} given by {N_cation} cations and {N_anion} anions")
        return c_salt_calculated

    def create_bond_in_espresso(self, bond_type, bond_parameters):
        '''
        Creates either a harmonic or a FENE bond in ESPREesSo

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.

        Note:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain:

                - k (`obj`)      : Magnitude of the bond. It should have units of energy/length**2 
                using the `pmb.units` UnitRegistry.
                - r_0 (`obj`)    : Equilibrium bond length. It should have units of length using 
                the `pmb.units` UnitRegistry.
           
            For a FENE bond the dictionay must additionally contain:
                
                - d_r_max (`obj`): Maximal stretching length for FENE. It should have 
                units of length using the `pmb.units` UnitRegistry. Default 'None'.

        Returns:
              bond_object (`obj`): a harmonic or a FENE bond object in ESPREesSo
        '''
        from espressomd import interactions

        valid_bond_types   = ["harmonic", "FENE"]
        
        if 'k' in bond_parameters:
            bond_magnitude     = bond_parameters['k'].to('reduced_energy / reduced_length**2')
        else:
            raise ValueError("Magnitud of the potential (k) is missing")
        
        if bond_type == 'harmonic':
            if 'r_0' in bond_parameters:
                bond_length        = bond_parameters['r_0'].to('reduced_length')
            else:
                raise ValueError("Equilibrium bond length (r_0) is missing")
            bond_object    = interactions.HarmonicBond(k   = bond_magnitude.magnitude,
                                                       r_0 = bond_length.magnitude)
        elif bond_type == 'FENE':
            if 'r_0' in bond_parameters:
                bond_length        = bond_parameters['r_0'].to('reduced_length').magnitude
            else:
                print("WARNING: No value provided for r_0. Defaulting to r_0 = 0")
                bond_length=0
            if 'd_r_max' in bond_parameters:
                max_bond_stret = bond_parameters['d_r_max'].to('reduced_length')
            else:
                raise ValueError("Maximal stretching length (d_r_max) is missing")
            bond_object    = interactions.FeneBond(r_0     = bond_length, 
                                                   k       = bond_magnitude.magnitude,
                                                   d_r_max = max_bond_stret.magnitude)
        else:
            raise ValueError(f"Bond type {bond_type} currently not implemented in pyMBE, valid types are {valid_bond_types}")

        return bond_object


    def create_counterions(self, object_name, cation_name, anion_name, espresso_system,verbose=True):
        """
        Creates particles of `cation_name` and `anion_name` in `espresso_system` to counter the net charge of `pmb_object`.
        
        Args:
            object_name(`str`): `name` of a pymbe object.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            verbose(`bool`): switch to activate/deactivate verbose. Defaults to True.

        Returns: 
            counterion_number(`dict`): {"name": number}
        """
        cation_charge = self.df.loc[self.df['name']==cation_name].state_one.charge.iloc[0]
        anion_charge = self.df.loc[self.df['name']==anion_name].state_one.charge.iloc[0]
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
            self.create_particle(espresso_system=espresso_system, name=cation_name, number_of_particles=counterion_number[cation_name])
        else:
            counterion_number[cation_name]=0
        if counterion_number[anion_name] > 0:
            self.create_particle(espresso_system=espresso_system, name=anion_name, number_of_particles=counterion_number[anion_name])
        else:
            counterion_number[anion_name] = 0
        if verbose:
            print('The following counter-ions have been created: ')
            for name in counterion_number.keys():
                print(f'Ion type: {name} created number: {counterion_number[name]}')
        return counterion_number
        
    def create_molecule(self, name, number_of_molecules, espresso_system, list_of_first_residue_positions=None, use_default_bond=False):
        """
        Creates `number_of_molecules` molecule of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the molecule type to be created. `name` must be defined in `pmb.df`
            espresso_system(`espressomd.system.System`): Instance of a system object from espressomd library.
            number_of_molecules(`int`): Number of molecules of type `name` to be created.
            list_of_first_residue_positions(`list`, optional): List of coordinates where the central bead of the first_residue_position will be created, random by default
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle with undefined bonds in `pymbe.df`

        Returns:
            molecules_info(`dict`):  {molecule_id: {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids": [particle_id1, ...]}}} 
        """
        if list_of_first_residue_positions is not None:
            for item in list_of_first_residue_positions:
                if not isinstance(item, list):
                    raise ValueError("The provided input position is not a nested list. Should be a nested list with elements of 3D lists, corresponding to xyz coord.")
                elif len(item) != 3:
                    raise ValueError("The provided input position is formatted wrong. The elements in the provided list does not have 3 coordinates, corresponding to xyz coord.")

            if len(list_of_first_residue_positions) != number_of_molecules:
                raise ValueError(f"Number of positions provided in {list_of_first_residue_positions} does not match number of molecules desired, {number_of_molecules}")
        if number_of_molecules <= 0:
            return 0
        self.check_if_name_is_defined_in_df(name=name,
                                            pmb_type_to_be_defined='molecule')
            
        first_residue = True
        molecules_info = {}
        residue_list = self.df[self.df['name']==name].residue_list.values [0]
        self.copy_df_entry(name=name,
                        column_name='molecule_id',
                        number_of_copies=number_of_molecules)

        molecules_index = np.where(self.df['name']==name)
        molecule_index_list =list(molecules_index[0])[-number_of_molecules:]
        used_molecules_id = self.df.molecule_id.dropna().drop_duplicates().tolist()
        pos_index = 0 
        for molecule_index in molecule_index_list:        
            molecule_id = self.assign_molecule_id(name=name,
                                                pmb_type='molecule',
                                                used_molecules_id=used_molecules_id,
                                                molecule_index=molecule_index)
            molecules_info[molecule_id] = {}
            for residue in residue_list:
                if first_residue:
                    if list_of_first_residue_positions is None:
                        residue_position = None
                    else:
                        for item in list_of_first_residue_positions:
                            residue_position = [np.array(list_of_first_residue_positions[pos_index])]
                    # Generate an arbitrary random unit vector
                    backbone_vector = self.generate_random_points_in_a_sphere(center=[0,0,0], 
                                                                radius=1, 
                                                                n_samples=1,
                                                                on_surface=True)[0]
                    residues_info = self.create_residue(name=residue,
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=residue_position,  
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    residue_id = next(iter(residues_info))
                    # Add the correct molecule_id to all particles in the residue
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        self.add_value_to_df(key=('molecule_id',''),
                                            index=int (index),
                                            new_value=molecule_id,
                                            overwrite=True,
                                            verbose=False)
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    previous_residue = residue
                    residue_position = espresso_system.part.by_id(central_bead_id).pos
                    previous_residue_id = central_bead_id
                    first_residue = False          
                else:                    
                    previous_central_bead_name=self.df[self.df['name']==previous_residue].central_bead.values[0]
                    new_central_bead_name=self.df[self.df['name']==residue].central_bead.values[0]       
                    bond = self.search_bond(particle_name1=previous_central_bead_name, 
                                            particle_name2=new_central_bead_name, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)                
                    l0 = self.get_bond_length(particle_name1=previous_central_bead_name, 
                                            particle_name2=new_central_bead_name, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)                
                    residue_position = residue_position+backbone_vector*l0
                    residues_info = self.create_residue(name=residue, 
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=[residue_position],
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    residue_id = next(iter(residues_info))      
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        self.add_value_to_df(key=('molecule_id',''),
                                            index=int (index),
                                            new_value=molecule_id,
                                            verbose=False,
                                            overwrite=True)            
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, previous_residue_id))
                    bond_index = self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=previous_residue_id,
                                        use_default_bond=use_default_bond) 
                    self.add_value_to_df(key=('molecule_id',''),
                                            index=int (bond_index),
                                            new_value=molecule_id,
                                            verbose=False,
                                            overwrite=True)           
                    previous_residue_id = central_bead_id
                    previous_residue = residue                    
                molecules_info[molecule_id][residue_id] = residues_info[residue_id]
            first_residue = True
            pos_index+=1
        
        return molecules_info
    
    def create_particle(self, name, espresso_system, number_of_particles, position=None, fix=False):
        """
        Creates `number_of_particles` particles of type `name` into `espresso_system` and bookkeeps them into `pymbe.df`.
        
        Args:
            name(`str`): Label of the particle type to be created. `name` must be a `particle` defined in `pmb_df`. 
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            number_of_particles(`int`): Number of particles to be created.
            position(list of [`float`,`float`,`float`], optional): Initial positions of the particles. If not given, particles are created in random positions. Defaults to None.
            fix(`bool`, optional): Controls if the particle motion is frozen in the integrator, it is used to create rigid objects. Defaults to False.
        Returns:
            created_pid_list(`list` of `float`): List with the ids of the particles created into `espresso_system`.
        """       
        if number_of_particles <=0:
            return 0
        self.check_if_name_is_defined_in_df(name=name,
                                       pmb_type_to_be_defined='particle')
        # Copy the data of the particle `number_of_particles` times in the `df`
        self.copy_df_entry(name=name,
                          column_name='particle_id',
                          number_of_copies=number_of_particles)
        # Get information from the particle type `name` from the df     
        q = self.df.loc[self.df['name']==name].state_one.charge.values[0]
        es_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
        # Get a list of the index in `df` corresponding to the new particles to be created
        index = np.where(self.df['name']==name)
        index_list =list(index[0])[-number_of_particles:]
        # Create the new particles into  `espresso_system`
        created_pid_list=[]
        for index in range (number_of_particles):
            df_index=int (index_list[index])
            self.clean_df_row(index=df_index)
            if position is None:
                particle_position = self.rng.random((1, 3))[0] *np.copy(espresso_system.box_l)
            else:
                particle_position = position[index]
            if len(espresso_system.part.all()) == 0:
                bead_id = 0
            else:
                bead_id = max (espresso_system.part.all().id) + 1
            created_pid_list.append(bead_id)
            
            if fix:
                espresso_system.part.add (id=bead_id, pos = particle_position, type = es_type, q = q,fix =[fix,fix,fix])        
            else:
                espresso_system.part.add (id=bead_id, pos = particle_position, type = es_type, q = q)        
            self.add_value_to_df(key=('particle_id',''),index=df_index,new_value=bead_id, verbose=False)                  
        return created_pid_list

    def create_pmb_object(self, name, number_of_objects, espresso_system, position=None, use_default_bond=False):
        """
        Creates all `particle`s associated to `pmb object` into  `espresso` a number of times equal to `number_of_objects`.
        
        Args:
            name(`str`): Unique label of the `pmb object` to be created. 
            number_of_objects(`int`): Number of `pmb object`s to be created.
            espresso_system(`espressomd.system.System`): Instance of an espresso system object from espressomd library.
            position(`list`): Coordinates where the particles should be created.
            use_default_bond(`bool`,optional): Controls if a `default` bond is used to bond particles with undefined bonds in `pmb.df`. Defaults to `False`.

        Note:
            - If no `position` is given, particles will be created in random positions. For bonded particles, they will be created at a distance equal to the bond length. 
        """
        allowed_objects=['particle','residue','molecule']
        pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]
        if pmb_type not in allowed_objects:
            raise ValueError('Object type not supported, supported types are ', allowed_objects)
        if pmb_type == 'particle':
            self.create_particle(name=name, 
                                number_of_particles=number_of_objects, 
                                espresso_system=espresso_system, 
                                position=position)
        elif pmb_type == 'residue':
            self.create_residue(name=name,  
                                espresso_system=espresso_system, 
                                central_bead_position=position,
                                use_default_bond=use_default_bond)
        elif pmb_type == 'molecule':
            self.create_molecule(name=name, 
                                number_of_molecules=number_of_objects, 
                                espresso_system=espresso_system, 
                                use_default_bond=use_default_bond, 
                                list_of_first_residue_positions=position)
        return

    def create_protein(self, name, number_of_proteins, espresso_system, topology_dict):
        """
        Creates `number_of_proteins` molecules of type `name` into `espresso_system` at the coordinates in `positions`

        Args:
            name(`str`): Label of the protein to be created. 
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            number_of_proteins(`int`): Number of proteins to be created.
            positions(`dict`): {'ResidueNumber': {'initial_pos': [], 'chain_id': ''}}
        """

        if number_of_proteins <=0:
            return
        self.check_if_name_is_defined_in_df(name=name,
                                        pmb_type_to_be_defined='protein')
        self.copy_df_entry(name=name,
                            column_name='molecule_id',
                            number_of_copies=number_of_proteins)
        protein_index = np.where(self.df['name']==name)
        protein_index_list =list(protein_index[0])[-number_of_proteins:]
        used_molecules_id = self.df.molecule_id.dropna().drop_duplicates().tolist()
        
        box_half=espresso_system.box_l[0]/2.0

        for molecule_index in protein_index_list:     

            molecule_id = self.assign_molecule_id (name=name,pmb_type='protein',used_molecules_id=used_molecules_id,molecule_index=molecule_index)

            protein_center = self.generate_coordinates_outside_sphere(radius = 1, 
                                                                        max_dist=box_half, 
                                                                        n_samples=1, 
                                                                        center=[box_half]*3)[0]
   
            for residue in topology_dict.keys():

                residue_name = re.split(r'\d+', residue)[0]
                residue_number = re.split(r'(\d+)', residue)[1]
                residue_position = topology_dict[residue]['initial_pos']
                position = residue_position + protein_center

                particle_id = self.create_particle(name=residue_name,
                                                            espresso_system=espresso_system,
                                                            number_of_particles=1,
                                                            position=[position], 
                                                            fix = True)
                
                index = self.df[self.df['particle_id']==particle_id[0]].index.values[0]
                
                self.add_value_to_df(key=('residue_id',''),
                                            index=int (index),
                                            new_value=residue_number,
                                            overwrite=True,
                                            verbose=False)

                self.add_value_to_df(key=('molecule_id',''),
                                        index=int (index),
                                        new_value=molecule_id,
                                        overwrite=True,
                                        verbose=False)

        return

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
            residues_info(`dict`): {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
        """
        self.check_if_name_is_defined_in_df(name=name,
                                            pmb_type_to_be_defined='residue')
        # Copy the data of a residue in the `df
        self.copy_df_entry(name=name,
                            column_name='residue_id',
                            number_of_copies=1)
        residues_index = np.where(self.df['name']==name)
        residue_index_list =list(residues_index[0])[-1:]
        for residue_index in residue_index_list:  
            side_chain_list = self.df.loc[self.df.index[residue_index]].side_chains.values[0]
            for side_chain_element in side_chain_list:
                if side_chain_element not in self.df.values:              
                    raise ValueError (side_chain_element +'is not defined')
        # Internal bookkepping of the residue info (important for side-chain residues)
        # Dict structure {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
        residues_info={}
        for residue_index in residue_index_list:     
            self.clean_df_row(index=int(residue_index))
            # Assign a residue_id
            if self.df['residue_id'].isnull().all():
                residue_id=0
            else:
                residue_id = self.df['residue_id'].max() + 1
            self.add_value_to_df(key=('residue_id',''),index=int (residue_index),new_value=residue_id, verbose=False)
            # create the principal bead   
            central_bead_name = self.df.loc[self.df['name']==name].central_bead.values[0]            
            central_bead_id = self.create_particle(name=central_bead_name,
                                                                espresso_system=espresso_system,
                                                                position=central_bead_position,
                                                                number_of_particles = 1)[0]
            central_bead_position=espresso_system.part.by_id(central_bead_id).pos
            #assigns same residue_id to the central_bead particle created.
            index = self.df[self.df['particle_id']==central_bead_id].index.values[0]
            self.df.at [index,'residue_id'] = residue_id
            # Internal bookkeeping of the central bead id
            residues_info[residue_id]={}
            residues_info[residue_id]['central_bead_id']=central_bead_id
            # create the lateral beads  
            side_chain_list = self.df.loc[self.df.index[residue_index]].side_chains.values[0]
            side_chain_beads_ids = []
            for side_chain_element in side_chain_list:
                
                pmb_type = self.df[self.df['name']==side_chain_element].pmb_type.values[0] 
                if pmb_type == 'particle':
                    bond = self.search_bond(particle_name1=central_bead_name, 
                                            particle_name2=side_chain_element, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)
                    l0 = self.get_bond_length(particle_name1=central_bead_name, 
                                              particle_name2=side_chain_element, 
                                              hard_check=True, 
                                              use_default_bond=use_default_bond)

                    if backbone_vector is None:
                        bead_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                 radius=l0, 
                                                                 n_samples=1,
                                                                 on_surface=True)[0]
                    else:
                        bead_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                                    magnitude=l0)
                    
                    side_bead_id = self.create_particle(name=side_chain_element, 
                                                                    espresso_system=espresso_system,
                                                                    position=[bead_position], 
                                                                    number_of_particles=1)[0]
                    index = self.df[self.df['particle_id']==side_bead_id].index.values[0]
                    self.add_value_to_df(key=('residue_id',''),
                                        index=int (index),
                                        new_value=residue_id, 
                                        verbose=False,
                                        overwrite=True)
                    side_chain_beads_ids.append(side_bead_id)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, side_bead_id))
                    index = self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=side_bead_id,
                                        use_default_bond=use_default_bond)
                    self.add_value_to_df(key=('residue_id',''),
                                        index=int (index),
                                        new_value=residue_id, 
                                        verbose=False,
                                        overwrite=True)

                elif pmb_type == 'residue':
                    central_bead_side_chain = self.df[self.df['name']==side_chain_element].central_bead.values[0]
                    bond = self.search_bond(particle_name1=central_bead_name, 
                                            particle_name2=central_bead_side_chain, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)
                    l0 = self.get_bond_length(particle_name1=central_bead_name, 
                                              particle_name2=central_bead_side_chain, 
                                              hard_check=True, 
                                              use_default_bond=use_default_bond)
                    if backbone_vector is None:
                        residue_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                    radius=l0, 
                                                                    n_samples=1,
                                                                    on_surface=True)[0]
                    else:
                        residue_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                                        magnitude=l0)
                    lateral_residue_info = self.create_residue(name=side_chain_element, 
                                                                espresso_system=espresso_system,
                                                                central_bead_position=[residue_position],
                                                                use_default_bond=use_default_bond)
                    lateral_residue_dict=list(lateral_residue_info.values())[0]
                    central_bead_side_chain_id=lateral_residue_dict['central_bead_id']
                    lateral_beads_side_chain_ids=lateral_residue_dict['side_chain_ids']
                    residue_id_side_chain=list(lateral_residue_info.keys())[0]
                    # Change the residue_id of the residue in the side chain to the one of the bigger residue
                    index = self.df[(self.df['residue_id']==residue_id_side_chain) & (self.df['pmb_type']=='residue') ].index.values[0]
                    self.add_value_to_df(key=('residue_id',''),
                                        index=int(index),
                                        new_value=residue_id, 
                                        verbose=False,
                                        overwrite=True)
                    # Change the residue_id of the particles in the residue in the side chain
                    side_chain_beads_ids+=[central_bead_side_chain_id]+lateral_beads_side_chain_ids
                    for particle_id in side_chain_beads_ids:
                        index = self.df[(self.df['particle_id']==particle_id) & (self.df['pmb_type']=='particle')].index.values[0]
                        self.add_value_to_df(key=('residue_id',''),
                                            index=int (index),
                                            new_value=residue_id, 
                                            verbose=False,
                                            overwrite=True)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, central_bead_side_chain_id))
                    index = self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=central_bead_side_chain_id,
                                        use_default_bond=use_default_bond)
                    self.add_value_to_df(key=('residue_id',''),
                                        index=int (index),
                                        new_value=residue_id, 
                                        verbose=False,
                                        overwrite=True)
                    # Change the residue_id of the bonds in the residues in the side chain to the one of the bigger residue
                    for index in self.df[(self.df['residue_id']==residue_id_side_chain) & (self.df['pmb_type']=='bond') ].index:        
                        self.add_value_to_df(key=('residue_id',''),
                                            index=int(index),
                                            new_value=residue_id, 
                                            verbose=False,
                                            overwrite=True)
            # Internal bookkeeping of the side chain beads ids
            residues_info[residue_id]['side_chain_ids']=side_chain_beads_ids
        return  residues_info

    def create_variable_with_units(self, variable):
        """
        Returns a pint object with the value and units defined in `variable`.

        Args:
            variable(`dict` or `str`): {'value': value, 'units': units}
        Returns:
            variable_with_units(`obj`): variable with units using the pyMBE UnitRegistry.
        """        
        
        if isinstance(variable, dict):

            value=variable.pop('value')
            units=variable.pop('units')

        elif isinstance(variable, str):

            value = float(re.split(r'\s+', variable)[0])
            units = re.split(r'\s+', variable)[1]
        
        variable_with_units=value*self.units(units)

        return variable_with_units 

    def define_AA_residues(self, sequence, model):
        """
        Defines in `pmb.df` all the different residues in `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein.
            model(`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.

        Returns:
            residue_list(`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.             
        """
        residue_list = []
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
            if residue_name not in residue_list:   
                self.define_residue(name = 'AA-'+residue_name, 
                                    central_bead = central_bead,
                                    side_chains = side_chains)              
            residue_list.append('AA-'+residue_name)
        return residue_list

    def define_bond(self, bond_type, bond_parameters, particle_pairs):
        
        '''
        Defines a pmb object of type `bond` in `pymbe.df`.

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.
            particle_pairs(`lst`): list of the `names` of the `particles` to be bonded.

        Note:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain the following parameters:

                - k (`obj`)      : Magnitude of the bond. It should have units of energy/length**2 
                using the `pmb.units` UnitRegistry.
                - r_0 (`obj`)    : Equilibrium bond length. It should have units of length using 
                the `pmb.units` UnitRegistry.
           
            For a FENE bond the dictionay must contain the same parameters as for a HARMONIC bond and:
                
                - d_r_max (`obj`): Maximal stretching length for FENE. It should have 
                units of length using the `pmb.units` UnitRegistry. Default 'None'.
        '''

        bond_object=self.create_bond_in_espresso(bond_type, bond_parameters)

        for particle_name1, particle_name2 in particle_pairs:

            lj_parameters=self.get_lj_parameters(particle_name1 = particle_name1,
                                                 particle_name2 = particle_name2,
                                                 combining_rule = 'Lorentz-Berthelot')

            l0 = self.calculate_initial_bond_length(bond_object = bond_object,
                                                    bond_type   = bond_type,
                                                    epsilon     = lj_parameters["epsilon"],
                                                    sigma       = lj_parameters["sigma"],
                                                    cutoff      = lj_parameters["cutoff"],
                                                    offset      = lj_parameters["offset"],)
            index = len(self.df)
            for label in [particle_name1+'-'+particle_name2,particle_name2+'-'+particle_name1]:
                self.check_if_name_is_defined_in_df(name=label, pmb_type_to_be_defined="bond")
            self.df.at [index,'name']= particle_name1+'-'+particle_name2
            self.df.at [index,'bond_object'] = bond_object
            self.df.at [index,'l0'] = l0
            self.add_value_to_df(index=index,
                                    key=('pmb_type',''),
                                    new_value='bond')
            self.add_value_to_df(index=index,
                                    key=('parameters_of_the_potential',''),
                                    new_value=bond_object.get_params(),
                                    non_standard_value=True)

        return

    
    def define_default_bond(self, bond_type, bond_parameters, epsilon=None, sigma=None, cutoff=None, offset=None):
        """
        Asigns `bond` in `pmb.df` as the default bond.
        The LJ parameters can be optionally provided to calculate the initial bond length

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.
            sigma(`float`, optional): LJ sigma of the interaction between the particles.
            epsilon(`float`, optional): LJ epsilon for the interaction between the particles.
            cutoff(`float`, optional): cutoff-radius of the LJ interaction.
            offset(`float`, optional): offset of the LJ interaction.

        Note:
            - Currently, only harmonic and FENE bonds are supported. 
        """
        
        bond_object=self.create_bond_in_espresso(bond_type, bond_parameters)

        if epsilon is None:
            epsilon=1*self.units('reduced_energy')
        if sigma is None:
            sigma=1*self.units('reduced_length')
        if cutoff is None:
            cutoff=2**(1.0/6.0)*self.units('reduced_length')
        if offset is None:
            offset=0*self.units('reduced_length')
        l0 = self.calculate_initial_bond_length(bond_object = bond_object, 
                                                bond_type   = bond_type, 
                                                epsilon     = epsilon,
                                                sigma       = sigma,
                                                cutoff      = cutoff,
                                                offset      = offset)

        if self.check_if_name_is_defined_in_df(name='default',pmb_type_to_be_defined='bond'):
            return
        if len(self.df.index) != 0:
            index = max(self.df.index)+1
        else:
            index = 0
        self.df.at [index,'name']        = 'default'
        self.df.at [index,'bond_object'] = bond_object
        self.df.at [index,'l0']          = l0
        self.add_value_to_df(index       = index,
                            key          = ('pmb_type',''),
                            new_value    = 'bond')
        self.add_value_to_df(index       = index,
                            key          = ('parameters_of_the_potential',''),
                            new_value    = bond_object.get_params(),
                            non_standard_value=True)
        return

    def define_molecule(self, name, residue_list):
        """
        Defines a pyMBE object of type `molecule` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `molecule`.
            residue_list(`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.  
        """
        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='molecule'):
            return
        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'molecule'
        self.df.at [index,('residue_list','')] = residue_list
        return

    def define_particle_entry_in_df(self,name):
        """
        Defines a particle entry in pmb.df.

        Args:
            name(`str`): Unique label that identifies this particle type.

        Returns:
            index(`int`): Index of the particle in pmb.df  
        """

        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='particle'):
            index = self.df[self.df['name']==name].index[0]                                   
        else:
            index = len(self.df)
            self.df.at [index, 'name'] = name
            self.df.at [index,'pmb_type'] = 'particle'
        return index

    def define_particle(self, name, q=0, acidity='inert', pka=None, sigma=None, epsilon=None, cutoff=None, offset=None,verbose=True,overwrite=False):
        """
        Defines the properties of a particle object.

        Args:
            name(`str`): Unique label that identifies this particle type.  
            q(`int`, optional): Permanent charge of this particle type. Defaults to 0.
            acidity(`str`, optional): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            pka(`float`, optional): If `particle` is an acid or a base, it defines its  pka-value. Defaults to None.
            sigma(`pint.Quantity`, optional): Sigma parameter used to set up Lennard-Jones interactions for this particle type. Defaults to None.
            cutoff(`pint.Quantity`, optional): Cutoff parameter used to set up Lennard-Jones interactions for this particle type. Defaults to None.
            offset(`pint.Quantity`, optional): Offset parameter used to set up Lennard-Jones interactions for this particle type. Defaults to None.
            epsilon(`pint.Quantity`, optional): Epsilon parameter used to setup Lennard-Jones interactions for this particle tipe. Defaults to None.
            verbose(`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False.

        Note:
            - `sigma`, `cutoff` and `offset` must have a dimensitonality of `[length]` and should be defined using pmb.units.
            - `epsilon` must have a dimensitonality of `[energy]` and should be defined using pmb.units.
            - `cutoff` defaults to `2**(1./6.) reduced_length`. 
            - `offset` defaults to 0.
            - The default setup corresponds to the Weeks−Chandler−Andersen (WCA) model, corresponding to purely steric interactions.
            - For more information on `sigma`, `epsilon`, `cutoff` and `offset` check `pmb.setup_lj_interactions()`.
        """ 
        index=self.define_particle_entry_in_df(name=name)
        
        # If `cutoff` and `offset` are not defined, default them to the following values
        if cutoff is None:
            cutoff=self.units.Quantity(2**(1./6.), "reduced_length")
        if offset is None:
            offset=self.units.Quantity(0, "reduced_length")

        # Define LJ parameters
        parameters_with_dimensionality={"sigma":{"value": sigma, "dimensionality": "[length]"},
                                        "cutoff":{"value": cutoff, "dimensionality": "[length]"},
                                        "offset":{"value": offset, "dimensionality": "[length]"},
                                        "epsilon":{"value": epsilon, "dimensionality": "[energy]"},}

        for parameter_key in parameters_with_dimensionality.keys():
            if parameters_with_dimensionality[parameter_key]["value"] is not None:
                self.check_dimensionality(variable=parameters_with_dimensionality[parameter_key]["value"], 
                                          expected_dimensionality=parameters_with_dimensionality[parameter_key]["dimensionality"])
                self.add_value_to_df(key=(parameter_key,''),
                                    index=index,
                                    new_value=parameters_with_dimensionality[parameter_key]["value"],
                                    verbose=verbose,
                                    overwrite=overwrite)

        # Define particle acid/base properties
        self.set_particle_acidity(name=name, 
                                acidity=acidity, 
                                default_charge=q, 
                                pka=pka,
                                verbose=verbose,
                                overwrite=overwrite)
        return 
    
    def define_particles(self, parameters, overwrite=False, verbose=True):
        '''
        Defines a particle object in pyMBE for each particle name in `particle_names`

        Args:
            parameters(`dict`):  dictionary with the particle parameters. 
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.

        Note:
            - parameters = {"particle_name1: {"sigma": sigma_value, "epsilon": epsilon_value, ...}, particle_name2: {...},}
        '''
        if not parameters:
            return 0
        for particle_name in parameters.keys():
            parameters[particle_name]["overwrite"]=overwrite
            parameters[particle_name]["verbose"]=verbose
            self.define_particle(**parameters[particle_name])
        return
    
    def define_peptide(self, name, sequence, model):
        """
        Defines a pyMBE object of type `peptide` in the `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the `peptide`.
            sequence (`string`): Sequence of the `peptide`.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
        """
        if not self.check_if_name_is_defined_in_df(name = name, pmb_type_to_be_defined='peptide'):
            valid_keys = ['1beadAA','2beadAA']
            if model not in valid_keys:
                raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
            clean_sequence = self.protein_sequence_parser(sequence=sequence)    
            residue_list = self.define_AA_residues(sequence=clean_sequence,
                                                    model=model)
            self.define_molecule(name = name, residue_list=residue_list)
            index = self.df.loc[self.df['name'] == name].index.item() 
            self.df.at [index,'model'] = model
            self.df.at [index,('sequence','')] = clean_sequence
        return
    
    def define_protein(self, name,model, topology_dict, lj_setup_mode="wca", overwrite=False, verbose=True):
        """
        Defines a globular protein pyMBE object  in `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the protein.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
            topology_dict (`dict`): {'initial_pos': coords_list, 'chain_id': id, 'radius': radius_value}
            lj_setup_mode(`str`): Key for the setup for the LJ potential. Defaults to "wca".
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.

        Note:
            - Currently, only `lj_setup_mode="wca"` is supported. This corresponds to setting up the WCA potential.
        """

        # Sanity checks
        valid_model_keys = ['1beadAA','2beadAA']
        valid_lj_setups = ["wca"]
        if model not in valid_model_keys:
            raise ValueError('Invalid key for the protein model, supported models are {valid_model_keys}')
        if lj_setup_mode not in valid_lj_setups:
            raise ValueError('Invalid key for the lj setup, supported setup modes are {valid_lj_setups}')
        if lj_setup_mode == "wca":
            sigma = 1*self.units.Quantity("reduced_length")
            epsilon = 1*self.units.Quantity("reduced_energy")
        part_dict={}
        sequence=[]
        metal_ions_charge_map=self.get_metal_ions_charge_map()
        for particle in topology_dict.keys():
            particle_name = re.split(r'\d+', particle)[0] 
            if particle_name not in part_dict.keys():
                if lj_setup_mode == "wca":
                    part_dict[particle_name]={"sigma": sigma,
                                        "offset": topology_dict[particle]['radius']*2-sigma,
                                        "epsilon": epsilon,
                                        "name": particle_name}
                if self.check_if_metal_ion(key=particle_name):
                    q=metal_ions_charge_map[particle_name]
                else:
                    q=0
                part_dict[particle_name]["q"]=q
            
            if self.check_aminoacid_key(key=particle_name):
                sequence.append(particle_name) 
            
        self.define_particles(parameters=part_dict,
                            overwrite=overwrite,  
                            verbose=verbose)
        residue_list = self.define_AA_residues(sequence=sequence, 
                                               model=model)
        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'protein'
        self.df.at [index,'model'] = model
        self.df.at [index,('sequence','')] = sequence  
        self.df.at [index,('residue_list','')] = residue_list    
        return 
    
    def define_residue(self, name, central_bead, side_chains):
        """
        Defines a pyMBE object of type `residue` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `residue`.
            central_bead(`str`): `name` of the `particle` to be placed as central_bead of the `residue`.
            side_chains(`list` of `str`): List of `name`s of the pmb_objects to be placed as side_chains of the `residue`. Currently, only pmb_objects of type `particle`s or `residue`s are supported.
        """
        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='residue'):
            return
        index = len(self.df)
        self.df.at [index, 'name'] = name
        self.df.at [index,'pmb_type'] = 'residue'
        self.df.at [index,'central_bead'] = central_bead
        self.df.at [index,('side_chains','')] = side_chains
        return 

    def destroy_pmb_object_in_system(self, name, espresso_system):
        """
        Destroys all particles associated with `name` in `espresso_system` amd removes the destroyed pmb_objects from `pmb.df` 

        Args:
            name(`str`): Label of the pmb object to be destroyed. The pmb object must be defined in `pymbe.df`.
            espresso_system(`espressomd.system.System`): Instance of a system class from espressomd library.

        Note:
            - If `name`  is a object_type=`particle`, only the matching particles that are not part of bigger objects (e.g. `residue`, `molecule`) will be destroyed. To destroy particles in such objects, destroy the bigger object instead.
        """
        allowed_objects = ['particle','residue','molecule']
        pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]
        if pmb_type not in allowed_objects:
            raise ValueError('Object type not supported, supported types are ', allowed_objects)
        if pmb_type == 'particle':
            particles_index = self.df.index[(self.df['name'] == name) & (self.df['residue_id'].isna()) 
                                                  & (self.df['molecule_id'].isna())]
            particle_ids_list= self.df.loc[(self.df['name'] == name) & (self.df['residue_id'].isna())
                                                 & (self.df['molecule_id'].isna())].particle_id.tolist()
            for particle_id in particle_ids_list:
                espresso_system.part.by_id(particle_id).remove()
            self.df = self.df.drop(index=particles_index)
        if pmb_type == 'residue':
            residues_id = self.df.loc[self.df['name']== name].residue_id.to_list()
            for residue_id in residues_id:
                molecule_name = self.df.loc[(self.df['residue_id']==molecule_id) & (self.df['pmb_type']=="residue")].name.values[0]
                particle_ids_list = self.get_particle_id_map(object_name=molecule_name)["all"]
                self.df = self.df.drop(self.df[self.df['residue_id'] == residue_id].index)
                for particle_id in particle_ids_list:
                    espresso_system.part.by_id(particle_id).remove()
                    self.df= self.df.drop(self.df[self.df['particle_id']==particle_id].index)    
        if pmb_type == 'molecule':
            molecules_id = self.df.loc[self.df['name']== name].molecule_id.to_list()
            for molecule_id in molecules_id:
                molecule_name = self.df.loc[(self.df['molecule_id']==molecule_id) & (self.df['pmb_type']=="molecule")].name.values[0]
                particle_ids_list = self.get_particle_id_map(object_name=molecule_name)["all"]
                self.df = self.df.drop(self.df[self.df['molecule_id'] == molecule_id].index)
                for particle_id in particle_ids_list:
                    espresso_system.part.by_id(particle_id).remove()   
                    self.df= self.df.drop(self.df[self.df['particle_id']==particle_id].index)             
        
        self.df.reset_index(drop=True,inplace=True)

        return

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
        print ('enable_motion_of_rigid_object requires that espressomd has the following features activated: ["VIRTUAL_SITES_RELATIVE", "MASS"]')
        pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]
        if pmb_type != 'protein':
            raise ValueError (f'The pmb_type: {pmb_type} is not currently supported. The supported pmb_type is: protein')
        molecule_ids_list = self.df.loc[self.df['name']==name].molecule_id.to_list()
        for molecule_id in molecule_ids_list:    
            particle_ids_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
            center_of_mass = self.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)
            rigid_object_center = espresso_system.part.add(pos=center_of_mass,
                                                           rotation=[True,True,True], 
                                                           type=self.propose_unused_type())
            for particle_id in particle_ids_list:
                pid = espresso_system.part.by_id(particle_id)
                pid.vs_auto_relate_to(rigid_object_center.id)
        return

    def filter_df(self, pmb_type):
        """
        Filters `pmb.df` and returns a sub-set of it containing only rows with pmb_object_type=`pmb_type` and non-NaN columns.
        
        Args:
            pmb_type(`str`): pmb_object_type to filter in `pmb.df`.

        Returns:
            pmb_type_df(`Pandas.Dataframe`): filtered `pmb.df`.
        """
        pmb_type_df = self.df.loc[self.df['pmb_type']== pmb_type]
        pmb_type_df = pmb_type_df.dropna( axis=1, thresh=1)
        return pmb_type_df

    def find_bond_key(self, particle_name1, particle_name2, use_default_bond=False):
        """
        Searches for the `name` of the bond between `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        
        Args:
            particle_name1(`str`): label of the type of the first particle type of the bonded particles.
            particle_name2(`str`): label of the type of the second particle type of the bonded particles.
            use_default_bond(`bool`, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to 'False'. 

        Returns:
            bond_key (str): `name` of the bond between `particle_name1` and `particle_name2` 

        Note:
            - If `use_default_bond`=`True`, it returns "default" if no key is found.
        """
        bond_keys = [particle_name1 +'-'+ particle_name2, particle_name2 +'-'+ particle_name1 ]
        bond_defined=False
        for bond_key in bond_keys:
            if bond_key in self.df.values:
                bond_defined=True
                correct_key=bond_key
                break
        if bond_defined:
            return correct_key
        elif use_default_bond:
            return 'default'
        else:
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
            index = np.where(self.df[(state, 'es_type')] == es_type)[0]
            if len(index) > 0:
                if column_name == 'label':
                    label = self.df.loc[idx[index[0]], idx[(state,column_name)]]
                    return label
                else: 
                    column_name_value = self.df.loc[idx[index[0]], idx[(column_name,'')]]
                    return column_name_value
        return None

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

        Note:
            - Algorithm from: https://baezortega.github.io/2018/10/14/hypersphere-sampling/
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
        np_vec /= np.linalg.norm(np_vec) 
        if np.all(np_vec == 0):
            raise ValueError('Zero vector')
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

    def get_bond_length(self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for bonds between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns the initial bond length.
        If `use_default_bond` is activated and a "default" bond is defined, returns the length of that default bond instead.
        If no bond is found, it prints a message and it does not return anything. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1(str): label of the type of the first particle type of the bonded particles.
            particle_name2(str): label of the type of the second particle type of the bonded particles.
            hard_check(bool, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond(bool, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            l0(`pint.Quantity`): bond length
        
        Note:
            - If `use_default_bond`=True and no bond is defined between `particle_name1` and `particle_name2`, it returns the default bond defined in `pmb.df`.
            - If `hard_check`=`True` stops the code when no bond is found.
        """
        bond_key = self.find_bond_key(particle_name1=particle_name1, 
                                    particle_name2=particle_name2, 
                                    use_default_bond=use_default_bond)
        if bond_key:
            return self.df[self.df['name']==bond_key].l0.values[0]
        else:
            print("Bond not defined between particles ", particle_name1, " and ", particle_name2)    
            if hard_check:
                sys.exit(1)
            else:
                return

    def get_charge_map(self):
        '''
        Gets the charge of each `espresso_type` in `pymbe.df`.
        
        Returns:
            charge_map(`dict`): {espresso_type: charge}.
        '''
        if self.df.state_one['es_type'].isnull().values.any():         
            df_state_one = self.df.state_one.dropna()     
            df_state_two = self.df.state_two.dropna()  
        else:    
            df_state_one = self.df.state_one
            if self.df.state_two['es_type'].isnull().values.any():
                df_state_two = self.df.state_two.dropna()   
            else:
                df_state_two = self.df.state_two
        state_one = pd.Series (df_state_one.charge.values,index=df_state_one.es_type.values)
        state_two = pd.Series (df_state_two.charge.values,index=df_state_two.es_type.values)
        charge_map  = pd.concat([state_one,state_two],axis=0).to_dict()
        return charge_map

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
        lj_parameters_keys=["sigma","epsilon","offset","cutoff"]
        if combining_rule not in supported_combining_rules:
            raise ValueError(f"Combining_rule {combining_rule} currently not implemented in pyMBE, valid keys are {supported_combining_rules}")
        lj_parameters={}
        for key in lj_parameters_keys:
            lj_parameters[key]=[]
        # Search the LJ parameters of the type pair
        for name in [particle_name1,particle_name2]:
            for key in lj_parameters_keys:
                lj_parameters[key].append(self.df[self.df.name == name][key].values[0])
        # If one of the particle has sigma=0, no LJ interations are set up between that particle type and the others    
        if not all(sigma_value.magnitude for sigma_value in lj_parameters["sigma"]):
            return {}
        # Apply combining rule
        if combining_rule == 'Lorentz-Berthelot':
            lj_parameters["sigma"]=(lj_parameters["sigma"][0]+lj_parameters["sigma"][1])/2
            lj_parameters["cutoff"]=(lj_parameters["cutoff"][0]+lj_parameters["cutoff"][1])/2
            lj_parameters["offset"]=(lj_parameters["offset"][0]+lj_parameters["offset"][1])/2
            lj_parameters["epsilon"]=np.sqrt(lj_parameters["epsilon"][0]*lj_parameters["epsilon"][1])
        return lj_parameters

    def get_metal_ions_charge_map(self):
        """
        Gets a map with the charge of all the metal ions supported.

        Returns:
            metal_charge_map(dict): Has the structure {"metal_name": metal_charge}

        """
        metal_charge_map = {"Ca": 2}
        return metal_charge_map

    def get_particle_id_map(self, object_name):
        '''
        Gets all the ids associated with the object with name `object_name` in `pmb.df`

        Args:
            object_name(`str`): name of the object
        
        Returns:
            id_map(`dict`): dict of the structure {"all": [all_ids_with_object_name], "residue_map": {res_id: [particle_ids_in_res_id]}, "molecule_map": {mol_id: [particle_ids_in_mol_id]}, }
        '''
        object_type = self.df.loc[self.df['name']== object_name].pmb_type.values[0]
        valid_types = ["particle", "molecule", "residue", "protein"]
        if object_type not in valid_types:
            raise ValueError(f"{object_name} is of pmb_type {object_type}, which is not supported by this function. Supported types are {valid_types}")
        id_list = []
        mol_map = {}
        res_map = {}
        def do_res_map(res_ids):
            for res_id in res_ids:
                res_list=self.df.loc[(self.df['residue_id']== res_id) & (self.df['pmb_type']== "particle")].particle_id.dropna().tolist()
                res_map[res_id]=res_list
            return res_map
        if object_type in ['molecule', 'protein']:
            mol_ids = self.df.loc[self.df['name']== object_name].molecule_id.dropna().tolist()
            for mol_id in mol_ids:
                res_ids = set(self.df.loc[(self.df['molecule_id']== mol_id) & (self.df['pmb_type']== "particle") ].residue_id.dropna().tolist())
                res_map=do_res_map(res_ids=res_ids)    
                mol_list=self.df.loc[(self.df['molecule_id']== mol_id) & (self.df['pmb_type']== "particle")].particle_id.dropna().tolist()
                id_list+=mol_list
                mol_map[mol_id]=mol_list
        elif object_type == 'residue':     
            res_ids = self.df.loc[self.df['name']== object_name].residue_id.dropna().tolist()
            res_map=do_res_map(res_ids=res_ids)
            id_list=[]
            for res_id_list in res_map.values():
                id_list+=res_id_list
        elif object_type == 'particle':
            id_list = self.df.loc[self.df['name']== object_name].particle_id.dropna().tolist()
        return {"all": id_list, "molecule_map": mol_map, "residue_map": res_map}

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
    
    def get_radius_map(self):
        '''
        Gets the effective radius of each `espresso type` in `pmb.df`. 
        
        Returns:
            radius_map(`dict`): {espresso_type: radius}.

        Note:
            The radius corresponds to (sigma+offset)/2
        '''
        df_state_one = self.df[[('sigma',''),('offset',''),('state_one','es_type')]].dropna().drop_duplicates()
        df_state_two = self.df[[('sigma',''),('offset',''),('state_two','es_type')]].dropna().drop_duplicates()
        state_one = pd.Series((df_state_one.sigma.values+df_state_one.offset.values)/2.0,index=df_state_one.state_one.es_type.values)
        state_two = pd.Series((df_state_two.sigma.values+df_state_two.offset.values)/2.0,index=df_state_two.state_two.es_type.values)
        radius_map  = pd.concat([state_one,state_two],axis=0).to_dict()  
        return radius_map

    def get_resource(self, path):
        '''
        Locate a file resource of the pyMBE package.

        Args:
            path(`str`): Relative path to the resource

        Returns:
            path(`str`): Absolute path to the resource

        '''
        import os
        return os.path.join(os.path.dirname(__file__), path)

    def get_type_map(self):
        """
        Gets all different espresso types assigned to particles  in `pmb.df`.
        
        Returns:
            type_map(`dict`): {"name": espresso_type}.
        """
        if self.df.state_one['es_type'].isnull().values.any():         
            df_state_one = self.df.state_one.dropna(how='all')     
            df_state_two = self.df.state_two.dropna(how='all')  
        else:    
            df_state_one = self.df.state_one
            if self.df.state_two['es_type'].isnull().values.any():
                df_state_two = self.df.state_two.dropna(how='all')   
            else:
                df_state_two = self.df.state_two
        state_one = pd.Series (df_state_one.es_type.values,index=df_state_one.label)
        state_two = pd.Series (df_state_two.es_type.values,index=df_state_two.label)
        type_map  = pd.concat([state_one,state_two],axis=0).to_dict()
        return type_map

    def load_interaction_parameters(self, filename, verbose=False, overwrite=False):
        """
        Loads the interaction parameters stored in `filename` into `pmb.df`
        
        Args:
            filename(`str`): name of the file to be read
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to False.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
        """
        without_units = ['q','es_type','acidity']
        with_units = ['sigma','epsilon','offset','cutoff']

        with open(filename, 'r') as f:
            interaction_data = json.load(f)
            interaction_parameter_set = interaction_data["data"]

        for key in interaction_parameter_set:
            param_dict=interaction_parameter_set[key]
            object_type=param_dict['object_type']
            if object_type == 'particle':
                not_requiered_attributes={}    
                for not_requiered_key in without_units+with_units:
                    if not_requiered_key in param_dict.keys():
                        if not_requiered_key in with_units:
                            not_requiered_attributes[not_requiered_key]=self.create_variable_with_units(variable=param_dict.pop(not_requiered_key))
                        elif not_requiered_key in without_units:
                            not_requiered_attributes[not_requiered_key]=param_dict.pop(not_requiered_key)
                    else:
                        if not_requiered_key == 'acidity':
                            not_requiered_attributes[not_requiered_key] = 'inert'
                        else:    
                            not_requiered_attributes[not_requiered_key]=None
                self.define_particle(name=param_dict.pop('name'),
                                q=not_requiered_attributes.pop('q'),
                                sigma=not_requiered_attributes.pop('sigma'),
                                offset=not_requiered_attributes.pop('offset'),
                                cutoff=not_requiered_attributes.pop('cutoff'),
                                acidity=not_requiered_attributes.pop('acidity'),
                                epsilon=not_requiered_attributes.pop('epsilon'),
                                verbose=verbose,
                                overwrite=overwrite)
            elif object_type == 'residue':
                self.define_residue (name = param_dict.pop('name'),
                                    central_bead = param_dict.pop('central_bead'),
                                    side_chains = param_dict.pop('side_chains'))
            elif object_type == 'molecule':
                self.define_molecule(name=param_dict.pop('name'),
                                    residue_list=param_dict.pop('residue_list'))
            elif object_type == 'peptide':
                self.define_peptide(name=param_dict.pop('name'),
                                    sequence=param_dict.pop('sequence'),
                                    model=param_dict.pop('model'))
            elif object_type == 'bond':
                particle_pairs = param_dict.pop('particle_pairs')
                bond_parameters = param_dict.pop('bond_parameters')
                bond_type = param_dict.pop('bond_type')
                if bond_type == 'harmonic':
                    k = self.create_variable_with_units(variable=bond_parameters.pop('k'))
                    r_0 = self.create_variable_with_units(variable=bond_parameters.pop('r_0'))
                    bond = {'r_0'    : r_0,
                            'k'      : k,
                            }

                elif bond_type == 'FENE':
                    k = self.create_variable_with_units(variable=bond_parameters.pop('k'))
                    r_0 = self.create_variable_with_units(variable=bond_parameters.pop('r_0'))
                    d_r_max = self.create_variable_with_units(variable=bond_parameters.pop('d_r_max'))
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
    
    def load_pka_set(self, filename, verbose=False, overwrite=True):
        """
        Loads the pka_set stored in `filename` into `pmb.df`.
        
        Args:
            filename(`str`): name of the file with the pka set to be loaded. Expected format is {name:{"acidity": acidity, "pka_value":pka_value}}.
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to False.
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
                                      verbose=verbose, 
                                      overwrite=overwrite)
        return

    def parse_sequence_from_file(self,sequence):
        """
        Parses the given sequence such that it can be used in pyMBE. This function has to be used if the df was read from a file.

        Args:
            sequence(`str`): sequence to be parsed

        Returns:
            sequence(`lst`): parsed sequence
        """
        sequence = sequence.replace(' ', '')
        sequence = sequence.replace("'", '')
        sequence = sequence.split(",")[1:-1]
        return sequence

    def print_reduced_units(self):
        """
        Prints the  current set of reduced units defined in pyMBE.units.
        """
        print("\nCurrent set of reduced units:")
        unit_length=self.units.Quantity(1,'reduced_length')
        unit_energy=self.units.Quantity(1,'reduced_energy')
        unit_charge=self.units.Quantity(1,'reduced_charge')
        print(f"{unit_length.to('nm'):.5g} = {unit_length}")
        print(f"{unit_energy.to('J'):.5g} = {unit_energy}")
        print(f"{unit_charge.to('C'):.5g} = {unit_charge}")
        print(f"Temperature: {(self.kT/self.Kb).to('K'):.5g}")
        print()

    def propose_unused_type(self):
        """
        Searches in `pmb.df` all the different particle types defined and returns a new one.

        Returns:
            unused_type(`int`): unused particle type
        """
        type_map=self.get_type_map()
        if type_map == {}:    
            unused_type = 0
        else:
            unused_type=max(type_map.values())+1    
        return unused_type

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
                        clean_sequence.append(keys[residue.upper()])
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
        columns_names = self.setup_df()
        
        multi_index = pd.MultiIndex.from_tuples(columns_names)
        df.columns = multi_index
        
        self.df = self.convert_columns_to_original_format(df)
        
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

        print (f'Loading protein coarse grain model file: {filename}')

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

    def search_bond(self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for bonds between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        If `use_default_bond` is activated and a "default" bond is defined, returns that default bond instead.
        If no bond is found, it prints a message and it does not return anything. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1(`str`): label of the type of the first particle type of the bonded particles.
            particle_name2(`str`): label of the type of the second particle type of the bonded particles.
            hard_check(`bool`, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond(`bool`, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            bond(`espressomd.interactions.BondedInteractions`): bond object from the espressomd library.
        
        Note:
            - If `use_default_bond`=True and no bond is defined between `particle_name1` and `particle_name2`, it returns the default bond defined in `pmb.df`.
            - If `hard_check`=`True` stops the code when no bond is found.
        """
        
        bond_key = self.find_bond_key(particle_name1=particle_name1, 
                                    particle_name2=particle_name2, 
                                    use_default_bond=use_default_bond)
        if use_default_bond:
            if not self.check_if_name_is_defined_in_df(name="default",pmb_type_to_be_defined='bond'):
                raise ValueError(f"use_default_bond is set to {use_default_bond} but no default bond has been defined. Please define a default bond with pmb.define_default_bond")
        if bond_key:
            return self.df[self.df['name']==bond_key].bond_object.values[0]
        else:
            print("Bond not defined between particles ", particle_name1, " and ", particle_name2)    
            if hard_check:
                sys.exit(1)
            else:
                return

    def search_particles_in_residue(self, residue_name):
        '''
        Searches for all particles in a given residue of name `residue_name`.

        Args:
            residue_name (`str`): name of the residue to be searched

        Returns:
            list_of_particles_in_residue (`lst`): list of the names of all particles in the residue

        Note:
            - The function returns a name per particle in residue, i.e. if there are multiple particles with the same type `list_of_particles_in_residue` will have repeated items.
 
        '''
        index_residue = self.df.loc[self.df['name'] == residue_name].index[0].item() 
        central_bead = self.df.at [index_residue, ('central_bead', '')]
        list_of_side_chains = self.df.at [index_residue, ('side_chains', '')]

        list_of_particles_in_residue = []
        list_of_particles_in_residue.append(central_bead)

        for side_chain in list_of_side_chains: 
            object_type = self.df[self.df['name']==side_chain].pmb_type.values[0]

            if object_type == "residue":
                list_of_particles_in_side_chain_residue = self.search_particles_in_residue(side_chain)
                list_of_particles_in_residue += list_of_particles_in_side_chain_residue
            elif object_type == "particle":
                list_of_particles_in_residue.append(side_chain)

        return list_of_particles_in_residue

    def set_particle_acidity(self, name, acidity='inert', default_charge=0, pka=None, verbose=True, overwrite=True):
        """
        Sets the particle acidity if it is acidic or basic, creates `state_one` and `state_two` with the protonated and 
        deprotonated states. In each state is set: `label`,`charge` and `es_type`. If it is inert, it will define only `state_one`.

        Args:
            name(`str`): Unique label that identifies the `particle`. 
            acidity(`str`): Identifies whether the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            default_charge (`int`): Charge of the particle. Defaults to 0.
            pka(`float`, optional):  If `particle` is an acid or a base, it defines its pka-value. Defaults to None.
            verbose(`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
        """
        acidity_valid_keys = ['inert','acidic', 'basic']
        if acidity not in acidity_valid_keys:
            raise ValueError(f"Acidity {acidity} provided for particle name  {name} is not supproted. Valid keys are: {acidity_valid_keys}")
        if acidity in ['acidic', 'basic'] and pka is None:
            raise ValueError(f"pKa not provided for particle with name {name} with acidity {acidity}. pKa must be provided for acidic or basic particles.")   
        
        self.define_particle_entry_in_df(name=name)
        
        for index in self.df[self.df['name']==name].index:       
            if pka:
                self.add_value_to_df(key=('pka',''),
                                    index=index,
                                    new_value=pka, 
                                    verbose=verbose,
                                    overwrite=overwrite)
           
            self.add_value_to_df(key=('acidity',''),
                                 index=index,
                                 new_value=acidity, 
                                 verbose=verbose,
                                 overwrite=overwrite) 
            if not self.check_if_df_cell_has_a_value(index=index,key=('state_one','es_type')):
                self.add_value_to_df(key=('state_one','es_type'),
                                     index=index,
                                     new_value=self.propose_unused_type(), 
                                     verbose=verbose,
                                     overwrite=overwrite)  
            if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'inert':
                self.add_value_to_df(key=('state_one','charge'),
                                     index=index,
                                     new_value=default_charge, 
                                     verbose=verbose,
                                     overwrite=overwrite)
                self.add_value_to_df(key=('state_one','label'),
                                     index=index,
                                     new_value=name, 
                                     verbose=verbose,
                                    overwrite=overwrite)
            else:
                protonated_label = f'{name}H'
                self.add_value_to_df(key=('state_one','label'),
                                     index=index,
                                     new_value=protonated_label, 
                                     verbose=verbose,
                                    overwrite=overwrite)
                self.add_value_to_df(key=('state_two','label'),
                                     index=index,
                                     new_value=name, 
                                     verbose=verbose,
                                    overwrite=overwrite)
                if not self.check_if_df_cell_has_a_value(index=index,key=('state_two','es_type')):
                    self.add_value_to_df(key=('state_two','es_type'),
                                         index=index,
                                         new_value=self.propose_unused_type(), 
                                         verbose=verbose,
                                         overwrite=overwrite)
                if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'acidic':        
                    self.add_value_to_df(key=('state_one','charge'),
                                         index=index,new_value=0, 
                                         verbose=verbose,
                                         overwrite=overwrite)
                    self.add_value_to_df(key=('state_two','charge'),
                                         index=index,
                                         new_value=-1, 
                                         verbose=verbose,
                                         overwrite=overwrite)
                elif self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'basic':
                    self.add_value_to_df(key=('state_one','charge'),
                                         index=index,new_value=+1, 
                                         verbose=verbose,
                                         overwrite=overwrite)
                    self.add_value_to_df(key=('state_two','charge'),
                                         index=index,
                                         new_value=0, 
                                         verbose=verbose,
                                         overwrite=overwrite)   
        return
    
    def set_reduced_units(self, unit_length=None, unit_charge=None, temperature=None, Kw=None, verbose=True):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args:
            unit_length(`pint.Quantity`,optional): Reduced unit of length defined using the `pmb.units` UnitRegistry. Defaults to None. 
            unit_charge(`pint.Quantity`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            temperature(`pint.Quantity`,optional): Temperature of the system, defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw(`pint.Quantity`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 
            verbose(`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.

        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
        """
        self.units=pint.UnitRegistry()
        if unit_length is None:
            unit_length=0.355*self.units.nm
        if temperature is None:
            temperature=298.15 * self.units.K
        if unit_charge is None:
            unit_charge=self.units.e
        if Kw is None:
            Kw = 1e-14
        self.N_A=6.02214076e23 / self.units.mol
        self.Kb=1.38064852e-23 * self.units.J / self.units.K
        self.e=1.60217662e-19 *self.units.C
        self.kT=temperature*self.Kb
        self.Kw=Kw*self.units.mol**2 / (self.units.l**2)
        self.units.define(f'reduced_energy = {self.kT} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')
        if verbose:        
            self.print_reduced_units()
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
                                                    exclusion_range=exclusion_range.magnitude, 
                                                    seed=self.SEED, 
                                                    constant_pH=constant_pH,
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )
        sucessfull_reactions_labels=[]
        charge_map = self.get_charge_map()
        for name in pka_set.keys():
            if self.check_if_name_is_defined_in_df(name,pmb_type_to_be_defined='particle') is False :
                print('WARNING: the acid-base reaction of ' + name +' has not been set up because its espresso type is not defined in sg.type_map')
                continue
            gamma=10**-pka_set[name]['pka_value']
            state_one_type   = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type   = self.df.loc[self.df['name']==name].state_two.es_type.values[0]
            counter_ion_type = self.df.loc[self.df['name']==counter_ion].state_one.es_type.values[0]
            RE.add_reaction(gamma=gamma,
                            reactant_types=[state_one_type],
                            product_types=[state_two_type, counter_ion_type],
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            counter_ion_type: charge_map[counter_ion_type]})
            sucessfull_reactions_labels.append(name)
        return RE, sucessfull_reactions_labels

    def setup_gcmc(self, c_salt_res, salt_cation_name, salt_anion_name, activity_coefficient=None, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up grand-canonical coupling to a reservoir of salt.
        For reactive systems coupled to a reservoir, the grand-reaction method has to be used instead.

        Args:
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
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
        
        #If no function for the activity coefficient is provided, the ideal case is assumed.
        if activity_coefficient is None:
            activity_coefficient = lambda x: 1.0
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range.magnitude, 
                                                    seed=self.SEED, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        determined_activity_coefficient = activity_coefficient(c_salt_res)
        K_salt = (c_salt_res.to('1/(N_A * reduced_length**3)')**2) * determined_activity_coefficient

        salt_cation_es_type = self.df.loc[self.df['name']==salt_cation_name].state_one.es_type.values[0]
        salt_anion_es_type = self.df.loc[self.df['name']==salt_anion_name].state_one.es_type.values[0]     

        salt_cation_charge = self.df.loc[self.df['name']==salt_cation_name].state_one.charge.values[0]
        salt_anion_charge = self.df.loc[self.df['name']==salt_anion_name].state_one.charge.values[0]     

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

    def setup_grxmc_reactions(self, pH_res, c_salt_res, proton_name, hydroxide_name, salt_cation_name, salt_anion_name, activity_coefficient=None, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
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
            activity_coefficient ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
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
        
        #If no function for the activity coefficient is provided, the ideal case is assumed.
        if activity_coefficient is None:
            activity_coefficient = lambda x: 1.0
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range.magnitude, 
                                                    seed=self.SEED, 
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

        proton_charge = self.df.loc[self.df['name']==proton_name].state_one.charge.values[0]
        hydroxide_charge = self.df.loc[self.df['name']==hydroxide_name].state_one.charge.values[0]     
        salt_cation_charge = self.df.loc[self.df['name']==salt_cation_name].state_one.charge.values[0]
        salt_anion_charge = self.df.loc[self.df['name']==salt_anion_name].state_one.charge.values[0]     

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
        charge_map = self.get_charge_map()
        for name in pka_set.keys():
            if self.check_if_name_is_defined_in_df(name,pmb_type_to_be_defined='particle') is False :
                print('WARNING: the acid-base reaction of ' + name +' has not been set up because its espresso type is not defined in the type map.')
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
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            proton_es_type: proton_charge})

            # Reaction in terms of salt cation: HA = A + Na+
            RE.add_reaction(gamma=(Ka * K_NACL / K_HCL).magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, salt_cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            salt_cation_es_type: salt_cation_charge})

            # Reaction in terms of hydroxide: OH- + HA = A
            RE.add_reaction(gamma=(Ka / K_W).magnitude,
                            reactant_types=[state_one_type, hydroxide_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            hydroxide_es_type: hydroxide_charge})

            # Reaction in terms of salt anion: Cl- + HA = A
            RE.add_reaction(gamma=(Ka / K_HCL).magnitude,
                            reactant_types=[state_one_type, salt_anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            salt_anion_es_type: salt_anion_charge})

            sucessful_reactions_labels.append(name)
        return RE, sucessful_reactions_labels, ionic_strength_res

    def setup_grxmc_unified (self, pH_res, c_salt_res, cation_name, anion_name, activity_coefficient=None, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
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
            activity_coefficient ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
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
        
        #If no function for the activity coefficient is provided, the ideal case is assumed.
        if activity_coefficient is None:
            activity_coefficient = lambda x: 1.0
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range.magnitude, 
                                                    seed=self.SEED, 
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
        cation_charge = self.df.loc[self.df['name']==cation_name].state_one.charge.values[0]
        anion_charge = self.df.loc[self.df['name']==anion_name].state_one.charge.values[0]     
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
        charge_map = self.get_charge_map()
        for name in pka_set.keys():
            if self.check_if_name_is_defined_in_df(name,pmb_type_to_be_defined='particle') is False :
                print('WARNING: the acid-base reaction of ' + name +' has not been set up because its espresso type is not defined in sg.type_map')
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
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            cation_es_type: cation_charge})

            # Reaction in terms of small anion: X- + HA = A
            RE.add_reaction(gamma=gamma_K_AX.magnitude / K_XX.magnitude,
                            reactant_types=[state_one_type, anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_map[state_one_type],
                            state_two_type: charge_map[state_two_type],
                            anion_es_type: anion_charge})

            sucessful_reactions_labels.append(name)
        return RE, sucessful_reactions_labels, ionic_strength_res

    def setup_df (self):
        """
        Sets up the pyMBE's dataframe `pymbe.df`.

        Returns:
            columns_names(`obj`): pandas multiindex object with the column names of the pyMBE's dataframe
        """
        
        columns_dtypes = {
            'name': {
                '': str},
            'pmb_type': {
                '': str},
            'particle_id': {
                '': object},
            'particle_id2':  {
                '': object},
            'residue_id':  {
                '': object},
            'molecule_id':  {
                '': object},
            'acidity':  {
                '': str},
            'pka':  {
                '': float},
            'central_bead':  {
                '': object},
            'side_chains': {
                '': object},
            'residue_list': {
                '': object},
            'model': {
                '': str},
            'sigma': {
                '': object},
            'cutoff': {
                '': object},
            'offset': {
                '': object},
            'epsilon': {
                '': object},
            'state_one': {
                'label': str,
                'es_type': object,
                'charge': object },
            'state_two': {
                'label': str,
                'es_type': object,
                'charge': object },
            'sequence': {
                '': object},
            'bond_object': {
                '': object},
            'parameters_of_the_potential':{
                '': object},
            'l0': {
                '': float},
        }
        
        self.df = pd.DataFrame(columns=pd.MultiIndex.from_tuples([(col_main, col_sub) for col_main, sub_cols in columns_dtypes.items() for col_sub in sub_cols.keys()]))
        
        for level1, sub_dtypes in columns_dtypes.items():
            for level2, dtype in sub_dtypes.items():
                self.df[level1, level2] = self.df[level1, level2].astype(dtype)
                
        columns_names = pd.MultiIndex.from_frame(self.df)
        columns_names = columns_names.names
                
        return columns_names

    def setup_lj_interactions(self, espresso_system, shift_potential=True, combining_rule='Lorentz-Berthelot', warnings=True):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle types with values for `sigma`, `offset`, and `epsilon` stored in `pymbe.df`.

        Args:
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            shift_potential(`bool`, optional): If True, a shift will be automatically computed such that the potential is continuous at the cutoff radius. Otherwise, no shift will be applied. Defaults to True.
            combining_rule(`string`, optional): combining rule used to calculate `sigma` and `epsilon` for the potential between a pair of particles. Defaults to 'Lorentz-Berthelot'.
            warning(`bool`, optional): switch to activate/deactivate warning messages. Defaults to True.

        Note:
            - LJ interactions will only be set up between particles with defined values of `sigma` and `epsilon` in the pmb.df. 
            - Currently, the only `combining_rule` supported is Lorentz-Berthelot.
            - Check the documentation of ESPResSo for more info about the potential https://espressomd.github.io/doc4.2.0/inter_non-bonded.html

        """
        from itertools import combinations_with_replacement
        implemented_combining_rules = ['Lorentz-Berthelot']
        compulsory_parameters_in_df = ['sigma','epsilon']
        # Sanity check
        if combining_rule not in implemented_combining_rules:
            raise ValueError('In the current version of pyMBE, the only combinatorial rules implemented are ', implemented_combining_rules)
        if shift_potential:
            shift="auto"
        else:
            shift=0
        # List which particles have sigma and epsilon values defined in pmb.df and which ones don't
        particles_types_with_LJ_parameters = []
        non_parametrized_labels= []
        for particle_type in self.get_type_map().values():
            check_list=[]
            for key in compulsory_parameters_in_df:
                value_in_df=self.find_value_from_es_type(es_type=particle_type,
                                                        column_name=key)
                check_list.append(np.isnan(value_in_df))
            if any(check_list):
                non_parametrized_labels.append(self.find_value_from_es_type(es_type=particle_type, 
                                                                            column_name='label'))
            else:
                particles_types_with_LJ_parameters.append(particle_type)
        # Set up LJ interactions between all particle types
        for type_pair in combinations_with_replacement(particles_types_with_LJ_parameters, 2): 
            particle_name1 = self.find_value_from_es_type(es_type=type_pair[0],
                                                        column_name="name")
            particle_name2 = self.find_value_from_es_type(es_type=type_pair[1],
                                                        column_name="name")
            lj_parameters= self.get_lj_parameters(particle_name1 = particle_name1,
                                                 particle_name2 = particle_name2,
                                                 combining_rule = combining_rule)
            
            # If one of the particle has sigma=0, no LJ interations are set up between that particle type and the others    
            if not lj_parameters:
                continue
            espresso_system.non_bonded_inter[type_pair[0],type_pair[1]].lennard_jones.set_params(epsilon = lj_parameters["epsilon"].to('reduced_energy').magnitude, 
                                                                                    sigma = lj_parameters["sigma"].to('reduced_length').magnitude, 
                                                                                    cutoff = lj_parameters["cutoff"].to('reduced_length').magnitude,
                                                                                    offset = lj_parameters["offset"].to("reduced_length").magnitude, 
                                                                                    shift = shift)                                                                                          
            index = len(self.df)
            label1 = self.find_value_from_es_type(es_type=type_pair[0], column_name="label")
            label2 = self.find_value_from_es_type(es_type=type_pair[1], column_name="label")
            self.df.at [index, 'name'] = f'LJ: {label1}-{label2}'
            lj_params=espresso_system.non_bonded_inter[type_pair[0], type_pair[1]].lennard_jones.get_params()

            self.add_value_to_df(index=index,
                                key=('pmb_type',''),
                                new_value='LennardJones')
            
            self.add_value_to_df(index=index,
                                key=('parameters_of_the_potential',''),
                                new_value=lj_params,
                                non_standard_value=True)
        if non_parametrized_labels and warnings:
            print(f'WARNING: The following particles do not have a defined value of sigma or epsilon in pmb.df: {non_parametrized_labels}. No LJ interaction has been added in ESPResSo for those particles.')
            
        return

    def write_pmb_df (self, filename):
        '''
        Writes the pyMBE dataframe into a csv file
        
        Args:
            filename(`str`): Path to the csv file 
        '''

        self.df.to_csv(filename)

        return
