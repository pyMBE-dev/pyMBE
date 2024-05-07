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
    import pint
    import re
    import numpy as np
    import pandas as pd 
    import json
    from scipy import optimize
    units = pint.UnitRegistry()
    N_A=6.02214076e23    / units.mol
    Kb=1.38064852e-23    * units.J / units.K
    e=1.60217662e-19 *units.C
    df=None
    kT=None
    Kw=None
    SEED=None
    rng=None

    def __init__(self, SEED, temperature=None, unit_length=None, unit_charge=None, Kw=None):
        """
        Initializes the pymbe_library by setting up the reduced unit system with `temperature` and `reduced_length` 
        and sets up  the `pmb.df` for bookkeeping.

        Args:
            temperature(`obj`,optional): Value of the temperature in the pyMBE UnitRegistry. Defaults to None.
            unit_length(`obj`, optional): Value of the unit of length in the pyMBE UnitRegistry. Defaults to None.
            unit_charge (`obj`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw (`obj`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 
        
        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
        """
        # Seed and RNG
        self.SEED=SEED
        self.rng = self.np.random.default_rng(SEED)
        # Default definitions of reduced units
        if temperature is None:
            temperature= 298.15 * self.units.K
        if unit_length is None:
            unit_length= 0.355 * self.units.nm
        if unit_charge is None:
            unit_charge=self.units.e
        if Kw is None:
            Kw = 1e-14
        self.kT=temperature*self.Kb
        self.Kw=Kw*self.units.mol**2 / (self.units.l**2)
        self.units.define(f'reduced_energy = {self.kT}')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define('reduced_charge = 1*e')
        self.setup_df()
        return
    
    def enable_motion_of_rigid_object(self, name, espresso_system):
        '''
        Enables the motion of the rigid object `name` in the `espresso_system`.

        Args:
            name(`str`): Label of the object.
            espresso_system(`obj`): Instance of a system object from the espressomd library.

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

    def add_bond_in_df(self, particle_id1, particle_id2, use_default_bond=False):
        """
        Adds a bond entry on the `pymbe.df` storing the particle_ids of the two bonded particles.

        Args:
            particle_id1 (`int`): particle_id of the type of the first particle type of the bonded particles
            particle_id2 (`int`): particle_id of the type of the second particle type of the bonded particles
            use_default_bond (`bool`, optional): Controls if a bond of type `default` is used to bond particle whose bond types are not defined in `pmb.df`. Defaults to False.
        """
        particle_name1 = self.df.loc[self.df['particle_id']==particle_id1].name.values[0]
        particle_name2 = self.df.loc[self.df['particle_id']==particle_id2].name.values[0]
        bond_key = self.find_bond_key(particle_name1=particle_name1,
                                    particle_name2=particle_name2, 
                                    use_default_bond=use_default_bond)
        if not bond_key:
            return
        self.copy_df_entry(name=bond_key,column_name='particle_id2',number_of_copies=1)
        indexs = self.np.where(self.df['name']==bond_key)
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
        return

    def add_bonds_to_espresso(self, espresso_system) :
        """
        Adds all bonds defined in `pmb.df` to `espresso_system`.

        Args:
            espresso_system (str): system object of espressomd library
        """

        if 'bond' in self.df.values:
            bond_df = self.df.loc[self.df ['pmb_type'] == 'bond']
            bond_list = bond_df.bond_object.values.tolist()
            for bond in bond_list:
                espresso_system.bonded_inter.add(bond)
        else:
            print ('WARNING: There are no bonds defined in pymbe.df')
        
        return

    def add_value_to_df(self,index,key,new_value, verbose=True):
        """
        Adds a value to a cell in the `pmb.df` DataFrame.

        Args:
            index(`int`): index of the row to add the value to.
            key(`str`): the column label to add the value to.
            verbose(`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.
        """
        # Make sure index is a scalar integer value
        index = int (index)
        assert isinstance(index, int), '`index` should be a scalar integer value.'
        idx = self.pd.IndexSlice
        if self.check_if_df_cell_has_a_value(index=index,key=key):
            old_value= self.df.loc[index,idx[key]]
            if verbose:
                if old_value != new_value:
                    name=self.df.loc[index,('name','')]
                    pmb_type=self.df.loc[index,('pmb_type','')]
                    print(f"WARNING: you are redefining the properties of {name} of pmb_type {pmb_type}")
                    print(f'WARNING: overwritting the value of the entry `{key}`: old_value = {old_value} new_value = {new_value}')
        self.df.loc[index,idx[key]] = new_value
        return
    
    def assign_molecule_id(self, name, molecule_index, pmb_type, used_molecules_id):
        """
        Assigns the `molecule_id` of the pmb object given by `pmb_type`
        
        Args:
            name(`str`): Label of the molecule type to be created. `name` must be defined in `pmb.df`
            pmb_type(`str`): pmb_object_type to assign the `molecule_id` 
            molecule_index (`int`): index of the current `pmb_object_type` to assign the `molecule_id`
            used_molecules_id (`lst`): list with the `molecule_id` values already used.
        
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
            espresso_system(`obj`): Instance of a system object from the espressomd library.
        
        Returns:
            center_of_mass(`lst`): Coordinates of the center of mass.
        """
        total_beads = 0
        center_of_mass = self.np.zeros(3)
        axis_list = [0,1,2]
        particle_id_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
        for pid in particle_id_list:
            total_beads +=1 
            for axis in axis_list:
                center_of_mass [axis] += espresso_system.part.by_id(pid).pos[axis]
        center_of_mass = center_of_mass /total_beads  
        return center_of_mass

    def calculate_HH(self, molecule_name, pH_list=None, pka_set=None):
        """
        Calculates the charge per molecule according to the ideal Henderson-Hasselbalch titration curve 
        for molecules with the name `molecule_name`.

        Args:
            molecule_name (`str`): name of the molecule to calculate the ideal charge for
            pH_list(`lst`): pH-values to calculate. 
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}

        Returns:
            Z_HH (`lst`): Henderson-Hasselbalch prediction of the charge of `sequence` in `pH_list`

        Note:
            - This function supports objects with pmb types: "molecule", "peptide" and "protein".
            - If no `pH_list` is given, 50 equispaced pH-values ranging from 2 to 12 are calculated
            - If no `pka_set` is given, the pKa values are taken from `pmb.df`
            - This function should only be used for single-phase systems. For two-phase systems `pmb.calculate_HH_Donnan`  should be used.
        """
        if pH_list is None:
            pH_list=self.np.linspace(2,12,50)    
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
            if self.np.any(self.pd.isnull(sequence)):
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
            c_macro ('dict'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
            c_salt ('float'): Salt concentration in the reservoir.
            pH_list ('lst'): List of pH-values in the reservoir. 
            pka_set ('dict'): {"name": {"pka_value": pka, "acidity": acidity}}.

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
            pH_list=self.np.linspace(2,12,50)    
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
                c_macro ('dic'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
                pH ('float'): pH-value that is used in the HH equation.

            Returns:
                charge ('dict'): {"molecule_name": charge}
            """
            charge = {}
            for name in c_macro:
                charge[name] = self.calculate_HH(name, [pH], pka_set)[0]
            return charge

        def calc_partition_coefficient(charge, c_macro):
            """
            Calculates the partition coefficients of positive ions according to the ideal Donnan theory.

            Args:
                charge ('dict'): {"molecule_name": charge}
                c_macro ('dic'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
            """
            nonlocal ionic_strength_res
            charge_density = 0.0
            for key in charge:
                charge_density += charge[key] * c_macro[key]
            return (-charge_density / (2 * ionic_strength_res) + self.np.sqrt((charge_density / (2 * ionic_strength_res))**2 + 1)).magnitude

        for pH_value in pH_list:    
            # calculate the ionic strength of the reservoir
            if pH_value <= 7.0:
                ionic_strength_res = 10 ** (-pH_value) * self.units.mol/self.units.l + c_salt 
            elif pH_value > 7.0:
                ionic_strength_res = 10 ** (-(14-pH_value)) * self.units.mol/self.units.l + c_salt

            #Determine the partition coefficient of positive ions by solving the system of nonlinear, coupled equations
            #consisting of the partition coefficient given by the ideal Donnan theory and the Henderson-Hasselbalch equation.
            #The nonlinear equation is formulated for log(xi) since log-operations are not supported for RootResult objects.
            equation = lambda logxi: logxi - self.np.log10(calc_partition_coefficient(calc_charges(c_macro, pH_value - logxi), c_macro))
            logxi = self.optimize.root_scalar(equation, bracket=[-1e2, 1e2], method="brentq")
            partition_coefficient = 10**logxi.root

            charges_temp = calc_charges(c_macro, pH_value-self.np.log10(partition_coefficient))
            for key in c_macro:
                Z_HH_Donnan[key].append(charges_temp[key])

            pH_system_list.append(pH_value - self.np.log10(partition_coefficient))
            partition_coefficients_list.append(partition_coefficient)

        return {"charges_dict": Z_HH_Donnan, "pH_system_list": pH_system_list, "partition_coefficients": partition_coefficients_list}

    def calculate_initial_bond_length(self, bond_object, bond_type, epsilon, sigma, cutoff):
        """
        Calculates the initial bond length that is used when setting up molecules,
        based on the minimum of the sum of bonded and short-range (LJ) interactions.
        
        Args:
            bond_object (`cls`): instance of a bond object from espressomd library
            bond_type (`str`): label identifying the used bonded potential
            epsilon (`float`): LJ epsilon of the interaction between the particles
            sigma (`float`): LJ sigma of the interaction between the particles
            cutoff (`float`): cutoff-radius of the LJ interaction 
        """    
        def truncated_lj_potential(x, epsilon, sigma, cutoff):
            if x>cutoff:
                return 0.0
            else:
                return 4*epsilon*((sigma/x)**12-(sigma/x)**6) - 4*epsilon*((sigma/cutoff)**12-(sigma/cutoff)**6)

        epsilon_red=epsilon.to('reduced_energy').magnitude
        sigma_red=sigma.to('reduced_length').magnitude
        cutoff_red=cutoff.to('reduced_length').magnitude

        if bond_type == "harmonic":
            r_0 = bond_object.params.get('r_0')
            k = bond_object.params.get('k')
            l0 = self.optimize.minimize(lambda x: 0.5*k*(x-r_0)**2 + truncated_lj_potential(x, epsilon_red, sigma_red, cutoff_red), x0=r_0).x
        elif bond_type == "FENE":
            r_0 = bond_object.params.get('r_0')
            k = bond_object.params.get('k')
            d_r_max = bond_object.params.get('d_r_max')
            l0 = self.optimize.minimize(lambda x: -0.5*k*(d_r_max**2)*self.np.log(1-((x-r_0)/d_r_max)**2) + truncated_lj_potential(x, epsilon_red, sigma_red, cutoff_red), x0=1.0).x
        return l0

    def calculate_net_charge (self, espresso_system, molecule_name):
        '''
        Calculates the net charge per molecule of molecules with `name` = molecule_name. 
        Returns the net charge per molecule and a maps with the net charge per residue and molecule.

        Args:
            espresso_system: system information 
            molecule_name (str): name of the molecule to calculate the net charge

        Returns:
            {"mean": mean_net_charge, "molecules": {mol_id: net_charge_of_mol_id, }, "residues": {res_id: net_charge_of_res_id, }}

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
        mean_charge=self.np.mean(self.np.array(list(net_charge_molecules.values())))
        return {"mean": mean_charge, "molecules": net_charge_molecules, "residues": net_charge_residues}

    def center_molecule_in_simulation_box(self, molecule_id, espresso_system):
        """
        Centers the pmb object matching `molecule_id` in the center of the simulation box in `espresso_md`.
        
        Args:
            molecule_id(`int`): Id of the molecule to be centered.
            espresso_system(`obj`): Instance of a system object from the espressomd library.
        """
        if len(self.df.loc[self.df['molecule_id']==molecule_id].pmb_type) == 0:
            raise ValueError("The provided molecule_id is not present in the pyMBE dataframe.")      
        center_of_mass = self.calculate_center_of_mass_of_molecule(molecule_id=molecule_id,espresso_system=espresso_system)
        box_center = [espresso_system.box_l[0]/2.0,
                      espresso_system.box_l[1]/2.0,
                      espresso_system.box_l[2]/2.0]
        particle_id_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
        for pid in particle_id_list:
            es_pos = espresso_system.part.by_id(pid).pos
            espresso_system.part.by_id(pid).pos = es_pos - center_of_mass + box_center
        return 

    def check_dimensionality(self, variable, expected_dimensionality):
        """
        Checks if the dimensionality of `variable` matches `expected_dimensionality`.

        Args:
            `variable`(`pint.Quantity`): Quantity to be checked.
            `expected_dimensionality(`str`): Expected dimension of the variable.

        Returns:
            `bool`: `True` if the variable if of the expected dimensionality, `False` otherwise.

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
        idx = self.pd.IndexSlice
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return not self.pd.isna(self.df.loc[index, idx[key]])

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
                raise ValueError ((f"The name {name} is already defined in the df with a pmb_type = {current_object_type}, pymMBE does not support objects with the same name but different pmb_types"))    
            return True            
        else:
            return False


    def check_pka_set(self, pka_set):
        """
        Checks that `pka_set` has the formatting expected by the pyMBE library.
       
        Args:
            pka_set (dict): {"name" : {"pka_value": pka, "acidity": acidity}}
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_entry in pka_set.values():
                if required_key not in pka_entry.keys():
                    raise ValueError('missing a requiered key ', required_keys, 'in the following entry of pka_set', pka_entry)
        return

    def clean_df_row(self, index, columns_keys_to_clean=("particle_id", "particle_id2", "residue_id", "molecule_id")):
        """
        Cleans the columns of `pmb.df` in `columns_keys_to_clean` of the row with index `index` by assigning them a np.nan value.

        Args:
            index(`int`): Index of the row to clean.
            columns_keys_to_clean(`list` of `str`, optional): List with the column keys to be cleaned. Defaults to [`particle_id`, `particle_id2`, `residue_id`, `molecule_id`].
        """   
        for column_key in columns_keys_to_clean:
            self.add_value_to_df(key=(column_key,''),index=index,new_value=self.np.nan, verbose=False)
        return

    def convert_columns_to_original_format (self, df):
        """
        Converts the columns of the Dataframe to the original format in pyMBE.
        
        Args:
            df(`DataFrame`): dataframe with pyMBE information as a string  
        
        """

        from ast import literal_eval

        columns_dtype_int = ['particle_id','particle_id2', 'residue_id','molecule_id', 'model',('state_one','es_type'),('state_two','es_type'),('state_one','charge'),('state_two','charge') ]  

        columns_with_units = ['sigma', 'epsilon', 'cutoff', 'offset']

        columns_with_list_or_dict = ['residue_list','side_chains', 'parameters_of_the_potential','sequence']

        for column_name in columns_dtype_int:
            df[column_name] = df[column_name].astype(object)
            
        for column_name in columns_with_list_or_dict:
            if df[column_name].isnull().all():
                df[column_name] = df[column_name].astype(object)
            else:
                df[column_name] = df[column_name].apply(lambda x: literal_eval(str(x)) if self.pd.notnull(x) else x)

        for column_name in columns_with_units:
            df[column_name] = df[column_name].apply(lambda x: self.create_variable_with_units(x) if self.pd.notnull(x) else x)

        df['bond_object'] = df['bond_object'].apply(lambda x: (self.convert_str_to_bond_object(x)) if self.pd.notnull(x) else x)

        return df
    
    def convert_str_to_bond_object (self, bond_str):
        
        """
        Convert a row read as a `str` to the corresponding bond object. There are two supported bonds: HarmonicBond and FeneBond

        Args:
            bond_str (`str`): string with the information of a bond object

        Returns:
            bond_object(`obj`): EsPRESSo bond object 
        """
        
        from ast import literal_eval
        from espressomd.interactions import HarmonicBond
        from espressomd.interactions import FeneBond

        supported_bonds = ['HarmonicBond', 'FeneBond']

        for bond in supported_bonds:

            variable = self.re.subn(f'{bond}', '', bond_str)

            if variable[1] == 1: 
            
                params = literal_eval(variable[0])

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
                df_by_name_repeated = self.pd.concat ([df_by_name]*(number_of_copies-1), ignore_index=True)
            else:
                df_by_name = df_by_name[df_by_name.index == df_by_name.index.min()] 
                df_by_name_repeated = self.pd.concat ([df_by_name]*(number_of_copies), ignore_index=True)
                df_by_name_repeated[column_name] =self.np.NaN
            # Concatenate the new particle rows to  `df`
            self.df = self.pd.concat ([self.df,df_by_name_repeated], ignore_index=True)
        else:
            if not df_by_name[column_name].isnull().values.any():     
                df_by_name = df_by_name[df_by_name.index == df_by_name.index.min()] 
                df_by_name_repeated = self.pd.concat ([df_by_name]*(number_of_copies), ignore_index=True)
                df_by_name_repeated[column_name] =self.np.NaN
                self.df = self.pd.concat ([self.df,df_by_name_repeated], ignore_index=True)        
        return

    def create_added_salt (self, espresso_system, cation_name, anion_name, c_salt, verbose=True):    
        """
        Creates a `c_salt` concentration of `cation_name` and `anion_name` ions into the `espresso_system`.

        Args:
            espresso_system (`obj`): instance of an espresso system object.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            c_salt(`float`): Salt concentration.
            verbose(`bool`): switch to activate/deactivate verbose. Defaults to True.
            
        Returns:
            c_salt_calculated (float): Calculated salt concentration added to `espresso_system`.
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
            N_ions= int((volume*c_salt.to('reduced_length**-3')*self.N_A).magnitude)
            c_salt_calculated=N_ions/volume
        else:
            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)
        N_cation = N_ions*abs(anion_name_charge)
        N_anion = N_ions*abs(cation_name_charge)
        self.create_particle(espresso_system=espresso_system, name=cation_name, number_of_particles=N_cation)
        self.create_particle(espresso_system=espresso_system, name=anion_name, number_of_particles=N_anion)
        if verbose:
            print(f"\n Added salt concentration of {c_salt_calculated.to('mol/L')} given by {N_cation} cations and {N_anion} anions")
        return c_salt_calculated

    def create_bond_in_espresso(self, bond_type, bond_parameters):
        '''
        Creates either a harmonic or a FENE bond in ESPREesSo

        Args:
            bond_type (`str`)        : label to identify the potential to model the bond.
            bond_parameters (`dict`) : parameters of the potential of the bond.

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
            espresso_system(`obj`): Instance of a system object from the espressomd library.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            verbose(`bool`): switch to activate/deactivate verbose. Defaults to True.

        Returns: 
            counterion_number (dict): {"name": number}
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
                object_charge['positive']+=1*(self.np.abs(espresso_system.part.by_id(id).q ))
            elif espresso_system.part.by_id(id).q < 0:
                object_charge['negative']+=1*(self.np.abs(espresso_system.part.by_id(id).q ))
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
            espresso_system(`obj`): Instance of a system object from espressomd library.
            number_of_molecules(`int`): Number of molecules of type `name` to be created.
            list_of_first_residue_positions(`list`, optional): List of coordinates where the central bead of the first_residue_position will be created, random by default
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle with undefined bonds in `pymbe.df`
        Returns:

            molecules_info (`dict`):  {molecule_id: {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids": [particle_id1, ...]}}} 
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
            return
        if not self.check_if_name_is_defined_in_df(name=name,
                                                    pmb_type_to_be_defined='molecule'):
            raise ValueError(f"{name} must correspond to a label of a pmb_type='molecule' defined on df")
        first_residue = True
        molecules_info = {}
        residue_list = self.df[self.df['name']==name].residue_list.values [0]

        self.copy_df_entry(name=name,
                        column_name='molecule_id',
                        number_of_copies=number_of_molecules)
        
        molecules_index = self.np.where(self.df['name']==name)
        molecule_index_list =list(molecules_index[0])[-number_of_molecules:]
        used_molecules_id = self.df.molecule_id.dropna().drop_duplicates().tolist()
        pos_index = 0 
        for molecule_index in molecule_index_list:        
            molecule_id = self.assign_molecule_id(name=name,pmb_type='molecule',used_molecules_id=used_molecules_id,molecule_index=molecule_index)
            molecules_info[molecule_id] = {}
            for residue in residue_list:
                if first_residue:
                    if list_of_first_residue_positions is None:
                        residue_position = None
                    else:
                        for item in list_of_first_residue_positions:
                            residue_position = [self.np.array(list_of_first_residue_positions[pos_index])]
                    # Generate an arbitrary random unit vector
                    backbone_vector = self.generate_random_points_in_a_sphere(center=[0,0,0], 
                                                                radius=1, 
                                                                n_samples=1,
                                                                on_surface=True)[0]
                    residues_info = self.create_residue(name=residue,
                                                        number_of_residues=1, 
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=residue_position,  
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    residue_id = next(iter(residues_info))
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        self.add_value_to_df(key=('molecule_id',''),
                                            index=int (index),
                                            new_value=molecule_id)
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
                                                        number_of_residues=1, 
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=[residue_position],
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    residue_id = next(iter(residues_info))      
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        if not self.check_if_df_cell_has_a_value(index=index,key=('molecule_id','')):
                            self.df.at[index,'molecule_id'] = molecule_id
                            self.add_value_to_df(key=('molecule_id',''),
                                                index=int (index),
                                                new_value=molecule_id,
                                                verbose=False)            
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, previous_residue_id))
                    self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=previous_residue_id,
                                        use_default_bond=use_default_bond)            
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
            name (`str`): Label of the particle type to be created. `name` must be a `particle` defined in `pmb_df`. 
            espresso_system (`cls`): Instance of a system object from the espressomd library.
            number_of_particles (`int`): Number of particles to be created.
            position (list of [`float`,`float`,`float`], optional): Initial positions of the particles. If not given, particles are created in random positions. Defaults to None.
            fix(`bool`, optional): Controls if the particle motion is frozen in the integrator, it is used to create rigid objects. Defaults to False.
        Returns:
            created_pid_list(`list` of `float`): List with the ids of the particles created into `espresso_system`.
        """       
        if number_of_particles <=0:
            return
        if not self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='particle'):
            raise ValueError(f"{name} must correspond to a label of a pmb_type='particle' defined on df")
        # Copy the data of the particle `number_of_particles` times in the `df`
        self.copy_df_entry(name=name,column_name='particle_id',number_of_copies=number_of_particles)
        # Get information from the particle type `name` from the df     
        q = self.df.loc[self.df['name']==name].state_one.charge.values[0]
        es_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
        # Get a list of the index in `df` corresponding to the new particles to be created
        index = self.np.where(self.df['name']==name)
        index_list =list(index[0])[-number_of_particles:]
        # Create the new particles into  `espresso_system`
        created_pid_list=[]
        for index in range (number_of_particles):
            df_index=int (index_list[index])
            self.clean_df_row(index=df_index)
            if position is None:
                particle_position = self.rng.random((1, 3))[0] *self.np.copy(espresso_system.box_l)
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
            espresso_system(`obj`): Instance of an espresso system object from espressomd library.
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
            self.create_particle(name=name, number_of_particles=number_of_objects, espresso_system=espresso_system, position=position)
        elif pmb_type == 'residue':
            self.create_residue(name=name, number_of_residues=number_of_objects, espresso_system=espresso_system, central_bead_position=position,use_default_bond=use_default_bond)
        elif pmb_type == 'molecule':
            self.create_molecule(name=name, number_of_molecules=number_of_objects, espresso_system=espresso_system, use_default_bond=use_default_bond, list_of_first_residue_positions=position)
        return

    def create_protein(self, name, number_of_proteins, espresso_system, topology_dict):
        """
        Creates `number_of_proteins` molecules of type `name` into `espresso_system` at the coordinates in `positions`

        Args:
            name(`str`): Label of the protein to be created. 
            espresso_system(`obj`): Instance of a system object from the espressomd library.
            number_of_proteins(`int`): Number of proteins to be created.
            positions(`dict`): {'ResidueNumber': {'initial_pos': [], 'chain_id': ''}}
        """

        if number_of_proteins <=0:
            return
        if not self.check_if_name_is_defined_in_df(name=name,
                                                  pmb_type_to_be_defined='protein'):
            raise ValueError(f"{name} must correspond to a name of a pmb_type='protein' defined on df")

        self.copy_df_entry(name=name,column_name='molecule_id',number_of_copies=number_of_proteins)

        protein_index = self.np.where(self.df['name']==name)
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

                residue_name = self.re.split(r'\d+', residue)[0]
                residue_number = self.re.split(r'(\d+)', residue)[1]
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
                                            new_value=residue_number)

                self.add_value_to_df(key=('molecule_id',''),
                                        index=int (index),
                                        new_value=molecule_id)

        return

    def create_residue(self, name, espresso_system, number_of_residues, central_bead_position=None,use_default_bond=False, backbone_vector=None):
        """
        Creates `number_of_residues` residues of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the residue type to be created. `name` must be defined in `pmb.df`
            espresso_system(`obj`): Instance of a system object from espressomd library.
            number_of_residue(`int`): Number of residues of type `name` to be created.
            central_bead_position(`list` of `float`): Position of the central bead.
            use_default_bond (`bool`): Switch to control if a bond of type `default` is used to bond a particle whose bonds types are not defined in `pmb.df`
            backbone_vector (`list` of `float`): Backbone vector of the molecule. All side chains are created perpendicularly to `backbone_vector`.

        Returns:
            residues_info (`dict`): {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
        """
        if number_of_residues <= 0:
            return
        if not self.check_if_name_is_defined_in_df(name=name,
                                                    pmb_type_to_be_defined='residue'):
            raise ValueError(f"{name} must correspond to a label of a pmb_type='residue' defined on df")
        # Copy the data of a residue a `number_of_residues` times in the `df
        self.copy_df_entry(name=name,
                            column_name='residue_id',
                            number_of_copies=number_of_residues)
        residues_index = self.np.where(self.df['name']==name)
        residue_index_list =list(residues_index[0])[-number_of_residues:]
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
            if self.df.loc[self.df['name']==name].central_bead.values[0] is self.np.NaN:
                raise ValueError("central_bead must contain a particle name")        
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
                if side_chain_element not in self.df.values:              
                    raise ValueError (side_chain_element +'is not defined')
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
                                        verbose=False)
                    side_chain_beads_ids.append(side_bead_id)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, side_bead_id))
                    self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=side_bead_id,
                                        use_default_bond=use_default_bond)
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
                    lateral_residue_info = self.create_residue(name=side_chain_element, espresso_system=espresso_system,
                        number_of_residues=1, central_bead_position=[residue_position],use_default_bond=use_default_bond)
                    lateral_residue_dict=list(lateral_residue_info.values())[0]
                    central_bead_side_chain_id=lateral_residue_dict['central_bead_id']
                    lateral_beads_side_chain_ids=lateral_residue_dict['side_chain_ids']
                    residue_id_side_chain=list(lateral_residue_info.keys())[0]
                    # Change the residue_id of the residue in the side chain to the one of the biger residue
                    index = self.df[(self.df['residue_id']==residue_id_side_chain) & (self.df['pmb_type']=='residue') ].index.values[0]
                    self.add_value_to_df(key=('residue_id',''),
                                        index=int(index),
                                        new_value=residue_id, 
                                        verbose=False)
                    # Change the residue_id of the particles in the residue in the side chain
                    side_chain_beads_ids+=[central_bead_side_chain_id]+lateral_beads_side_chain_ids
                    for particle_id in side_chain_beads_ids:
                        index = self.df[(self.df['particle_id']==particle_id) & (self.df['pmb_type']=='particle')].index.values[0]
                        self.add_value_to_df(key=('residue_id',''),
                                            index=int (index),
                                            new_value=residue_id, 
                                            verbose=False)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, central_bead_side_chain_id))
                    self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=central_bead_side_chain_id,
                                        use_default_bond=use_default_bond)
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

            value = float(self.re.split(r'\s+', variable)[0])
            units = self.re.split(r'\s+', variable)[1]
        
        variable_with_units=value*self.units(units)

        return variable_with_units

    def define_AA_particles_in_sequence(self, sequence, sigma_dict=None):
        '''
        Defines in `pmb.df` all the different particles in `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein. 

        Note:
            - It assumes that the names in `sequence` correspond to amino acid names using the standard one letter code.
        '''

        already_defined_AA=[]
        
        for residue_name in sequence:
            if residue_name in already_defined_AA:
                continue
            self.define_particle (name=residue_name, q=0)
                
        if sigma_dict:
            self.define_particles_parameter_from_dict(param_dict = sigma_dict, 
                                                      param_name = 'sigma')
        return 

    def define_AA_residues(self, sequence, model):
        """
        Defines in `pmb.df` all the different residues in `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.

        Returns:
            residue_list (`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.             
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
            bond_type (`str`)        : label to identify the potential to model the bond.
            bond_parameters (`dict`) : parameters of the potential of the bond.
            particle_pairs (`lst`)   : list of the `names` of the `particles` to be bonded.

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

            cutoff = 2**(1.0/6.0) * lj_parameters["sigma"]

            l0 = self.calculate_initial_bond_length(bond_object = bond_object,
                                                    bond_type   = bond_type,
                                                    epsilon     = lj_parameters["epsilon"],
                                                    sigma       = lj_parameters["sigma"],
                                                    cutoff      = cutoff)
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
                                    new_value=self.json.dumps(bond_object.get_params()))

        self.df['parameters_of_the_potential'] = self.df['parameters_of_the_potential'].apply (lambda x: eval(str(x)) if self.pd.notnull(x) else x) 
        return

    
    def define_default_bond(self, bond_type, bond_parameters, epsilon=None, sigma=None, cutoff=None):
        """
        Asigns `bond` in `pmb.df` as the default bond.
        The LJ parameters can be optionally provided to calculate the initial bond length

        Args:
            bond_type (`str`)          : label to identify the potential to model the bond.
            bond_parameters (`dict`)   : parameters of the potential of the bond.
            sigma (`float`, optional)  : LJ sigma of the interaction between the particles
            epsilon (`float`, optional): LJ epsilon for the interaction between the particles
            cutoff (`float`, optional) : cutoff-radius of the LJ interaction

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
        l0 = self.calculate_initial_bond_length(bond_object = bond_object, 
                                                bond_type   = bond_type, 
                                                epsilon     = epsilon,
                                                sigma       = sigma,
                                                cutoff      = cutoff)

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
                            new_value    = self.json.dumps(bond_object.get_params()))
        self.df['parameters_of_the_potential'] = self.df['parameters_of_the_potential'].apply (lambda x: eval(str(x)) if self.pd.notnull(x) else x) 
        return

    def define_epsilon_value_of_particles(self, eps_dict):
        '''
        Defines the epsilon value of the particles in `eps_dict` into the `pmb.df`.

        Args:
            eps_dict(`dict`):  {'name': epsilon}
        '''
        for residue in eps_dict.keys():
            label_list = self.df[self.df['name'] == residue].index.tolist()
            for index in label_list:
                epsilon = eps_dict[residue]
                self.add_value_to_df(key= ('epsilon',''),
                        index=int (index),
                        new_value=epsilon)
        return 

    def define_molecule(self, name, residue_list):
        """
        Defines a pyMBE object of type `molecule` in `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the `molecule`.
            residue_list (`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.  
        """
        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='molecule'):
            return
        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'molecule'
        self.df.at [index,('residue_list','')] = residue_list
        return

    def define_particle(self, name, q=0, acidity='inert', pka=None, sigma=None, epsilon=None, cutoff=None, offset=None,verbose=True):
        """
        Defines the properties of a particle object.

        Args:
            name (`str`): Unique label that identifies this particle type.  
            q (`int`, optional): Permanent charge of this particle type. Defaults to 0.
            acidity (`str`, optional): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            pka (`float`, optional): If `particle` is an acid or a base, it defines its  pka-value. Defaults to None.
            sigma (`pint.Quantity`, optional): Sigma parameter used to set up Lennard-Jones interactions for this particle type. Defaults to None.
            cutoff (`pint.Quantity`, optional): Cutoff parameter used to set up Lennard-Jones interactions for this particle type. Defaults to None.
            offset (`pint.Quantity`, optional): Offset parameter used to set up Lennard-Jones interactions for this particle type. Defaults to None.
            epsilon (`pint.Quantity`, optional): Epsilon parameter used to setup Lennard-Jones interactions for this particle tipe. Defaults to None.
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.

        Note:
            - `sigma`, `cutoff` and `offset` must have a dimensitonality of `[length]` and should be defined using pmb.units.
            - `epsilon` must have a dimensitonality of `[energy]` and should be defined using pmb.units.
            - `cutoff` defaults to `2**(1./6.) reduced_length`. 
            - `offset` defaults to 0.
            - The default setup corresponds to the Weeks−Chandler−Andersen (WCA) model, corresponding to purely steric interactions.
            - For more information on `sigma`, `epsilon`, `cutoff` and `offset` check `pmb.setup_lj_interactions()`.
        """ 

        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='particle'):
            index = self.df[self.df['name']==name].index[0]                                   
        else:
            index = len(self.df)
            self.df.at [index, 'name'] = name
            self.df.at [index,'pmb_type'] = 'particle'
        
        # If `cutoff` and `offset` are not defined, default them to 0

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
                                    verbose=verbose)

        # Define particle acid/base properties
        self.set_particle_acidity (name=name, 
                                acidity = acidity , 
                                default_charge=q, 
                                pka=pka,
                                verbose=verbose)
        return 
    
    def define_particles_parameter_from_dict(self, param_dict, param_name):
        '''
        Defines the `param_name` value of the particles in `param_dict` into the `pmb.df`.

        Args:
            param_dict(`dict`):  {'name': `param_name` value}
        '''
        for residue in param_dict.keys():
            label_list = self.df[self.df['name'] == residue].index.tolist()
            for index in label_list:
                value = param_dict[residue]
                self.add_value_to_df(key= (param_name,''),
                        index=int (index),
                        new_value=value)
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
            residue_list = self.define_AA_residues(sequence=clean_sequence,model=model)

            self.define_molecule(name = name, residue_list=residue_list)
            index = self.df.loc[self.df['name'] == name].index.item() 
            self.df.at [index,'model'] = model
            self.df.at [index,('sequence','')] = clean_sequence
        return
    
    def define_protein(self, name,model, topology_dict):
        """
        Defines a pyMBE object of type `protein` in `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the `protein`.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
            topology_dict (`dict`): {'initial_pos': coords_list, 'chain_id': id, 'sigma': sigma_value}
        """

        protein_seq_list = []     
    
        valid_keys = ['1beadAA','2beadAA']

        if model not in valid_keys:
            raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
        
        if model == '2beadAA':
            self.define_particle(name='CA')

        sigma_dict = {}
        
        for residue in topology_dict.keys():
            residue_name = self.re.split(r'\d+', residue)[0]
            
            if residue_name not in sigma_dict.keys():
                sigma_dict [residue_name] =  topology_dict[residue]['sigma']
            if residue_name not in ('CA', 'Ca'):
                protein_seq_list.append(residue_name)                  

        protein_sequence = ''.join(protein_seq_list)
        clean_sequence = self.protein_sequence_parser(sequence=protein_sequence)

        self.define_AA_particles_in_sequence (sequence=clean_sequence, sigma_dict = sigma_dict)
        residue_list = self.define_AA_residues(sequence=clean_sequence, model=model)

        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'protein'
        self.df.at [index,'model'] = model
        self.df.at [index,('sequence','')] = clean_sequence  
        self.df.at [index,('residue_list','')] = residue_list    

        return 
    
    def define_residue(self, name, central_bead, side_chains):
        """
        Defines a pyMBE object of type `residue` in `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the `residue`.
            central_bead (`str`): `name` of the `particle` to be placed as central_bead of the `residue`.
            side_chains (`list` of `str`): List of `name`s of the pmb_objects to be placed as side_chains of the `residue`. Currently, only pmb_objects of type `particle`s or `residue`s are supported.
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
            name (str): Label of the pmb object to be destroyed. The pmb object must be defined in `pymbe.df`.
            espresso_system (cls): Instance of a system class from espressomd library.

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
                particle_ids_list = self.df.loc[self.df['residue_id']==residue_id].particle_id.dropna().to_list()
                self.df = self.df.drop(self.df[self.df['residue_id'] == residue_id].index)
                for particle_id in particle_ids_list:
                    espresso_system.part.by_id(particle_id).remove()
                    self.df= self.df.drop(self.df[self.df['particle_id']==particle_id].index)    
        if pmb_type == 'molecule':
            molecules_id = self.df.loc[self.df['name']== name].molecule_id.to_list()
            for molecule_id in molecules_id:
                particle_ids_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
                
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
            pH_res ('float'): pH-value in the reservoir.
            c_salt_res ('float'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            activity_coefficient_monovalent_pair ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.

        Returns:
            cH_res ('float'): Concentration of H+ ions.
            cOH_res ('float'): Concentration of OH- ions.
            cNa_res ('float'): Concentration of Na+ ions.
            cCl_res ('float'): Concentration of Cl- ions.
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
        determined_pH = -self.np.log10(cH_res.to('mol/L').magnitude * self.np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))

        while abs(determined_pH-pH_res)>1e-6:
            if determined_pH > pH_res:
                cH_res *= 1.005
            else:
                cH_res /= 1.003
            cH_res, cOH_res, cNa_res, cCl_res = determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res)
            ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
            determined_pH = -self.np.log10(cH_res.to('mol/L').magnitude * self.np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))
            self_consistent_run=0

        return cH_res, cOH_res, cNa_res, cCl_res

    def filter_df(self, pmb_type):
        """
        Filters `pmb.df` and returns a sub-set of it containing only rows with pmb_object_type=`pmb_type` and non-NaN columns.
        
        Args:
            pmb_type(`str`): pmb_object_type to filter in `pmb.df`.

        Returns:
            pmb_type_df(`obj`): filtered `pmb.df`.
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
        idx = self.pd.IndexSlice
        for state in ['state_one', 'state_two']:
            index = self.np.where(self.df[(state, 'es_type')] == es_type)[0]      
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
            if self.np.linalg.norm(coord-self.np.asarray(center))>=radius:
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
        center=self.np.array(center)
        d = center.shape[0]
        # sample n_samples points in d dimensions from a standard normal distribution
        samples = self.rng.normal(size=(n_samples, d))
        # make the samples lie on the surface of the unit hypersphere
        normalize_radii = self.np.linalg.norm(samples, axis=1)[:, self.np.newaxis]
        samples /= normalize_radii
        if not on_surface:
            # make the samples lie inside the hypersphere with the correct density
            uniform_points = self.rng.uniform(size=n_samples)[:, self.np.newaxis]
            new_radii = self.np.power(uniform_points, 1/d)
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
        np_vec = self.np.array(vector) 
        np_vec /= self.np.linalg.norm(np_vec) 
        if self.np.all(np_vec == 0):
            raise ValueError('Zero vector')
        # Generate a random vector 
        random_vector = self.generate_random_points_in_a_sphere(radius=1, 
                                                                center=[0,0,0],
                                                                n_samples=1, 
                                                                on_surface=True)[0]
        # Project the random vector onto the input vector and subtract the projection
        projection = self.np.dot(random_vector, np_vec) * np_vec
        perpendicular_vector = random_vector - projection
        # Normalize the perpendicular vector to have the same magnitude as the input vector
        perpendicular_vector /= self.np.linalg.norm(perpendicular_vector) 
        return perpendicular_vector*magnitude

    def get_bond_length(self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for bonds between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns the initial bond length.
        If `use_default_bond` is activated and a "default" bond is defined, returns the length of that default bond instead.
        If no bond is found, it prints a message and it does not return anything. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1 (str): label of the type of the first particle type of the bonded particles.
            particle_name2 (str): label of the type of the second particle type of the bonded particles.
            hard_check (bool, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond (bool, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            l0 (`float`): bond length
        
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
                raise KeyboardInterrupt
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
        state_one = self.pd.Series (df_state_one.charge.values,index=df_state_one.es_type.values)
        state_two = self.pd.Series (df_state_two.charge.values,index=df_state_two.es_type.values)
        charge_map  = self.pd.concat([state_one,state_two],axis=0).to_dict()
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
            {"epsilon": LJ epsilon, "sigma": LJ sigma}

        Note:
            Currently, the only `combining_rule` supported is Lorentz-Berthelot.
        """
        supported_combining_rules=["Lorentz-Berthelot"]

        if combining_rule not in supported_combining_rules:
            raise ValueError(f"Combining_rule {combining_rule} currently not implemented in pyMBE, valid keys are {supported_combining_rules}")

        eps1 = self.df.loc[self.df["name"]==particle_name1].epsilon.iloc[0]
        sigma1 = self.df.loc[self.df["name"]==particle_name1].sigma.iloc[0]
        eps2 = self.df.loc[self.df["name"]==particle_name2].epsilon.iloc[0]
        sigma2 = self.df.loc[self.df["name"]==particle_name2].sigma.iloc[0]
        
        return {"epsilon": self.np.sqrt(eps1*eps2), "sigma": (sigma1+sigma2)/2.0}

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
                res_list=self.df.loc[self.df['residue_id']== res_id].particle_id.dropna().tolist()
                res_map[res_id]=res_list
            return res_map
        if object_type in ['molecule', 'protein']:
            mol_ids = self.df.loc[self.df['name']== object_name].molecule_id.dropna().tolist()
            for mol_id in mol_ids:
                res_ids = set(self.df.loc[self.df['molecule_id']== mol_id].residue_id.dropna().tolist())
                res_map=do_res_map(res_ids=res_ids)    
                mol_list=self.df.loc[self.df['molecule_id']== mol_id].particle_id.dropna().tolist()
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
        '''

        df_state_one = self.df[[('sigma',''),('state_one','es_type')]].dropna().drop_duplicates()
        df_state_two = self.df[[('sigma',''),('state_two','es_type')]].dropna().drop_duplicates()

        state_one = self.pd.Series(df_state_one.sigma.values/2.0,index=df_state_one.state_one.es_type.values)
        state_two = self.pd.Series(df_state_two.sigma.values/2.0,index=df_state_two.state_two.es_type.values)
        radius_map  = self.pd.concat([state_one,state_two],axis=0).to_dict()
        
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
        state_one = self.pd.Series (df_state_one.es_type.values,index=df_state_one.label)
        state_two = self.pd.Series (df_state_two.es_type.values,index=df_state_two.label)
        type_map  = self.pd.concat([state_one,state_two],axis=0).to_dict()
        return type_map

    def load_interaction_parameters(self, filename, verbose=True):
        """
        Loads the interaction parameters stored in `filename` into `pmb.df`
        
        Args:
            filename(`str`): name of the file to be read
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.
        """
        without_units = ['q','es_type','acidity']
        with_units = ['sigma','epsilon']

        with open(filename, 'r') as f:
            interaction_data = self.json.load(f)
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
                                acidity=not_requiered_attributes.pop('acidity'),
                                epsilon=not_requiered_attributes.pop('epsilon'),
                                verbose=verbose)
            elif object_type == 'residue':
                self.define_residue (name = param_dict.pop('name'),
                                    central_bead = param_dict.pop('central_bead_name'),
                                    side_chains = param_dict.pop('side_chains_names'))
            elif object_type == 'molecule':
                self.define_molecule(name=param_dict.pop('name'),
                                    residue_list=param_dict.pop('residue_name_list'))
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
                self.define_bond(bond_type=bond_type, bond_parameters=bond, particle_pairs=particle_pairs)
            else:
                raise ValueError(object_type+' is not a known pmb object type')
            
        return
    
    def load_pka_set(self, filename, verbose=True):
        """
        Loads the pka_set stored in `filename` into `pmb.df`.
        
        Args:
            filename(`str`): name of the file with the pka set to be loaded. Expected format is {name:{"acidity": acidity, "pka_value":pka_value}}.
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.

        Note:
            - If `name` is already defined in the `pymbe.df`, it prints a warning.
        """
        with open(filename, 'r') as f:
            pka_data = self.json.load(f)
            pka_set = pka_data["data"]

        self.check_pka_set(pka_set=pka_set)

        for key in pka_set:
            acidity = pka_set[key]['acidity']
            pka_value = pka_set[key]['pka_value']
            self.define_particle(
                    name=key,
                    acidity=acidity, 
                    pka=pka_value,
                    verbose=verbose)
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
        print(unit_length.to('nm'), "=", unit_length)
        print(unit_energy.to('J'), "=", unit_energy)
        print('Temperature:', (self.kT/self.Kb).to("K"))
        print(unit_charge.to('C'), "=", unit_charge)
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
            filename(str): path to file.

        Note:
            This function only accepts files with CSV format. 
        """
        
        if filename.rsplit(".", 1)[1] != "csv":
            raise ValueError("Only files with CSV format are supported")
        df = self.pd.read_csv (filename,header=[0, 1], index_col=0)
        columns_names = self.setup_df()
        
        multi_index = self.pd.MultiIndex.from_tuples(columns_names)
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
                        sigma = 2*radius
                        particles_dict [int(atom_id)] = [atom_name , atom_resname, chain_id, sigma]
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
                                                    'sigma':numbered_label[i][2] }

        return topology_dict

    def search_bond(self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for bonds between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        If `use_default_bond` is activated and a "default" bond is defined, returns that default bond instead.
        If no bond is found, it prints a message and it does not return anything. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1 (str): label of the type of the first particle type of the bonded particles.
            particle_name2 (str): label of the type of the second particle type of the bonded particles.
            hard_check (bool, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond (bool, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            bond (cls): bond object from the espressomd library.
        
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
                exit()
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

    def set_particle_acidity(self, name, acidity='inert', default_charge=0, pka=None, verbose=True):
        """
        Sets the particle acidity if it is acidic or basic, creates `state_one` and `state_two` with the protonated and 
        deprotonated states. In each state is set: `label`,`charge` and `es_type`. If it is inert, it will define only `state_one`.

        Args:
            name (`str`): Unique label that identifies the `particle`. 
            acidity (`str`): Identifies whether the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            default_charge (`int`): Charge of the particle. Defaults to 0.
            pka (`float`, optional):  If `particle` is an acid or a base, it defines its pka-value. Defaults to None.
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.
        """
        acidity_valid_keys = ['inert','acidic', 'basic']
        if acidity not in acidity_valid_keys:
            raise ValueError(f"Acidity {acidity} provided for particle name  {name} is not supproted. Valid keys are: {acidity_valid_keys}")
        if acidity in ['acidic', 'basic'] and pka is None:
            raise ValueError(f"pKa not provided for particle with name {name} with acidity {acidity}. pKa must be provided for acidic or basic particles.")   
        for index in self.df[self.df['name']==name].index:       
            if pka:
                self.add_value_to_df(key=('pka',''),
                                    index=index,
                                    new_value=pka, 
                                    verbose=verbose)
            self.add_value_to_df(key=('acidity',''),
                                 index=index,
                                 new_value=acidity, 
                                 verbose=verbose) 
            if not self.check_if_df_cell_has_a_value(index=index,key=('state_one','es_type')):
                self.add_value_to_df(key=('state_one','es_type'),
                                     index=index,
                                     new_value=self.propose_unused_type(), 
                                     verbose=verbose)  
            if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'inert':
                self.add_value_to_df(key=('state_one','charge'),
                                     index=index,
                                     new_value=default_charge, 
                                     verbose=verbose)
                self.add_value_to_df(key=('state_one','label'),
                                     index=index,
                                     new_value=name, 
                                     verbose=verbose)
            else:
                protonated_label = f'{name}H'
                self.add_value_to_df(key=('state_one','label'),
                                     index=index,
                                     new_value=protonated_label, 
                                     verbose=verbose)
                self.add_value_to_df(key=('state_two','label'),
                                     index=index,
                                     new_value=name, 
                                     verbose=verbose)
                if not self.check_if_df_cell_has_a_value(index=index,key=('state_two','es_type')):
                    self.add_value_to_df(key=('state_two','es_type'),
                                         index=index,
                                         new_value=self.propose_unused_type(), 
                                         verbose=verbose)
                if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'acidic':        
                    self.add_value_to_df(key=('state_one','charge'),
                                         index=index,new_value=0, 
                                         verbose=verbose)
                    self.add_value_to_df(key=('state_two','charge'),
                                         index=index,
                                         new_value=-1, 
                                         verbose=verbose)
                elif self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'basic':
                    self.add_value_to_df(key=('state_one','charge'),
                                         index=index,new_value=+1, 
                                         verbose=verbose)
                    self.add_value_to_df(key=('state_two','charge'),
                                         index=index,
                                         new_value=0, 
                                         verbose=verbose)   
        return
    
    def set_reduced_units(self, unit_length=None, unit_charge=None, temperature=None, Kw=None, verbose=True):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args:
            unit_length (`obj`,optional): Reduced unit of length defined using the `pmb.units` UnitRegistry. Defaults to None. 
            unit_charge (`obj`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            temperature (`obj`,optional): Temperature of the system, defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw (`obj`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 
            verbose (`bool`, optional): Switch to activate/deactivate verbose. Defaults to True.

        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
        """
        self.units=self.pint.UnitRegistry()
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
            exclusion_range(`float`, optional): Below this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.

        Returns:
            RE (`obj`): Instance of a reaction_ensemble.ConstantpHEnsemble object from the espressomd library.
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
            c_salt_res ('float'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`float`, optional): For distances shorter than this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.

        Returns:
            RE (`obj`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
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
            c_salt_res ('float'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            proton_name ('str'): Name of the proton (H+) particle.
            hydroxide_name ('str'): Name of the hydroxide (OH-) particle.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`float`, optional): For distances shorter than this value, no particles will be inserted.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.

        Returns:
            RE (`obj`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
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
            c_salt_res ('float'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            cation_name ('str'): Name of the cationic particle.
            anion_name ('str'): Name of the anionic particle.
            activity_coefficient ('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`float`, optional): Below this value, no particles will be inserted.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius per each espresso_type. Defaults to `False`.

        Returns:
            RE (`obj`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
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
        a_cation = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * self.np.sqrt(determined_activity_coefficient)
        a_anion = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * self.np.sqrt(determined_activity_coefficient)
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
        
        self.df = self.pd.DataFrame(columns=self.pd.MultiIndex.from_tuples([(col_main, col_sub) for col_main, sub_cols in columns_dtypes.items() for col_sub in sub_cols.keys()]))
        
        for level1, sub_dtypes in columns_dtypes.items():
            for level2, dtype in sub_dtypes.items():
                self.df[level1, level2] = self.df[level1, level2].astype(dtype)
                
        columns_names = self.pd.MultiIndex.from_frame(self.df)
        columns_names = columns_names.names
                
        return columns_names

    def setup_lj_interactions(self, espresso_system, shift_potential=True, combining_rule='Lorentz-Berthelot', warnings=True):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle types with values for `sigma`, `offset`, and `epsilon` stored in `pymbe.df`.

        Args:
            espresso_system(`obj`): Instance of a system object from the espressomd library.
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
                check_list.append(self.np.isnan(value_in_df))
            if any(check_list):
                non_parametrized_labels.append(self.find_value_from_es_type(es_type=particle_type, 
                                                                            column_name='label'))
            else:
                particles_types_with_LJ_parameters.append(particle_type)
        lj_parameters_keys=["label","sigma","epsilon","offset","cutoff"]
        # Set up LJ interactions between all particle types
        for type_pair in combinations_with_replacement(particles_types_with_LJ_parameters, 2):
            lj_parameters={}
            for key in lj_parameters_keys:
                lj_parameters[key]=[]
            # Search the LJ parameters of the type pair
            for ptype in type_pair:
                for key in lj_parameters_keys:
                    lj_parameters[key].append(self.find_value_from_es_type(es_type=ptype, column_name=key))
            # If one of the particle has sigma=0, no LJ interations are set up between that particle type and the others    
            if not all(sigma_value.magnitude for sigma_value in lj_parameters["sigma"]):
                continue
            # Apply combining rule
            if combining_rule == 'Lorentz-Berthelot':
                sigma=(lj_parameters["sigma"][0]+lj_parameters["sigma"][1])/2
                cutoff=(lj_parameters["cutoff"][0]+lj_parameters["cutoff"][1])/2
                offset=(lj_parameters["offset"][0]+lj_parameters["offset"][1])/2
                epsilon=self.np.sqrt(lj_parameters["epsilon"][0]*lj_parameters["epsilon"][1])
            espresso_system.non_bonded_inter[type_pair[0],type_pair[1]].lennard_jones.set_params(epsilon = epsilon.to('reduced_energy').magnitude, 
                                                                                    sigma = sigma.to('reduced_length').magnitude, 
                                                                                    cutoff = cutoff.to('reduced_length').magnitude,
                                                                                    offset = offset.to("reduced_length").magnitude, 
                                                                                    shift = shift)                                                                                          
            index = len(self.df)
            self.df.at [index, 'name'] = f'LJ: {lj_parameters["label"][0]}-{lj_parameters["label"][1]}'
            lj_params=espresso_system.non_bonded_inter[type_pair[0], type_pair[1]].lennard_jones.get_params()

            self.add_value_to_df(index=index,
                                key=('pmb_type',''),
                                new_value='LennardJones')
            
            self.add_value_to_df(index=index,
                                key=('parameters_of_the_potential',''),
                                new_value=self.json.dumps(lj_params))                      
        if non_parametrized_labels and warnings:
            print(f'WARNING: The following particles do not have a defined value of sigma or epsilon in pmb.df: {non_parametrized_labels}. No LJ interaction has been added in ESPResSo for those particles.')
            
        self.df['parameters_of_the_potential'] = self.df['parameters_of_the_potential'].apply (lambda x: eval(str(x)) if self.pd.notnull(x) else x) 
        return

    def setup_particle_sigma(self, topology_dict):
        '''
        Sets up sigma of the particles in `topology_dict`.

        Args:
            topology_dict(`dict`): {'initial_pos': coords_list, 'chain_id': id, 'sigma': sigma_value}
        '''
        for residue in topology_dict.keys():
            residue_name = self.re.split(r'\d+', residue)[0]
            residue_number = self.re.split(r'(\d+)', residue)[1]
            residue_sigma  = topology_dict[residue]['sigma']
            sigma = residue_sigma*self.units.nm
            index = self.df[(self.df['residue_id']==residue_number) & (self.df['name']==residue_name) ].index.values[0]
            self.add_value_to_df(key= ('sigma',''),
                        index=int (index),
                        new_value=sigma)           
        return 

    def write_pmb_df (self, filename):
        '''
        Writes the pyMBE dataframe into a csv file
        
        Args:
            filename (`str`): Path to the csv file 
        '''

        self.df.to_csv(filename)

        return

    def write_output_vtf_file(self, espresso_system, filename):
        '''
        Writes a snapshot of `espresso_system` on the vtf file `filename`.

        Args:
            espresso_system(`obj`): Instance of a system object from the espressomd library.
            filename(`str`): Path to the vtf file.

        '''
        box = espresso_system.box_l[0]
        with open(filename, mode='w+t') as coordinates:
            coordinates.write (f'unitcell {box} {box} {box} \n')
            for particle in espresso_system.part: 
                type_label = self.find_value_from_es_type(es_type=particle.type, column_name='label')
                coordinates.write (f'atom {particle.id} radius 1 name {type_label} type {type_label}\n' )
            coordinates.write ('timestep indexed\n')
            for particle in espresso_system.part:
                coordinates.write (f'{particle.id} \t {particle.pos[0]} \t {particle.pos[1]} \t {particle.pos[2]}\n')
        return 
