class pymbe_library():
    """
    The library for the Molecular Brewer for ESPResSo (pyMBE)

    Attributes:
        N_A(`obj`): Avogadro number using the `pmb.units` UnitRegistry.
        Kb(`obj`): Boltzmann constant using the `pmb.units` UnitRegistry.
        e(`obj`): Elemental charge using the `pmb.units` UnitRegistry.
        df(`obj`): PandasDataframe used to bookkeep all the information stored in pyMBE. Typically refered as `pmb.df`. 
        kT(`obj`): Thermal energy using the `pmb.units` UnitRegistry. it is used as the unit of reduced energy.
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

    def __init__(self, temperature=None, unit_length=None, unit_charge=None):
        """
        Initializes the pymbe_library by setting up the reduced unit system with `temperature` and `reduced_length` 
        and sets up  the `pmb.df` for bookkepping.

        Args:
            temperature(`obj`,optional): Value of the temperature in the pyMBE UnitRegistry. Defaults to None.
            unit_length(`obj`, optional): Value of the unit of length in the pyMBE UnitRegistry. Defaults to None.
            unit_charge (`obj`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
        
        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
        """
        # Default definitions of reduced units
        if temperature is None:
            temperature= 298.15 * self.units.K
        if unit_length is None:
            unit_length= 0.355 * self.units.nm
        if unit_charge is None:
            unit_charge=self.units.e
        self.kT=temperature*self.Kb
        self.units.define(f'reduced_energy = {self.kT}')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = 1*e')
        self.setup_df()
        return
    
    def activate_motion_of_rigid_object (self, name, espresso_system):
        '''
        Activates the motion of the rigid object `name` in the `espresso_system`.

        Args:
            name(`str`): Label of the object.
            espresso_system(`obj`): Instance of a system object from the espressomd library.

        Note:
            - It requires that espressodmd has the following feautures activated: ["VIRTUAL_SITES_RELATIVE", "MASS"].
        '''
        print ('activate_motion_of_rigid_object requires that espressodmd has the following feautures activated: ["VIRTUAL_SITES_RELATIVE", "MASS"]')
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

    def add_bond_in_df (self,particle_id1,particle_id2, use_default_bond=False):
        """
        Adds a bond entry on the `pymbe.df` storing the particle_ids of the two bonded.

        Args:
            particle_id1 (`int`): particle_id of the type of the first particle type of the bonded particles
            particle_id2 (`int`): particle_id of the type of the second particle type of the bonded particles
            use_default_bond (`bool`, optional): Controls if a bond of type `default` is used to bond particle whose bonds types are not defined in `pmb.df`. Defaults to False.
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
        if len(used_bond_df.index) > 1: 
            used_bond_df = used_bond_df.drop_duplicates(keep='first')
        used_bond_index = used_bond_df.index.to_list()
        for index in index_list:
            if index not in used_bond_index:
                self.clean_df_row(index=int(index))
                self.df.at[index,'particle_id'] = particle_id1
                self.df.at[index,'particle_id2'] = particle_id2
                break
        return

    def add_bonds_to_espresso (self, espresso_system) :
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

    def add_value_to_df(self,index,key,new_value, warning=True):
        """
        Adds a value to a cell in the `pmb.df` DataFrame.

        Args:

            index(`int`): index of the row to add the value to.
            key(`str`): the column label to add the value to.
            warning(`bool`, optional): If true, prints a warning if a value is being overwritten in `pmb.df`. Defaults to true.
        """
        # Make sure index is a scalar integer value
        index = int (index)
        assert isinstance(index, int), '`index` should be a scalar integer value.'
        idx = self.pd.IndexSlice
        if self.check_if_df_cell_has_a_value(index=index,key=key):
            old_value= self.df.loc[index,idx[key]]
            if warning:
                if old_value != new_value:
                    name=self.df.loc[index,('name','')]
                    pmb_type=self.df.loc[index,('pmb_type','')]
                    print(f"WARNING: you are redefining the properties of {name} of pmb_type {pmb_type}")
                    print(f'WARNING: overwritting the value of the entry `{key}`: old_value = {old_value} new_value = {new_value}')
        self.df.loc[index,idx[key]] = new_value
        return
    
    def assign_molecule_id (self, name, molecule_index, pmb_type, used_molecules_id):
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
            pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]                
            if not check_residue_name.empty and pmb_type == pmb_type :              
                for value in check_residue_name.index.to_list():                  
                    if value not in used_molecules_id:                              
                        molecule_id = self.df.loc[value].molecule_id.values[0]                    
                        break
            else:
                molecule_id = self.df['molecule_id'].max() +1

        self.add_value_to_df (key=('molecule_id',''),
                                index=int(molecule_index),
                                new_value=molecule_id, 
                                warning=False)

        return molecule_id
    
    def calculate_center_of_mass_of_molecule (self,molecule_id, espresso_system):
        """
        Calculates the center of mass of type `name`.

        Args:
            molecule_id(`int`): Id of the molecule to be centered.
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

    def calculate_HH(self, sequence, pH_list=None, pka_set=None):
        """
        Calculates the ideal Henderson-Hasselbach titration curve of `sequence`.

        Args:
            sequence(`lst`): List of `name`s of particles.
            pH_list(`lst`): pH-values to calculate. 
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}

        Returns:
            Z_HH (`lst`): Henderson-Hasselbach prediction of the charge of `sequence` in `pH_list`

        Note:
            - If no `pH_list` is given, 50 equispaced pH-values ranging from 2 to 12 are calculated
            - If no `pka_set` is given, the pKa values are taken from `pmb.df`
        """
        if pH_list is None:
            pH_list=self.np.linspace(2,12,50)    
        if pka_set is None:
            pka_set=self.get_pka_set() 
        self.check_pka_set(pka_set=pka_set)
        Z_HH=[]
        for pH_value in pH_list:    
            Z=0
            for name in sequence: 
                if name in pka_set.keys():
                    if pka_set[name]['acidity'] == 'acidic':
                        psi=-1
                    elif pka_set[name]['acidity']== 'basic':
                        psi=+1
                    else:
                        psi=0
                    Z+=psi/(1+10**(psi*(pH_value-pka_set[name]['pka_value'])))                      
            Z_HH.append(Z)
        return Z_HH

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

    def calculate_net_charge_in_molecules(self,espresso_system, object_name):
        '''
        Calculates the net charge per molecule of molecules  with `name` = object_name. 
        Returns the net charge per molecule and a maps with the net charge per residue and molecule.

        Args:
            espresso_system: system information 
            object_name (str): name of the molecule to calculate de net charge

        Returns:
            {"mean": mean_net_charge, "molecules": {mol_id: net_charge_of_mol_id, }, "residues": {res_id: net_charge_of_res_id, }}

        Note:
            - The net charge of the molecule is averaged over all molecules of type `name` 
            - The net charge of each particle type is averaged over all particle of the same type in all molecules of type `name`
        '''        
        valid_pmb_types = ["molecule", "protein"]
        pmb_type=self.df.loc[self.df['name']==object_name].pmb_type.values[0]
        if pmb_type not in valid_pmb_types:
            raise ValueError("The pyMBE object with name {object_name} has a pmb_type {pmb_type}. This function only supports pyMBE types {valid_pmb_types}")      
        charge_in_residues = {}
        id_map = self.get_particle_id_map(object_name=object_name)
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

    def center_molecule_in_simulation_box (self, molecule_id, espresso_system):
        """
        Centers the pmb object matching `molecule_id` in the center of the simulation box in `espresso_md`.
        
        Args:
            molecule_id(`int`): Id of the molecule to be centered.
            espresso_system(`obj`): Instance of a system object from the espressomd library.
        """
        center_of_mass = self.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)
        box_center = [espresso_system.box_l[0]/2.0]*3
        pmb_type = self.df.loc[self.df['molecule_id']==molecule_id].pmb_type.values[0]
        pmb_objects = ['protein','molecule','peptide']
        if pmb_type in pmb_objects:
            particle_id_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
            for pid in particle_id_list:
                es_pos = espresso_system.part.by_id(pid).pos
                espresso_system.part.by_id(pid).pos = es_pos - center_of_mass + box_center
        return 

    def check_if_df_cell_has_a_value(self,index,key):
        """
        Checks if a cell in the `pmb.df` at the specified index and column has a value.

        Args:
            index(`int`): Index of the row to check.
            key(`str`): Column label to check.

        Returns:
            `bool`: `True` if the cell has a value, `False` otherwise.
        """
        idx = self.pd.IndexSlice
        return not self.pd.isna(self.df.loc[index, idx[key]])

    def check_if_name_is_defined_in_df (self, name, pmb_type_to_be_defined):
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
            if  current_object_type != pmb_type_to_be_defined:
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

    def clean_df_row(self, index, columns_keys_to_clean=["particle_id", "particle_id2", "residue_id", "molecule_id"]):
        """
        Cleans the columns of `pmb.df` in `columns_keys_to_clean` of the row with index `index` by asigning them a np.nan value.

        Args:
            index(`int`): Index of the row to clean.
            columns_keys_to_clean(`list` of `str`, optional): List with the column keys to be cleaned. Defaults to [`particle_id`, `particle_id2`, `residue_id`, `molecule_id`].
        """   
        for column_key in columns_keys_to_clean:
            self.add_value_to_df(key=(column_key,''),index=index,new_value=self.np.nan, warning=False)
        return

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

    def create_added_salt_in_espresso(self, espresso_system, cation_name, anion_name, c_salt):    
        """
        Creates a `c_salt` concentration of `cation_name` and `anion_name` ions into the `espresso_system`.

        Args:
            espresso_system (`obj`): instance of a espresso system object.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            c_salt (float): Salt concentration. 
            
        Returns:
            c_salt_calculated (float): Calculated salt concentration added to `espresso_system`.
        """
        cation_name_charge = self.df.loc[self.df['name']==cation_name].state_one.charge.values[0]
        anion_name_charge = self.df.loc[self.df['name']==anion_name].state_one.charge.values[0]     
        if cation_name_charge <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ',cation_name_charge)
        if anion_name_charge >= 0:
            raise ValueError('ERROR anion charge must be positive, charge ', anion_name_charge)
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
        self.create_particle_in_espresso(espresso_system=espresso_system, name=cation_name, number_of_particles=N_cation)
        self.create_particle_in_espresso(espresso_system=espresso_system, name=anion_name, number_of_particles=N_anion)
        print('\n Added salt concentration of ', c_salt_calculated.to('mol/L'), 'given by ', N_cation, ' cations and ', N_anion, ' anions')
        return c_salt_calculated

    def create_counterions_in_espresso(self, object_name, cation_name, anion_name, espresso_system):
        """
        Creates particles of `cation_name` and `anion_name` in `espresso_system` to counter the net charge of `pmb_object`.
        
        Args:
            object_name(`str`): `name` of a pymbe object.
            espresso_system(`obj`): Instance of a system object from the espressomd library.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.

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
        if (object_charge['positive'] % abs(anion_charge) == 0):
            counterion_number[anion_name]=int(object_charge['positive']/abs(anion_charge))
        else:
            raise ValueError('The number of positive charges in the pmb_object must be divisible by the  charge of the anion')
        if (object_charge['negative'] % abs(cation_charge) == 0):
            counterion_number[cation_name]=int(object_charge['negative']/cation_charge)
        else:
            raise ValueError('The number of negative charges in the pmb_object must be divisible by the  charge of the cation')
        if counterion_number[cation_name] > 0: 
            self.create_particle_in_espresso(espresso_system=espresso_system, name=cation_name, number_of_particles=counterion_number[cation_name])
        else:
            counterion_number[cation_name]=0
        if counterion_number[anion_name] > 0:
            self.create_particle_in_espresso(espresso_system=espresso_system, name=anion_name, number_of_particles=counterion_number[anion_name])
        else:
            counterion_number[anion_name] = 0
        print('The following counter-ions have been created: ')
        for name in counterion_number.keys():
            print(f'Ion type: {name} created number: {counterion_number[name]}')
        return counterion_number


    def create_molecule_in_espresso(self, name, number_of_molecules, espresso_system, first_residue_position=None, use_default_bond=False):
        """
        Creates `number_of_molecules` molecule of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the molecule type to be created. `name` must be defined in `pmb.df`
            espresso_system(`obj`): Instance of a system object from espressomd library.
            number_of_molecules(`int`): Number of molecules of type `name` to be created.
            first_residue_position(`list`, optional): coordinates where the first_residue_position will be created, random by default
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle with undefined bonds in `pymbe.df`
        Returns:
            molecules_info (`dict`):  {molecule_id: {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids": [particle_id1, ...]}}} 
        """
        if number_of_molecules <= 0:
            return
        if not self.check_if_name_is_defined_in_df(name=name,
                                                    pmb_type_to_be_defined='molecule'):
            raise ValueError(f"{name} must correspond to a label of a pmb_type='molecule' defined on df")
        first_residue = True
        molecule_info = []
        residue_list = self.df[self.df['name']==name].residue_list.values [0]

        self.copy_df_entry(name=name,
                        column_name='molecule_id',
                        number_of_copies=number_of_molecules)
        
        molecules_index = self.np.where(self.df['name']==name)
        molecule_index_list =list(molecules_index[0])[-number_of_molecules:]
        used_molecules_id = self.df.molecule_id.dropna().drop_duplicates().tolist()

        for molecule_index in molecule_index_list:         

            molecule_id = self.assign_molecule_id (name=name,pmb_type='molecule',used_molecules_id=used_molecules_id,molecule_index=molecule_index)

            for residue in residue_list:
                if first_residue:
                    residue_position = first_residue_position
                    backbone_vector = self.generate_trialvectors(center=[0,0,0], 
                                                                 radius=1, 
                                                                 n_samples=1,
                                                                 on_surface=True)[0]
                    residues_info = self.create_residue_in_espresso(name=residue,
                                                                    number_of_residues=1, 
                                                                    espresso_system=espresso_system, 
                                                                    central_bead_position=first_residue_position,  
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
                    residues_info = self.create_residue_in_espresso(name=residue, 
                                                                    number_of_residues=1, 
                                                                    espresso_system=espresso_system, 
                                                                    central_bead_position=[residue_position],
                                                                    use_default_bond= use_default_bond)
                    residue_id = next(iter(residues_info))      
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        if not self.check_if_df_cell_has_a_value(index=index,key=('molecule_id','')):
                            self.df.at[index,'molecule_id'] = molecule_id
                            self.add_value_to_df(key=('molecule_id',''),
                                                index=int (index),
                                                new_value=molecule_id,
                                                warning=False)            
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, previous_residue_id))
                    self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=previous_residue_id,
                                        use_default_bond=use_default_bond)            
                    previous_residue_id = central_bead_id
                    previous_residue = residue                    
                molecule_info.append(residues_info)
            first_residue = True
        return molecule_info
    
    def create_particle_in_espresso(self, name, espresso_system, number_of_particles, position=None, fix=False):
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
                particle_position = self.np.random.random((1, 3))[0] *self.np.copy(espresso_system.box_l)
            else:
                particle_position = position[index]
            if len(espresso_system.part.all()) == 0:
                bead_id = 0
            else:
                bead_id = max (espresso_system.part.all().id) + 1
            created_pid_list.append(bead_id)
            
            espresso_system.part.add (id=bead_id, pos = particle_position, type = es_type, q = q,fix =[fix,fix,fix])        
            self.add_value_to_df(key=('particle_id',''),index=df_index,new_value=bead_id, warning=False)                  
        return created_pid_list

    def create_pmb_object_in_espresso(self, name, number_of_objects, espresso_system, position=None, use_default_bond=False):
        """
        Creates all `particle`s associated to `pmb object` into  `espresso` a number of times equal to `number_of_objects`.
        
        Args:
            name(`str`): Unique label of the `pmb object` to be created. 
            number_of_objects(`int`): Number of `pmb object`s to be created.
            espresso_system(`obj`): Instance of a espresso system object from espressomd library.
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
            self.create_particle_in_espresso(name=name, number_of_particles=number_of_objects, espresso_system=espresso_system, position=position)
        elif pmb_type == 'residue':
            self.create_residue_in_espresso(name=name,number_of_residues=number_of_objects, espresso_system=espresso_system, central_bead_position=position,use_default_bond=use_default_bond)
        elif pmb_type == 'molecule':
            self.create_molecule_in_espresso(name=name,number_of_molecules=number_of_objects, espresso_system=espresso_system, use_default_bond=use_default_bond, first_residue_position=position)
        return

    def create_protein_in_espresso(self, name, number_of_proteins, espresso_system, positions):
        """
        Creates `number_of_proteins` molecule of type `name` into `espresso_system` in the coordinates in `positions`

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
   
            for residue in positions.keys():

                residue_name = self.re.split(r'\d+', residue)[0]
                residue_number = self.re.split(r'(\d+)', residue)[1]
                residue_position = positions[residue]['initial_pos']
                position = residue_position + protein_center

                particle_id = self.create_particle_in_espresso(name=residue_name,
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

    def create_residue_in_espresso (self, name, espresso_system, number_of_residues, central_bead_position=None,use_default_bond=False, backbone_vector=None):
        """
        Creates `number_of_residues` residues of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the residue type to be created. `name` must be defined in `pmb.df`
            espresso_system(`obj`): Instance of a system object from espressomd library.
            number_of_residue(`int`): Number of residues of type `name` to be created.
            central_bead_position(`list` of `float`): Position of the central bead.
            use_default_bond (`bool`): Switch to control if a bond of type `default` is used to bond particle whose bonds types are not defined in `pmb.df`
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
            self.add_value_to_df(key=('residue_id',''),index=int (residue_index),new_value=residue_id, warning=False)
            # create the principal bead
            if self.df.loc[self.df['name']==name].central_bead.values[0] is self.np.NaN:
                raise ValueError("central_bead must contain a particle name")        
            central_bead_name = self.df.loc[self.df['name']==name].central_bead.values[0]            
            central_bead_id = self.create_particle_in_espresso(name=central_bead_name,
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
                        bead_position=self.generate_trialvectors(center=central_bead_position, 
                                                                 radius=l0, 
                                                                 n_samples=1,
                                                                 on_surface=True)[0]
                    else:
                        bead_position=self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                            center=central_bead_position, 
                                                                            radius=l0)
                        
                    side_bead_id = self.create_particle_in_espresso(name=side_chain_element, 
                                                                    espresso_system=espresso_system,
                                                                    position=[bead_position], 
                                                                    number_of_particles=1)[0]
                    index = self.df[self.df['particle_id']==side_bead_id].index.values[0]
                    self.add_value_to_df(key=('residue_id',''),
                                        index=int (index),
                                        new_value=residue_id, 
                                        warning=False)
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
                                              particle_name2=side_chain_element, 
                                              hard_check=True, 
                                              use_default_bond=use_default_bond)
                    if backbone_vector is None:
                        residue_position=self.generate_trialvectors(center=central_bead_position, 
                                                                    radius=l0, 
                                                                    n_samples=1,
                                                                    on_surface=True)[0]
                    else:
                        residue_position=self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                center=central_bead_position, 
                                                                                radius=l0)
                    lateral_residue_info = self.create_residue_in_espresso(name=side_chain_element, espresso_system=espresso_system,
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
                                        warning=False)
                    # Change the residue_id of the particles in the residue in the side chain
                    side_chain_beads_ids+=[central_bead_side_chain_id]+lateral_beads_side_chain_ids
                    for particle_id in side_chain_beads_ids:
                        index = self.df[(self.df['particle_id']==particle_id) & (self.df['pmb_type']=='particle')].index.values[0]
                        self.add_value_to_df(key=('residue_id',''),
                                            index=int (index),
                                            new_value=residue_id, 
                                            warning=False)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, central_bead_side_chain_id))
                    self.add_bond_in_df(particle_id1=central_bead_id,
                                        particle_id2=central_bead_side_chain_id,
                                        use_default_bond=use_default_bond)
            # Internal bookkeeping of the side chain beads ids
            residues_info[residue_id]['side_chain_ids']=side_chain_beads_ids
        return  residues_info

    def create_variable_with_units(self, variable_dict):
        """
        Returns a pint object with the value and units defined in `variable_dict`.

        Args:
            variable_dict(`dict`): {'value': value, 'units': units}
        Returns:
            variable_with_units(`obj`): variable with units using the pyMBE UnitRegistry.
        """        
        value=variable_dict.pop('value')
        units=variable_dict.pop('units')
        variable_with_units=value*self.units(units)

        return variable_with_units

    def define_AA_particles_in_sequence (self, sequence):
        '''
        Defines in `pmb.df` all the different particles in `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein. 

        Note:
            - It assumes that the names in `sequence` correspond to aminoacid names using the standard  one letter code.
        '''

        already_defined_AA=[]
        acidic_aminoacids = ['c','E','D','Y','C']
        basic_aminoacids  = ['R','n','K','H']

        for residue_name in sequence:
            if residue_name in already_defined_AA:
                continue
            if residue_name in acidic_aminoacids:
                self.define_particle (name=residue_name, acidity='acidic')
            elif residue_name in basic_aminoacids:
                self.define_particle (name=residue_name, acidity='basic')
            else:
                self.define_particle (name=residue_name, q=0)
        return 

    def define_AA_residues (self,sequence, model):
        """
        Defines in `pmb.df` all the different residues in `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per aminoacid are supported.
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
 
    def define_bond(self, bond_object, bond_type, particle_name1, particle_name2):
        """
        Defines a pmb object of type `bond` in `pymbe.df`
        
        Args:
            bond_object (`cls`): instance of a bond object from espressomd library
            bond_type (`str`): label identifying the used bonded potential
            particle_name1 (`str`): `name` of the first `particle` to be bonded.
            particle_name2 (`str`): `name` of the second `particle` to be bonded.
        
        Note:
            - Currently, only harmonic and FENE bonds are supported.
        """    
        
        valid_bond_types=["harmonic", "FENE"]
        if bond_type not in valid_bond_types:
                       raise ValueError(f"Bond type {bond_type} currently not implemented in pyMBE, valid types are {valid_bond_types}")

        lj_parameters=self.get_lj_parameters(particle_name1=particle_name1, 
                                        particle_name2=particle_name2, 
                                        combining_rule='Lorentz-Berthelot')

        cutoff = 2**(1.0/6.0) * lj_parameters["sigma"]
        l0 = self.calculate_initial_bond_length(bond_object=bond_object, 
                                                bond_type=bond_type, 
                                                epsilon=lj_parameters["epsilon"],
                                                sigma=lj_parameters["sigma"],
                                                cutoff=cutoff)

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
        return

    def define_default_bond(self, bond_object, bond_type, epsilon=None, sigma=None, cutoff=None):
        """
        Asigns `bond` in `pmb.df` as the default bond.
        The LJ parameters can be optionally provided to calculate the initial bond length

        Args:
            bond (`obj`): instance of a bond object from the espressomd library.
            bond_type (`str`): label identifying the used bonded potential
            sigma (`float`, optional): LJ sigma of the interaction between the particles
            epsilon (`float`, optional): LJ epsilon for the interaction between the particles
            cutoff (`float`, optional): cutoff-radius of the LJ interaction

        Note:
            - Currently, only harmonic and FENE bonds are supported. 
        """
        valid_bond_types=["harmonic", "FENE"]
        if bond_type not in valid_bond_types:
            raise ValueError(f"Bond type {bond_type} currently not implemented in pyMBE, valid types are {valid_bond_types}")
        if epsilon is None:
            epsilon=1*self.units('reduced_energy')
        if sigma is None:
            sigma=1*self.units('reduced_length')
        if cutoff is None:
            cutoff=2**(1.0/6.0)*self.units('reduced_length')
        l0 = self.calculate_initial_bond_length(bond_object=bond_object, 
                                                bond_type=bond_type, 
                                                epsilon=epsilon,
                                                sigma=sigma,
                                                cutoff=cutoff)

        if self.check_if_name_is_defined_in_df(name='default',pmb_type_to_be_defined='bond'):
            return
        if len(self.df.index) != 0:
            index = max(self.df.index)+1
        else:
            index = 0
        self.df.at [index,'name'] = 'default'
        self.df.at [index,'bond_object'] = bond_object
        self.df.at [index,'l0'] = l0
        self.add_value_to_df(index=index,
                            key=('pmb_type',''),
                            new_value='bond')
        self.add_value_to_df(index=index,
                            key=('parameters_of_the_potential',''),
                            new_value=self.json.dumps(bond_object.get_params()))
        return

    def define_epsilon_value_of_particles (self, eps_dict):
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

    def define_particle(self, name, q=0, diameter=None, acidity='inert', epsilon=None, pka=None):
        """
        Defines a pyMBE object of type `particle` in  `pymbe.df`

        Args:
            name (`str`): Unique label that identifies the `particle`.  
            q (`int`, optional): Charge of the `particle`. Defaults to 0.
            diameter (`obj`, optional): Diameter used to setup Lennard-Jones interactions for the `particle`, should have units of length using the `pmb.units` UnitRegistry. Defaults to None.
            epsilon (`obj`, optional): Epsilon parameter used to setup Lennard-Jones interactions for the `particle`, should have units of energy using the `pmb.units` UnitRegistry.. Defaults to None
            acidity (`str`, optional): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            pka (`float`, optional): If `particle` is an acid or a base, it defines its  pka-value. Defaults to None.
        """ 
        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='particle'):
            index = self.df[self.df['name']==name].index[0]                                   
        else:
            index = len(self.df)
            self.df.at [index, 'name'] = name
            self.df.at [index,'pmb_type'] = 'particle'
            
        if diameter:
            self.add_value_to_df(key=('diameter',''),index=index,new_value=diameter)
        if epsilon:
            self.add_value_to_df(key=('epsilon',''),index=index,new_value=epsilon)        
        
        self.set_particle_acidity ( name=name, acidity = acidity , default_charge=q, pka=pka)
        return 
    
    def define_particles_parameter_from_dict (self, param_dict, param_name):
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
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per aminoacid are supported.
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
    
    def define_protein (self, name,model, topology_dict ):
        """
        Defines a pyMBE object of type `protein` in `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the `protein`.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per aminoacid are supported.
            topology_dict (`dict`): {'initial_pos': coords_list, 'chain_id': id, 'diameter': diameter_value}
        """

        protein_seq_list = []     
    
        valid_keys = ['1beadAA','2beadAA']

        if model not in valid_keys:
            raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
        
        if model == '2beadAA':
            self.define_particle(name='CA')

        for residue in topology_dict.keys():
            residue_name = self.re.split(r'\d+', residue)[0]
            if residue_name != 'CA' and residue_name != 'Ca':
                protein_seq_list.append(residue_name)                  

        protein_sequence = ''.join(protein_seq_list)
        clean_sequence = self.protein_sequence_parser(sequence=protein_sequence)

        self.define_AA_particles_in_sequence (sequence=clean_sequence)
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
            - If `name`  is a object_type=`particle`, only the mathcing particles that are not part of bigger objects (e.g. `residue`, `molecule`) will be destroyed. To destroy particles in such objects, destroy the bigger object instead.
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

    def filter_df (self, pmb_type):
        """
        Filters `pmb.df` and returns a sub-set of it containing only rows with pmb_object_type=`pmb_type` and non-Nan columns.
        
        Args:
            pmb_type(`str`): pmb_object_type to filter in `pmb.df`.
        Returns:
            pmb_type_df(`obj`): filtered `pmb.df`.
        """
        pmb_type_df = self.df.loc[self.df['pmb_type']== pmb_type]
        pmb_type_df = pmb_type_df.dropna( axis=1, thresh=1)
        return pmb_type_df

    def find_bond_key(self,particle_name1, particle_name2, use_default_bond=False):
        """
        Searches for the `name`  of the bond between `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        
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

    def generate_coordinates_outside_sphere (self, center, radius, max_dist, n_samples):

        """
        Generates coordinates outside a sphere centered.

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
            coord = self.generate_trialvectors(center=center, radius=max_dist,n_samples=1)[0]
            if self.np.linalg.norm(coord-self.np.asarray(center))>=radius:
                coord_list.append (coord)
                counter += 1
        return coord_list

    def generate_trial_perpendicular_vector(self, vector, center, radius):
        """
        Generates a random vector perpendicular to `vector`.
        
        Args:
            vector(`lst`): Coordinates of the vector.
            magnitude(`float`): magnitude of the vector perpendicular to `vector`.

        Returns:
            perpendicular_vector(`np.array`): random vector perpendicular to `vector`.
        """
        np_vec=self.np.array(vector)
        if np_vec[1] == 0 and np_vec[2] == 0 and np_vec[0] == 1:
            raise ValueError('zero vector')
        perp_vec = self.np.cross(np_vec, self.generate_trialvectors(center=center, radius=radius, n_samples=1, on_surface=True)[0])
        norm_perp_vec = perp_vec/self.np.linalg.norm(perp_vec)
        return center+norm_perp_vec*radius  
    
    def generate_trialvectors(self,center, radius, n_samples, seed=None, on_surface=False):
        """
        Uniformly samples points from a hypersphere. If on_surface is set to True, the points are
        uniformly sampled from the surface of the hypersphere.
        
        Args:
            center(`lst`): Array with the coordinates of the center of the spheres.
            radius(`float`): Radius of the sphere.
            n_samples(`int`): Number of sample points to generate inside the sphere.
            seed (`int`, optional): Seed for the random number generator
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
        rng = self.np.random.default_rng(seed)
        samples = rng.normal(size=(n_samples, d))
        # make the samples lie on the surface of the unit hypersphere
        normalize_radii = self.np.linalg.norm(samples, axis=1)[:, self.np.newaxis]
        samples /= normalize_radii

        if not on_surface:
            # make the samples lie inside the hypersphere with the correct density
            uniform_points = rng.uniform(size=n_samples)[:, self.np.newaxis]
            new_radii = self.np.power(uniform_points, 1/d)
            samples *= new_radii

        # scale the points to have the correct radius and center
        samples = samples * radius + center
        return samples 

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
                exit()
            else:
                return

    def get_charge_map (self):
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
        sigma1 = self.df.loc[self.df["name"]==particle_name1].diameter.iloc[0]
        eps2 = self.df.loc[self.df["name"]==particle_name2].epsilon.iloc[0]
        sigma2 = self.df.loc[self.df["name"]==particle_name2].diameter.iloc[0]
        
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
        Gets the diameter of each `espresso type` in `pmb.df`. 
        
        Returns:
            radius_map(`dict`): {espresso_type: radius}.

        '''

        df_state_one = self.df[[('diameter',''),('state_one','es_type')]].dropna().drop_duplicates()
        df_state_two = self.df[[('diameter',''),('state_two','es_type')]].dropna().drop_duplicates()

        state_one = self.pd.Series (df_state_one.diameter.values/2.0,index=df_state_one.state_one.es_type.values)
        state_two = self.pd.Series (df_state_two.diameter.values/2.0,index=df_state_two.state_two.es_type.values)
        radius_map  = self.pd.concat([state_one,state_two],axis=0).to_dict()
        
        return radius_map

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

    def load_interaction_parameters (self, filename, verbose=False):
        """
        Loads the interaction parameters stored in filename into `pmb.df`
        
        Args:
            filename(`str`): name of the file to be read
            verbose(`bool`): switches on/off the reading prints
        """
        from espressomd import interactions
        without_units = ['q','es_type','acidity']
        with_units = ['diameter','epsilon']
        with open(filename) as f:
            for line in f:
                if line[0] == '#':
                    continue
                param_dict=self.json.loads(line)
                object_type=param_dict['object_type']
                if object_type == 'particle':
                    not_requiered_attributes={}    
                    for not_requiered_key in without_units+with_units:
                        if not_requiered_key in param_dict.keys():
                            if not_requiered_key in with_units:
                                not_requiered_attributes[not_requiered_key]=self.create_variable_with_units(variable_dict=param_dict.pop(not_requiered_key))
                            elif not_requiered_key in without_units:
                                not_requiered_attributes[not_requiered_key]=param_dict.pop(not_requiered_key)
                        else:
                            if not_requiered_key == 'acidity':
                                not_requiered_attributes[not_requiered_key] = 'inert'
                            else:    
                                not_requiered_attributes[not_requiered_key]=None
                    self.define_particle(name=param_dict.pop('name'),
                                    q=not_requiered_attributes.pop('q'),
                                    diameter=not_requiered_attributes.pop('diameter'),
                                    acidity=not_requiered_attributes.pop('acidity'),
                                    epsilon=not_requiered_attributes.pop('epsilon'))
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
                    name1 = param_dict.pop('name1')
                    name2 = param_dict.pop('name2')
                    bond_type = param_dict.pop('bond_type')
                    if bond_type == 'harmonic':
                        k = self.create_variable_with_units(variable_dict=param_dict.pop('k'))
                        r_0 = self.create_variable_with_units(variable_dict=param_dict.pop('r_0'))
                        bond = interactions.HarmonicBond(k=k.to('reduced_energy / reduced_length**2').magnitude, r_0=r_0.to('reduced_length').magnitude)
                    elif bond_type == 'FENE':
                        k = self.create_variable_with_units(variable_dict=param_dict.pop('k'))
                        r_0 = self.create_variable_with_units(variable_dict=param_dict.pop('r_0'))
                        d_r_max = self.create_variable_with_units(variable_dict=param_dict.pop('d_r_max'))
                        bond = interactions.FeneBond(k=k.to('reduced_energy / reduced_length**2').magnitude, r_0=r_0.to('reduced_length').magnitude, d_r_max=d_r_max.to('reduced_length').magnitude)
                    else:
                        raise ValueError("Current implementation of pyMBE only supports harmonic and FENE bonds")
                    self.define_bond(bond_object=bond, bond_type=bond_type, particle_name1=name1, particle_name2=name2)
                else:
                    raise ValueError(object_type+' is not a known pmb object type')
                if verbose:
                    print('Added: '+line)
        return

    
    
    def load_pka_set(self,filename, verbose=False):
        """
        Loads the pka_set stored in `filename` into `pmb.df`.
        
        Args:
            filename(`str`): name of the file with the pka set to be loaded. Expected format is {name:{"acidity": acidity, "pka_value":pka_value}}.
            verbose(`bool`, optional): If activated, the function reports each pKa value loaded. Defaults to False

        Note:
            - If `name` is already defined in the `pymbe.df`, it prints a warning.
        """
        pKa_list=[]
        with open(filename) as f:
            for line in f:
                if line[0] == '#':
                    continue
                param_dict=self.json.loads(line)
                pKa_list.append(param_dict)    
                if verbose:
                    print('Added: '+line)
        for pka_set in pKa_list:            
            self.check_pka_set(pka_set=pka_set)
            for pka_key in pka_set: 
                acidity =   pka_set[pka_key]['acidity']        
                pka_value = pka_set[pka_key]['pka_value']     
                self.define_particle(name=pka_key,acidity=acidity, pka=pka_value)
        return

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
        Parses `sequence` to the one letter code for aminoacids.
        
        Args:
            sequence(`str` or `lst`): Sequence of the aminoacid. 

        Returns:
            clean_sequence(`list`): `sequence` using the one letter code.
        
        Note:
            - Accepted formats for `sequence` are:
                - `lst` with one letter or three letter code of each aminoacid in each element
                - `str` with the sequence using the one letter code
                - `str` with the squence using the three letter code, each aminoacid must be separated by a hyphon "-"
        
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
            if (sequence.find("-") != -1):
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
                            if (residue.upper() in keys.keys()):
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
        if filename[-3:] != "csv":
            raise ValueError("Only files with CSV format are supported")
        df = self.pd.read_csv (filename,header=[0, 1], index_col=0)
        columns_names = self.setup_df()
        df.columns = columns_names
        self.df=df
        return df
    
    def read_protein_vtf_in_df (self,filename,unit_length=None):
        """
        Loads a coarse-grained protein model in a vtf file `filename` into the `pmb.df` and it labels it with `name`.

        Args:
            filename(`str`): Path to the vtf file with the coarse-grained model.
            unit_length(`obj`): unit of lenght of the the coordinates in `filename` using the pyMBE UnitRegistry. Defaults to None.

        Returns:
            topology_dict(`dict`): {'initial_pos': coords_list, 'chain_id': id, 'diameter': diameter_value}

        Note:
            - If no `unit_length` is provided, it is assumed that the coordinates are in Angstrom.
        """

        print (f'Loading protein coarse grain model file: {filename}')

        coord_list = []
        particles_dict = {}

        if unit_length == None:
            unit_length = 1 * self.units.angstrom 

        with open (filename,'r') as protein_model:

                for line in protein_model :
                    line_split = line.split ()
    
                    if line_split : 
                        line_header = line_split [0]

                        if line_header == 'atom':

                            atom_id  = line_split [1]
                            atom_name = line_split [3]
                            atom_resname = line_split [5]
                            chain_id = line_split [9]

                            particles_dict [int(atom_id)] = [atom_name , atom_resname, chain_id]
         
                        elif line_header.isnumeric (): 

                            atom_coord = line_split [1:] 
                            atom_coord = [(float(i)*unit_length).to('reduced_length').magnitude for i in atom_coord]
                            coord_list.append (atom_coord)

        numbered_label = []
        i = 0   
        
        for atom_id in particles_dict.keys():
    
            if atom_id == 1:
                atom_name = particles_dict[atom_id][0]
                numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2]]
                numbered_label.append(numbered_name)

            elif atom_id != 1: 
            
                if particles_dict[atom_id-1][1] != particles_dict[atom_id][1]:
                    i += 1                    
                    count = 1
                    atom_name = particles_dict[atom_id][0]
                    numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2]]
                    numbered_label.append(numbered_name)
                    
                elif particles_dict[atom_id-1][1] == particles_dict[atom_id][1]:
                    if count == 2 or particles_dict[atom_id][1] == 'GLY':
                        i +=1  
                        count = 0
                    atom_name = particles_dict[atom_id][0]
                    numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2]]
                    numbered_label.append(numbered_name)
                    count +=1

        topology_dict = {}

        for i in range (0, len(numbered_label)):   
            topology_dict [numbered_label[i][0]] = {'initial_pos': coord_list[i] ,'chain_id':numbered_label[i][1] }

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
        if bond_key:
            return self.df[self.df['name']==bond_key].bond_object.values[0]
        else:
            print("Bond not defined between particles ", particle_name1, " and ", particle_name2)    
            if hard_check:
                exit()
            else:
                return

    


    def set_particle_acidity(self, name, acidity='inert', default_charge=0, pka=None):
        """
        Sets the particle acidity if its acidic or basic creates `state_one` and `state_two` with the protonated and 
        deprotonated states. In each state is set: `label`,`charge` and `es_type`. If its inert it will define only `state_one`

        Args:
            name (`str`): Unique label that identifies the `particle`. 
            acidity (`str`): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            default_charge (`int`): Charge of the particle. Defaults to 0.
            pka (`float`, optional):  If `particle` is an acid or a base, it defines its  pka-value. Defaults to None.
        """
        acidity_valid_keys = ['inert','acidic', 'basic']
        if acidity not in acidity_valid_keys:
            raise ValueError(self.df[name].name +' provided acidity not supported, valid keys are ', acidity_valid_keys)   
        for index in self.df[self.df['name']==name].index:       
            if pka:
                self.add_value_to_df(key=('pka',''),index=index,new_value=pka)
            self.add_value_to_df(key=('acidity',''),index=index,new_value=acidity) 
            if not self.check_if_df_cell_has_a_value(index=index,key=('state_one','es_type')):
                self.add_value_to_df(key=('state_one','es_type'),index=index,new_value=self.propose_unused_type())  
            if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'inert':
                self.add_value_to_df(key=('state_one','charge'),index=index,new_value=default_charge)
                self.add_value_to_df(key=('state_one','label'),index=index,new_value=name)
            else:
                protonated_label = f'{name}H'
                self.add_value_to_df(key=('state_one','label'),index=index,new_value=protonated_label)
                self.add_value_to_df(key=('state_two','label'),index=index,new_value=name)
                if not self.check_if_df_cell_has_a_value(index=index,key=('state_two','es_type')):
                    self.add_value_to_df(key=('state_two','es_type'),index=index,new_value=self.propose_unused_type())
                if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'acidic':        
                    self.add_value_to_df(key=('state_one','charge'),index=index,new_value=0)
                    self.add_value_to_df(key=('state_two','charge'),index=index,new_value=-1)
                elif self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'basic':
                    self.add_value_to_df(key=('state_one','charge'),index=index,new_value=+1)
                    self.add_value_to_df(key=('state_two','charge'),index=index,new_value=0)   
        return

    
    
    def set_reduced_units(self, unit_length=None,unit_charge=None,temperature=None):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args:
            unit_length (`obj`,optional): Reduced unit of length defined using the `pmb.units` UnitRegistry. Defaults to None. 
            unit_charge (`obj`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            temperature (`obj`,optional): Temperature of the system, defined using the `pmb.units` UnitRegistry. Defaults to None. 

        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
        """
        self.units=self.pint.UnitRegistry()
        if unit_length is None:
            unit_length=0.355*self.units.nm
        if temperature is None:
            temperature=298.15 * self.units.K
        if unit_charge is None:
            unit_charge=self.units.e
        self.N_A=6.02214076e23 / self.units.mol
        self.Kb=1.38064852e-23 * self.units.J / self.units.K
        self.e=1.60217662e-19 *self.units.C
        self.kT=temperature*self.Kb
        self.units.define(f'reduced_energy = {self.kT} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')        
        self.print_reduced_units()
        return

    def setup_constantpH_reactions_in_espresso(self, counter_ion, constant_pH, SEED, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up the Acid/Base reactions for acidic/basic `particles` defined in `pmb.df` to be sampled in the constant pH ensamble. 

        Args:
            counter_ion(`str`): `name` of the counter_ion `particle`.
            constant_pH(`float`): pH-value.
            SEED(`int`): Seed for the random number generator.
            exclusion_range(`float`, optional): Bellow this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius per each espresso_type. Defaults to `False`.
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
                                                    seed=SEED, 
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

    def setup_df (self):
        """
        Sets up the pyMBE's dataframe `pymbe.df`.

        Returns:
            columns_names(`obj`): pandas multiindex object with the column names of the pyMBE's dataframe
        """

        columns_names = self.pd.MultiIndex.from_tuples ([
                        ('name',''),
                        ('pmb_type',''),
                        ('particle_id',''), 
                        ('particle_id2',''),
                        ('residue_id',''),
                        ('molecule_id',''),
                        ('acidity',''),
                        ('pka',''),
                        ('central_bead',''),
                        ('side_chains',''),
                        ('residue_list',''),
                        ('model',''),
                        ('diameter',''),
                        ('epsilon',''),
                        ('state_one','label'),
                        ('state_one','es_type'),
                        ('state_one','charge'),
                        ('state_two','label'),
                        ('state_two','es_type'),
                        ('state_two','charge'),
                        ('sequence',''),
                        ('bond_object',''),
                        ('parameters_of_the_potential','')
                        ])

        self.df = self.pd.DataFrame (columns = columns_names)
      
        return columns_names


    def setup_lj_interactions_in_espresso (self, espresso_system, cutoff=None, shift='auto', combining_rule='Lorentz-Berthelot'):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle types with values for `diameter` and `epsilon` stored in `pymbe.df`.
        Stores the parameters loaded into  ESPResSo for each type pair in `pymbe.df`.
        Assumes that the LJ interactions between all particles will have the same `cutoff`, `shift` and `combining_rule`.
        Check the documentation of ESPResSo for more info about the potential https://espressomd.github.io/doc4.2.0/inter_non-bonded.html

        Args:
            espresso_system(`obj`): Instance of a system object from the espressomd library.
            cutoff(`float`, optional): cut-off length of the LJ potential. Defaults to None.
            shift (`string`, optional): If set to `auto` shifts the potential to be continous at `cutoff`. Defaults to `auto`.
            combining_rule (`string`, optional): combining rule used to calculate `sigma` and `epsilon` for the potential betwen a pair of particles. Defaults to 'Lorentz-Berthelot'.

        Note:
            If no `cutoff`  is given, its value is set to 2**(1./6.) in reduced_lenght units, corresponding to a purely steric potential.

        Note:
            Currently, the only `combining_rule` supported is Lorentz-Berthelot.
        """
        from itertools import combinations_with_replacement
        sigma=1*self.units('reduced_length')
        implemented_combinatiorial_rules = ['Lorentz-Berthelot']
        compulsory_parameters_in_df = ['diameter','epsilon']
        if combining_rule not in implemented_combinatiorial_rules:
            raise ValueError('In the current version of pyMBE, the only combinatorial rules implemented are ', implemented_combinatiorial_rules)
        if cutoff is None:
            cutoff=2**(1./6.)*self.units('reduced_length')
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
        for type_pair in combinations_with_replacement(particles_types_with_LJ_parameters, 2):
            diameter_list=[]
            epsilon_list=[]
            label_list=[]
            Non_interacting_pair=False
            # Search the LJ parameters of the type pair
            for ptype in type_pair:
                diameter=self.find_value_from_es_type(es_type=ptype, column_name='diameter').to('reduced_length').magnitude
                label=self.find_value_from_es_type(es_type=ptype, column_name='label')
                epsilon=self.find_value_from_es_type(es_type=ptype, column_name='epsilon').to('reduced_energy').magnitude
                if diameter == 0:
                    Non_interacting_pair=True
                    break
                elif diameter < 0:
                    raise ValueError(f"Particle {label} has a negative diameter = {diameter}, check your pmb.df")
                else:   
                    label_list.append(label)
                    diameter_list.append(diameter)
                    epsilon_list.append(epsilon)
            # If the diameter of one of the particle types is 0, no LJ interaction is setup
            if Non_interacting_pair:
                pass
            if combining_rule == 'Lorentz-Berthelot':
                diameter_sum=self.np.sum(diameter_list)
                epsilon=self.np.sqrt(self.np.prod(epsilon_list))
                offset=diameter_sum/2.-sigma.to('reduced_length').magnitude
            espresso_system.non_bonded_inter[type_pair[0],type_pair[1]].lennard_jones.set_params(epsilon = epsilon, 
                                                                                    sigma = sigma.to('reduced_length').magnitude, 
                                                                                    cutoff = cutoff.to('reduced_length').magnitude,
                                                                                    offset = offset, 
                                                                                    shift = shift)                                                                                          
            index = len(self.df)
            self.df.at [index, 'name'] = f'LJ: {label_list[0]}-{label_list[1]}'
            lj_params=espresso_system.non_bonded_inter[type_pair[0], type_pair[1]].lennard_jones.get_params()
            self.add_value_to_df(index=index,
                                key=('pmb_type',''),
                                new_value='LennardJones')
            self.add_value_to_df(index=index,
                                key=('parameters_of_the_potential',''),
                                new_value=self.json.dumps(lj_params))                
        if non_parametrized_labels:
            print(f'WARNING: No LJ interaction has been added in ESPResSo for particles with labels: {non_parametrized_labels}. Please, check your pmb.df to ensure if this is your desired setup.')
        return

    def setup_particle_diameter (self,topology_dict):
        '''
        Sets up the diameter of the particles in `positions`.

        Args:
            topology_dict(`dict`): {'initial_pos': coords_list, 'chain_id': id, 'diameter': diameter_value}
        '''
        for residue in topology_dict.keys():
            residue_name = self.re.split(r'\d+', residue)[0]
            residue_number = self.re.split(r'(\d+)', residue)[1]
            residue_diameter  = topology_dict[residue]['diameter']
            diameter = residue_diameter*self.units.nm
            index = self.df[(self.df['residue_id']==residue_number) & (self.df['name']==residue_name) ].index.values[0]

            self.add_value_to_df(key= ('diameter',''),
                        index=int (index),
                        new_value=diameter)
            
        return 

    def write_output_vtf_file (self, espresso_system, filename):
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
            coordinates.write (f'timestep indexed\n')
            for particle in espresso_system.part:
                coordinates.write (f'{particle.id} \t {particle.pos[0]} \t {particle.pos[1]} \t {particle.pos[2]}\n')
        return 
    
