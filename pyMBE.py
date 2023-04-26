class pymbe_library():
    """
    The library for the Molecular Brewer for ESPResSo (pyMBE)
    """
    import pint
    units = pint.UnitRegistry()
    import numpy as np
    import pandas as pd 
    import json
    import math
    # Default values    
    TEMPERATURE = 298.15 * units.K
    PARTICLE_SIZE = 0.355 * units.nm
    N_A=6.02214076e23    / units.mol
    Kb=1.38064852e-23    * units.J / units.K
    e=1.60217662e-19 *units.C
    pi=3.14159265359
    df=None
    
    def __init__(self):
        """
        Initializes the pymbe_library object with default settings and sets up `pymbe.df` for storing the data.
        """
        # Default definitions of reduced units
        self.units.define(f'reduced_energy = {self.TEMPERATURE * self.Kb}')
        self.units.define(f'reduced_length = {self.PARTICLE_SIZE}')
        self.units.define(f'reduced_charge = 1*e')
        self.kT=self.TEMPERATURE*self.Kb
        self.df = None
        self.setup_df()
        return

    def setup_df (self):
        """
        Sets up `pymbe.df`
        """

        columns_name = self.pd.MultiIndex.from_tuples ([
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

        self.df = self.pd.DataFrame (columns = columns_name)
      
        return 

    def print_reduced_units(self):
        """
        Prints the set of reduced units defined in pyMBE.units.
        """
        print("\nCurrent set of reduced units:")
        unit_length=self.units.Quantity(1,'reduced_length')
        unit_energy=self.units.Quantity(1,'reduced_energy')
        unit_charge=self.units.Quantity(1,'reduced_charge')
        print(unit_length.to('nm'), "=", unit_length)
        print(unit_energy.to('J'), "=", unit_energy)
        print('Temperature:', self.TEMPERATURE)
        print(unit_charge.to('C'), "=", unit_charge)
        print()
        
    def set_reduced_units(self, unit_length=0.355*units.nm,  unit_charge=units.e, temperature=298.15 * units.K):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args: 
            unit_length (cls): self.units object with units of length which will be used as reduced unit of lenght
            unit_charge (cls): self.units object with units of charge which will be used as reduced unit of charge
            temperature (cls): self.units object with units of temperature which will be used to calculate the reduced unit of energy = k_BT
        """
        self.units=self.pint.UnitRegistry()
        self.N_A=6.02214076e23 / self.units.mol
        self.Kb=1.38064852e-23 * self.units.J / self.units.K
        self.e=1.60217662e-19 *self.units.C
        self.TEMPERATURE=temperature.to('K').magnitude*self.units.K
        unit_energy=self.TEMPERATURE*self.Kb
        self.units.define(f'reduced_energy = {unit_energy} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')        
        self.kT=self.TEMPERATURE*self.Kb

        self.print_reduced_units()

    def define_particle(self, name, q=0, diameter=None, acidity='inert', epsilon=None, pka=None):
        """
        Defines a pyMBE object of type `particle` in  `pymbe.df`

        Args:
            name (str): Unique label to identify the `particle` type 
            q (int, optional): Charge of the `particle`. Defaults to 0.
            diameter (cls, optional): Diameter used to setup Lennard-Jones interactions for the `particle`, should be given as a pyMBE.units object with units of length. Defaults to None
            epsilon (cls, optional): Epsilon parameter used to setup Lennard-Jones interactions for the `particle`, should be given as a pyMBE.units object with units of energy. Defaults to None
            acidity (str, optional): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            pka (float, optional): If `particle` is an acid or a base, it its  pka-value.
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

    def set_particle_acidity(self, name, acidity='inert', default_charge=0, pka=None):
        """
        Sets the particle acidity if its acidic or basic creates `state_one` and `state_two` with the protonated and 
        deprotonated states. In each state is set: `label`,`charge` and `es_type`. If its inert it will define only `state_one`

        Args:
            name (str): Unique label to identify the `particle` type 
            acidity (str): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to `inert`.
            default_charge (int): Charge of the particle if its not set up the default is 0
            pka (float, optional): If `particle` is an acid or a base, it its  pka-value.
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

    def define_residue(self, name, central_bead, side_chains):
        """
        Defines a pyMBE object of type `residue` in the `pymbe.df`

        Args:
            name (str): Unique label to identify the `residue` type 
            central_bead (str): Label of `name` of the `particle` pmb_type to be placed as central_bead of the residue
            side_chains (list of str): List of `name`s corresponding to the `pmb_type`s to be placed as side_chains of the residue. Currently, only `pmb_types` of `particle` or `residue` are supported.
        """
        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='residue'):
            return
        index = len(self.df)
        self.df.at [index, 'name'] = name
        self.df.at [index,'pmb_type'] = 'residue'
        self.df.at [index,'central_bead'] = central_bead
        self.df.at [index,('side_chains','')] = side_chains
        return 

    def define_molecule(self, name, residue_list):
        """
        Defines a pyMBE object of type `molecule` in the `pymbe.df`

        Args:
            name (str): Unique label to identify the `molecule` type
            residue_list (list of str): List of `name`s corresponding to `residue`s  of the `molecule`.  
        """
        if self.check_if_name_is_defined_in_df(name=name,pmb_type_to_be_defined='molecule'):
            return
        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'molecule'
        self.df.at [index,('residue_list','')] = residue_list
        return 

    def define_peptide(self, name, sequence, model):
        """
        Returns a pyMBE object of type `molecule` in the `pymbe.df`

        Args:
            name (str): Unique label to identify the `peptide` type
            sequence (string): string with the peptide sequence
            model (string): string with the model name: '1beadAA' or '2beadAA'
        """
        if not self.check_if_name_is_defined_in_df(name = name, pmb_type_to_be_defined='peptide'):
            valid_keys = ['1beadAA','2beadAA']
            if model not in valid_keys:
                raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
            clean_sequence = self.protein_sequence_parser(sequence=sequence)
            residue_list = []
            for residue_name in clean_sequence:
                if model == '1beadAA':
                    central_bead = residue_name
                    side_chains = []
                elif model == '2beadAA':
                    # terminal groups + glycine are represented with only 1 bead
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
            self.define_molecule(name = name, residue_list=residue_list)
            index = self.df.loc[self.df['name'] == name].index.item() 
            self.df.at [index,'model'] = model
            self.df.at [index,('sequence','')] = clean_sequence
        return 

    def define_bond(self, bond_object, particle_name1, particle_name2):
        """
        Defines a pmb object of type `bond` in `pymbe.df`
        
        Args:
            bond_object (cls): instance of a bond object from espressomd library
            particle_name1 (str): label of the type of the first particle type of the bonded particles.
            particle_name2 (str): label of the type of the second particle type of the bonded particles.
        """    
        index = len(self.df)
        for label in [particle_name1+'-'+particle_name2,particle_name2+'-'+particle_name1]:
            self.check_if_name_is_defined_in_df(name=label, pmb_type_to_be_defined="bond")
        self.df.at [index,'name']= particle_name1+'-'+particle_name2
        self.df.at [index,'bond_object'] = bond_object
        self.add_value_to_df(index=index,
                                key=('pmb_type',''),
                                new_value='bond')
        self.add_value_to_df(index=index,
                                key=('parameters_of_the_potential',''),
                                new_value=self.json.dumps(bond_object.get_params()))
        return

    def add_bond_in_df (self,particle_id1,particle_id2, use_default_bond=False):
        """
        Creates a bond entry on the `pymbe.df` storing the particle_ids of the two particles in the bond

        Args: 
            particle_id1 (int): particle_id of the type of the first particle type of the bonded particles
            particle_id2 (int): particle_id of the type of the second particle type of the bonded particles
            use_default_bond (bool, optional): Switch to control if a bond of type `default` is used to bond particle whose bonds types are not defined in `df`. Defaults to False.
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

    def create_particle_in_espresso(self, name, espresso_system, number_of_particles, position=None, fix=False):
        """
        Creates `number_of_particles` particles of type `name` into `espresso_system`.
        Bookeeps the particle ids of the particle created into `pymbe.df`.
        Returns a list with the ids of the particles created.
        `name` must be defined into `pymbe.df` as a pmb_object of type particle.

        Args:
            name (str): Label of the particle type to be created. The particle type must be defined in `df
            espresso_system (cls): Instance of a system class from espressomd library.
            number_of_particles (int): Number of particles of type `name` to be created.
            position (array of array, optional): Initial positions of the particles. If not given, particles are created in random positions
        Returns:
            created_pid_list (:obj:`list` of :obj:`int`): List with the ids of the particles created into `espresso_system`.
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

    def create_residue_in_espresso (self, name, espresso_system, number_of_residues, central_bead_position=None,use_default_bond=False, backbone_vector=None):
        """
        Creates `number_of_residues` residues of type `name` into `espresso_system`.
        Bookeeps the residue ids of the residues created into `pymbe.df` in the correponding residues and particle entries.
        Bookeeps the particles ids of the particles created into `pymbe.df` in the correponding particle entries.
        Returns a dict with the ids of the particles created {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
        `name` must be defined into `pymbe.df` as a pmb_object of type residue.

        Args:
            name (str): Label of the residue type to be created. The residue type must be defined in `df`
            espresso_system (cls): Instance of a system class from espressomd library.
            number_of_residue (int): Number of residues of type `name` to be created.
            central_bead_position (:obj:`list` of :obj:`float`): List with the coordinates of the position of the central bead
            use_default_bond (bool): Switch to control if a bond of type `default` is used to bond particle whose bonds types are not defined in `df`
            backbone_vector (:obj:`list` of :obj:`float`): List with the coordinates of the backbone vector of the molecule. All side chains will be created in perpendicular positions to `backbone_vector`
        Returns:
            residues_info (dict of int: dict): dict with the ids of the particles created {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
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
                            hard_check=True, use_default_bond=use_default_bond)
                    if backbone_vector is None:
                        bead_position=self.generate_trialvectors(center=central_bead_position, 
                                                              radius=bond.params.get('r_0'), 
                                                            n_samples=1)[0]
                    else:
                        bead_position=self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                            center=central_bead_position, 
                                                                            radius=bond.params.get('r_0'))
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
                    if backbone_vector is None:
                        residue_position=self.generate_trialvectors(center=central_bead_position, 
                                                                radius=bond.params.get('r_0'), 
                                                                n_samples=1)[0]
                    else:
                        residue_position=self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                center=central_bead_position, 
                                                                                radius=bond.params.get('r_0'))
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

    def create_molecule_in_espresso(self, name, number_of_molecules,espresso_system, first_residue_position=None, use_default_bond=False):
        """
        Creates a molecule in the espresso system

        Creates `number_of_molecules` molecule of type `name` into `espresso_system`.
        Bookeeps the molecule ids of the residues created into `pymbe.df` in the correponding residues and particle entries.
        Bookeeps the particles ids of the particles created into `pymbe.df` in the correponding particle entries.
        Returns a dict with the ids of the particles created {molecule_id: {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":
        [particle_id1, ...]}}} 
        `name` must be defined into `pymbe.df` as a pmb_object of type molecule.

        Args:
            name (str): Label of the molecule type to be created. The molecule type must be defined in `pymbe.df`
            espresso_system (cls): Instance of a system class from espressomd library.
            number_of_molecules (int): Number of molecules of type `name` to be created.
            first_residue_position (list, optional): coordinates where the first_residue_position will be created, random by default
            use_default_bond (bool, optional): Switch to control if a bond of type `default` is used to bond particle whose bonds 
            types are not defined in `pymbe.df`
        Returns:
            molecules_info (dict of int: dict): dict with the ids of the particles created {molecule_id: {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":
        [particle_id1, ...]}}} 
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
            self.clean_df_row(index=int(molecule_index))
            if self.df['molecule_id'].isnull().values.all():
                molecule_id = 0        
            else:
                # check if a residue is part of another molecule
                check_residue_name = self.df[self.df['residue_list'].astype(str).str.contains(name)]
                pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]                
                if not check_residue_name.empty and pmb_type == 'molecule' :              
                    for value in check_residue_name.index.to_list():                  
                        if value not in used_molecules_id:                              
                            molecule_id = self.df.loc[value].molecule_id.values[0]                    
                            break
                else:
                    molecule_id = self.df['molecule_id'].max() +1
            #assigns molecule_id to the residue defined       
            self.add_value_to_df (key=('molecule_id',''),
                                index=int (molecule_index),
                                new_value=molecule_id, 
                                warning=False)
            for residue in residue_list:
                if first_residue:
                    residue_position = first_residue_position
                    backbone_vector = self.generate_trialvectors(center=[0,0,0], 
                                                                radius=1, 
                                                                n_samples=1)[0]
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
                    residue_position = residue_position+backbone_vector*bond.params.get('r_0')  
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
    
    def check_if_df_cell_has_a_value(self,index,key):
        """
        Checks if a cell in the `pymbe.df` DataFrame at the specified index and column has a value.

        Args:
            index (int): the index of the row to check.
            key (str): the column label to check.

        Returns:
            (bool): True if the cell has a value, False otherwise.
        """
        idx = self.pd.IndexSlice
        return not self.pd.isna(self.df.loc[index, idx[key]])

    def clean_df_row(self, index, columns_keys_to_clean=["particle_id", "particle_id2", "residue_id", "molecule_id"]):
        """
        Cleans the columns of `pymbe.df` in `columns_keys_to_clean` of the row with index `index` by asigning them a np.nan value.

        Args:
            index (int): the index of the row to clean.
            columns_keys_to_clean (list of str, optional): list with the column keys to be cleaned. 
            Defaults to `particle_id`, `particle_id2`, `residue_id`, `molecule_id`.
        """   
        for column_key in columns_keys_to_clean:
            self.add_value_to_df(key=(column_key,''),index=index,new_value=self.np.nan, warning=False)
        return
    def add_value_to_df(self,index,key,new_value, warning=True):
        """
        Adds a value to a cell in the `pymbe.df` DataFrame.

        Args:
            index (int): the index of the row to add the value to.
            key (str): the column label to add the value to.
            warning (bool, optional): If true, prints a warning if a value is being overwritten in `pymbe.df`. Defaults to true.
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

    def check_if_name_is_defined_in_df (self, name, pmb_type_to_be_defined):
        """
        Checks if `name` is defined in `pymbe.df`.
        If it is defined and it coresponds to a pmb_type different than `pmb_type_to_be_defined`, raises a ValueError.

        Args:
            name (str): label to check if defined in `pymbe.df`.
            pmb_type_to_be_defined (str): pmb object type corresponding to `name`.

        Returns:
            bool: True for success, False otherwise.
        """

        if name in self.df['name'].unique():
            current_object_type = self.df[self.df['name']==name].pmb_type.values[0]
            if  current_object_type != pmb_type_to_be_defined:
                raise ValueError ((f"The name {name} is already defined in the df with a pmb_type = {current_object_type}, pymMBE does not support objects with the same name but different pmb_types"))    
            return True            
        else:
            return False
  
    def copy_df_entry(self, name, column_name, number_of_copies):
        '''
        Creates 'number_of_copies' of a given 'name' into `pymbe.df`.
        Returns the updated `pymbe.df`

        Args:
            name (str): Label of the particle/residue/molecule type to be created. 
                        must be defined in `pymbe.df`
            column_name (str): Column name to use as a filter should be: particle_id/residue_id/molecule_id 
            number_of_copies (int): number of copies of `name` to be created.
        '''
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

    def propose_unused_type(self):
        """
        Searches in `pymbe.df` all particle types defined and returns a new type

        Returns:
            unused_type: (int) unused particle type
        """
        type_map=self.get_type_map()
        if type_map == {}:    
            unused_type = 0
        else:
            unused_type=max(type_map.values())+1    
        return unused_type
        
    def search_bond (self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        If `use_default_bond` is activated and a "default" bond is defined, returns that default bond instead.
        If no bond is found, it prints a message and it does not return nothing. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1 (str): label of the type of the first particle type of the bonded particles.
            particle_name2 (str): label of the type of the second particle type of the bonded particles.
            hard_check (bool, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond (bool, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            bond (cls): bond object from the espressomd library. 
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
            
    def find_bond_key(self,particle_name1, particle_name2, use_default_bond=False):
        """
        Searches for the `name`  of the bond between `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        If `use_default_bond` is activated and a "default" bond is defined, returns the default bond if no key is found.
        Otherwise it does not return anything.

        Args:
            particle_name1 (str): label of the type of the first particle type of the bonded particles.
            particle_name2 (str): label of the type of the second particle type of the bonded particles.
            use_default_bond (bool, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            bond_key (str): `name` of the bond between `particle_name1` and `particle_name2` 
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
        
    def add_bonds_to_espresso (self, espresso_system) :
        """
        Adds all bonds defined in `pymbe.df` to `espresso_system`.
        If no bonds are defined in `pymbe.df` prints a warning message.

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
    
    def load_interaction_parameters (self, filename, verbose=False):
        """
        Loads the interaction parameters stored in filename into `pymbe.df`
        
        Args:
            filename (string): name of the file to be read
            verbose (boolean): switches on/off the reading prints
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
                    else:
                        raise ValueError("current implementation of pyMBE only supports harmonic bonds")
                    self.define_bond(bond_object=bond, particle_name1=name1, particle_name2=name2)
                else:
                    raise ValueError(object_type+' is not a known pmb object type')
                if verbose:
                    print('Added: '+line)
        return
    
    def load_pka_set(self,filename, verbose=False):
        """
        Loads the parameters stored in `filename` into `pymbe.df`
        Checks if `name` is defined in the `pymbe.df`, if it is prints a warning

        Args:
            filename (str): name of the file with the pka set to be loaded. Expected format is {name:{"acidity": acidity, "pka_value":pka_value}}.
            verbose (bool, optional): If activated, the function reports each pKa value loaded. Defaults to False
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

    def define_default_bond(self, bond_object):
        """
        Asigns `bond` in `pymbe.df`  as the default bond type to be used by pyMBE.

        Args:
            bond (cls): bond object from espressomd library
        """
        if self.check_if_name_is_defined_in_df(name='default',pmb_type_to_be_defined='bond'):
            return
        if len(self.df.index) != 0:
            index = max(self.df.index)+1
        else:
            index = 0
        self.df.at [index,'name'] = 'default'
        self.df.at [index,'bond_object'] = bond_object
        self.add_value_to_df(index=index,
                            key=('pmb_type',''),
                            new_value='bond')
        self.add_value_to_df(index=index,
                            key=('parameters_of_the_potential',''),
                            new_value=self.json.dumps(bond_object.get_params()))
        return

    def create_variable_with_units(self, variable_dict):
        """
        Returns a pint object with the value and units defined in variable_dict

        Args:
            variable_dict(dict): contains the value and the units of the variable
        Returns:
            variable_with_units(pint object): pint object containing the value and the desired units
        """        
        value=variable_dict.pop('value')
        units=variable_dict.pop('units')
        variable_with_units=value*self.units(units)

        return variable_with_units

    def filter_df (self, pmb_type):
        """
        Filters `pymbe.df` and returns a sub-set of it containing only 
        rows with pmb object type=`pmb_type` and non-Nan columns
        
        Args:
            pmb_type (str): pmb object type to filter in `pymbe.df`
        Returns:
            pmb_type_df (cls): pandas df filtered for rows with pmb object type=`pmb_type` and non-Nan columns
        """
        pmb_type_df = self.df.loc[self.df['pmb_type']== pmb_type]
        pmb_type_df = pmb_type_df.dropna( axis=1, thresh=1)
        return pmb_type_df


    def get_type_map(self):
        """
        Searches all different espresso types assigned to particles  in `pymbe.df` and returns them in a dictionary {"name": espresso_type}
        
        Returns:
            type_map(dict): Dictionary with the type map {"name": espresso_type}.
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

    def get_charge_map (self):
        '''
        Maps the charge associated to each `espresso_type` in `pymbe.df` and 
        returns them in a dictionary {espresso_type: charge}
        
        Returns:
            charge_map(dict): Dictionary with the charge map {espresso_type: charge}.
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
    
    def get_radius_map(self):
        '''
        Searches the diameter associate to each `espresso type` in `pymbe.df` and 
        returns the radius in a dictionary {espresso_type: radius}
        
        Returns:
            radius_map(dict): Dictionary with the radius map {espresso_type: radius}.

        '''

        df_state_one = self.df[[('diameter',''),('state_one','es_type')]].dropna().drop_duplicates()
        df_state_two = self.df[[('diameter',''),('state_two','es_type')]].dropna().drop_duplicates()

        state_one = self.pd.Series (df_state_one.diameter.values/2.0,index=df_state_one.state_one.es_type.values)
        state_two = self.pd.Series (df_state_two.diameter.values/2.0,index=df_state_two.state_two.es_type.values)
        radius_map  = self.pd.concat([state_one,state_two],axis=0).to_dict()
        
        return radius_map
    
    def get_pka_set(self):
        '''
        Searches particles with `acidity` equal to "acidic" or "basic"  in `pymbe.df` and 
        returns a dict of dict with the associated pka-value and acidity
        
        Returns:
            pka_set (dict): Dictionary with the pka_set {"name" : {"pka_value": pka, "acidity": acidity}}
        '''
        titratables_AA_df = self.df[[('name',''),('pka',''),('acidity','')]].drop_duplicates().dropna()
        pka_set = {}
        for index in titratables_AA_df.name.keys():
            name = titratables_AA_df.name[index]
            pka_value = titratables_AA_df.pka[index]
            acidity = titratables_AA_df.acidity[index]   
            pka_set[name] = {'pka_value':pka_value,'acidity':acidity}
        return pka_set

    def protein_sequence_parser(self, sequence):
        '''
        Parses `sequence` to complay to the one letter code used for aminoacids.
        Accepted formats are:
        1) List with one letter or three letter code of each aminoacid in each element
        2) String with the sequence using the one letter code
        3) String with the squence using the three letter code, each aminoacid must be separated by a hyphon "-"
        Prints an error if finds a key unknown to the parser.
        
        Args:
            sequence (string or list): sequence of the aminoacid. 

        Output:
            clean_sequence (list): list with each aminoacid in `sequence` in the one letter code
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

    def setup_constantpH_reactions_in_espresso(self, counter_ion, constant_pH, SEED, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Set up the Acid/Base reactions for acidic/basidic residues in mol. The reaction steps are done following the constant pH ensamble procedure. 

        Args:      
            counter_ion (str): name of counter_ion particle defined in pyMBE library
            constant_pH (float): pH value to set up the reaction
            SEED (int): Seed for the random number generator.
            exclusion_radius (float, optional): exclusion radius for the constant pH ensamble
            use_exclusion_radius_per_type (bool,optional): switch to use an specific exclusion radius per espresso type, by default is set to False.
            pka_set (dict,optional): dictionary with the desired pka_set, by default the one stored in pyMBE will be used
        Returns:
            RE (cls): instance of the espresso class reaction_ensemble.ConstantpHEnsemble
            sucessfull_reactions_labels (list): list with the labels of the reactions that has been set up by pyMBE 
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

    def check_pka_set(self, pka_set):
        """"
        Checks that `pka_set` has the formatting expected by the pyMBE library.
        Raises a ValueError if `pka_set` does not comply with the expected format. 
        Args:
            pka_set (dict): dictionary to be checked. Should be structured following the format:
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_entry in pka_set.values():
                if required_key not in pka_entry.keys():
                    raise ValueError('missing a requiered key ', required_keys, 'in the following entry of pka_set', pka_entry)
        return

    def generate_trialvectors (self,center, radius, n_samples, seed=None):
        """
        Uniformly samples points from a hypersphere.
        Algorithm from
        https://baezortega.github.io/2018/10/14/hypersphere-sampling/

        Args:
            center (array): array with the coordinates of the center of the spheres.
            radius (float): magnitude of the radius of the sphere
            n_samples (int): number of sample points to generate inside the sphere
            seed (int, opt): seed for the random number generator

        Returns:
            samples (list): list with the generated points
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
        # make the samples lie inside the hypersphere with the correct density
        uniform_points = rng.uniform(size=n_samples)[:, self.np.newaxis]
        new_radii = self.np.power(uniform_points, 1/d)
        samples *= new_radii
        # scale the points to have the correct radius and center
        samples = samples * radius + center
        return samples
     
    def generate_trial_perpendicular_vector(self, vector, center, radius):
        """
        Generates a random vector perpendicular to `vector`.
        
        Args:
            vector (list): array with the coordinates of the vector.
            magnitude (float): magnitude desired of the vector perpendicular to `vector`.

        Returns:
            perpendicular_vector (numpy array):  random vector perpendicular to the input vector.
        """
        np_vec=self.np.array(vector)
        if np_vec[1] == 0 and np_vec[2] == 0 and np_vec[0] == 1:
            raise ValueError('zero vector')
        perp_vec = self.np.cross(np_vec, self.generate_trialvectors(center=center, radius=radius, n_samples=1)[0])
        norm_perp_vec = perp_vec/self.np.linalg.norm(perp_vec)
        return center+norm_perp_vec*radius

    def calculate_HH(self, sequence, pH=None, pka_set=None):
        """
        Calculates the ideal Henderson-Hasselbach titration curve in the given pH range

        Args:
            pH (list of floats): list of pH values
            sequence (list of strings or string): list of aminoacids 

        Returns:
            Z_HH (list of floats): Henderson-Hasselbach prediction of the protein charge in the given pH
        """
        if pH is None:
            pH=self.np.linspace(2,12,50)    

        if pka_set is None:
            pka_set=self.get_pka_set() 
        self.check_pka_set(pka_set=pka_set)

        Z_HH=[]

        for pH_value in pH:    
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

    def create_added_salt_in_espresso(self, espresso_system, cation_name, anion_name, c_salt):    
        """
        Adds a c_salt concentration of cations and anions to the espresso_system

        Args:
            espresso_system (class): espresso class object with all system variables.
            cation_name (str): particle defined in `pymbe.df` with `name`: cation_name that has a positive charge
            anion_name (str): particle defined in `pymbe.df` with `name`: anion_name that has a negative charge
            c_salt (float): Added salt concentration
            
        Returns:
            c_salt_calculated (float): Calculated added salt concentration added to solution    
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

    def create_pmb_object_in_espresso(self, name, number_of_objects, espresso_system, position=None, use_default_bond=False):
        """
        Creates all particles contained in the pmb object in the system of espresso
        
        Args:
            name (str): Label of the pmb object to be created. The pmb object must be defined in `pymbe.df`
            number_of_objects (int): Number of pmb object of type `name` to be created.
            espresso_system (cls): Instance of a system class from espressomd library.
            use_default_bond (bool): Switch to control if a bond of type `default` is used to bond particle whose bonds types are not defined in `df`
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
    
    def destroy_pmb_object_in_system(self, name, espresso_system):
        """
        Destroys the object labeled `name` and all particles associated with it from the espresso system.
        Removes the destroyed objects from `df`.
        NOTE: if `name` corresponds to a object_type `particle`, only the particles that are not part of bigger objects
        (i.e. `residue`, `molecule`) will be destroyed. To destroy particles in such objects, destroy the bigger object instead.

        Args:
            name (str): Label of the pmb object to be destroyed. The pmb object must be defined in `pymbe.df`.
            espresso_system (cls): Instance of a system class from espressomd library.
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

    def create_counterions_in_espresso(self, pmb_object, cation_name, anion_name, espresso_system):
        """
        Creates cation and anion particles in the espresso_system to counter the charge of the particles contained in object
        
        Args:
            pmb_object (str): particle, residue or molecule/peptide object
            espresso_system (class): espresso class object with all system variables.
            cation_name (str): particle defined in `pymbe.df` with `name`: cation_name that has a positive charge
            anion_name (str): particle defined in `pymbe.df` with `name`: anion_name that has a negative charge

        Returns: 
            counterion_number (dict): dictionary with the number of counterions created {"name": number}
        """
        cation_charge = self.df.loc[self.df['name']==cation_name].state_one.charge.iloc[0]
        anion_charge = self.df.loc[self.df['name']==anion_name].state_one.charge.iloc[0]
        object_ids = self.df.loc [self.df['pmb_type']== pmb_object].particle_id.dropna().tolist()
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
    
    def find_value_from_es_type(self, es_type, column_name):
        """
        Finds a value in `pymbe.df`  for a `column_name` and `es_type` pair.

        Args:
            es_type (int): value of the espresso type
            column_name (str): name of the column in `pymbe.df`

        Returns:
            column_name_value (any): value in `pymbe.df` matching  `column_name` and `es_type`
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
    
    def setup_lj_interactions_in_espresso (self, espresso_system, cutoff=None, shift='auto', combining_rule='Lorentz-Berthelot'):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle types with values for `diameter` and `epsilon` stored in `pymbe.df`.
        Stores the parameters loaded into  ESPResSo for each type pair in `pymbe.df`.
        Assumes that the LJ interactions between all particles will have the same `cutoff`, `shift` and `combining_rule`.
        Check the documentation of ESPResSo for more info about the potential https://espressomd.github.io/doc4.2.0/inter_non-bonded.html

        Args:
            espresso_system (class): espresso class object with all system variables.
            cutoff (float, optional): cut-off length of the potential. Defaults to 2**(1./6.) in reduced_lenght units, corresponding to a purely steric potential
            shift (string, optional): if set to `auto` shifts the potential to be continous at `cutoff`. Defaults to `auto`.
            combining_rule (string, optional): combining rule used to calculate `sigma` and `epsilon` for the potential betwen a pair of particles. Defaults to 'Lorentz-Berthelot'.
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
        particles_types_without_LJ_parameters= []
        for particle_type in self.get_type_map().values():
            check_list=[]
            for key in compulsory_parameters_in_df:
                value_in_df=self.find_value_from_es_type(es_type=particle_type,
                                                        column_name=key)
                check_list.append(self.np.isnan(value_in_df))
            if any(check_list):
                particles_types_without_LJ_parameters.append(particle_type)
            else:
                particles_types_with_LJ_parameters.append(particle_type)
        for type_pair in combinations_with_replacement(particles_types_with_LJ_parameters, 2):
            pair_parameters={}
            for ptype in type_pair:
                pair_parameters[ptype]={}
                pair_parameters[ptype]["label"]=self.find_value_from_es_type(es_type=ptype, column_name='label')
                pair_parameters[ptype]["diameter"]=self.find_value_from_es_type(es_type=ptype, column_name='diameter')
                pair_parameters[ptype]["epsilon"]=self.find_value_from_es_type(es_type=ptype, column_name='epsilon')
            if combining_rule == 'Lorentz-Berthelot':
                combined_sigma=0
                epsilon=1
                for ptype in type_pair:
                    combined_sigma+=pair_parameters[ptype]["diameter"]
                    epsilon*=pair_parameters[ptype]["epsilon"]
                offset=combined_sigma/2.-sigma
                epsilon=self.np.sqrt(epsilon)
            else:
                raise ValueError('Undefined combination rule, please report this bug to the developers')
            espresso_system.non_bonded_inter[type_pair[0],type_pair[1]].lennard_jones.set_params(epsilon = epsilon.to('reduced_energy').magnitude, 
                                                                                    sigma = sigma.to('reduced_length').magnitude, 
                                                                                    cutoff = cutoff.to('reduced_length').magnitude,
                                                                                    offset = offset.to('reduced_length').magnitude, 
                                                                                    shift = shift)                                                                                          
            index = len(self.df)
            self.df.at [index, 'name'] = f'LJ: {pair_parameters[type_pair[0]]["label"]}-{pair_parameters[type_pair[1]]["label"]}'
            lj_params=espresso_system.non_bonded_inter[type_pair[0], type_pair[1]].lennard_jones.get_params()
            self.add_value_to_df(index=index,
                                key=('pmb_type',''),
                                new_value='LennardJones')
            self.add_value_to_df(index=index,
                                key=('parameters_of_the_potential',''),
                                new_value=self.json.dumps(lj_params))                
        if particles_types_without_LJ_parameters:
            non_parametrized_labels=[]
            for ptype in particles_types_without_LJ_parameters:
                non_parametrized_labels.append(self.find_value_from_es_type(es_type=ptype, 
                                                                            column_name='label')) 
            print('WARNING: missing LJ parameters for the following particle labels:')
            print(non_parametrized_labels)
            print('no interaction has been added for those particles in ESPResSo')
        return

    def load_protein_vtf_in_df (self, name, filename,unit_length=None):
        """
        Reads the input VTF file of the protein model

        Args:
            filename: path of the protein file
        """

        print (f'Loading protein coarse grain model file: {filename}')

        protein_seq_list = []
        coord_list = []

        label_dict = {}
        label_CA = False 

        if unit_length == None:
            unit_length = 1 * self.units.angstrom 

        with open (filename,'r') as protein_model:

                for line in protein_model :
                    line_split = line.split ()
    
                    if line_split : 
                        line_header = line_split [0]

                        if line_header == 'unitcell':    
                            box = line_split[1:4]

                        elif line_header == 'atom':

                            atom_id  = line_split [1]
                            atom_name = line_split [3]
                            atom_resname = line_split [5]
                            chain_id = line_split [7]
                            atom_radius = line_split [9]

                            diameter = ((float(atom_radius)*2.0)*unit_length).to('reduced_length').magnitude 

                            label_dict [int(atom_id)] = [atom_name , atom_resname, chain_id, diameter]

                            if atom_name == 'CA' or label_CA == True:
                                label_CA = True

                            if atom_name != 'CA' and atom_name != 'Ca':
                                protein_seq_list.append(atom_name)                        

                        elif line_header == 'bond' :
                            atom_bond = line_split [1] 

                        elif line_header.isnumeric (): 
                            particle_id = line_split [0]
                            particle_coord = line_split [1:] 
                            particle_coord = [float(i) for i in particle_coord]

                            coord_list.append (particle_coord)

        axes_list = [0,1,2]
        updated_coordinates_list = []



        for pos in coord_list:
            updated_pos = self.np.zeros(3,)
            for axis in axes_list:
                updated_pos[axis] = (pos[axis]*unit_length).to('reduced_length').magnitude 
            updated_coordinates_list.append (updated_pos.tolist())

        protein_sequence = ''.join(protein_seq_list)

        print (f'Protein Sequence: {protein_sequence}')

        numbered_resname = []
        i = 0   
        
        for atom_id in label_dict.keys():
    
            if atom_id == 1:
                atom_name = label_dict[atom_id][0]
                updated_name = [f'{atom_name}{i}',label_dict[atom_id][2],label_dict[atom_id][3]]
                numbered_resname.append(updated_name)

            elif atom_id != 1: 
            
                if label_dict[atom_id-1][1] != label_dict[atom_id][1]:
                    i += 1                    
                    count = 1
                    atom_name = label_dict[atom_id][0]
                    updated_name = [f'{atom_name}{i}',label_dict[atom_id][2],label_dict[atom_id][3]]
                    numbered_resname.append(updated_name)
                    
                elif label_dict[atom_id-1][1] == label_dict[atom_id][1]:
                    if count == 2 or label_dict[atom_id][1] == 'GLY':
                        i +=1  
                        count = 0
                    atom_name = label_dict[atom_id][0]
                    updated_name = [f'{atom_name}{i}',label_dict[atom_id][2],label_dict[atom_id][3]]
                    numbered_resname.append(updated_name)
                    count +=1

        protein_coordinates = {}

        for i in range (0, len(numbered_resname)):   
            protein_coordinates [numbered_resname[i][0]] = {'initial_pos': updated_coordinates_list[i] ,'chain_id':numbered_resname[i][1], 'diameter':numbered_resname[i][2] }

        clean_sequence = self.protein_sequence_parser(sequence=protein_sequence)

        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'protein'
        self.df.at [index,'model'] = '2beadAA' 

        if label_CA:
            self.define_particle(name='CA')

        self.define_AA_particles_in_sequence (clean_sequence=clean_sequence)

        residue_list = []

        for residue_name in clean_sequence:

            if residue_name not in residue_list:   

                if residue_name in ['c','n', 'G']: 
                    central_bead = residue_name
                    side_chains = []
                else:
                    central_bead = 'CA'              
                    side_chains = [residue_name]

                self.define_residue(name = 'AA-'+residue_name, 
                                    central_bead = central_bead,
                                    side_chains = side_chains)
            residue_list.append('AA-'+residue_name)

        self.df.at [index,('sequence','')] = protein_sequence  
        self.df.at [index,('residue_list','')] = residue_list    

        return protein_coordinates
    
    def define_AA_particles_in_sequence (self, clean_sequence):

        '''
        Defines in the df the aminoacids presente in a 'clean_sequence' from a peptide/molecule/protein

        Args:
            clean_sequence: 

        '''

        already_defined_AA=[]

        acidic_aminoacids = ['c','E','D','Y','C']
        basic_aminoacids  = ['R','n','K','H']

        for residue_name in clean_sequence:
            if residue_name in already_defined_AA:
                continue
            if residue_name in acidic_aminoacids:
                self.define_particle (name=residue_name, acidity='acidic')
            elif residue_name in basic_aminoacids:
                self.define_particle (name=residue_name, acidity='basic')
            else:
                self.define_particle (name=residue_name, q=0)

        return 

    def create_protein_in_espresso(self, name, number_of_proteins, espresso_system, positions):

        """
        Creates a protein in the espresso: system
        
        Creates `number_of_proteins` molecule of type `name` into `espresso_system` in the coordinates obtained from the coarse grained model loaded

        `name` must be defined into `pymbe.df` as a pmb_object of type molecule.
        Args:
            name (str): Label of the protein type to be created. The protein type must be defined in `pymbe.df`
            espresso_system (cls): Instance of a system class from espressomd library.
            number_of_proteins (int): Number of proteins of type `name` to be created.
            positions (dict): dictionary with the coordinates that has the format {'ResidueNumber': {'initial_pos': [], 'chain_id': ''}}

        """

        import re 

        if number_of_proteins <=0:
            return
        if not self.check_if_name_is_defined_in_df(name=name,
                                                    pmb_type_to_be_defined='protein'):
            raise ValueError(f"{name} must correspond to a label of a pmb_type='protein' defined on df")

        self.copy_df_entry(name=name,column_name='molecule_id',number_of_copies=number_of_proteins)

        protein_index = self.np.where(self.df['name']==name)
        protein_index_list =list(protein_index[0])[-number_of_proteins:]
        used_molecules_id = self.df.molecule_id.dropna().drop_duplicates().tolist()

        for molecule_index in protein_index_list:          

            self.clean_df_row(index=int(molecule_index))
            if self.df['molecule_id'].isnull().values.all():
                molecule_id = 0        
            else:
                # check if a residue is part of another molecule
                check_residue_name = self.df[self.df['residue_list'].astype(str).str.contains(name)]
                pmb_type = self.df.loc[self.df['name']==name].pmb_type.values[0]                
                if not check_residue_name.empty and pmb_type == 'protein' :              
                    for value in check_residue_name.index.to_list():                  
                        if value not in used_molecules_id:                              
                            molecule_id = self.df.loc[value].molecule_id.values[0]                    
                            break
                else:
                    molecule_id = self.df['molecule_id'].max() +1

            #assigns molecule_id to the protein defined       
            self.add_value_to_df (key=('molecule_id',''),
                                index=int (molecule_index),
                                new_value=molecule_id, 
                                warning=False)
            
            # protein_center = self.np.random.random((1, 3))[0] *self.np.copy(espresso_system.box_l)

            protein_center = self.generate_coordinates_outside_sphere(espresso_system = espresso_system, min_dist = 1, max_dist=espresso_system.box_l[0]/2.0 , n_samples=1, center=[0,0,0])[0]
   
            for residue in positions.keys():

                residue_name = re.split(r'\d+', residue)[0]
                residue_number = re.split(r'(\d+)', residue)[1]
                residue_position = positions[residue]['initial_pos']
                
                residue_diameter  = positions[residue]['diameter']
                diameter = residue_diameter #*self.units.nm

                position = residue_position + protein_center
   
                particle_id = self.create_particle_in_espresso(name=residue_name,espresso_system=espresso_system,number_of_particles=1,position=[position], fix = True)

                index = self.df[self.df['particle_id']==particle_id[0]].index.values[0]

                #NOTE there is still a bug trying to add the diameter in each AA.
            
                # self.add_value_to_df(key=('diameter',''),
                #                             index=int (index),
                #                             new_value=diameter)

                self.add_value_to_df(key=('residue_id',''),
                                            index=int (index),
                                            new_value=residue_number)

                self.add_value_to_df(key=('molecule_id',''),
                                        index=int (index),
                                        new_value=molecule_id)
        return
    
    def activate_motion_of_rigid_object (self, name, espresso_system):
        '''
        Activates the motion of rigid object using the features of Virtual Sites from EsPRessoMD

        Args:
            name (str): Label of the protein type to be created. The protein type must be defined in `pymbe.df`
            espresso_system (cls): Instance of a system class from espressomd library.
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

    def center_molecule_in_simulation_box (self, molecule_id, espresso_system):

        """
        Centers the pmb object of type `name` in the center of the simulation box. Returns the position updated on the system.
        
        Args:
            molecule_id (int): ID of the molecue type to be centered. The molecule id must be defined in `df`
            espresso_system (cls): Instance of a system class from espressomd library.

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

    def calculate_center_of_mass_of_molecule (self,molecule_id, espresso_system):
        
        """
        Calculates the center of mass of type `name`
        Args:
            molecule_id (int): ID of the molecue type to be centered. The molecule id must be defined in `df`
            espresso_system (cls): Instance of a system class from espressomd library.
        Return:
            center_of_mass (lst): a list with the coordinates of the center of mass [ X, Y, Z]
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

    def generate_coordinates_outside_sphere (self, espresso_system, center, min_dist, max_dist, n_samples):

        """
        Generates coordinates outside a sphere 

        Args:
            espresso_system (cls): Instance of a system class from espressomd library.
            center (array): array with the coordinates of the center of the spheres.
            min_dist (int): minimun distance from the center to generate coordinates
            max_dist (int): maximum distance from the center to generate coordinates
            n_samples (int): number of sample points to generate inside the sphere

        """

        coord_list = []
        box_l = espresso_system.box_l[0]

        if not min_dist > 0: 
            raise ValueError (f'The value of {min_dist} must be a positive value')
        if not min_dist < max_dist:
            raise ValueError(f'The min_dist ({min_dist} must be lower than the max_dist ({max_dist}))')
        if not min_dist <= box_l/2.0 and max_dist <= box_l/2.0:
            raise ValueError(f'The min_dist and max_dist parameter should have values between 0 and {box_l/2.0} (box_l/2)')

        for _ in range(n_samples):

            rand = self.np.random.random() 
            rad = min_dist + rand*(max_dist-min_dist)

            coord = self.generate_trialvectors(center=center, radius=rad,n_samples=1)[0]
            coord_list.append (coord)

        return coord_list
    
    def write_output_vtf_file (self, espresso_system, n_frame):
    
        '''
        Writes the system information on a vtf file for each frame with additional information the format it follows is:
        'atom {particle_id} radius {} name {label} type {label}'

        Args: 
            espresso_system (cls): Instance of a system class from espressomd library.
            n_frame (int): number of current frame 

        '''
        number = 1000 + n_frame
        box = espresso_system.box_l[0]

        with open(f'frames/trajectory_{number}.vtf', mode='w+t') as coordinates:

            coordinates.write (f'unitcell {box} {box} {box} \n')
            
            for particle in espresso_system.part: 
                type_label = self.find_value_from_es_type(es_type=particle.type, column_name='label')
                coordinates.write (f'atom {particle.id} radius 1 name {type_label} type {type_label}\n' )

            coordinates.write (f'#timestep indexed\n')

            for particle in espresso_system.part:
                coordinates.write (f'{particle.id} \t {particle.pos[0]} \t {particle.pos[1]} \t {particle.pos[2]}\n')

        return 
    