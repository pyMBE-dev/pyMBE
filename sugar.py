class sugar_library():

    import pint
    units = pint.UnitRegistry()
    import numpy as np
    import random as rn
    import time

    # Default values    

    TEMPERATURE = 298.15 * units.K
    PARTICLE_SIZE = 0.355 * units.nm
    SEED=None
    N_A=6.02214076e23    / units.mol
    Kb=1.38064852e-23    * units.J / units.K
    e=1.60217662e-19 *units.C
    pi=3.14159265359
    initial_simulation_time=None
    stored_sugar_objects={}
    id_map={}
    type_map={}
    pka_set={}

    # Library output

    filename_parameters = 'simulation_parameters.txt'
     
    def __init__(self):

        # Default definitions of reduced units

        self.units.define(f'reduced_energy = {self.TEMPERATURE * self.Kb}')
        self.units.define(f'reduced_length = {self.PARTICLE_SIZE}')
        self.units.define(f'reduced_charge = 1*e')
        self.kT=self.TEMPERATURE*self.Kb
        self.initial_simulation_time=self.time.time()*self.units.s

    def print_reduced_units(self):

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

        # Change the unit registry of the initial simulation time

        self.initial_simulation_time=self.initial_simulation_time.magnitude*self.units.s

        self.print_reduced_units()

    def particle(self, name, es_type=None, q=0, diameter=None, acidity=None, epsilon=None):
        """
        Returns a sugar particle object. 
        Checks if the user has created another particle object with a shared type.
        Updates the sugar list of particles.
        Inputs:
        name: (string) label of the particle
        es_type: (int) particle type for espresso
        q: (int) particle charge
        Returns:
        particle (class) sugar particle object
        """
        class particle:

            pass
        
        particle.name=name
        particle.sugar_object_type='particle'
        particle.diameter=diameter
        particle.epsilon=epsilon
        particle.q=q

        if es_type is None:
            es_type=self.propose_unused_type()
        
        particle.es_type=es_type

        if self.check_particle_type_exists(particle=particle):
            raise ValueError("ERROR you have already created a particle object with the same type as ", particle.name)
        
        if acidity is None:
            acidity='inert'
        
        self.set_particle_acidity(particle=particle, acidity=acidity)
        self.store_sugar_object(sugar_object=particle, verbose=True)
        self.type_map[particle.name]=particle.es_type

        return particle

    def set_particle_acidity(self, particle, acidity):
        """
        Sets up a particle to become an acidic or basic particle
        Inputs:
        particle: (class) sugar particle object
        acidity: (string) defines the particle acidic properties can be 'inert', 'acidic' or 'basic'
        """

        acidity_valid_keys=['inert','acidic', 'basic']
        
        if acidity not in acidity_valid_keys:
            raise ValueError(particle.name+' provided acidity not supported, valid keys are ', acidity_valid_keys)

        particle.acidity=acidity

        if particle.acidity != 'inert':

            type_labels=['protonated','unprotonated']
            if particle.es_type is None:
                particle.es_type={}
                for label in type_labels:
                    particle.es_type[label]=self.propose_unused_type()
                    self.type_map[particle.name]=particle.es_type

            elif isinstance(particle.es_type,int):
                particle.es_type={'protonated': particle.es_type}
                self.type_map[particle.name]=particle.es_type
                particle.es_type['unprotonated']= self.propose_unused_type()
            elif isinstance(particle.es_type,dict):
                for key in particle.es_type.keys():
                    if key not in type_labels:
                        raise ValueError('Label not supported for particle type, valid options are ', type_labels, ' given ', particle.keys())
            else:
                raise ValueError('Unsopported variable type given in particle.type: ', particle.es_type)
        
        if particle.acidity == 'acidic':

            particle.q={'unprotonated': -1, 'protonated': 0}

        elif particle.acidity == 'basic':

            particle.q={'protonated': 1,  'unprotonated': 0}
        
        return

    def residue(self, name, central_bead, side_chains):
        """
        Returns a sugar residue object. 
        Updates the sugar list of residues.
        Inputs:
        name: (string) label of the residue
        central_bead: (instance of particle class of sugar) central bead of the residue
        side_chains: (list) list of particle or residue objects
        Returns:
        residue (class) sugar residue object
        """
        class residue:

            pass

        residue.name=name
        residue.central_bead=central_bead
        residue.side_chains=side_chains
        residue.sugar_object_type='residue'
        self.store_sugar_object(sugar_object=residue, verbose=True)

        return residue

    def molecule(self, name, residue_list):
        """
        Returns a sugar molecule object. 
        Updates the sugar list of molecules.
        Inputs:
        name: (string) label of the molecule
        sequence: (list) ordered list of residues of the molecule 
        Returns:
        molecule (class) sugar molecule object
        """
        class molecule:

            pass

        molecule.residue_list=residue_list
        molecule.name=name
        molecule.sugar_object_type='molecule'
        self.store_sugar_object(sugar_object=molecule, verbose=True)

        return molecule

    def peptide(self, name, sequence, model):
        """
        Returns a sugar molecule object, . 
        Updates the sugar list of molecules.
        Inputs:
        name: (string) label of the peptide
        sequence: (string) string with the peptide sequence
        Returns:
        peptide (class) sugar molecule object
        """

        valid_keys=['1beadAA','2beadAA']

        if model not in valid_keys:

            raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')

        clean_sequence=self.protein_sequence_parser(sequence=sequence)
        residue_list=[]

        if model == '2beadAA':
            C_particle=self.stored_sugar_objects['particle']['CA']
        
        if 'residue' not in self.stored_sugar_objects.keys():
            self.stored_sugar_objects['residue']={}

        for residue_name in clean_sequence:
            
            if residue_name not in self.stored_sugar_objects['residue'].keys():
                
                AA_particle=self.stored_sugar_objects['particle'][residue_name]

                if model == '1beadAA':

                    central_bead=AA_particle
                    side_chains=[]

                elif model == '2beadAA':
                    
                    # terminal groups + glycine are represented with only 1 bead
                    
                    if residue_name in ['c','n']: 
                        central_bead=AA_particle
                        side_chains=[]
                    if residue_name == 'G':
                        central_bead=C_particle
                        side_chains=[]
                    else:
                        central_bead=C_particle
                        side_chains=[AA_particle]

                residue=self.residue(name=residue_name, central_bead=central_bead,side_chains=side_chains)
            else:
                residue=self.stored_sugar_objects['residue'][residue_name]
            
            residue_list.append(residue)
        peptide=self.molecule(name=name, residue_list=residue_list)
        peptide.model=model
        peptide.sequence=clean_sequence
        return peptide

    def store_sugar_object(self,sugar_object, verbose=False):
        """
        Stores object in sugar for bookeeping
        Input:
        object: (class) particle object from sugar library
        verbose: (boolean) prints a warning if the object is already stored
        """

        if not self.check_sugar_object_stored(sugar_object=sugar_object, verbose=verbose):
            if sugar_object.sugar_object_type not in self.stored_sugar_objects:
                self.stored_sugar_objects[sugar_object.sugar_object_type]={}
            self.stored_sugar_objects[sugar_object.sugar_object_type][sugar_object.name]=sugar_object
        return    
        
    def check_sugar_object_stored(self, sugar_object, verbose=False):
        """
        Checks if the user has already created a  object with the same name
        Input:
        object: (class) particle object from sugar library
        verbose: (boolean) prints a warning if the object is already stored
        """

        if sugar_object.sugar_object_type not in self.stored_sugar_objects.keys():

            return False

        if sugar_object.name not in self.stored_sugar_objects[sugar_object.sugar_object_type].keys():

            return False

        else:

            if verbose:
                print('WARNING, you have already created a particle object with the same name')
            return True

    def define_bond(self, bond, particle1, particle2):
        """
        Adds a bond to the bond storage of sugar
        Inputs:
        bond: (class) instance of a bond object from espressomd library
        particle1: (class) instance of a particle object from sugar library
        particle2: (class) instance of a particle object from sugar library
        """


        bond_key=frozenset([particle1.name, particle2.name])
        exists_bond=self.check_bond_defined(particle1=particle1, particle2=particle2)
        
        if exists_bond:

            print('WARNING, you have already defined a bond between particle1 and particle2')
            print('The previously defined bond have been overwritten')

        self.stored_sugar_objects['bonds'][bond_key]=bond

        return

    def create_particle_in_espresso(self, particle, espresso_system, position=None, state=None, id=None):
        '''
        Creates a particle in the espresso system.
        particle = instance of a particle object from sugar library
        espresso_system = instance of a system object from espressomd library
        position = array with the cartesian coordinates to create the particle, by default random
        state = for particles with more than one possible state, defines in which state are created
        id = desired particle id
        '''

        if particle.sugar_object_type != 'particle':

            raise ValueError("particle must be an instance of a sugar particle object")

        if position is None:
        
        # By default, the particle is placed in a random position
            
            position=self.np.random.random((1, 3))[0] *self.np.copy(espresso_system.box_l)
                
        if state is None:

            if particle.acidity in ['acidic','basic']:

                state='protonated'
        
        if isinstance(particle.q, int) or isinstance(particle.q, float) :
            
            q=particle.q
        
        elif isinstance(particle.q, dict):

            if state is None:

                state=self.rn.choice(q.keys())

            q=particle.q[state]

        else:

            raise ValueError("Unvalid charge for bead: ", particle.name, " given ", particle.q)

        if isinstance(particle.es_type, int):
            
            es_type=particle.es_type

        elif isinstance(particle.es_type, dict):

            if state is None:

                state=self.rn.choice(q.keys())

            es_type=particle.es_type[state]

        else:

            raise ValueError("Unvalid type for bead: ", particle.name, " given ", particle.es_type)

        
        if id is None:

            if len(espresso_system.part[:].id) == 0:
                
                bead_id=0

            else:
                
                bead_id=max(espresso_system.part[:].id)+1
        espresso_system.part.add(id=[bead_id], pos=[position], type=[es_type], q=[q])

        # particle id bookeeping

        if 'particle' not in self.id_map.keys():
            self.id_map['particle']={} 

        if particle.name in self.id_map['particle'].keys():
            self.id_map['particle'][particle.name].append(bead_id)
        else:
            self.id_map['particle'][particle.name]=[bead_id]
        
        return  bead_id
        
    def create_residue_in_espresso(self, residue, espresso_system, central_bead_position=None, use_default_bond=False, backbone_vector=None):
        """
        Creates a residue in the espresso system
        Inputs:
        residue: instance of a residue object as defined in sugar library
        espresso_system: instance of a system object of espresso library
        position: (list) coordinates where the central bead will be created, random by default
        use_default_bond: (boolean, default=False) returns the default bond if no bond has been defined between  particle1 and particle2
        Returns:
        particle_ids: (list) list with the ids of the particles created
        central_bead_id: (id) id of the central bead of the residue
        """
        
        if residue.sugar_object_type != 'residue':

            raise ValueError("residue must be an instance of a sugar residue object")

        residue_ids_dict={}        

        # create the principal bead
        if self.check_sugar_object_stored(sugar_object=residue.central_bead):
            
            central_bead_id=self.create_particle_in_espresso(particle=residue.central_bead, espresso_system=espresso_system, position=central_bead_position)
            residue_ids_dict['central-'+residue.central_bead.name]=[central_bead_id]

        else:

            raise ValueError("Residue.central_bead must contain a particle object.")

        central_bead_position=espresso_system.part[central_bead_id].pos
                
        # create the lateral beads        
        
        for sg_object in residue.side_chains:

            if sg_object.sugar_object_type == 'particle': # Only one bead
                
                bond=self.search_bond(particle1=residue.central_bead, particle2=sg_object, hard_check=True, use_default_bond=use_default_bond)
                if backbone_vector is None:
                    bond_vector=self.generate_trialvectors(bond.params.get('r_0'))
                else:
                    bond_vector=self.generate_trial_perpendicular_vector(vector=backbone_vector,magnitude=bond.params.get('r_0'))
                bead_position=central_bead_position+bond_vector
                side_bead_id=self.create_particle_in_espresso(particle=sg_object, espresso_system=espresso_system, position=bead_position)
                espresso_system.part[central_bead_id].add_bond((bond, side_bead_id))
                if 'side-'+sg_object.name in residue_ids_dict.keys():
                    residue_ids_dict['side-'+sg_object.name].append(side_bead_id)
                else:
                    residue_ids_dict['side-'+sg_object.name]=[side_bead_id]

            elif sg_object.sugar_object_type == 'residue': # More than one bead
                
                bond=self.search_bond(particle1=residue.central_bead, particle2=sg_object.central_bead, hard_check=True, use_default_bond=use_default_bond)
                if backbone_vector is None:
                    bond_vector=self.generate_trialvectors(bond.params.get('r_0'))
                else:
                    bond_vector=self.generate_trial_perpendicular_vector(vector=backbone_vector,magnitude=bond.params.get('r_0'))
                residue_position=central_bead_position+bond_vector
                lateral_residue_ids_dict=self.create_residue_in_espresso(residue=sg_object, espresso_system=espresso_system, central_bead_position=residue_position)
                lateral_residue_central_bead_id=next(value for key,value in lateral_residue_ids_dict.items() if 'central-' in key)[0]
                espresso_system.part[central_bead_id].add_bond((bond, lateral_residue_central_bead_id))
                
                for key in lateral_residue_ids_dict:
                    
                    if 'central-' in key:
                        clean_key='side-'+key.replace('central-','')
                    elif 'side-' in key:
                        clean_key=key

                    for id in lateral_residue_ids_dict[key]:

                        if clean_key in residue_ids_dict.keys():
                            residue_ids_dict[clean_key].append(id)
                        else:
                            residue_ids_dict[clean_key]=[id]

            else:

                raise ValueError("Residue.side_chains must contain a list of particle or residue objects.")

        # for particle id bookeeping

        if 'residue' not in self.id_map.keys():
            self.id_map['residue']={} 

        if residue.name not in self.id_map['residue'].keys():
            self.id_map['residue'][residue.name]=[residue_ids_dict]
        else:
            self.id_map['residue'][residue.name].append(residue_ids_dict)

        return  residue_ids_dict

    def create_molecule_in_espresso(self, molecule, espresso_system, first_residue_position=None, use_default_bond=False):
        """
        Creates a molecule in the espresso system
        Inputs:
        molecule: instance of a molecule object as defined in sugar library
        espresso_system: instance of a system object of espresso library
        first_residue_position: (list) coordinates where the first_residue_position will be created, random by default
        use_default_bond: (boolean, default=False) returns the default bond if no bond has been defined between  particle1 and particle2
        Returns:
        particle_ids: (list) list with the ids of the particles created
        """
        
        if molecule.sugar_object_type != 'molecule':

            raise ValueError("molecule must be an instance of a sugar residue object")

        first_residue=True
        molecule_dicts=[]
        
        for residue in molecule.residue_list:

            if first_residue:

                residue_position=first_residue_position
                backbone_vector=self.generate_trialvectors(1)
                residue_ids_dict=self.create_residue_in_espresso(residue=residue, espresso_system=espresso_system, central_bead_position=first_residue_position,  use_default_bond= use_default_bond, backbone_vector=backbone_vector)
                central_bead_id=next(value for key,value in residue_ids_dict.items() if 'central-' in key)[0]
                previous_residue=residue
                residue_position=espresso_system.part[central_bead_id].pos
                previous_residue_id=central_bead_id
                first_residue=False
                
                
            else:

                bond=self.search_bond(particle1=residue.central_bead, particle2=previous_residue.central_bead, hard_check=True, use_default_bond=use_default_bond)                
                residue_position=residue_position+backbone_vector*bond.params.get('r_0')
                residue_ids_dict=self.create_residue_in_espresso(residue=residue, espresso_system=espresso_system, central_bead_position=residue_position,use_default_bond= use_default_bond)
                central_bead_id=next(value for key,value in residue_ids_dict.items() if 'central-' in key)[0]
                espresso_system.part[central_bead_id].add_bond((bond, previous_residue_id))
                previous_residue_id=central_bead_id
                previous_residue=residue

            molecule_dicts.append(residue_ids_dict)

        # for particle id bookeeping

        if 'molecule' not in self.id_map.keys():
            self.id_map['molecule']={} 

        if molecule.name not in self.id_map['molecule'].keys():
            self.id_map['molecule'][molecule.name]=[molecule_dicts]
        else:    
            self.id_map['molecule'][molecule.name].append(molecule_dicts)
        return molecule_dicts

    def propose_unused_type(self):
        """
        Searches in sugar database the currently used particle types and returns an unused type
        Returns:
        unused_type: (int) unused particle type
        """
        type_list=[]
        if not self.type_map:
            
            unused_type=0

        else:

            for value in self.type_map.values():
                if isinstance(value,int):
                    type_list.append(value)
                elif isinstance(value,dict):
                    for es_type in value.values():
                        type_list.append(es_type)
                
            unused_type=max(type_list)+1

        return unused_type

    def check_bond_defined(self, particle1, particle2):
        """
        Checks if the user has already defined a bond between particle1 and particle2 in sugar
        Input:
        particle1 = instance of a particle object from sugar library
        particle2 = instance of a particle object from sugar library
        """

        if 'bonds' not in self.stored_sugar_objects.keys():
            self.stored_sugar_objects['bonds']={}
            return False

        bond_key=frozenset([particle1.name, particle2.name])

        if bond_key not in self.stored_sugar_objects['bonds'].keys():

            return False

        else:

            return True

    def check_particle_type_exists(self, particle):
        """
        Checks if the user has already created a particle object with the same name as particle
        Input:
        particle: sugar particle object
        """

        if 'particle' not in self.stored_sugar_objects.keys():
            return False

        for stored_particle in self.stored_sugar_objects['particle'].values():

            if particle.es_type in self.get_particle_types(particle=stored_particle):

                return True

        return False

    def search_bond(self, particle1, particle2, hard_check=False, use_default_bond=False):
        """
        Searches for a bond between particle1 and particle2 defined in sugar library
        Input:
        particle1 = instance of a particle object from sugar library
        particle2 = instance of a particle object from sugar library
        hard_check: (boolean, default=False) break if the bond is not defined
        use_default_bond: (boolean, default=False) returns the default bond if no bond has been defined between  particle1 and particle2
        Returns:
        bond: (class) instance of a bond object from espressomd library
        """

                
        if self.check_bond_defined(particle1=particle1, particle2=particle2):

            bond_key=frozenset([particle1.name, particle2.name])
            bond=self.stored_sugar_objects['bonds'][bond_key]
            return bond
        
        else:

            if  use_default_bond:
            
                if 'default' not in self.stored_sugar_objects['bonds'].keys():

                    raise ValueError('Default bond is not defined')

                else:

                    return self.stored_sugar_objects['bonds']['default']

            else:
                
                print("Bond not defined between particles ", particle1.name, " and ", particle2.name)
                
                if hard_check:

                    exit()
                else:

                    return

    def add_bonds_to_espresso(self, espresso_system):
        """
        Adds all the bonds stored in sugar to the espressomd system, 
        including the default bond if has been defined
        Inputs:
        espresso_system: instance of a system object of espresso library
        """

        for bond in self.stored_sugar_objects['bonds'].values():
            espresso_system.bonded_inter.add(bond)

        return
    
    def load_interaction_parameters(self,filename, verbose=False):
        """
        Loads the interaction parameters stored in filename into sugar
        Inputs:
        filename: (string) name of the file to be read
        verbose: (boolean) switches on/off the reading prints
        """
        import json
        from espressomd import interactions

        particle_param_list=[]
        residue_param_list=[]
        molecule_param_list=[]
        peptide_param_list=[]
        bond_param_list=[]

        with open(filename) as f:
            for line in f:
                if line[0] == '#':
                    continue
                param_dict=json.loads(line)
                object_type=param_dict['object_type']
                if object_type == 'particle':
                    particle_param_list.append(param_dict)
                elif object_type == 'residue':
                    residue_param_list.append(param_dict)
                elif object_type == 'molecule':
                    molecule_param_list.append(param_dict)
                elif object_type == 'peptide':
                    peptide_param_list.append(param_dict)
                elif object_type == 'bond':
                    bond_param_list.append(param_dict)
                else:
                    raise ValueError(object_type+' is not a known sugar object type')
                if verbose:
                    print('Added: '+line)

        without_units=['q','es_type','acidity']
        with_units=['diameter','epsilon']

        for particle_param in particle_param_list:
            not_requiered_attributes={}    
            for not_requiered_key in without_units+with_units:
                if not_requiered_key in particle_param.keys():
                    if not_requiered_key in with_units:
                        not_requiered_attributes[not_requiered_key]=self.create_variable_with_units(variable_dict=particle_param.pop(not_requiered_key))
                    elif not_requiered_key in without_units:
                        not_requiered_attributes[not_requiered_key]=particle_param.pop(not_requiered_key)
                else:
                    not_requiered_attributes[not_requiered_key]=None

            self.particle(name=particle_param.pop('name'),
                            es_type=not_requiered_attributes.pop('es_type'),
                            q=not_requiered_attributes.pop('q'),
                            acidity=not_requiered_attributes.pop('acidity'),
                            diameter=not_requiered_attributes.pop('diameter'),
                            epsilon=not_requiered_attributes.pop('epsilon'))

        for residue_param in residue_param_list:
            central_bead=self.stored_sugar_objects['particle'][residue_param.pop('central_bead_name')]
            side_chains=[]
            for side_chain_name in self.residue_param.pop('side_chains_names'):

                if side_chain_name in self.stored_sugar_objects['particle'].keys() and side_chain_name in self.stored_objects['residue'].keys():
                    raise ValueError(side_chain_name+ 'is defined both as a particle name and a residue name, please rename to avoid the ambiguity')

                elif side_chain_name in self.stored_sugar_objects['particle'].keys():
                    side_chains.append(self.stored_sugar_objects['particle'][side_chain_name])
                elif side_chain_name in self.stored_sugar_objects['residue'].keys():
                    side_chains.append(self.stored_sugar_objects['residue'][side_chain_name])
                else:
                    raise ValueError('Objects in residue side chains must be either particles or residues')
                
            self.residue(name=residue_param.pop('name'),
            central_bead=central_bead,
            side_chains=side_chains)

        for molecule_param in molecule_param_list:
            residue_name_list=molecule_param.pop('residue_name_list')
            residue_list=[]
            for residue_name in residue_name_list:
                residue_list.append(self.stored_sugar_objects['residue'][residue_name])

            self.molecule(name=molecule_param.pop('name'),
                        residue_list=residue_list)
        
        for peptide_param in peptide_param_list:

            self.peptide(name=peptide_param.pop('name'),
                        sequence=peptide_param.pop('sequence'),
                        model=peptide_param.pop('model'))

        for bond_param in bond_param_list:
            
            name1=bond_param.pop('name1')
            name2=bond_param.pop('name2')
            particle1=self.stored_sugar_objects['particle'][name1]
            particle2=self.stored_sugar_objects['particle'][name2]
            bond_type=bond_param.pop('bond_type')

            if bond_type == 'harmonic':

                k=self.create_variable_with_units(variable_dict=bond_param.pop('k'))
                r_0=self.create_variable_with_units(variable_dict=bond_param.pop('r_0'))
                bond = interactions.HarmonicBond(k=k.to('reduced_energy / reduced_length**2').magnitude, r_0=r_0.to('reduced_length').magnitude)

            else:

                raise ValueError("current implementation of sugar only supports harmonic bonds")

            self.define_bond(bond=bond, particle1=particle1, particle2=particle2)
            
        return
    
    def load_pka_set(self,filename, verbose=False):
        """
        Loads the parameters stored in filename into sugar
        Inputs:
        filename: (string) name of the file to be read
        verbose: (boolean) switches on/off the reading prints
        """
        import json
        from espressomd import interactions

        pKa_list=[]

        with open(filename) as f:
            for line in f:
                if line[0] == '#':
                    continue
                param_dict=json.loads(line)
                pKa_list.append(param_dict)    
                if verbose:
                    print('Added: '+line)

        if 'particle' not in self.stored_sugar_objects.keys():
            self.stored_sugar_objects['particle']={}

        for pka_set in pKa_list:
            self.check_pka_set(pka_set=pka_set)
            for pka_key in pka_set: 
                if pka_key in self.pka_set.keys():

                    print("WARNING overwritting stored pKa value for ", pka_key)

                self.pka_set[pka_key]=pka_set[pka_key]
                if pka_key in self.stored_sugar_objects['particle'].keys():
                    self.set_particle_acidity(particle=self.stored_sugar_objects['particle'][pka_key],acidity=pka_set[pka_key]['acidity'])
            
        return

    def define_default_bond(self, bond):
        """
        Defines the default bond, used when a bonded interaction is not parametrized
        Input:
        bond: (class) instance of a bond object from espressomd library
        """
        if 'bonds' not in self.stored_sugar_objects.keys():
            self.stored_sugar_objects['bonds']={}
        self.stored_sugar_objects['bonds']['default']=bond

        return

    def define_default_lennard_jones(self, sigma, epsilon, cutoff, offset, shift):
        """
        Defines the default lennard jones interaction, used when the interaction between two particle types is not properly parametrized
        Input:
        See section 6.1.2 of the user guide from espresso for a complete description of the parameters
        """

        default_param={}
        default_param={'epsilon' : epsilon, 
                            "sigma" : sigma, 
                            'cutoff' : cutoff,
                            'offset' : offset, 
                            'shift' : shift}
        if 'LennardJones' not in self.stored_sugar_objects.keys():
            self.stored_sugar_objects['LennardJones']={}

        self.stored_sugar_objects['LennardJones']['default']=default_param

    def create_variable_with_units(self, variable_dict):
        """
        Returns a pint object with the value and units defined in variable_dict
        Inputs:
        variable_dict:(dict) contains the value and the units of the variable
        Returns:
        variable_with_units:(pint object) pint object containing the value and the desired units
        """
        
        value=variable_dict.pop('value')
        units=variable_dict.pop('units')
        variable_with_units=value*self.units(units)

        return variable_with_units
    
    def get_particle_types(self, particle):
        """
        Returns all the particle types stored in particle
        Inputs:
        particle: sugar particle object
        Returns:
        type_list: (list) list of types in particle
        """
        type_list=[]
        if isinstance(particle.es_type, int):
            
            type_list.append(particle.es_type)

        elif isinstance(particle.es_type, dict):

            for es_type in particle.es_type.values():

                type_list.append(es_type)

        else:

            raise ValueError("Unvalid type for bead: ", particle.name, " given ", particle.es_type)

        return type_list

    def protein_sequence_parser(self, sequence):
        '''
        Reads the input residue sequence and
        1) Checks that all residues are in the parameters key
        2) Transforms the aminoacids on the sequence from three letter to one letter format

        Input:
        sequence: (string or list) aminoacid sequence

        Output:
        clean_sequence: (string) 
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

    def setup_constantpH_reactions_in_espresso(self, counter_ion, exclusion_radius=None, pka_set=None):
        """
        Set up the Acid/Base reactions for acidic/basidic residues in mol. The reaction steps are done following the constant pH ensamble procedure. 

        Inputs:
        
        counter_ion:(class) particle class object as defined in sugar library
        exclusion_radius:(float, optional) exclusion radius for the constant pH ensamble
        pka_set:(dict,optional) dictionary with the desired pka_set, by default the one stored in sugar will be used
        Output:
        RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
        sucessfull_reactions_labels: (list) list with the labels of the reactions that has been set up by sugar 

        """

        from espressomd import reaction_ensemble

        if exclusion_radius is None:

            exclusion_radius=self.search_largest_particle_diameter_in_sugar()

        if self.SEED is None:

            self.create_random_seed()

        if pka_set is None:
            pka_set=self.pka_set    

        self.check_pka_set(pka_set=pka_set)

        RE = reaction_ensemble.ConstantpHEnsemble(temperature=self.kT.to('reduced_energy').magnitude, exclusion_radius=exclusion_radius.magnitude, seed=self.SEED)
        sucessfull_reactions_labels=[]

        for name in pka_set.keys():

            if name not in self.type_map.keys():

                print('WARNING: the acid-base reaction of ' + name +' has not been set up because its espresso type is not defined in sg.type_map')
                continue

            if name in self.stored_sugar_objects['particle'].keys():
                particle=self.stored_sugar_objects['particle'][name]
                self.set_particle_acidity(particle=particle,acidity=pka_set[name]['acidity'])

            gamma=10**-pka_set[name]['pka_value']

            if pka_set[name]['acidity'] == 'acidic':
                charge_protonated=0
                charge_unprotonated=-1
            elif pka_set[name]['acidity'] == 'basic':
                charge_protonated=1
                charge_unprotonated=0

            RE.add_reaction(gamma=gamma,
                                    reactant_types=[self.type_map[name]["protonated"]],
                                    reactant_coefficients=[1],
                                    product_types=[self.type_map[name]["unprotonated"], counter_ion.es_type],
                                    product_coefficients=[1,1],
                                    default_charges={self.type_map[name]["unprotonated"]: charge_unprotonated,
                                    self.type_map[name]["protonated"]: charge_protonated,
                                    counter_ion.es_type: counter_ion.q})
            sucessfull_reactions_labels.append(name)
            
        return RE, sucessfull_reactions_labels

    def check_pka_set(self, pka_set):
        """"
        Checks that the given pka_set has the formatting expected by sugar
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_entry in pka_set.values():
                if required_key not in pka_entry.keys():
                    raise ValueError('missing a requiered key ', required_keys, 'in the following entry of pka_set', pka_entry)

    def generate_trialvectors(self,mag):
        """
        Generates a random 3D unit vector (direction) with a uniform spherical distribution
        Algorithm from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
        Input:
        
        mag: magnitude of the vector
        
        Output:
        
        vec: random vector of dimension 3
        """

        phi = self.np.random.uniform(0,self.np.pi*2)
        costheta = self.np.random.uniform(-1,1)

        theta = self.np.arccos( costheta )
        x = self.np.sin( theta) * self.np.cos( phi )
        y = self.np.sin( theta) * self.np.sin( phi )
        z = self.np.cos( theta )

        vec=self.np.array([x,y,z])
        vec=vec*mag

        return vec

    def generate_trial_perpendicular_vector(self, vector, magnitude=1):
        """
        Generates a random vector perpendicular to the input vector
        Inputs:
        vector: (list) array with the componentes of the input vector
        magnitude: (float) magnitude desired for the output vector
        Returns:
        perpendicular_vector: (numpy array) random vector perpendicular to the input vector
        """

        np_vec=self.np.array(vector)
        if np_vec[1] == 0 and np_vec[2] == 0 and np_vec[0] == 1:
            raise ValueError('zero vector')
            
        return self.np.cross(np_vec, self.generate_trialvectors(1))*magnitude

    def setup_electrostatic_interactions_in_espresso(self, espresso_system, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3):
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

        #Initial checks

        valid_methods_list=['p3m', 'DH']

        if method not in valid_methods_list:

            raise ValueError('provided an unknown label for method, valid values are', valid_methods_list)

        if c_salt is None and method == 'DH':

            raise ValueError('Please provide the added salt concentration c_salt to settup the Debye-Huckel potential')
            

        BJERRUM_LENGTH = self.e.to('reduced_charge')**2 / (4 * self.units.pi * self.units.eps0 * solvent_permittivity * self.kT.to('reduced_energy'))

        print('\n Bjerrum length ', BJERRUM_LENGTH.to('nm'), '=', BJERRUM_LENGTH.to('reduced_length'))

        COULOMB_PREFACTOR=BJERRUM_LENGTH.to('reduced_length') * self.kT.to('reduced_energy') 
        
        if c_salt is not None:

            if c_salt.check('[substance] [length]**-3'):

                KAPPA=1./self.np.sqrt(8*self.units.pi*BJERRUM_LENGTH*self.N_A*c_salt)

            elif c_salt.check('[length]**-3'):
                
                KAPPA=1./self.np.sqrt(8*self.units.pi*BJERRUM_LENGTH*c_salt)

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

    def minimize_espresso_system_energy(self, espresso_system, skin=1, gamma=1, Nsteps=10000, time_step=1e-5, max_displacement=0.1, verbose=True, reset=True):
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

    def setup_langevin_dynamics_in_espresso(self,espresso_system, time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
        """
        Sets up Langevin Dynamics in espressomd.
        espresso_system: instance of espressmd system class
        time_step: time s

        """        

        kT=self.TEMPERATURE*self.units.k

        if self.SEED is None:

            # Take the random seed from the system time
            self.create_random_seed()
            
        espresso_system.time_step=time_step
        espresso_system.integrator.set_vv()
        espresso_system.thermostat.set_langevin(kT=kT.to('reduced_energy').magnitude, gamma=gamma, seed=self.SEED)

        # Optimize the value of skin

        if tune_skin:

            print("\n*** Optimizing skin ... ***")

            if max_skin is None:

                max_skin=espresso_system.box_l[0]/2

            espresso_system.cell_system.tune_skin(min_skin=min_skin, max_skin=max_skin, tol=tolerance, int_steps=int_steps, adjust_max_skin=adjust_max_skin)

            print("Optimized skin value: ", espresso_system.cell_system.skin, "\n")

        return

    def create_random_seed(self):
        """
        Creates the seed for the random number generator from the system hour
        """

        import time 

        SEED=int(time.time())
        print('\n The chosen seed for the random number generator is ', SEED)
        self.SEED=SEED

    def block_analyze(self, input_data, n_blocks=16):
        '''         
        Performs a binning analysis of input_data. 
        Divides the samples in ``n_blocks`` equispaced blocks
        and returns the mean, its uncertainty, the correlation time 
        and the block size        
        '''

        data = self.np.asarray(input_data)
        block = 0
        # this number of blocks is recommended by Janke as a reasonable compromise
        # between the conflicting requirements on block size and number of blocks
        block_size = int(data.shape[1] / n_blocks)
        print(f"block_size: {block_size}")
        # initialize the array of per-block averages
        block_average = self.np.zeros((n_blocks, data.shape[0]))
        # calculate averages per each block
        for block in range(n_blocks):
            block_average[block] = self.np.average(data[:, block * block_size: (block + 1) * block_size], axis=1)
        # calculate the average and average of the square
        av_data = self.np.average(data, axis=1)
        av2_data = self.np.average(data * data, axis=1)
        # calculate the variance of the block averages
        block_var = self.np.var(block_average, axis=0)
        # calculate standard error of the mean
        err_data = self.np.sqrt(block_var / (n_blocks - 1))
        # estimate autocorrelation time using the formula given by Janke
        # this assumes that the errors have been correctly estimated
        tau_data = self.np.zeros(av_data.shape)
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
        if self.np.any(block_size / tau_data < threshold):
            print("\nWarning: some blocks may contain less than ", threshold, "uncorrelated samples."
          "\nYour error estimated may be unreliable."
          "\nPlease, check them using a more sophisticated method or run a longer simulation.")
            print("? block_size/tau > threshold ? :", block_size/tau_data > threshold)
        else:
            print("\nAll blocks seem to contain more than ", threshold, "uncorrelated samples.\
            Error estimates should be OK.")

        return av_data, err_data, tau_data, block_size

    def write_progress(self, step, total_steps):
        """
            Writes the progress of the loop and estimates the time for its completion
            
            Inputs:
            step: (int) actual step of the loop
            total_steps: (int) total number of loop steps

            Assumptions:
            It assumes that the simulation starts with step = 0
        """
        
        time_act=self.time.time()*self.units.s
        perc_sim=100 *(step+1) / (total_steps)
        time_per_step= (time_act - self.initial_simulation_time)/(step+1)
        remaining_time=(total_steps - step +1)*time_per_step
        elapsed_time=time_act-self.initial_simulation_time

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

    def search_particles(self, sugar_object, sort_different_particles=False):
        """
        Searches for all particles in object 
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        sort_different_particles:(boolean) filters the particle_list to return only different particles
        Returns:
        particle_list:(list) list of all particles in object
        """
        particle_list=[]
        
        if sugar_object.sugar_object_type == 'particle':
            particle_list.append(object)
        elif sugar_object.sugar_object_type == 'residue':
            particle_list.append(sugar_object.central_bead)
            for side_object in sugar_object.side_chains:
                if side_object.sugar_object_type=='particle': 
                    side_particle_list=[side_object]
                elif side_object.sugar_object_type=='residue':
                    side_particle_list=self.search_particles(sugar_object=side_object)
                particle_list+=side_particle_list
            
        elif sugar_object.sugar_object_type == 'molecule':
            for residue in sugar_object.residue_list:
                particle_list+=self.search_particles(sugar_object=residue)
        if sort_different_particles:
            sorted_particle_list=[]
            for particle_name in [particle.name for particle in particle_list]:
                if particle_name not in [particle.name for particle in sorted_particle_list]:
                    sorted_particle_list.append(self.sugar_stored_objects['particle'][particle_name])
                
            return sorted_particle_list
        else:

            return particle_list

    def count_particles(self, sugar_object):
        """
        Counts how many particle of each type are in object and stores them in a dictionary
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        Returns:
        particle_number:(dict) number of particles of each type in object
        """

        particle_list=self.search_particles(sugar_object=sugar_object,sort_different_particles=False)

        particle_number={}

        for particle in particle_list:

            if particle.name in particle_number.keys():
                particle_number[particle.name]+=1
            else:
                particle_number[particle.name]=1

        return particle_number

    def count_titrable_particles(self, sugar_object):
        """
        Counts how many titrable particle of each type are in object and stores them in a dictionary
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        Returns:
        titrable_particle_number:(dict) number of particles of each type in object
        """

        particle_dict=self.count_particles(sugar_object=sugar_object)
        titrable_particle_number={}

        for particle_name in particle_dict.keys():

            if particle_name in self.pka_set.keys():
                titrable_particle_number[particle_name]=particle_dict[particle_name]

        return titrable_particle_number

    def calculate_HH(self, sequence, pH=None, pka_set=None):
        """
        Calculates the ideal Henderson-Hasselbach titration curve in the given pH range

        Inputs:
        pH: (list of floats) list of pH values
        sequence: (list of strings or string)

        Outputs:
        Z_HH: (list of floats) Henderson-Hasselbach prediction of the protein charge in the given pH
        """

        if pH is None:

            pH=self.np.linspace(2,12,50)    

        if pka_set is None:
            pka_set=self.pka_set
        
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

    def create_added_salt_in_espresso(self, espresso_system, cation, anion, c_salt):
        
        """
        Adds a c_salt concentration of cations and anions to the espresso_system

        Inputs:
        espresso_system: (class)espresso class object with all system variables.
        cation: (class) particle class object, as defined in sugar, with a positive charge
        anion: (class) particle class object, as defined in sugar library, with a negative charge
        c_salt: (float) Added salt concentration
        
        Output:
        c_salt_calculated: (float) Calculated added salt concentration added to solution    
        """

        if cation.q <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ',cation.q)
        if anion.q >= 0:
            raise ValueError('ERROR anion charge must be positive, charge ', anion.q)

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

        N_cation=N_ions*abs(anion.q)
        N_anion=N_ions*abs(cation.q)
        
        for _ in range(N_cation):
            self.create_particle_in_espresso(espresso_system=espresso_system, particle=cation)
        
        for _ in range(N_anion):
            self.create_particle_in_espresso(espresso_system=espresso_system, particle=anion)

        print('\n Added salt concentration of ', c_salt_calculated.to('mol/L'), 'given by ', N_cation, ' cations and ', N_anion, ' anions')
        
        return c_salt_calculated

    def setup_espresso_to_track_ionization(self, espresso_system):
        """
        Sets up espresso to track the average number of particles of the acid/base particles 
        
        Inputs:
        espresso_system: espresso class object with all espresso_system variables
        """

        acidbase_types_list=[]

        for particle in self.stored_sugar_objects['particle'].values():

            if particle.name in self.pka_set.keys():

                for acidbase_type in particle.es_type.values():

                    acidbase_types_list.append(acidbase_type)

        espresso_system.setup_type_map(acidbase_types_list)

        return

    def get_ids_from_sugar(self,sugar_object):
        """
        Returns all particles ids stored in sugar that match the sugar_object.
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        Returns:
        id_list:(list) list with corresponding ids.
        """

        self.check_sugar_object_stored(sugar_object=sugar_object)

        if sugar_object.sugar_object_type == 'particle':
            id_list=[]
            for id in self.id_map['particle'][sugar_object.name]:
                id_list.append([id]) 
            return id_list

        elif sugar_object.sugar_object_type == 'residue':
            all_residue_id_list=[]
            for residue_dict in self.id_map['residue'][sugar_object.name]:
                residue_id_list=[]
                for id_list in residue_dict.values():
                    for id in id_list:
                        residue_id_list.append(id)
                all_residue_id_list.append(residue_id_list)
            return all_residue_id_list

        elif sugar_object.sugar_object_type == 'molecule':
            all_molecule_id_list=[]
            for molecule in self.id_map['molecule'][sugar_object.name]:
                molecule_id_list=[]
                for residue_dict in molecule:
                    for id_list in residue_dict.values():    
                        for id in id_list:
                            molecule_id_list.append(id)
                all_molecule_id_list.append(molecule_id_list)

            return all_molecule_id_list
        return

    def create_sugar_object_in_espresso(self, sugar_object, espresso_system, position=None, use_default_bond=False):
        """
        Creates all particles contained in the sugar object in the system of espresso
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        espresso_system: (class)espresso class object with all system variables.
        use_default_bond: (boolean, default=False) returns the default bond if no bond has been defined between  particle1 and particle2
        """

        allowed_objects=['particle','residue','molecule']
        if sugar_object.sugar_object_type not in allowed_objects:
            raise ValueError('Object type not supported, supported types are ', allowed_objects)

        if sugar_object.sugar_object_type == 'particle':
            self.create_particle_in_espresso(particle=sugar_object, espresso_system=espresso_system, position=position)

        elif sugar_object.sugar_object_type == 'residue':
            self.create_residue_in_espresso(residue=sugar_object, espresso_system=espresso_system, central_bead_position=position,use_default_bond=use_default_bond)

        elif sugar_object.sugar_object_type == 'molecule':
            self.create_molecule_in_espresso(molecule=sugar_object, espresso_system=espresso_system, use_default_bond=use_default_bond, first_residue_position=position)

        return
    
    def destroy_sugar_object_in_system(self, sugar_object, espresso_system):
        """
        Destroys all particles contained in the object from the espresso system 
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        espresso_system: (class)espresso class object with all espresso_system variables.
        """
        allowed_objects=['particle','residue','molecule']
        if sugar_object.sugar_object_type not in allowed_objects:
            raise ValueError('Object type not supported, supported types are ', allowed_objects)
        
        ids_lists_in_object=self.get_ids_from_sugar(sugar_object=sugar_object)

        for id_list in ids_lists_in_object:
            for id in id_list:
                espresso_system.part[id].remove()

        # Update sugar storage of ids

        if sugar_object.sugar_object_type == 'particle':
            for id_list in ids_lists_in_object:
                for id in id_list:
                    self.id_map['particle'][sugar_object.name].remove(id)

        elif sugar_object.sugar_object_type == 'residue':
            self.id_map['residue'].pop(object.name)

        elif sugar_object.sugar_object_type == 'molecule':
            self.id_map['molecule'].pop(sugar_object.name)


        return

    def create_counterions_in_espresso(self, sugar_object, cation, anion, espresso_system):
        """
        Creates cation and anion particles in the espresso_system to counter the charge of the particles contained in object
        Inputs:
        object:(class) particle, residue or molecule/peptide object
        espresso_system: (class)espresso class object with all system variables.
        cation: (class) particle class object, as defined in sugar, with a positive charge
        anion: (class) particle class object, as defined in sugar library, with a negative charge
        """
        
        object_ids=self.get_ids_from_sugar(sugar_object=sugar_object)

        N_pos=0
        N_neg=0

        for id_list in object_ids:
            for id in id_list:
                if espresso_system.part[id].q > 0:
                    N_pos+=1
                elif espresso_system.part[id].q < 0:
                    N_neg+=1

        if (N_pos % abs(anion.q) == 0):

            N_anion=int(N_pos/abs(anion.q))

        else:

            raise ValueError('The number of positive charges in the sugar_object must be divisible by the  charge of the anion')

        if (N_neg % abs(cation.q) == 0):

            N_cation=int(N_neg/cation.q)

        else:

            raise ValueError('The number of negative charges in the sugar_object must be divisible by the  charge of the cation')
        
        for _ in range(N_anion):
            self.create_particle_in_espresso(particle=anion,espresso_system=espresso_system)

        for _ in range(N_cation):
            self.create_particle_in_espresso(particle=cation,espresso_system=espresso_system)

        print('Created ', N_cation, ' cations and ', N_anion, 'anions as counterions')

        return N_pos, N_neg

    def get_net_charge_from_espresso(self, espresso_system, sugar_object):
        """ 
        Calculates the net charge of all the objects in espresso compatibles with sugar_object

        Inputs:
        sugar_object:(class) particle, residue or molecule/peptide object
        espresso_system: (class)espresso class object with all system variables.

        Returns:
        Z_list: (list) list with the net charge of all the objects in espresso compatibles with sugar_object
        """

        ids_lists_in_object=self.get_ids_from_sugar(sugar_object=sugar_object)
        Z_list=[]
        
        for id_list in ids_lists_in_object:
            z_one_object=0
            for id in id_list:
                z_one_object+=espresso_system.part[id].q
            Z_list.append(z_one_object)

        return Z_list
        
    def get_charge_in_residues(self, espresso_system, molecule):
        """
        Returns a list with the charge in each residue of molecule stored in  dictionaries
        Inputs:
        sugar_object:(class) particle, residue or molecule/peptide object
        espresso_system: (class) espresso class object with all system variables.
        Returns:
        charge_in_residues: (list) list with the charge in each residue of molecule stored in dictionaries
        """

        charge_in_residues=[]

        for molecule in self.id_map['molecule'][molecule.name]:
            molecule_list=[]
            for residue_dict in molecule:
                charge_dict={}
                for key in residue_dict.keys():
                    charge_key=[]                      
                    for id in residue_dict[key]:                         
                        charge_key.append(espresso_system.part[id].q)
                    
                    if 'side-' in key:
                        new_key=key.replace('side-', '')
                    elif 'central-' in key:
                        new_key=key.replace('central-', '')
                    else:
                        new_key=key
                    charge_dict[new_key]=charge_key
                molecule_list.append(charge_dict)
            charge_in_residues.append(molecule_list)
        
        return charge_in_residues

    def setup_lj_interactions_in_espresso(self, espresso_system, sigma=None, cutoff=None, shift='auto', use_default_values=False, combining_rule='Lorentz-Berthelot'):
        """
        Setup lennard-jones interactions for all particles types present in the  espresso_system. 

        Inputs:
        espresso_system: (class) espresso class object with all system variables.
        See section 6.1.2 of the user guide from espresso for a complete description of the rest of parameters
        """
    
        from itertools import combinations_with_replacement

        implemented_combinatiorial_rules=['Lorentz-Berthelot']

        if combining_rule not in implemented_combinatiorial_rules:
            raise ValueError('In the current version of sugar, the only combinatorial rules implemented are ', implemented_combinatiorial_rules)

        particle_name_list=self.stored_sugar_objects['particle'].keys()

        if sigma is None:

            sigma=1*self.units('reduced_length')

        if cutoff is None:

            cutoff=2**(1./6.)*self.units('reduced_length')

        # For bookeeping
        if 'LennardJones' not in self.stored_sugar_objects.keys():
            self.stored_sugar_objects['LennardJones']={}

        all_types=[]
        type_dict={}
        for particle_name in particle_name_list:
            particle_type_list=self.get_particle_types(particle=self.stored_sugar_objects['particle'][particle_name])
            all_types+=particle_type_list
            for es_type in particle_type_list:
                type_dict[es_type]=particle_name

        for type_pair in combinations_with_replacement(all_types, 2):

            type1=type_pair[0]
            type2=type_pair[1]
            name1=type_dict[type1]
            name2=type_dict[type2]
 
            particle1=self.stored_sugar_objects['particle'][name1]
            particle2=self.stored_sugar_objects['particle'][name2]

            pair_interaction_defined=True
            interaction_key=frozenset(type_pair)

            if particle1.diameter is not None and particle2.diameter is not None:

                if combining_rule == 'Lorentz-Berthelot':
                    
                    combined_sigma=(particle1.diameter.to("reduced_length")+particle2.diameter.to("reduced_length"))/2.
                    offset=combined_sigma.to("reduced_length")-sigma.to("reduced_length")
                    

            elif use_default_values and 'default' in self.stored_sugar_objects['LennardJones'].keys():

                sigma=self.stored_sugar_objects['LennardJones']['default']['sigma']
                offset=self.stored_sugar_objects['LennardJones']['default']['offset']

            else:
            
                pair_interaction_defined=False                

            if particle1.epsilon is not None and particle2.epsilon is not None:

                if combining_rule == 'Lorentz-Berthelot':

                    epsilon=self.np.sqrt(particle1.epsilon*particle2.epsilon)

            elif use_default_values and 'default' in self.stored_sugar_objects['LennardJones'].keys():

                epsilon=self.stored_sugar_objects['LennardJones']['default']['epsilon']
                
            else:
            
                pair_interaction_defined=False

            if interaction_key in self.stored_sugar_objects['LennardJones'].keys():
                pair_interaction_defined=False
                print('Lennard Jones interaction is already defined for ', interaction_key)
                continue          

            if pair_interaction_defined:
                
                espresso_system.non_bonded_inter[type1, type2].lennard_jones.set_params(epsilon = epsilon.to('reduced_energy').magnitude, 
                                                                                sigma = sigma.to('reduced_length').magnitude, 
                                                                                cutoff = cutoff.to('reduced_length').magnitude,
                                                                                offset = offset.to('reduced_length').magnitude, 
                                                                                shift = shift)
                param_dict={'epsilon' : epsilon.to('reduced_energy').magnitude, 
                            "sigma" : sigma.to('reduced_length').magnitude, 
                            'cutoff' : cutoff.to('reduced_length').magnitude,
                            'offset' : offset.to('reduced_length').magnitude, 
                            'shift' : shift}
                self.stored_sugar_objects['LennardJones'][interaction_key]=param_dict
            else:

                print('WARNING missing parameters for Lennard Jones interaction between ', [name1,name2], ', no interaction has been defined')
                
        return

        
    def get_all_stored_types(self):
        """
        Returns a dictionary with all particle types stored in sugar
        Returns:
        type_dict: (dict)  dictionary with all particle types stored and its keys
        """

        particle_name_list=self.stored_sugar_objects['particle'].keys()

        type_dict={}
        for particle_name in particle_name_list:
            particle_type_list=self.get_particle_types(particle=self.stored_sugar_objects['particle'][particle_name])
            for es_type in particle_type_list:
                type_dict[es_type]=particle_name

        return type_dict

    def search_largest_particle_diameter_in_sugar(self):
        """
        Returns the largest particle diameter stored in sugar
        Returns:
        largest_diameter: (float) largest particle diameter in stored in sugar
        """

        diameter_list=[]

        if 'particle' not in self.id_map.keys():
            return 0*self.units('reduced_length')
        
        for particle_name in self.id_map['particle'].keys():
            particle=self.stored_sugar_objects['particle'][particle_name]
            if particle.diameter is not None:
                diameter_list.append(particle.diameter)

        if len(diameter_list) >0:
            largest_diameter=max(diameter_list)
            return largest_diameter
        else:
            return 0*self.units('reduced_length')

    def visualize_espresso_system(self, espresso_system):
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

    def do_snapshot_espresso_system(self,espresso_system, filename):
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
