class sugar_library(object):

    import pint
    units = pint.UnitRegistry()
    import numpy as np
    import random as rn
    import scipy.constants

    # Default values    

    TEMPERATURE = 298.15 * units.K
    PARTICLE_SIZE = 0.355 * units.nm
    SEED=None
    N_A=scipy.constants.Avogadro / units.mol
    Kb=scipy.constants.Boltzmann * units.J / units.K

    # Library output

    filename_parameters = 'simulation_parameters.txt'
    
    # Aminoacid key

    aminoacid_key={"ALA": "A",
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

    def __init__(self):

        # Default definitions of reduced units

        self.units.define(f'reduced_energy = {self.TEMPERATURE * self.Kb}')
        self.units.define(f'reduced_length = {self.PARTICLE_SIZE}')
        self.units.define(f'reduced_charge = 1*e')
        self.print_reduced_units()
        self.kT=self.TEMPERATURE*self.Kb

        # Load parameters
        
        self.param=self.parameters(units=self.units, particle=self.particle)

        # Create molecule object

        self.create_molecule_object()

        # Open the parameters file

        with open(self.filename_parameters, 'w') as par_file:
            par_file.write('')

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
        self.N_A=self.scipy.constants.Avogadro / self.units.mol
        self.Kb=self.scipy.constants.Boltzmann * self.units.J / self.units.K
        self.TEMPERATURE=temperature.to('K').magnitude*self.units.K
        unit_energy=self.TEMPERATURE*self.Kb
        self.units.define(f'reduced_energy = {unit_energy} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')
        
        self.kT=self.TEMPERATURE*self.Kb

        self.print_reduced_units()

    class particle:

        radius=None
        type=None
        q=None
        ids=[]
        name=None
        N=None
        acidity=None
        state=None
        pKa=None
        model=None

    class residue:

        def __init__(self, name):

            self.name=name

        principal_bead=None
        lateral_beads=[]
        bonds=[]
        beads=[]
        ids=[]
        N=None
        model=None

    def create_molecule_object(sugar_self):

        class molecule:

            N = None
            Nm = None
            ids = []
            model = None
            
            def __init__(self, sequence, model=None, param_custom=None, pKa_set=None, pKa_custom=None):


                model_param=None
                model_dict={}
            
                if model is not None:

                    model_list=sugar_self.get_subclasses(sugar_self.param)

                    # Get the model parameters

                    for param_set in model_list:

                        if param_set.name is not None:
                            
                            if model.lower() == param_set.name.lower():

                                model_param=param_set

                    if model_param is None:

                        model_names=sugar_self.get_modelnames(cls=sugar_self)
                        raise ValueError("Unknown model chosen: ", model, ". Valid options are ", model_names)
                    
                    # Get a list of model parameters
                    
                    list_model=sugar_self.get_attributes(model_param)

                    for par in list_model:

                        model_dict[par[0]]=par[1]

                    # Get valid keys for the model

                    model_actors=sugar_self.get_particles(model_param)
                    keys=[]

                    for actor in model_actors:
            
                        keys.append(actor.name)
                        
                else:
                    
                    keys=list(sugar_self.aminoacid_key.values())
                    
                # Load Custom parameters if any

                # Check that the custom parameters are in the proper format:

                names_custom=[]
                
                if  param_custom is not None:

                    if isinstance(param_custom, sugar_self.param.custom_parameters):
                        
                        

                        if model_param is not None:
                            
                        # Overwrite model parameters by custom ones

                            list_custom=sugar_self.get_attributes(param_custom)
                            custom_dict={}
                            
                            for par in list_custom:

                                custom_dict[par[0]]=par[1]

                            for atr in custom_dict.keys():

                                if custom_dict[atr] is not None:

                                    if isinstance(custom_dict[atr], sugar_self.param.particle):
                                            
                                        if atr in model_dict.keys():

                                            part_custom=sugar_self.get_attributes(custom_dict[atr])

                                            for par in part_custom:

                                                if par[1] is not None:

                                                    setattr(model_dict[atr],par[0],par[1])
                            
                                        else:
                                            
                                            setattr(model_param,atr,custom_dict[atr])

                                    if isinstance(custom_dict[atr], dict):

                                        if atr in model_dict.keys():
                                            
                                            atr_model_dict=model_dict[atr]
                                            atr_custom_dict=custom_dict[atr]

                                            for key in atr_custom_dict.keys():

                                                atr_model_dict[key]=atr_custom_dict[key]

                                        else:

                                            setattr(model_param,atr,custom_dict[atr])

                                    else:
                                        
                                        setattr(model_param,atr,custom_dict[atr])
                                        
                        else:

                            model_param=param_custom

                    else:

                        raise ValueError("Unrecognized format for the custom parameters. Please, use the library function 'setup_custom_parameters' to define your custom parameters")
                
                if model_param is None:

                    raise ValueError("A model is needed to construct a molecule, please either chose one of the default ones or provide a custom one")

                model_lj=sugar_self.get_lj(model_param)
                
                if (len(model_lj) == 0): # By default, WCA lennard jones interactions are included in the model
                    
                    setattr(model_param,sugar_self.param.default.lj_WCA.name, sugar_self.param.default.lj_WCA)

                else:

                    WCA_in_model=False

                    for lj_par in model_lj:

                        if lj_par.name == "WCA":

                            WCA_in_model=True

                    if not WCA_in_model:

                        setattr(model_param,sugar_self.param.default.lj_WCA.name,sugar_self.param.default.lj_WCA)

                # Store the model in the molecule object

                self.model=model_param

                if pKa_custom is not None:

                    for key in pKa_custom.keys():

                        keys.append(key)
                
                if (model == '1beadpeptide' or model == '2beadpeptide'):

                    clean_sequence=sugar_self.sequence_parser(sequence, keys)

                else:

                    clean_sequence=sequence
                
                # Import the pKa values

                param_pKa_set=sugar_self.get_subclasses(sugar_self.param.pKa_set)
                names_pKa_set=[]

                for sets in param_pKa_set:
                    
                    names_pKa_set.append(sets.name)

                if isinstance(pKa_set, str):

                    for sets in param_pKa_set:

                        if (pKa_set.lower() == sets.name.lower()):

                            pKa_set=sets.pKa
                            break

                    if not isinstance(pKa_set,dict):

                        raise ValueError("Unknown key provided for a pKa set, valid options are " , names_pKa_set)
                    
                elif pKa_set is not None:
                

                    raise ValueError("The desired pKa set must be given as a string, valid options are ", names_pKa_set)

                if pKa_custom is not None:

                    if isinstance(pKa_custom,dict):

                        for key in pKa_custom.keys():

                            pKa_set[key]=pKa_custom[key]

                    else:

                        raise ValueError("The custom pKa-values must be provided in a dictionary")

                # Create an object residue per each residue in sequence

                self.sequence=clean_sequence
                self.residues=[]
                self.Nm=len(clean_sequence)

                for res in clean_sequence:

                    monomer=sugar_self.residue(name=res)
                    self.residues.append([monomer])

                # Set-up the model

                if model_param is None:

                    # If no model nor custom parameters are given, create one bead particle per residue

                    for r_chain in self.residues:
                        
                        for res in r_chain:
                            bead=sugar_self.particle()
                            bead.name=res.name
                            bead_list=[]
                        
                        # If the residue has a pKa value list in the pKa_set put in the bead
                            if pKa_set is not None:
                                if res.name in pKa_set.keys():

                                    bead.pKa=pKa_set[res.name]

                                bead_list.append(bead)
                                res.beads=bead_list
                
                else:

                    param_part=sugar_self.get_particles(model_param)
                    param_bonds=sugar_self.get_bonds(model_param)
                    model_keys=[]
                    
                    for p_set in param_part:
                        
                        model_keys.append(p_set.name)
                    
                    for r_chain in self.residues:

                        for res in r_chain:

                            bead_list=[]
                            lateral_beads=[]

                        # Bead of the principal chain

                            p_bead_name=None
                            if model_param.principal_chain is not None:

                                if res.name in model_param.principal_chain.keys():

                                    p_bead_name=model_param.principal_chain[res.name]

                                elif "default" in model_param.principal_chain.keys():

                                    p_bead_name=model_param.principal_chain["default"]
                            
                            if model_param.principal_chain is None or p_bead_name ==  'sequence':
      
                                bead_unparameterized=True

                                for p_set in param_part:

                                    if (p_set.name == res.name):

                                        bead=p_set
                                        bead_unparameterized=False
                                            
                                if bead_unparameterized:

                                    bead=sugar_self.particle()
                                    bead.name=res.name

                                if pKa_set is not None:
                                    if res.name in pKa_set.keys():

                                        bead.pKa=pKa_set[res.name]

                                            
                            elif p_bead_name in model_keys:
                                
                                for p_set in param_part:

                                    if (p_set.name == p_bead_name):

                                        bead=p_set

                                        if p_bead_name in pKa_set.keys():

                                            bead.pKa=pKa_set[res.name]
                            
                            else:
                                print(p_bead_name,model_keys)
                                raise ValueError("Unknown key for the principal chain: ", p_bead_name)
                                
                            res.principal_bead=bead.name
                            bead_list.append(bead)

                        # Beads on the lateral chain

                            if model_param.side_chain is not None:

                                side_dict=model_param.side_chain.copy()

                                if res.name  in side_dict.keys():

                                    chains=side_dict[res.name]

                                else:

                                    if "default"  in side_dict.keys():
                                    
                                        chains=[side_dict["default"]]

                                    else:

                                        raise ValueError("Unknown lateral chain for :", res.name)
                                    
                            
                                for chain in chains:
                                    
                                    for part in chain:
                                        
                                        if part == "sequence":
                                    
                                            name= res.name
                                            for p_set in param_part:

                                                if (p_set.name == name):

                                                    bead=p_set
                                            if pKa_set is not None:
                                                if name in pKa_set.keys():

                                                    bead.pKa=pKa_set[name]

                                            bead_list.append(bead)
                                            lateral_beads.append(bead.name)
                                        

                                        else:     
                                                            
                                            if part in model_keys:
                                    
                                                name=part

                                                for p_set in param_part:

                                                    if (p_set.name == name):

                                                        bead=p_set
                                                if pKa_set is not None:
                                                    if name in pKa_set.keys():

                                                        bead.pKa=pKa_set[name]
                                
                                            else:

                                                raise ValueError("Unknown key for the side chain: ", part)

                                            bead_list.append(bead)
                                            lateral_beads.append(bead.name)
                        
                            res.lateral_beads=lateral_beads
                            res.beads=[bead_list]
                    
                        res_bond_list=[]
                                            
                        # All residues contain at least the default bond

                        for bond in param_bonds:
                                    
                            if bond.actors[0] == "default":

                                if bond not in res_bond_list:
                                    
                                    res_bond_list.append(bond)

                        for chain in lateral_beads:
                            
                            bead_name=chain[0]
                            actors=[res.principal_bead,bead_name]

                        # Check if there is a specific parametrization in model for the bond between  res.principal_bead and bead_name

                            for bond in param_bonds:

                                if actors == bond.actors or actors[::-1] == bond.actors:
                                    
                                    bond_assigned=True
                                    if bond not in res_bond_list:
                                        res_bond_list.append(bond)
                                    break

                        res.bonds=res_bond_list
                    
        sugar_self.molecule = molecule
                
    def create_custom_model(self, principal_chain=None, side_chain=None, custom_particles=None):
        '''
        Helps the user to setup custom parameters for a model or a full custom model

        Inputs:
        beads_per_residue(int) = number of beads per residue
        principal_chain(string) = key of the particle type in the main chain of the molecule. 
        side_chain(list) = keys of the particles in the side chains of the molecule
        For principal_chain and side_chain the keyword "sequence"
        custom_particles (dict) = custom parameters for existing/new particle types. The dictionary must have the following structure:

            dict = {'particle1_name': {'propertie1_name': propertie1_value,
                                        'propertie2_name': propertie2_value
                                        }
                    'particle2_name': {'propertie1_name': propertie1_value,
                                        'propertie2_name': propertie2_value
                                        }
                                        }
        custom_bonds (dict) = custom parameters for existing/new bond types. The dictionary must have the following structure:

            dict = {'bond1_name': {'actors': [bonded_name1, bonded_name2] ,
                                'type': bond_type,
                                'bondl': value,
                                'k': value
                                        }
                    {'bond2_name': {'actors': [bonded_name1, bonded_name2] ,
                                'type': bond_type,
                                'bondl': value,
                                'k': value
                                        }

        '''
        
        custom_param=self.param.custom_parameters()
        
        if principal_chain is not None:

            if isinstance(principal_chain, dict):

                custom_param.principal_chain=principal_chain

            else:

                raise ValueError("principal_chain must contain a string with the key of the beads desired in the principal chain of the molecule. Provided:", principal_chain)


        if side_chain is not None:

            if isinstance(side_chain, dict):

                custom_param.side_chain=side_chain

            else:

                raise ValueError("side_chain must contain a list with the keys of the beads desired in the side_chain of the molecule. Provided:", side_chain)

        if custom_particles is not None:

            if isinstance(custom_particles, dict):

                for key in custom_particles.keys():

                    if not isinstance(key,str):

                        raise ValueError("Particle keys must be given as strings. Key given: ", key)
                    
                    elif not isinstance(custom_particles[key], dict):

                        raise ValueError("Particle properties must be given as a dictionary. Properties given: ", custom_particles[key])

                    particle_dict=custom_particles[key]
                    custom_part=self.particle()
                    properties=self.get_attributes(custom_part)
                    properties_names=[]

                    for propertie in properties:

                        properties_names.append(propertie[0])

                    custom_part.name=key
                    
                    for propertie in particle_dict.keys():

                        if (propertie in properties_names):

                            if propertie == 'radius':

                                if hasattr(particle_dict[propertie],'check'):
                                
                                    if particle_dict[propertie].check('[length]'):

                                        setattr(custom_part, propertie, particle_dict[propertie])

                                    else:

                                        raise ValueError("The bead radius must have units of length. Given: ", particle_dict[propertie]) 

                                else:
                                    
                                    raise ValueError("Bead radius must be given as a pint quantity. Given:", particle_dict[propertie])
                            
                            else:    

                                setattr(custom_part, propertie, particle_dict[propertie])

                        else:

                            raise ValueError("Unknown particle propertie given: ", propertie, " valid properties names are ", properties_names)

                    setattr(custom_param, key, custom_part)

            else:

                raise ValueError("Custom particle properties must be given in a dictionary. Please refer to the library documentation for more information on how to set-up custom particle properties. Input given: ", custom_particles)


        return custom_param

    def sequence_parser(self, sequence, keys=aminoacid_key.values()):
        '''
        Reads the input residue sequence and
        1) Checks that all residues are in the parameters key
        2) Transforms the aminoacids on the sequence from three letter to one letter format

        Input:
        sequence: (string or list) aminoacid sequence

        Output:
        clean_sequence: (string) 

        '''

        clean_sequence=[]
        
        if isinstance(sequence, str):
                    
            if (sequence.find("-") != -1):

                splited_sequence=sequence.split("-")

                for residue in splited_sequence:

                    if len(residue) == 1:

                        if residue in keys:

                            residue_ok=residue

                        else:

                            if residue.upper() in keys:

                                residue_ok=residue.upper()

                            else:

                                raise ValueError("Unknown one letter code for a residue given: ", residue, " please review the input sequence")

                        clean_sequence.append(residue_ok)
                    
                    else:

                        if residue in keys:

                            clean_sequence.append(residue)

                        else:

                            if (residue.upper() in self.aminoacid_key()):

                                clean_sequence.append(self.aminoacid_key[residue.upper()])

                            else:

                                raise ValueError("Unknown  code for a residue: ", residue, " please review the input sequence")

            else:

                for residue in sequence:

                    if residue in keys:

                        residue_ok=residue

                    else:

                        if residue.upper() in keys:

                            residue_ok=residue.upper()

                        else:

                            raise ValueError("Unknown one letter code for a residue: ", residue, " please review the input sequence")

                    clean_sequence.append(residue_ok)

        if isinstance(sequence, list):

            for residue in sequence:
                
                    if residue in keys:

                        residue_ok=residue

                    else:

                        if residue.upper() in keys:

                            residue_ok=residue.upper()

                        elif (residue.upper() in self.aminoacid_key.keys()):

                            clean_sequence.append(self.aminoacid_key[residue.upper()])

                        else:

                            raise ValueError("Unknown code for a residue: ", residue, " please review the input sequence")

                    clean_sequence.append(residue_ok)

        return clean_sequence

    def write_parameters(self, mol):

        if isinstance(mol,self.molecule):

            print('molecule parameters')
            for atr in self.get_attributes(mol):
                print(atr)

            for chain in mol.residues:
                
                for res in chain:

                    print("\t residue ", res.name, 'parameters')
            
                    for atr in self.get_attributes(res):

                        print('\t' ,atr)


                    for chain_bead in res.beads:

                        for bead in chain_bead:

                            print("\t \t Particle", bead.name)
                
                            for atr in self.get_attributes(bead):

                                print('\t \t ', atr)

        if isinstance(mol,self.residue):
            
            print("residue ", mol.name, 'parameters')

            for atr in self.get_attributes(mol):
                
                print(atr)

            for chain_bead in res.beads:

                for bead in chain_bead:

                    print("\t \t Particle", bead.name)
                
                    for atr in self.get_attributes(bead):

                        print('\t \t ', atr)

        if isinstance(mol, self.particle): 

            print("Particle", mol.name)
                
            for atr in self.get_attributes(mol):

                print(atr)   

    def create_particle(self, particle, system, position=None, state=None, id=None):
        '''
        Creates a bead particle in the espresso system with the properties stored in bead.
        part = instance of a particle object from sugar library
        system = instance of a system object from espressomd library
        position = array of particle positions of size part.N
        state = for particles with more than one possible state, defines in which state are created
        id = array of particle ids of size part.N
        '''

        if not isinstance(particle,self.particle):

            raise ValueError("Please provide an instance of a particle object from sugar library. Given:", part)

        if position is None:
        
        # By default, the beads are placed in random positions
            
            pos=self.np.random.random((1, 3))[0] *self.np.copy(system.box_l)

        else:

            if (len(position) == particle.N):

                pos=position[0]
            
            else:
                
                raise ValueError("Please provide ", particle.N ," position arrays")
                
        if state is not None:
            
            if isinstance(particle.q,dict):
                
                if state in particle.q.keys():

                    particle.state=state

                else:

                    raise ValueError("Unknown state for bead: ", particle.name, " please review the input state:", state, ' valid keys are', particle.q.keys())
        else:

            if isinstance(particle.type, dict):
                
                # Inicialice beads with more than one state in a random state

                state=self.rn.choice(list(particle.type.keys()))
                particle.state=state

        
        if isinstance(particle.q, int) or isinstance(particle.q, float) :
            
            q=particle.q
        
        elif isinstance(particle.q, dict):

            q=particle.q[state]

        else:

            raise ValueError("Unvalid charge for bead: ", particle.name, " given ", particle.q)

        if isinstance(particle.type, int):
            
            type=particle.type

        elif isinstance(particle.type, dict):

            type=particle.type[state]

        else:

            raise ValueError("Unvalid type for bead: ", particle.name, " given ", particle.type)

        
        if id is None:

            if len(system.part[:].id) == 0:
                
                bead_id=0

            else:
                
                bead_id=max(system.part[:].id)+1
            

        else:
            if isinstance(id,list):
                
                if len(id) == particle.N:

                    bead_id=id[0]

                else:

                    raise ValueError('Please provide ', particle.N, 'ids. Provided:', id)

            else:

                raise ValueError('particle ids must be provided as a list. Given:', id)

        bead_id_lst=[]

        # By default only one bead is created

        if not isinstance(particle.N, int): 

            particle.N=1

        for bead_i in range(particle.N):
                
            system.part.add(id=[bead_id], pos=[pos], type=[type], q=[q])
            bead_id_lst.append(bead_id)
            
            if id is None:

                bead_id+=1
            
            else:

                id.remove(bead_id)

                if len(id) != 0:
                    
                    bead_id=id[0]
            
            if position is not None:
                    
                position.remove(pos)
                    
                if len(position) != 0:
                        
                    pos=position[0]
                
            else:

                pos=self.np.random.random((1, 3))[0] *self.np.copy(system.box_l)

        
        particle.ids=bead_id_lst   

        return

    def create_residue(self, res, system, position=None, state=None, id=None):
        """
        Creates a residue in the espresso given with the properties stored in res
        res: instance of a residue object as defined in sugar library
        system: instance of a system object of espresso library
        position: array of position of each of the particles of the residue. 
        The positions are asigned in the same order than the particles stored in residue (first the principal bead and then the lateral ones)
        state = for particles with more than one possible state, defines in which state are created
        id = array of particle ids of size part.N
        """

        bead_list=[]
        residue_ids=[]

        if id is not None:

            if isinstance(id,list):

                bead_id=[id[0]]

            else:

                raise ValueError("Please provide the inicialization state as an array, provided:", id)

        else:

            bead_id=None

        if state is not None:

            if isinstance(state,list):

                bead_state=state[0]

            else:

                raise ValueError("Please provide the inicialization state as an array, provided:", state)

        else:

            bead_state=None

        if position is None:

            pos=position

        else:

            if isinstance(position,list):

                pos=[position[0]]

            else:

                raise ValueError("Please provide the desired position for the particle as a list, given:", position)


        # create the principal bead
        
        if not isinstance(res.N, int): 

            res.N=1

        for n_res in range(res.N):
        
            for b_chain in res.beads:
                
                for bead in b_chain:

                    bead.N=1
                    if bead.name == res.principal_bead:
                
                # Create an individual particle object

                        new_bead=self.particle()
                        bd_atrb=self.get_attributes(bead)

                        for atr in bd_atrb:

                            setattr(new_bead,atr[0],atr[1])
                    
                        new_bead.N=1
                
                        self.create_particle(particle=new_bead, system=system, position=pos, state=bead_state, id=bead_id)                        

                        if state is not None:
                    
                            state.remove(bead_state)

                            if len(state) > 0:

                                bead_state=state[0]
                    
                            else:

                                state=None

                
                        if id is not None:
                    
                            id.remove(new_bead.ids[0])

                            if len(id) > 0:

                                bead_id=[id[0]]
                    
                            else:

                                id=None

                        if position is not None:
                    
                            position.remove(position[0])

                        central_bead=new_bead
                        central_id=new_bead.ids
                        central_name=new_bead.name
                        central_pos=system.part[central_id].pos
                        central_pos=central_pos[0]
                        bead_list.append(new_bead)
                        residue_ids.append(central_id[0])
                        break
        
                    # create the lateral beads

            for chain in res.lateral_beads:

                start_position=True
                
                for bead_name in chain:
                
                    if start_position:

                        for b_chain in res.beads:
                            
                            for bead in b_chain:

                                if bead.name == bead_name:

                        # Create an individual particle object
                                    new_bead=self.particle()
                                    bd_atrb=self.get_attributes(bead)

                                    central_bond=self.search_bond(bead1=central_bead,bead2=bead, res=res)

                                    bond_vector=self.generate_trialvectors(central_bond.bondl.to('nm').magnitude)

                                    if position is not None and len(position) > 0:
                    
                                        prebead_position=position[0]

                                    else:
                                        
                                        prebead_position=central_pos+bond_vector

                                    for atr in bd_atrb:

                                        setattr(new_bead,atr[0],atr[1])

                                    new_bead.N=1
                            
                                    self.create_particle(particle=new_bead, system=system, position=[prebead_position], state=bead_state, id=bead_id)

                                    if state is not None:
                                
                                        state.remove(bead_state)

                                        if len(state) > 0:

                                            bead_state=state[0]
                    
                                        else:

                                            state=None
                            
                                    if id is not None:
                                
                                        id.remove(new_bead.ids[0])

                                        if len(id) > 0:

                                            bead_id=[id[0]]
                    
                                        else:

                                            id=None

                                    if position is not None and len(position) > 0:
                    
                                        position.remove(position[0])                        

                                # Create the bond   
                                    
                                    self.create_bond(system,bond=central_bond,id1=new_bead.ids[0], id2=central_id[0])

                                    prebead=new_bead
                                    prebead_id=new_bead.ids[0]
                                    prebead_name=new_bead.name
                                    bead_list.append(new_bead)
                                    residue_ids.append(new_bead.ids[0])
                                    break

                            start_position=False

                    else:

                    
                        if position is not None and len(position) > 0:
                    
                            prebead_position=position[0]

                        else:
                        
                            prebead_position=prebead_position+bond_vector
                    
                        for b_chain in res.beads:
                        
                            for bead in b_chain:

                                if bead.name == bead_name:

                                # Create an individual particle object
                            
                                    new_bead=self.particle()
                                    bd_atrb=self.get_attributes(bead)


                                    for atr in bd_atrb:

                                        setattr(new_bead,atr[0],atr[1])
                            
                                    new_bead.N=1

                                    lateral_bond=self.search_bond(bead1=prebead, bead2=new_bead, res=res)
                            
                                    self.create_particle(particle=new_bead, system=system, position=[prebead_position], state=bead_state, id=bead_id)

                                    if state is not None:
                                
                                        state.remove(bead_state)
                                
                                        if len(state) > 0:

                                            bead_state=state[0]
                    
                                        else:

                                            state=None
                            
                                    if id is not None:
                                
                                        id.remove(new_bead.ids[0])
                                
                                        if len(id) > 0:

                                            bead_id=[id[0]]
                    
                                        else:

                                            id=None                        
                            
                                    if position is not None and len(position) > 0:
                    
                                        position.remove(position[0])


                                    self.create_bond(system,bond=lateral_bond,id1=prebead_id,id2=new_bead.ids[0])

                                    prebead=new_bead
                                    prebead_name=new_bead.name
                                    prebead_id=new_bead.ids[0]
                                    bead_list.append(new_bead)
                                    residue_ids.append(new_bead.ids[0])

                                    break

        # Store the data in res

            res.ids=[residue_ids]
            res.beads=[bead_list]    

        return

    def create_molecule(self, mol,system):
        '''
        Creates a molecules in the espresso given with the properties stored in mol
        mol: instance of a molecule object as defined in sugar library
        system: instance of a system object of espresso library
        '''
        
        if mol.N is None:

            mol.N=1

        elif not isinstance(mol.N,int):

            raise ValueError("The number of molecules must be an integer number, given: ", mol.N)
        
        print('Parameters used to create ' + ''.join(mol.sequence) +' stored in' + self.filename_parameters)
        with open(self.filename_parameters, 'a') as par_file:
            par_file.write('\n Created molecule ' + ''.join(mol.sequence)+ ' with bonds:\n')

        for mol_i in range(mol.N):

            first_res_inexistent=True
            bond_vector=self.generate_trialvectors(1)
            id_list=[]

            for r_chain in mol.residues:

                for res in r_chain:

                    if (first_res_inexistent):

                        self.create_residue(res, system)
                        first_res_inexistent=False
                        pre_backbone_bead=res.beads[0]
                        pre_backbone_bead=pre_backbone_bead[0]
                        pre_bead_id=pre_backbone_bead.ids[0]
                        pre_bead_pos=system.part[pre_bead_id].pos
                        pre_bead_name=pre_backbone_bead.name
                        ids_res=res.ids[0]
                        id_list+=ids_res

                    else:

                        new_backbone_bead=res.beads[0]
                        new_backbone_bead=new_backbone_bead[0]
                        new_bead_name=new_backbone_bead.name

                        actors=[pre_bead_name,new_bead_name]

                        backbone_bond=self.search_bond(bead1=pre_backbone_bead, bead2=new_backbone_bead,res=res)
        
                        new_bead_pos=pre_bead_pos+bond_vector*backbone_bond.bondl.to('nm').magnitude
                        self.create_residue(res, system, position=[new_bead_pos])
                        
                        # Create the harmonic bond for the principal chain

                        ids_res=res.ids[0]
                        id_ppal_bead=ids_res[0]
                        
                        self.create_bond(system,bond=backbone_bond,id1=pre_bead_id,id2=id_ppal_bead)

                        # Update lists and variables
                        
                        id_list+=ids_res
                        pre_bead_id=id_ppal_bead
                        pre_backbone_bead=new_backbone_bead
                        pre_bead_name=new_bead_name
                        pre_bead_pos=new_bead_pos
            
            mol.ids.append(id_list)

        return 

    def count_titrable_groups(self, mol):
        """
        Counts the number of titrable groups in the protein object

        Input:
        mol: instance of a molecule class object as defined in sugar library.

        Output:
        N_ti: (int) number of titrable groups in the sequence
        """

        N_ti=0
        
        for chain in mol.residues:
                
            for res in chain:

                for chain_bead in res.beads:

                    for bead in chain_bead:

                        if bead.pKa is not None:

                            N_ti+=1

        N_ti*=mol.N

        return N_ti

    def setup_acidbase_reactions(self, mol, counter_ion, method='constant_pH', exclusion_radius=None):
        """
        Set up the Acid/Base reactions for acidic/basidic residues in mol. The reaction steps are done following the constant pH ensamble procedure. 

        Inputs:
        
        mol: molecule class object as defined in sugar library
        cation: particle class object as defined in sugar library

        Output:
        RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble

        """

        import espressomd

        valid_method_list=['constant_pH']
        
        if method in valid_method_list:

            if exclusion_radius is None:

                exclusion_radius= 1*self.units('reduced_length')

            if self.SEED is None:

                self.create_random_seed()

            RE = espressomd.reaction_ensemble.ConstantpHEnsemble(temperature=self.kT.to('reduced_energy').magnitude, exclusion_radius=exclusion_radius.magnitude, seed=self.SEED)

        else:

            raise ValueError('Only the implementation for the constant pH ensamble is available, method =', method)

        reaction_absent={}

        if isinstance(mol,self.molecule):

            for chain in mol.residues:
                
                for res in chain:

                    for chain_bead in res.beads:

                        for bead in chain_bead:

                            if bead.pKa is not None and  bead.acidity in ['acid','basic']:

                                if bead.name not in reaction_absent.keys():

                                    reaction_absent[bead.name]=True

                                if (reaction_absent[bead.name]):

                                    self.setup_bead_acidbase_reaction(RE, bead, counter_ion)
                                    reaction_absent[bead.name]=False

        if isinstance(mol,self.residue):
            
            for chain_bead in res.beads:

                for bead in chain_bead:

                    if bead.pKa is not None and  bead.acidity in ['acid','basic']:

                        if bead.name not in reaction_absent.keys():

                            reaction_absent[bead.name]=True

                        if (reaction_absent[bead.name]):

                            self.setup_bead_acidbase_reaction(RE, bead, counter_ion)
                            reaction_absent[bead.name]=False

        if isinstance(mol, self.particle): 

            if mol.pKa is not None and  mol.acidity in ['acid','basic']:

                self.setup_bead_acidbase_reaction(RE, bead, counter_ion)
                reaction_absent[bead.name]=False

        return RE

    def setup_bead_acidbase_reaction(self, RE, part, counter_ion):
        """
        Set up the Acid/Base reactions for acidic/basidic residues in protein. The reaction steps are done following the constant pH ensamble procedure. 

        Inputs:
        RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
        mol: particle/residue/molecule class object as defined in sugar library
        cation: particle class object as defined in sugar library
        """

        if not isinstance(part, self.particle) or not isinstance(counter_ion, self.particle):

            raise ValueError("part and cation must be instances of a particle object from sugar library")

        if part.acidity not in ['acid','basic']:

            print("WARNING, the added particle does not have its acidity defined, please define part.acidity to 'acid' or 'basic'. No reaction has been added")
        
        if  isinstance(part.type, dict):

            if 'protonated' not in part.type.keys() or 'unprotonated' not in part.type.keys():

                raise ValueError('part.type must contain as keys "protonated" and "unprotoanted". Given: ',  part.type.keys())

        else:

            raise ValueError("Particle type must store the tipe of the protonated and unprotonated species so that part.type['protonated'] returns the type of the protonated specie and part.type['unprotonated'] returns the type of the unprotonated specie. Given: ", part.type)


        if part.pKa is not None and  part.acidity in ['acid','basic']:

            if (part.acidity == 'basic') : # Basic residue
                        
                RE.add_reaction(gamma=10**-part.pKa,
                                reactant_types=[part.type["protonated"]],
                                reactant_coefficients=[1],
                                product_types=[part.type["unprotonated"], counter_ion.type],
                                product_coefficients=[1,1],
                                default_charges={part.type["protonated"]: 1,
                                part.type["unprotonated"]: 0,
                                counter_ion.type: 1})

            elif (part.acidity == 'acid') : # Acid residue

                RE.add_reaction(gamma=10**-part.pKa,
                                reactant_types=[part.type["protonated"]],
                                reactant_coefficients=[1],
                                product_types=[part.type["unprotonated"], counter_ion.type],
                                product_coefficients=[1, 1],
                                default_charges={part.type["protonated"]: 0,
                                                part.type["unprotonated"]: -1,
                                                counter_ion.type: 1})

        return

    def setup_lj_interactions(self, mol_list, system):
        """
        Setup lennard-jones interactions for all the molecules in molecule list, using the parametrization defined in its model. 
        If no specific parametrization is included, a purely repulsive LJ is assumed. 

        Inputs:
        molecule_list(list): list of instances of molecule/residue/bead objects, as defined in sugar library
        system: instance of a system object of espresso library

        """

        if not isinstance(mol_list,list):

            raise ValueError("molecule list must be a list of instances of molecule/residue/bead objects, as defined in sugar library. Given: ", mol_list)

        # Search for all possible particle types and for all the lj parameters

        name_dict={}
        radi_dict={}
        lj_list=[]

        for mol in mol_list:

            if isinstance(mol,self.molecule):

                for chain in mol.residues:
                
                    for res in chain:

                        for chain_bead in res.beads:

                            for bead in chain_bead:
                                    
                                if isinstance (bead.type,dict):

                                    for type in bead.type.values():
                                            
                                        if type not in name_dict.keys():

                                            name_dict[type]=bead.name
                                            radi_dict[type]=bead.radius

                                if isinstance (bead.type,int):
                                    
                                    if bead.type not in name_dict.keys():

                                        name_dict[bead.type]=bead.name
                                        radi_dict[bead.type]=bead.radius

            elif isinstance(mol,self.residue):
            
                for chain_bead in res.beads:

                    for bead in chain_bead:
                                
                        if isinstance (bead.type,dict):

                            for type in bead.type.values():
                                            
                                if type not in name_dict.keys():

                                    name_dict[type]=bead.name
                                    radi_dict[type]=bead.radius

                        if isinstance (bead.type,int):

                            if bead.type not in name_dict.keys():

                                name_dict[bead.type]=bead.name
                                radi_dict[bead.type]=bead.radius

            elif isinstance(mol, self.particle): 
                            
                if isinstance (bead.type,dict):

                    for type in bead.type.values():
                                            
                        if type not in name_dict.keys():

                            name_dict[type]=bead.name
                            radi_dict[type]=bead.radius

                if isinstance(bead.type,int):
                            
                    if bead.type not in name_dict.keys():

                        name_dict[bead.type]=bead.name
                        radi_dict[bead.type]=bead.radius

            else:

                raise ValueError("molecule list must contain a list of instances of molecule/residue/bead objects, as defined in sugar library. Given: ", mol)

        # By default, Lennard-Jones parameters with WCA potential are assumed

        lj_WCA=self.param.default.lj_WCA

        if mol.model is not None:

            lj_parameters=self.get_lj(mol.model)

            for lj_param in lj_parameters:

                if lj_param not in lj_list:

                    lj_list.append(lj_param)

        type_list=[]

        for type in name_dict.keys():

            if type not in type_list:

                type_list.append(type)

        print('\n Parameters used for the Lennard Jones potential stored in ' + self.filename_parameters)
        with open(self.filename_parameters, 'a') as par_file:
            par_file.write('\n \n ***** Lennard Jones parameters ***** \n')


            for type1 in type_list:

                index_type=type_list.index(type1)

                for type2 in type_list[index_type:]:

                    name1=name_dict[type1]
                    name2=name_dict[type2]

                    radius1=radi_dict[type1]
                    radius2=radi_dict[type2]

                    specific_lj=False

                    for lj_param in lj_list:

                        if name1 in lj_param.actors and name2 in lj_param.actors:

                            specific_lj=True
                            
                            if lj_param.sigma is not None:                    
                                
                                sigma=lj_param.sigma

                            else:

                                if radius1 is None:
                                    
                                    radius1=lj_WCA.sigma/2.

                                if radius2 is None:

                                    radius2=lj_WCA.sigma/2.

                                sigma=(radius1.to('nm').magnitude+radius2.to('nm').magnitude)*self.units.nm

                            epsilon=lj_param.epsilon
                            cutoff=lj_param.cutoff
                            shift=lj_param.shift
                            name=lj_param.name

                    if not specific_lj:

                        if radius1 is None:

                            radius1=lj_WCA.sigma/2

                        if radius2 is None:

                            radius2=lj_WCA.sigma/2

                        sigma=radius1 + radius2
                        epsilon=lj_WCA.epsilon
                        cutoff=sigma.to('reduced_length')*2**(1./6.)
                        shift=lj_WCA.shift
                        name=lj_WCA.name

                    pair_lj=self.param.custom_lj()
                    pair_lj.sigma=sigma
                    pair_lj.epsilon=epsilon
                    pair_lj.cutoff=cutoff
                    pair_lj.shift=shift
                    pair_lj.name=name

                    lj_parameters=self.get_attributes(pair_lj)
                    lj_parameters.remove(('actors', None))

                    par_file.write("\n The LJ interactions between "+name1+ '(type = ' + str(type1) + ') and ' + str(name2) + '(type = ' + str(type2) + ") is created using the parameters" + str(lj_parameters))
                    self.setup_lj_pair(type1,type2,pair_lj, system)
                
        return

    def setup_lj_pair(self,type1,type2,lj_param, system):
        """
        Creates a Lennard-Jones interaction between particle of types 'type1' and 'type2' with the parameters contained in lj_param

        Inputs:

        type1(int): first particle type
        type2(int): second particle type 
        lj_param: instance of a lennard-jones object, as defined in sugar library
        system: instance of a system object of espresso library
        """

        if not isinstance(type1, int) or not isinstance(type2, int): 

            raise ValueError("type1 and type2 must be integers. Given type1 =  ", type1, " and type2 = ", type2)

        if not isinstance(lj_param, self.param.custom_lj):

            raise ValueError("lj_param must be an instance of a lennard-jones object, as defined in sugar library. lj_param  =  ", lj_param)


        

        system.non_bonded_inter[type1, type2].lennard_jones.set_params(
                            epsilon = lj_param.epsilon.to('reduced_energy').magnitude,
                            sigma = lj_param.sigma.to('reduced_length').magnitude,
                            cutoff = lj_param.cutoff.to('reduced_length').magnitude,
                            shift = lj_param.shift)

        return

    def calculate_molecule_charge(self,system, mol):
        """ 
        Calculates the charge of the protein and its square

        Inputs:
        system: espresso class object with all system variables
        mol: particle/residue/molecule class object as defined in sugar library


        Outputs:
        Z_prot: (float) total charge of the protein
        Z2: (float) square of Z_prot
        """

        if not isinstance(mol.N,int):

            raise ValueError("The number of objects must be given in mol.N, given:", mol.N)

        if (mol.N == 0):

            raise ValueError("The number of objects must not be zero, given: ", mol.N)

        Z=0

        if isinstance(mol,self.molecule) or isinstance(mol,self.residue):

            for chain in mol.ids:
                
                for id in chain:

                    Z+=system.part[id].q

        if isinstance(mol, self.particle): 

            for id in mol.ids:

                Z+=system.part[id].q

        Z=Z/mol.N
        Z2=Z**2

        return Z, Z2

    def track_ionization(self,system, mol):
        """
        Sets up espresso to track the average number of particles of the types contained in mol
        
        Inputs:
        system: espresso class object with all system variables
        mol: particle/residue/molecule class object as defined in sugar library

        """

        types=[]

        if isinstance(mol,self.molecule) or isinstance(mol,self.residue):

            for chain in mol.ids:
                
                for id in chain:

                    b_type=system.part[id].type
                    
                    if b_type not in types:

                        types.append(b_type)


        elif isinstance(mol, self.particle): 

            for id in mol.ids:

                b_type=system.part[id].type
                    
                if b_type not in types:

                    types.append(b_type)

        else:
            
            raise ValueError("mol must be a particle/residue/molecule class object as defined in sugar library given: ", mol)

        system.setup_type_map(types)

        return

    def calculate_HH(self, mol,pH=None):
        """
        Calculates the ideal Henderson-Hassebach titration curve in the given pH range

        Inputs:
        pH: (list of floats) list of pH values
        mol: particle/residue/molecule class object as defined in sugar library

        Outputs:
        Z_HH: (list of floats) Henderson-Hasselnach prediction of the protein charge in the given pH
        """

        if pH is None:

            pH=self.np.linspace(2,12,50)
        
        elif not isinstance(pH,list):

            raise ValueError("pH must contain a list with the pH-values where the Henderson-Hassebach titration curve will be calculated. Given: ", pH)        

        Z_HH=[]

        for pH_value in pH:
            
            Z=0

            if isinstance(mol,self.molecule):

                for chain in mol.residues:
                
                    for res in chain:

                        for chain_bead in res.beads:

                            for bead in chain_bead:

                                Z+=self.calculate_HH_part(pH_value, bead)


            if isinstance(mol,self.residue):
            
                for chain_bead in res.beads:

                    for bead in chain_bead:

                        Z+=self.calculate_HH_part(pH_value, bead)                    

            if isinstance(mol, self.particle): 

                if mol.pKa is not None and  mol.acidity in ['acid','basic']:

                    Z+=self.calculate_HH_part(pH_value, bead)
                            
            Z_HH.append(Z)

        return Z_HH

    def calculate_HH_part(self, pH, part):
        """
        Calculates the ideal Henderson-Hassebach titration curve of part at one pH-value

        Inputs:
        pH: (float) pH value
        part: particle class object as defined in sugar library

        Outputs:
        z: (float) Henderson-Hasselnach prediction of charge of part at the given pH
        """  

        if not isinstance(part, self.particle): 

            raise ValueError("part must an instance of a particle object as defined in sugar library. Given: ", part)

        if not isinstance(pH,float) or   isinstance(pH,int):

            raise ValueError("pH  must contain a float or integer number with the pH-value, given: ", pH)

        if part.pKa is not None and  part.acidity in ['acid','basic']:

            if part.acidity == 'acid':

                q=-1

            elif part.acidity == 'basic':

                q=+1

            else:

                raise ValueError("Unkown particle acidity, known options are 'acid' or 'basic'. Given:  ", part.acidity)

            z=q/(1+10**(q*(pH-part.pKa)))

        else:
            
            z=0

        return z

    def create_counterions(self,system,mol, cation=None, anion=None):
        """
        Adds one monovalent counter-ion (cation or anion) of opposite charge per each charge in mol

        Inputs:

        system: espresso class object with all system variables.
        mol: molecule/residue/particle class object as defined in this library

        In/Out:
        cation: particle class object as defined in sugar library
        anion: particle class object as defined in sugar library


        """

        # Load the parameters

        if cation is None:

            cation=self.param.small_ions.cation
            cation.N=0

        if anion is None:

            anion=self.param.small_ions.anion
            anion.N=0

        
        # The ids of the counter ions are created to be larger then the ids in the system

        if  len(system.part[:].id) != 0:

            I_id=max(system.part[:].id)

        total_q=0

        # Create a counter-ion for each charged group in the peptide

        cation_id=[]
        anion_id=[]

        if isinstance(mol,self.molecule) or isinstance(mol,self.residue):

            for chain in mol.ids:
                
                for id in chain:

                    q=system.part[id].q
                    
                    if q == +1:

                        I_id+=1
                        I_pos=self.np.random.random((1, 3)) * system.box_l
                        system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[anion.q])
                        anion.N+=1
                        anion_id.append(I_id)
                        total_q+=q+anion.q

                    elif q == -1:

                        I_id+=1
                        I_pos=self.np.random.random((1, 3)) * system.box_l
                        system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[cation.q])
                        cation.N+=1
                        cation_id.append(I_id)
                        total_q+=q+cation.q


        if isinstance(mol, self.particle): 

            for id in mol.ids:

                q=system.part[id].q
                    
                if q == +1:

                    I_id+=1
                    I_pos=self.np.random.random((1, 3)) * system.box_l
                    system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[q])
                    anion.N+=1
                    anion_id.append(I_id)
                    total_q+=q+anion.q

                elif q == -1:

                    I_id+=1
                    I_pos=self.np.random.random((1, 3)) * system.box_l
                    system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[q])
                    cation.N+=1
                    cation_id.append(I_id)
                    total_q+=q+cation.q
        
        cation.ids=cation_id
        anion.ids=anion_id

        if total_q != 0:

            raise ValueError("Subrutine attempted to create a not electroneutral system, please review your peptide setup")


        return cation, anion

    def create_added_salt(self, system, cation, anion, c_salt):
        
        """
        Adds extra cations/anions to the system

        Inputs:

        system: espresso class object with all system variables.
        cation/anion:    particle class object as defined in sugar library
        c_salt: Added salt concentration
        
        Output:
        c_salt_calculated: Calculated added salt concentration added to solution    
        """

        # The current implementation is only valid for 1:1 salt

        if abs(cation.q) != abs(anion.q):

            raise ValueError('The current implementation is only valid for 1:1 salt, charge cation', cation.q ,'charge anion', anion.q) 

        # Calculate the number of ions in the simulation box

        volume=self.units.Quantity(system.volume(), 'reduced_length**3')

        if c_salt.check('[substance] [length]**-3'):
            
            N_ions= int((volume*c_salt.to('mol/reduced_length**3')*self.N_A).magnitude)
            c_salt_calculated=N_ions/(volume*self.N_A)
            
        elif c_salt.check('[length]**-3'):
            
            N_ions= int((volume*c_salt.to('reduced_length**-3')*self.N_A).magnitude)
            c_salt_calculated=N_ions/volume

        else:

            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)

        cation.N=N_ions
        anion.N=N_ions

        self.create_particle(system=system, particle=cation)
        self.create_particle(system=system, particle=anion)

        print('\n Added an added salt concentration of ', c_salt_calculated.to('mol/L'), 'given by ', N_ions, 'cations/anions')
        
        return c_salt_calculated

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

    def get_subclasses(self, cls):

        import inspect
        model_list=[]

        for attribute in cls.__dict__.items():

            name_attribute=attribute[0]

            if name_attribute[0:2] != "__" and inspect.isclass(attribute[1]):

                model_list.append(attribute[1])

        return model_list

    def get_particles(self,cls):

        particle_list=[]

        for attribute in cls.__dict__.items():

            if isinstance(attribute[1], self.param.particle):

                particle_list.append(attribute[1])

        return particle_list

    def get_bonds(self,cls):

        bond_list=[]

        for attribute in cls.__dict__.items():

            if isinstance(attribute[1], self.param.bond):

                bond_list.append(attribute[1])

        return bond_list

    def get_lj(self,cls):

        lj_list=[]

        for attribute in cls.__dict__.items():

            if isinstance(attribute[1], self.param.custom_lj):

                lj_list.append(attribute[1])

        return lj_list

    def get_attributes(self,cls):
        import inspect

        attribute_list=[]

        for attribute in  inspect.getmembers(cls, lambda a:not(inspect.isroutine(a))):

            name_attribute=attribute[0]

            if name_attribute[0:2] != "__":

                attribute_list.append(attribute)

        return attribute_list

    def newAttr(self, attr):

        setattr(self, attr.name, attr)

    def get_modelnames(self, cls):

        model_names=[]
        model_list=self.get_subclasses(cls=cls.param)

        for param_set in model_list:
            
            if param_set.name is not None:
                        
                model_names.append(param_set.name)

        return model_names

    def search_bond(self,bead1,bead2,res):
        """
        Search for a bond in res for joining bead1 and bead2. If there are specific parametrization for such a bond in res, adds the default bond as defined in param.default 

        Inputs:
        bead1: particle object as defined in sugar library
        bead2: particle object as defined in sugar library
        res: residue object as defined in sugar library

        Outputs:
        asigned_bond: bond object as defined in sugar library

        """

        unasigned_bond=True
        actors=[bead1.name, bead2.name]

        for bond in res.bonds:

            if actors == bond.actors or actors[::-1] == bond.actors:

                asigned_bond=bond
                unasigned_bond=False
                break

        if unasigned_bond:

            asigned_bond=self.param.default.default_bond
            asigned_bond.bondl=bead1.radius+bead2.radius

        with open(self.filename_parameters, 'a') as par_file:
            bond_parameters=self.get_attributes(asigned_bond)
            bond_parameters.remove(('actors', ['default']))
            par_file.write("\nThe bond between " +  ' '.join(actors) + " is created using the parameters" + str(bond_parameters))

        return asigned_bond

    def create_bond(self, system,bond,id1,id2):
        """
        Creates a bond between id1 and id2 in system

        Inputs:
        system: instance of the espresso library system class
        bond: instance of a bond object as defined in sugar library
        id1:(int) id number of first particle
        id2: (int) id number of the second particle
        
        """                         

        import espressomd
        from espressomd import interactions

        if not isinstance(bond, self.param.bond):

            raise ValueError("bond should be a sugar bond object, given: ", bond )

        if (bond.type == "harmonic"):

            bond_potential = interactions.HarmonicBond(k=bond.k.to('reduced_energy / reduced_length**2').magnitude, r_0=bond.bondl.to('reduced_length').magnitude)
                                

        else:

            raise ValueError("Unknown bond type", bond.type)

        
        system.bonded_inter.add(bond_potential)
        system.part[id1].add_bond((bond_potential, id2))

        return

    def setup_electrostatic_interactions(self, system, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3):
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
            
        kT=self.TEMPERATURE*self.units.k

        BJERRUM_LENGTH = self.units.e**2 / (4 * self.units.pi * self.units.eps0 * solvent_permittivity * kT)

        print('\n Bjerrum length ', BJERRUM_LENGTH.to('nm'), '=', BJERRUM_LENGTH.to('reduced_length'))

        COULOMB_PREFACTOR=BJERRUM_LENGTH.to('reduced_length') * kT.to('reduced_energy') 
        
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
                system.time_step=0.01
                system.actors.add(coulomb)

                # save the optimal parameters and add them by hand

                p3m_params = coulomb.get_params()
                system.actors.remove(coulomb)
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

        
        system.actors.add(coulomb)
        print("\n Electrostatics successfully added to the system \n")

        return

    def minimize_system_energy(self, system, force_change=1e-2, gamma=1, max_steps=10000, time_step=1e-2, max_displacement=0.1, steps_per_iteration=10):
        """
        Does a steppest descent minimization to relax the system energy

        Inputs:
        system: instance of espressmd system class
        force_change: relative force change ratio threshoold. (Default=1e-2)
        gamma: dammping constant (Default=1 reduced length)
        max_steps: maximum number of steps of the minimization (Default=10000)
        time_step: Time step used for the energy minimization (Default=1e-2)
        max_displacement: maximum particle displacement allowed (Default=0.1 reduced unit)
        steps_per_iteration: Number integration steps per minimization iteration
        """

        print("\n*** Minimazing system energy... ***\n")

        # Set up steepest descent integration

        system.integrator.set_steepest_descent(f_max=0, gamma=gamma, max_displacement=max_displacement)

        # Initialize integrator to obtain initial forces
        
        system.time_step=time_step
        system.integrator.run(0)
        old_force = self.np.max(self.np.linalg.norm(system.part[:].f, axis=1))

        # Minimize the energy

        while system.time / system.time_step < max_steps:
            system.integrator.run(steps_per_iteration)
            force = self.np.max(self.np.linalg.norm(system.part[:].f, axis=1)) # Calculate the norm of the F acting in each particle and select the larger one
            rel_force = self.np.abs((force - old_force) / old_force)
            print(f'rel. force change:{rel_force:.2e}')
            if rel_force < force_change:
                break
            old_force = force

        # Reset the time of the system to 0

        system.time = 0.

        print("\n Minimization finished \n")

        return

    def setup_langevin_dynamics(self,system, time_step=1e-3, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
        """
        Sets up Langevin Dynamics in espressomd.
        system: instance of espressmd system class
        time_step: time s

        """
        

        kT=self.TEMPERATURE*self.units.k

        if self.SEED is None:

            # Take the random seed from the system time
            self.create_random_seed()
            
        system.time_step=time_step
        system.integrator.set_vv()
        system.thermostat.set_langevin(kT=kT.to('reduced_energy').magnitude, gamma=gamma, seed=self.SEED)

        # Optimize the value of skin

        if tune_skin:

            print("\n*** Optimizing skin ... ***")

            if max_skin is None:

                max_skin=system.box_l[0]/2

            system.cell_system.tune_skin(min_skin=min_skin, max_skin=max_skin, tol=tolerance, int_steps=int_steps, adjust_max_skin=adjust_max_skin)

            print("Optimized skin value: ", system.cell_system.skin, "\n")

        return

    def create_random_seed(self):
        """
        Creates the seed for the random number generator from the system hour
        """

        import time 

        SEED=int(time.time())
        print('\n The chosen seed for the random number generator is ', SEED)
        self.SEED=SEED

    class parameters:

        def __init__(self, units, particle):

            self.units=units
            self.particle=particle
            self.create_default_parameter_sets()

        class bond:
            name=None
            actors=None
            type=None
            bondl=None
            k=None

        class custom_lj:
            name=None
            actors=None
            epsilon = None
            sigma = None
            cutoff = None
            shift = None

        class custom_parameters:
            name="custom"
            beads_per_residue=None
            principal_chain=None
            side_chain=None

            def newAttr(self, attr):

                setattr(self, attr.name, attr)

        def create_default_parameter_sets(parameters_self):

            class default:
                
                name=None
                lj_WCA = parameters_self.custom_lj()
                lj_WCA.name='WCA'
                lj_WCA.actors=['default']
                lj_WCA.epsilon= 1 * parameters_self.units('reduced_energy')
                lj_WCA.sigma=1 * parameters_self.units('reduced_length')
                lj_WCA.cutoff=2**(1.0 / 6)*parameters_self.units('reduced_length')
                lj_WCA.shift='auto'

                default_bond=parameters_self.bond()
                default_bond.actors=['default']
                default_bond.type='harmonic'
                default_bond.bondl=1 * parameters_self.units('reduced_length')
                default_bond.k=300 * parameters_self.units('reduced_energy / reduced_length**2')
            
            parameters_self.default=default 
            
            class small_ions:
                name='small_ions'
                beads_per_residue=1
                principal_chain='sequence'
                side_chain=None
                

                cation=parameters_self.particle()
                cation.name="cation"
                cation.type=18
                cation.radius=0.5 * parameters_self.units('reduced_length')
                cation.q=1

                anion=parameters_self.particle()
                anion.name="anion"
                anion.type=19
                anion.radius=0.5 * parameters_self.units('reduced_length')
                anion.q=-1

            parameters_self.small_ions=small_ions

            class pKa_set:
                name=None

                    # Values for the phosphorilated aminoacids J U and Z are always taken from Bienkiewicz & K.J. Lumb, J Biomol NMR 15: 203-206 (1999).

                class Hass:
                    # Values from Hass MA, Mulder FAA. Contemporary NMR Studies of Protein Electrostatics. Annu Rev Biophys. 2015;44:53-75.
                    name="hass"
                    pKa= { "D" : 4.0,
                                        "E" : 4.4,
                                        "H" : 6.8,
                                        "Y" : 9.6,
                                        "K" : 10.4,
                                        "R" : 13.5,
                                        "C" : 8.3,
                                        "J" : 5.96,
                                        "U" : 6.30,
                                        "Z" : 5.96,
                                        "n" : 8.0,
                                        "c" : 3.6
                                }


                class Platzer:
                    # Platzer G, Okon M, McIntosh LP. 2014. pH-dependent random coil 1 H, 13 C, and 15 N chemical shifts of the ionizable amino acids: a guide for protein pK a measurements. J. Biomol. NMR 60:109–29
                    name ='platzer'        
                    pKa = { "D" : 3.86,
                                        "E" : 4.34,
                                        "H" : 6.45,
                                        "Y" : 9.76,
                                        "K" : 10.34,
                                        "R" : 13.9,
                                        "C" : 8.49,
                                        "J" : 5.96,
                                        "U" : 6.30,
                                        "Z" : 5.96,
                                        "n" : 8.23,
                                        "c" : 3.55
                                }

                class CRC:
                    # Values from Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
                    name='crc'
                    pKa = { "D" : 3.65,
                                        "E" : 4.25,
                                        "H" : 6.00,
                                        "Y" : 10.07,
                                        "K" : 10.54,
                                        "R" : 12.48,
                                        "C" : 8.18,
                                        "J" : 5.96,
                                        "U" : 6.30,
                                        "Z" : 5.96,
                                        "n" : 8.0,
                                        "c" : 3.6
                                }

                class Nozaki:
                    # Y. Nozaki and C. Tanford, Methods Enzymol., 1967, 11, 715–734.
                    name = 'nozaki'
                    pKa  = { "D" : 4.00,
                                        "E" : 4.40,
                                        "H" : 6.30,
                                        "Y" : 9.6,
                                        "K" : 10.4,
                                        "R" : 12.0,
                                        "C" : 9.5,
                                        "J" : 5.96,
                                        "U" : 6.30,
                                        "Z" : 5.96,
                                        "n" : 7.5,
                                        "c" : 3.8
                                    }

            parameters_self.pKa_set=pKa_set

            class one_bead_peptide:

                name='1beadpeptide'
                principal_chain={"default": 'sequence'}
                side_chain=None
                
                A=parameters_self.particle()
                A.name="A"
                A.acidity="inert"
                A.radius=0.5 * parameters_self.units('reduced_length')
                A.type=21
                A.q=0

                N=parameters_self.particle()
                N.name="N"
                N.radius=0.5* parameters_self.units('reduced_length')
                N.type=22
                N.q=0
                N.acidity="inert"
                
                Q=parameters_self.particle()
                Q.name="Q"
                Q.radius=0.5 * parameters_self.units('reduced_length')
                Q.type=23
                Q.q=0
                Q.acidity="inert"
                
                G=parameters_self.particle()
                G.name="G"
                G.radius=0.5 * parameters_self.units('reduced_length')
                G.type=24
                G.q=0
                G.acidity="inert"
                
                I=parameters_self.particle()
                I.name="I"
                I.radius=0.5 * parameters_self.units('reduced_length')
                I.type=25
                I.q=0
                I.acidity="inert"
                
                L=parameters_self.particle()
                L.name="L"
                L.radius=0.5 * parameters_self.units('reduced_length')
                L.type=26
                L.q=0
                L.acidity="inert"
                
                M=parameters_self.particle()
                M.name="M"
                M.radius=0.5 * parameters_self.units('reduced_length')
                M.type=27
                M.q=0
                M.acidity="inert"
                
                F=parameters_self.particle()
                F.name="F"
                F.radius=0.5* parameters_self.units('reduced_length')
                F.type=28
                F.q=0
                F.acidity="inert"
                
                P=parameters_self.particle()
                P.name="P"
                P.radius=0.5* parameters_self.units('reduced_length')
                P.type=29
                P.q=0
                P.acidity="inert"
                
                S=parameters_self.particle()
                S.name="S"
                S.radius=0.5* parameters_self.units('reduced_length')
                S.type=30
                S.q=0
                S.acidity="inert"
                
                T=parameters_self.particle()
                T.name="T"
                T.radius=0.5 * parameters_self.units('reduced_length')
                T.type=31
                T.q=0
                T.acidity="inert"
                
                W=parameters_self.particle()
                W.name="W"
                W.radius=0.5 * parameters_self.units('reduced_length')
                W.type=32
                W.q=0
                W.acidity="inert"
                
                Y=parameters_self.particle()
                Y.name="Y"
                Y.radius=0.5* parameters_self.units('reduced_length')
                Y.type=33
                Y.q=0
                Y.acidity="inert"
                
                V=parameters_self.particle()
                V.name="V"
                V.radius=0.5* parameters_self.units('reduced_length')
                V.type=34
                V.q=0
                V.acidity="inert"
                
                K=parameters_self.particle()
                K.name="K"
                K.radius=0.5 * parameters_self.units('reduced_length')
                K.type={"unprotonated": 35, "protonated": 36}
                K.q={"unprotonated": 0, "protonated": 1}
                K.acidity="basic"
                
                H=parameters_self.particle()
                H.name="H"
                H.radius=0.5 * parameters_self.units('reduced_length')
                H.type={"unprotonated": 37, "protonated": 38}
                H.q={"unprotonated": 0, "protonated": 1}
                H.acidity="basic"
                
                C=parameters_self.particle()
                C.name="C"
                C.radius=0.5 * parameters_self.units('reduced_length')
                C.type={"unprotonated": 39, "protonated": 40}
                C.q={"unprotonated": 0, "protonated": 1}
                C.acidity="basic"
                
                R=parameters_self.particle()
                R.name="R"
                R.radius=0.5* parameters_self.units('reduced_length')
                R.type={"unprotonated": 41, "protonated": 42}
                R.q={"unprotonated": 0, "protonated": 1}
                R.acidity="basic"
                
                n=parameters_self.particle()
                n.name="n"
                n.radius=0.5 * parameters_self.units('reduced_length')
                n.type={"unprotonated": 43, "protonated": 44}
                n.q={"unprotonated": 0, "protonated": 1}
                n.acidity="basic"

                c=parameters_self.particle()
                c.name="c"
                c.radius=0.5* parameters_self.units('reduced_length')
                c.type={"unprotonated": 45, "protonated": 46}
                c.q={"protonated": 0, "unprotonated": -1}
                c.acidity="acid"

                D=parameters_self.particle()
                D.name="D"
                D.radius=0.5 * parameters_self.units('reduced_length')
                D.type={"unprotonated": 47, "protonated": 48}
                D.q={"protonated": 0, "unprotonated": -1}
                D.acidity="acid"

                E=parameters_self.particle()
                E.name="E"
                E.radius=0.5 * parameters_self.units('reduced_length')
                E.type={"unprotonated": 49, "protonated": 50}
                E.q={"protonated": 0, "unprotonated": -1}
                E.acidity="acid"
                
                J=parameters_self.particle()
                J.name="J"
                J.radius=0.5 * parameters_self.units('reduced_length')
                J.type={"unprotonated": 51, "protonated": 52}
                J.q={"protonated": -1, "unprotonated": -2}
                J.acidity="acid"
                
                U=parameters_self.particle()
                U.name="U"
                U.radius=0.5 * parameters_self.units('reduced_length')
                U.type={"unprotonated": 53, "protonated": 54}
                U.q={"protonated": -1, "unprotonated": -2}
                U.acidity="acid"
                
                Z=parameters_self.particle()
                Z.name="Z"
                Z.radius=0.5* parameters_self.units('reduced_length')
                Z.type={"unprotonated": 55, "protonated": 56}
                Z.q={"protonated": -1, "unprotonated": -2}
                Z.acidity="acid"

            parameters_self.one_bead_peptide=one_bead_peptide

            class two_bead_peptide:

                name='2beadpeptide'
                principal_chain={"c": "sequence", "n": "sequence", "default": "C_alpha"}
                side_chain={"c": [], "n": [], "default": ["sequence"]}
                
                principal_bond=parameters_self.bond()
                principal_bond.actors=['C_alpha','C_alpha']
                principal_bond.type='harmonic'
                principal_bond.bondl=0.388*parameters_self.units.nm
                principal_bond.k=2000*parameters_self.units.kJ / parameters_self.units.nm**2 / parameters_self.units.mol

                C_alpha=parameters_self.particle()
                C_alpha.name="C_alpha"
                C_alpha.acidity="inert"
                C_alpha.radius=0.5 * parameters_self.units('reduced_length')
                C_alpha.type=20
                C_alpha.q=0
                
                A=parameters_self.particle()
                A.name="A"
                A.acidity="inert"
                A.radius=0.5 * parameters_self.units('reduced_length')
                A.type=21
                A.q=0

                N=parameters_self.particle()
                N.name="N"
                N.radius=0.5 * parameters_self.units('reduced_length')
                N.type=22
                N.q=0
                N.acidity="inert"
                
                Q=parameters_self.particle()
                Q.name="Q"
                Q.radius=0.5 * parameters_self.units('reduced_length')
                Q.type=23
                Q.q=0
                Q.acidity="inert"
                
                G=parameters_self.particle()
                G.name="G"
                G.radius=0.5 * parameters_self.units('reduced_length')
                G.type=24
                G.q=0
                G.acidity="inert"
                
                I=parameters_self.particle()
                I.name="I"
                I.radius=0.5 * parameters_self.units('reduced_length')
                I.type=25
                I.q=0
                I.acidity="inert"
                
                L=parameters_self.particle()
                L.name="L"
                L.radius=0.5 * parameters_self.units('reduced_length')
                L.type=26
                L.q=0
                L.acidity="inert"
                
                M=parameters_self.particle()
                M.name="M"
                M.radius=0.5 * parameters_self.units('reduced_length')
                M.type=27
                M.q=0
                M.acidity="inert"
                
                F=parameters_self.particle()
                F.name="F"
                F.radius=0.5 * parameters_self.units('reduced_length')
                F.type=28
                F.q=0
                F.acidity="inert"
                
                P=parameters_self.particle()
                P.name="P"
                P.radius=0.5 * parameters_self.units('reduced_length')
                P.type=29
                P.q=0
                P.acidity="inert"
                
                S=parameters_self.particle()
                S.name="S"
                S.radius=0.5* parameters_self.units('reduced_length')
                S.type=30
                S.q=0
                S.acidity="inert"
                
                T=parameters_self.particle()
                T.name="T"
                T.radius=0.5* parameters_self.units('reduced_length')
                T.type=31
                T.q=0
                T.acidity="inert"
                
                W=parameters_self.particle()
                W.name="W"
                W.radius=0.5 * parameters_self.units('reduced_length')
                W.type=32
                W.q=0
                W.acidity="inert"
                
                Y=parameters_self.particle()
                Y.name="Y"
                Y.radius=0.5* parameters_self.units('reduced_length')
                Y.type=33
                Y.q=0
                Y.acidity="inert"
                
                V=parameters_self.particle()
                V.name="V"
                V.radius=0.5* parameters_self.units('reduced_length')
                V.type=34
                V.q=0
                V.acidity="inert"
                
                K=parameters_self.particle()
                K.name="K"
                K.radius=0.5 * parameters_self.units('reduced_length')
                K.type={"unprotonated": 35, "protonated": 36}
                K.q={"unprotonated": 0, "protonated": 1}
                K.acidity="basic"
                
                H=parameters_self.particle()
                H.name="H"
                H.radius=0.5* parameters_self.units('reduced_length')
                H.type={"unprotonated": 37, "protonated": 38}
                H.q={"unprotonated": 0, "protonated": 1}
                H.acidity="basic"
                
                C=parameters_self.particle()
                C.name="C"
                C.radius=0.5* parameters_self.units('reduced_length')
                C.type={"unprotonated": 39, "protonated": 40}
                C.q={"unprotonated": 0, "protonated": 1}
                C.acidity="basic"
                
                R=parameters_self.particle()
                R.name="R"
                R.radius=0.5 * parameters_self.units('reduced_length')
                R.type={"unprotonated": 41, "protonated": 42}
                R.q={"unprotonated": 0, "protonated": 1}
                R.acidity="basic"
                
                n=parameters_self.particle()
                n.name="n"
                n.radius=0.5* parameters_self.units('reduced_length')
                n.type={"unprotonated": 43, "protonated": 44}
                n.q={"unprotonated": 0, "protonated": 1}
                n.acidity="basic"

                c=parameters_self.particle()
                c.name="c"
                c.radius=0.5 * parameters_self.units('reduced_length')
                c.type={"unprotonated": 45, "protonated": 46}
                c.q={"protonated": 0, "unprotonated": -1}
                c.acidity="acid"

                D=parameters_self.particle()
                D.name="D"
                D.radius=0.5 * parameters_self.units('reduced_length')
                D.type={"unprotonated": 47, "protonated": 48}
                D.q={"protonated": 0, "unprotonated": -1}
                D.acidity="acid"

                E=parameters_self.particle()
                E.name="E"
                E.radius=0.5 * parameters_self.units('reduced_length')
                E.type={"unprotonated": 49, "protonated": 50}
                E.q={"protonated": 0, "unprotonated": -1}
                E.acidity="acid"
                
                J=parameters_self.particle()
                J.name="J"
                J.radius=0.5 * parameters_self.units('reduced_length')
                J.type={"unprotonated": 51, "protonated": 52}
                J.q={"protonated": -1, "unprotonated": -2}
                J.acidity="acid"
                
                U=parameters_self.particle()
                U.name="U"
                U.radius=0.5 * parameters_self.units('reduced_length')
                U.type={"unprotonated": 53, "protonated": 54}
                U.q={"protonated": -1, "unprotonated": -2}
                U.acidity="acid"
                
                Z=parameters_self.particle()
                Z.name="Z"
                Z.radius=0.5 * parameters_self.units('reduced_length')
                Z.type={"unprotonated": 55, "protonated": 56}
                Z.q={"protonated": -1, "unprotonated": -2}
                Z.acidity="acid"

            parameters_self.two_bead_peptide=two_bead_peptide
