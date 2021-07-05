#from operator import ne
import numpy as np
import math as mt
import parameters as param
import random as rn
from parameters import particle, bond

class residue:

    def __init__(self, name):

        self.name=name

    principal_bead=None
    lateral_beads=[]
    bonds=[]
    beads=[]
    ids=[]
    N=None

class molecule:

    N = None
    Nm = None
    ids = []
    model = None
     

    def __init__(self, sequence, model=None, param_custom=None, pKa_set=None, pKa_custom=None):


        model_param=None
        model_dict={}
    
        if model is not None:

            model_list=get_subclasses(param)

            # Get the model parameters

            for param_set in model_list:

                if param_set.name is not None:
                    
                    if model.lower() == param_set.name.lower():

                        model_param=param_set

            if model_param is None:

                model_names=get_modelnames()
                raise ValueError("Unknown model chosen: ", model, ". Valid options are ", model_names)
            
            # Get a list of model parameters
            
            list_model=get_attributes(model_param)

            for par in list_model:

                model_dict[par[0]]=par[1]

            # Get valid keys for the model

            model_actors=get_particles(model_param)
            keys=[]

            for actor in model_actors:
       
                keys.append(actor.name)
                
        else:
            
            keys=list(param.general.aminoacid.values())
               
        # Load Custom parameters if any

        # Check that the custom parameters are in the proper format:

        names_custom=[]
        
        if  param_custom is not None:

            if isinstance(param_custom, param.custom_parameters):

                if model_param is not None:
                    
                # Overwrite model parameters by custom ones

                    list_custom=get_attributes(param_custom)
                    custom_dict={}
            
                    for par in list_custom:

                        custom_dict[par[0]]=par[1]

                    for atr in custom_dict.keys():

                        if custom_dict[atr] is not None:

                            if isinstance(custom_dict[atr], param.particle):
                                    
                                if atr in model_dict.keys():

                                    part_custom=get_attributes(custom_dict[atr])

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
        
        # Store the model in the molecule object

        self.model=model_param

        if pKa_custom is not None:

            for key in pKa_custom.keys():

                keys.append(key)
        
        clean_sequence=sequence_parser(sequence, keys)
        
        # Import the pKa values

        param_pKa_set=get_subclasses(param.pKa_set)
        names_pKa_set=[]

        for sets in param_pKa_set:
            
            names_pKa_set.append(sets.name)

        if pKa_set is None:

            pKa_set=param.pKa_set.Hass.pKa

        elif isinstance(pKa_set, str):

            for sets in param_pKa_set:

                if (pKa_set.lower() == sets.name.lower()):

                    pKa_set=sets.pKa
                    break

            if not isinstance(pKa_set,dict):

                raise ValueError("Unknown key provided for a pKa set, valid options are " , names_pKa_set)
            
        else:
        

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

            monomer=residue(name=res)
            self.residues.append([monomer])

        # Set-up the model

        if model_param is None:

            # If no model nor custom parameters are given, create one bead particle per residue

            for r_chain in self.residues:
                
                for res in r_chain:
                    bead=particle()
                    bead.name=res.name
                    bead_list=[]
                
                # If the residue has a pKa value list in the pKa_set put in the bead

                    if res.name in pKa_set.keys():

                        bead.pKa=pKa_set[res.name]

                    bead_list.append(bead)
                    res.beads=bead_list
        
        else:

            param_part=get_particles(model_param)
            param_bonds=get_bonds(model_param)
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

                            bead=particle()
                            bead.name=res.name

                        if res.name in pKa_set.keys():

                            bead.pKa=pKa_set[res.name]

                                    
                    elif p_bead_name in model_keys:

                        for p_set in param_part:

                            if (p_set.name == p_bead_name):

                                bead=p_set

                                if p_bead_name in pKa_set.keys():

                                    bead.pKa=pKa_set[res.name]

                        
                    else:
                    
                        raise ValueError("Unknown key for the principal chain: ", model_param.principal_chain)
                        
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
            
def create_custom_model(beads_per_residue=None, principal_chain=None, side_chain=None, custom_particles=None, units=None):
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
    
    custom_param=param.custom_parameters()

    
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
                custom_part=param.particle()
                properties=get_attributes(custom_part)
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

def sequence_parser(sequence, keys=param.general.aminoacid.values()):
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

                        if (residue.upper() in param.general.aminoacid.keys()):

                            clean_sequence.append(param.general.aminoacid[residue.upper()])

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

                    elif (residue.upper() in param.general.aminoacid.keys()):

                        clean_sequence.append(param.general.aminoacid[residue.upper()])

                    else:

                        raise ValueError("Unknown code for a residue: ", residue, " please review the input sequence")

                clean_sequence.append(residue_ok)

    return clean_sequence

def write_parameters(mol):

    if isinstance(mol,molecule):

        print('molecule parameters')
        for atr in get_attributes(mol):
            print(atr)

        for chain in mol.residues:
            
            for res in chain:

                print("\t residue ", res.name, 'parameters')
        
                for atr in get_attributes(res):

                    print('\t' ,atr)


                for chain_bead in res.beads:

                    for bead in chain_bead:

                        print("\t \t Particle", bead.name)
            
                        for atr in get_attributes(bead):

                            print('\t \t ', atr)

    if isinstance(mol,residue):
        
        print("residue ", mol.name, 'parameters')

        for atr in get_attributes(mol):
            
            print(atr)

        for chain_bead in res.beads:

            for bead in chain_bead:

                print("\t \t Particle", bead.name)
            
                for atr in get_attributes(bead):

                    print('\t \t ', atr)

    if isinstance(mol, particle): 

        print("Particle", mol.name)
            
        for atr in get_attributes(mol):

            print(atr)   

def create_particle(part, system, position=None, state=None, id=None):
    '''
    Creates a bead particle in the espresso system with the properties stored in bead.
    part = instance of a particle object from sugar library
    system = instance of a system object from espressomd library
    position = array of particle positions of size part.N
    state = for particles with more than one possible state, defines in which state are created
    id = array of particle ids of size part.N
    '''

    if not isinstance(part,particle):

        raise ValueError("Please provide an instance of a particle object from sugar library. Given:", part)

    if position is None:
    
    # By default, the beads are placed in random positions
        
        pos=np.random.random((1, 3))[0] *np.copy(system.box_l)

    else:

        if (len(position) == part.N):

            pos=position[0]
        
        else:
            
            raise ValueError("Please provide ", part.N ," position arrays")
            
    if state is not None:
        
        if isinstance(part.q,dict):
            
            if state in part.q.keys():

                part.state=state

            else:

                raise ValueError("Unknown state for bead: ", part.name, " please review the input state:", state, ' valid keys are', part.q.keys())
    else:

        if isinstance(part.type, dict):
            
            # Inicialice beads with more than one state in a random state

            state=rn.choice(list(part.type.keys()))
            part.state=state

    
    if isinstance(part.q, int) or isinstance(part.q, float) :
        
        q=part.q
    
    elif isinstance(part.q, dict):

        q=part.q[state]

    else:

        raise ValueError("Unvalid charge for bead: ", part.name, " given ", part.q)

    if isinstance(part.type, int):
        
        type=part.type

    elif isinstance(part.type, dict):

        type=part.type[state]

    else:

        raise ValueError("Unvalid type for bead: ", part.name, " given ", part.type)

    
    if id is None:

        if len(system.part[:].id) == 0:
            
            bead_id=0

        else:
            
            bead_id=max(system.part[:].id)+1
           

    else:
        if isinstance(id,list):
            
            if len(id) == part.N:

                bead_id=id[0]

            else:

                raise ValueError('Please provide ', part.N, 'ids. Provided:', id)

        else:

            raise ValueError('particle ids must be provided as a list. Given:', id)

    bead_id_lst=[]

    # By default only one bead is created

    if not isinstance(part.N, int): 

        part.N=1

    for bead_i in range(part.N):
            
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

            pos=np.random.random((1, 3))[0] *np.copy(system.box_l)

    
    part.ids=bead_id_lst   

    return

def create_residue(res, system, position=None, state=None, id=None):
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
                    new_bead=particle()
                    bd_atrb=get_attributes(bead)

                    for atr in bd_atrb:

                        setattr(new_bead,atr[0],atr[1])
                
                    new_bead.N=1
            
                    create_particle(part=new_bead, system=system, position=pos, state=bead_state, id=bead_id)                        

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
                
                    actors=[central_name,bead_name]
                    central_bond=search_bond(actors,res)

                # Generate the position 

                    bond_vector=generate_trialvectors(central_bond.bondl.to('nm').magnitude)

                    if position is not None and len(position) > 0:
                
                        prebead_position=position[0]

                    else:
                        prebead_position=central_pos+bond_vector

                    for b_chain in res.beads:
                        
                        for bead in b_chain:

                            if bead.name == bead_name:

                     # Create an individual particle object
                                new_bead=particle()
                                bd_atrb=get_attributes(bead)

                                for atr in bd_atrb:

                                    setattr(new_bead,atr[0],atr[1])
                                new_bead.N=1
                        
                                create_particle(part=new_bead, system=system, position=[prebead_position], state=bead_state, id=bead_id)

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

                                create_bond(system,bond=central_bond,id1=new_bead.ids[0], id2=central_id[0])

                                prebead_id=new_bead.ids[0]
                                prebead_name=new_bead.name
                                bead_list.append(new_bead)
                                residue_ids.append(new_bead.ids[0])
                                break

                        start_position=False

                else:

                    actors=[prebead_name,bead_name]
                    lateral_bond=search_bond(actors,res)
                
                    if position is not None and len(position) > 0:
                
                        prebead_position=position[0]

                    else:
                    
                        prebead_position=prebead_position+bond_vector
                
                    for b_chain in res.beads:
                    
                        for bead in b_chain:

                            if bead.name == bead_name:

                            # Create an individual particle object
                        
                                new_bead=particle()
                                bd_atrb=get_attributes(bead)

                                for atr in bd_atrb:

                                    setattr(new_bead,atr[0],atr[1])
                        
                                new_bead.N=1

                        
                                create_particle(part=new_bead, system=system, position=[prebead_position], state=bead_state, id=bead_id)

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


                                create_bond(system,bond=lateral_bond,id1=prebead_id,id2=new_bead.ids[0])

                                prebead_name=new_bead.name
                                prebead_id=new_bead.ids[0]
                                bead_list.append(new_bead)
                                residue_ids.append(new_bead.ids[0])

                                break

    # Store the data in res

        res.ids=[residue_ids]
        res.beads=[bead_list]    

    return

def create_molecule(mol,system):
    '''
    Creates a molecules in the espresso given with the properties stored in mol
    mol: instance of a molecule object as defined in sugar library
    system: instance of a system object of espresso library
    position: list of arrays with the position of each residue of the molecule. 
    The positions are asigned in the same order than the particles stored in residue (first the principal bead and then the lateral ones)
    state = for particles with more than one possible state, defines in which state are created
    id = list of arrays with the ids of each residue of the molecule. 
    '''
    
    if mol.N is None:

        mol.N=1

    elif not isinstance(mol.N,int):

        raise ValueError("The number of molecules must be an integer number, given: ", mol.N)



    for mol_i in range(mol.N):

        first_res_inexistent=True
        bond_vector=generate_trialvectors(1)
        id_list=[]

        for r_chain in mol.residues:

            for res in r_chain:

                if (first_res_inexistent):

                    create_residue(res, system)
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
                    backbone_bond=search_bond(actors,res)
       
                    new_bead_pos=pre_bead_pos+bond_vector*backbone_bond.bondl.to('nm').magnitude
                    create_residue(res, system, position=[new_bead_pos])
                    
                    # Create the harmonic bond for the principal chain

                    ids_res=res.ids[0]
                    id_ppal_bead=ids_res[0]
                    
                    create_bond(system,bond=backbone_bond,id1=pre_bead_id,id2=id_ppal_bead)

                    # Update lists and variables
                    
                    id_list+=ids_res
                    pre_bead_id=res.ids[0]
                    pre_bead_name=new_bead_name
                    pre_bead_pos=new_bead_pos
        
        mol.ids.append(id_list)

    return 

def count_titrable_groups(mol):
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

def setup_acidbase_reactions(RE, mol, cation):
    """
    Set up the Acid/Base reactions for acidic/basidic residues in mol. The reaction steps are done following the constant pH ensamble procedure. 

    Inputs:
    RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
    mol: molecule class object as defined in sugar library
    cation: particle class object as defined in sugar library
    """

    reaction_absent={}

    if isinstance(mol,molecule):

        for chain in mol.residues:
            
            for res in chain:

                for chain_bead in res.beads:

                    for bead in chain_bead:

                        if bead.pKa is not None and  bead.acidity in ['acid','basic']:

                            if bead.name not in reaction_absent.keys():

                                reaction_absent[bead.name]=True

                            if (reaction_absent[bead.name]):

                                setup_bead_acidbase_reaction(RE, bead, cation)
                                reaction_absent[bead.name]=False

    if isinstance(mol,residue):
        
        for chain_bead in res.beads:

            for bead in chain_bead:

                if bead.pKa is not None and  bead.acidity in ['acid','basic']:

                    if bead.name not in reaction_absent.keys():

                        reaction_absent[bead.name]=True

                    if (reaction_absent[bead.name]):

                        setup_bead_acidbase_reaction(RE, bead, cation)
                        reaction_absent[bead.name]=False

    if isinstance(mol, particle): 

        if mol.pKa is not None and  mol.acidity in ['acid','basic']:

            setup_bead_acidbase_reaction(RE, bead, cation)
            reaction_absent[bead.name]=False

    return

def setup_bead_acidbase_reaction(RE, part, cation):
    """
    Set up the Acid/Base reactions for acidic/basidic residues in protein. The reaction steps are done following the constant pH ensamble procedure. 

    Inputs:
    RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
    mol: particle/residue/molecule class object as defined in sugar library
    cation: particle class object as defined in sugar library
    """

    if not isinstance(part, particle) or not isinstance(cation, particle):

        raise ValueError("part and cation must be instances of a particle object from sugar library")

    if part.acidity not in ['acid','basic']:

        print("WARNING, the added particle does not have its acidity defined, please define part.acidity to 'acid' or 'basic'. No reaction has been added")
    
    if cation.q != 1:

        print("the subroutine is concived for monovalent cations. The charge of the cation is set to 1, given cation.q = ", cation.q)

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
                            product_types=[part.type["unprotonated"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={part.type["protonated"]: 1,
                            part.type["unprotonated"]: 0,
                            cation.type: 1})

        elif (part.acidity == 'acid') : # Acid residue

            RE.add_reaction(gamma=10**-part.pKa,
                            reactant_types=[part.type["protonated"]],
                            reactant_coefficients=[1],
                            product_types=[part.type["unprotonated"], cation.type],
                            product_coefficients=[1, 1],
                            default_charges={part.type["protonated"]: 0,
                                            part.type["unprotonated"]: -1,
                                            cation.type: 1})

    return

def calculate_molecule_charge(system, mol):
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

    if isinstance(mol,molecule) or isinstance(mol,residue):

        for chain in mol.ids:
            
            for id in chain:

                Z+=system.part[id].q

    if isinstance(mol, particle): 

        for id in mol.ids:

            Z+=system.part[id].q

    Z=Z/mol.N
    Z2=Z**2

    return Z, Z2

def track_ionization(system, mol):
    """
    Sets up espresso to track the average number of particles of the types contained in mol
    
    Inputs:
    system: espresso class object with all system variables
    mol: particle/residue/molecule class object as defined in sugar library

    """

    types=[]

    if isinstance(mol,molecule) or isinstance(mol,residue):

        for chain in mol.ids:
            
            for id in chain:

                b_type=system.part[id].type
                
                if b_type not in types:

                    types.append(b_type)


    elif isinstance(mol, particle): 

        for id in mol.ids:

            b_type=system.part[id].type
                
            if b_type not in types:

                types.append(b_type)

    else:
        
        raise ValueError("mol must be a particle/residue/molecule class object as defined in sugar library given: ", mol)

    system.setup_type_map(types)

    return

def calculate_HH( mol,pH=None):
    """
    Calculates the ideal Henderson-Hassebach titration curve in the given pH range

    Inputs:
    pH: (list of floats) list of pH values
    mol: particle/residue/molecule class object as defined in sugar library

    Outputs:
    Z_HH: (list of floats) Henderson-Hasselnach prediction of the protein charge in the given pH
    """

    if pH is None:

        pH=np.linspace(2,12,50)
    
    elif not isinstance(pH,list):

        raise ValueError("pH must contain a list with the pH-values where the Henderson-Hassebach titration curve will be calculated. Given: ", pH)        

    Z_HH=[]

    for pH_value in pH:
        
        Z=0

        if isinstance(mol,molecule):

            for chain in mol.residues:
            
                for res in chain:

                    for chain_bead in res.beads:

                        for bead in chain_bead:

                            Z+=calculate_HH_part(pH_value, bead)


        if isinstance(mol,residue):
        
            for chain_bead in res.beads:

                for bead in chain_bead:

                    Z+=calculate_HH_part(pH_value, bead)                    

        if isinstance(mol, particle): 

            if mol.pKa is not None and  mol.acidity in ['acid','basic']:

                Z+=calculate_HH_part(pH_value, bead)
                        
        Z_HH.append(Z)

    return Z_HH

def calculate_HH_part(pH, part):
    """
    Calculates the ideal Henderson-Hassebach titration curve of part at one pH-value

    Inputs:
    pH: (float) pH value
    part: particle class object as defined in sugar library

    Outputs:
    z: (float) Henderson-Hasselnach prediction of charge of part at the given pH
    """  

    if not isinstance(part, particle): 

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

def create_counterions(system,mol, cation=None, anion=None):
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

        cation=param.small_ions.cation
        cation.N=0

    if anion is None:

        anion=param.small_ions.anion
        anion.N=0

    
    # The ids of the counter ions are created to be larger then the ids in the system

    if  len(system.part[:].id) != 0:

        I_id=max(system.part[:].id)

    total_q=0

    # Create a counter-ion for each charged group in the peptide

    cation_id=[]
    anion_id=[]

    if isinstance(mol,molecule) or isinstance(mol,residue):

        for chain in mol.ids:
            
            for id in chain:

                q=system.part[id].q
                
                if q == +1:

                    I_id+=1
                    I_pos=np.random.random((1, 3)) * system.box_l
                    system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[anion.q])
                    anion.N+=1
                    anion_id.append(I_id)
                    total_q+=q+anion.q

                elif q == -1:

                    I_id+=1
                    I_pos=np.random.random((1, 3)) * system.box_l
                    system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[cation.q])
                    cation.N+=1
                    cation_id.append(I_id)
                    total_q+=q+cation.q


    if isinstance(mol, particle): 

        for id in mol.ids:

            q=system.part[id].q
                
            if q == +1:

                I_id+=1
                I_pos=np.random.random((1, 3)) * system.box_l
                system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[q])
                anion.N+=1
                anion_id.append(I_id)
                total_q+=q+anion.q

            elif q == -1:

                I_id+=1
                I_pos=np.random.random((1, 3)) * system.box_l
                system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[q])
                cation.N+=1
                cation_id.append(I_id)
                total_q+=q+cation.q
    
    cation.ids=cation_id
    anion.ids=anion_id

    if total_q != 0:

        raise ValueError("Subrutine attempted to create a not electroneutral system, please review your peptide setup")


    return cation, anion

def create_added_salt(system, cation, anion, N_ions):
    
    """
    Adds extra cations/anions to the system

    Inputs:

    system: espresso class object with all system variables.
    cation/anion:    particle class object as defined in sugar library
    N_ions: number of ions to be added 
    
    Assumptions:

    -The number of anions and cations provided must ensure electroneutrality
  
    """

        # The ids of the counter ions are created to be larger then the ids in the system

    if  len(system.part[:].id) != 0:

        I_id=max(system.part[:].id)

    total_q=0

    N_anion=0
    N_cation=0

    while N_anion < N_ions:

        I_id+=1
        I_pos=np.random.random((1, 3)) * system.box_l
        system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[anion.q])
        anion.N+=1
        anion.ids.append(I_id)
        total_q+=anion.q
        N_anion+=1

    while N_cation < N_ions:

        I_id+=1
        I_pos=np.random.random((1, 3)) * system.box_l
        system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[cation.q])
        cation.N+=1
        cation.ids.append(I_id)
        total_q+=cation.q
        N_cation+=1

    if total_q != 0:

        raise ValueError("Subrutine attempted to create a not electroneutral system, please review your cation/anion setup")


    return 

def generate_trialvectors(mag):
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algorithm from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    Input:
    
    mag: magnitude of the vector
    
    Output:
    
    vec: random vector of dimension 3
    """

    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )

    vec=np.array([x,y,z])
    vec=vec*mag

    return vec

def get_subclasses(cls):

    import inspect
    model_list=[]

    for attribute in cls.__dict__.items():

        name_attribute=attribute[0]

        if name_attribute[0:2] != "__" and inspect.isclass(attribute[1]):

            model_list.append(attribute[1])

    return model_list

def get_particles(cls):

    particle_list=[]

    for attribute in cls.__dict__.items():

        if isinstance(attribute[1], param.particle):

            particle_list.append(attribute[1])

    return particle_list

def get_bonds(cls):

    bond_list=[]

    for attribute in cls.__dict__.items():

        if isinstance(attribute[1], param.bond):

            bond_list.append(attribute[1])

    return bond_list

def get_attributes(cls):
    import inspect

    attribute_list=[]

    for attribute in  inspect.getmembers(cls, lambda a:not(inspect.isroutine(a))):

        name_attribute=attribute[0]

        if name_attribute[0:2] != "__":

            attribute_list.append(attribute)

    return attribute_list

def newAttr(self, attr):

    setattr(self, attr.name, attr)

def get_modelnames():

    model_names=[]
    model_list=get_subclasses(param)

    for param_set in model_list:

        if param_set.name is not None:
                    
            model_names.append(param_set.name)

    return model_names

def search_bond(actors,res):
    """
    Search for a bond in res that joins the beads contained in actors

    Inputs:
    actors:list containing the bead objects
    res: residue object as defined in sugar library

    Outputs:
    asigned_bond: bond object as defined in sugar library

    """

    unasigned_bond=True

    for bond in res.bonds:

        if actors == bond.actors or actors[::-1] == bond.actors:

            asigned_bond=bond
            unasigned_bond=False
            break

    if unasigned_bond:

        for bond in res.bonds:

            if bond.actors[0] == "default":

                asigned_bond=bond
                unasigned_bond=False
                break

    if unasigned_bond:
        write_parameters(res)
        raise ValueError("The bond between ", actors, "is not defined in the residue")

    return asigned_bond

def create_bond(system,bond,id1,id2):
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


    if not isinstance(id1,int) and isinstance(id2,int):

        raise ValueError("id1 and id2 should be integer numbers, given id1 =", id1, "and id2 = ", id2)

    if not isinstance(bond, param.bond):

        raise ValueError("bond should be a sugar bond object, given: ", bond )

    if (bond.type == "harmonic"):

        bond_potential = interactions.HarmonicBond(k=bond.k.to('kJ / mol / nm**2').magnitude, r_0=bond.bondl.to('nm').magnitude)
                            

    else:

        raise ValueError("Unknown bond type", bond.type)


    system.bonded_inter.add(bond_potential)
    system.part[id1].add_bond((bond_potential, id2))

    return