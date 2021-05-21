import numpy as np
import math as mt
import parameters as param

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
    
class residue:

    def __init__(self, name):

        self.name=name

    bondl=None
    k=None
    part=[]
    ids=[]


class molecule:

    N = None
    Nm = None
    ids = []
    bondl = None
    k = None

    def __init__(self, sequence, model=None, param_custom=None, pKa_set=None, pKa_custom=None):


        model_param=None

        if model is not None:

            model_names=[]
 
            model_list=param.get_subclasses(param)

            for param_set in model_list:

                if param_set.name is not None:

                    if model.lower() == param_set.name.lower():

                        model_param=param_set

            if model_param is None:

                raise ValueError("Unknown model chosen: ", model, ". Valid options are ", model_names)

            model_actors=param.get_subclasses(model_param)
            keys=[]

            for actor in model_actors:
                
                keys.append(actor.name)

            if param_custom is not None:

                for key in param_custom.keys():

                    keys.append(key)

            clean_sequence=sequence_parser(sequence, keys)
            
        else:

            if param_custom is None:
            
                clean_sequence=sequence_parser(sequence)

            else:
                
                keys=list(param.general.aminoacid.values())
                
                for key in param_custom.keys():
                    
                    keys.append(key)

                clean_sequence=sequence_parser(sequence, keys)

        print(clean_sequence)
        exit()

        for residue in clean_sequence:

            monomer=residue(name=residue)
            monomer.bondl=model_param.bondl[residue]
            monomer.k=model_param.k[residue]
            
            bead_list=[]

            if residue in pKa.keys():

                monomer.pKa=pKa[residue]
            
            # Create an alpha carbon bead


            if beads_per_monomer == 2 and residue != "c" and residue != "n":

                alpha_carbon=particle()
                alpha_carbon.name="C_alpha"
                alpha_carbon.q={"neutral": 0}
                alpha_carbon.type={"neutral": 20}
                bead_list.append(alpha_carbon)

            # Create the lateral chain bead

            lateral_chain=particle()
            lateral_chain.name=residue
            lateral_chain.q=model_param.q[residue]
            lateral_chain.type=model_param.type[residue]
            lateral_chain.radi=model_param.radi[residue]
          
            bead_list.append(lateral_chain)
            monomer.part=bead_list

            self.sequence.append(monomer)

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

        for item in sequence:
               
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

def create_protein(system, protein, initial_state='neutral'):
    '''
    Creates a protein chain in the system with the sequence stored in protein, joined by harmonic_bond)
    Inputs:
    system: espresso class object with all system variables
    protein: class object as defined in this library
    '''
    
    import espressomd
    from espressomd import interactions

    param=parameters()

    if initial_state.lower() not in ['neutral', 'charged']:

        raise ValueError("Unvalid key for the initial state of the residue, valid options are 'neutral' or 'charged'. Key given: ",  initial_state)

    if protein.N is None:

        protein.N=1

    # Create a director vector for the principal chain

    backbone_vec=generate_trialvectors(protein.bondl.to('sigma').magnitude)

    # Create the harmonic bond for the principal chain

    backbone_bond = interactions.HarmonicBond(k=protein.k.to('sigma**-2').magnitude, r_0=protein.bondl.to('sigma').magnitude)
    
    system.bonded_inter.add(backbone_bond)
    # If there are other particles in system, the first id = max(current ids) +1 else first id = 0

    if  not system.part[:].id:

        bead_id=0

    else:

        bead_id=max(system.part[:].id)

    # create the protein chains in the espresso system

    for pep_chain in range(protein.N):

        first_monomer=True
        ids_chain=[]

        for residue in protein.sequence:

            n_bead=0
            list_ids_aa=residue.ids.copy()
            ids_aa=[]

            for bead in residue.part:

                bead_id+=1
                list_ids_bead=bead.ids.copy()

                if initial_state.lower() == 'charged':

                    if initial_state.lower() in bead.q.keys():

                        state=initial_state

                    else:
                        
                        state='neutral'
                else:

                    state=initial_state

                if first_monomer: # The first monomer is placed in a random place of the simulation box

                    first_monomer=False
                    pos_backbone=np.random.random((1, 3)) *np.copy(system.box_l)
                    id_backbone=bead_id
                    name_backbone=residue.name
                    system.part.add(id=[bead_id], pos=pos_backbone, type=[bead.type[state]], q=[bead.q[state]])
                    bead.state=state
                    bead.N+=1
                    list_ids_bead.append([bead_id])
                    bead.ids=list_ids_bead.copy()
                    ids_aa.append(bead_id)
                    ids_chain.append(bead_id)

                else:

                    if n_bead ==0:

                        pos_backbone+=backbone_vec
                        system.part.add(id=[bead_id], pos=pos_backbone, type=[bead.type[state]], q=[bead.q[state]])
                        
                        if (protein.beads_per_monomer == 1):

                            sidechain_bond = interactions.HarmonicBond(k=residue.k.to('sigma**-2').magnitude, r_0=residue.bondl.to('sigma').magnitude)
                            system.bonded_inter.add(sidechain_bond)
                            system.part[bead_id].add_bond((sidechain_bond, id_backbone))

                        else:


                            if (name_backbone =="n") or ( residue.name == "n"):

                                sidechain_bond = interactions.HarmonicBond(k=model_param.k["n"].to('sigma**-2').magnitude, r_0=model_param.bondl["n"].to('sigma').magnitude)
                                system.bonded_inter.add(sidechain_bond)
                                system.part[bead_id].add_bond((sidechain_bond, id_backbone))

                            elif (name_backbone =="c") or ( residue.name == "c"):
                                sidechain_bond = interactions.HarmonicBond(k=model_param.k["c"].to('sigma**-2').magnitude, r_0=model_param.bondl["c"].to('sigma').magnitude)
                                system.bonded_inter.add(sidechain_bond)
                                system.part[bead_id].add_bond((sidechain_bond, id_backbone))

                            else:

                                system.part[bead_id].add_bond((backbone_bond, id_backbone))
                        
                        list_ids_bead.append([bead_id])
                        bead.ids=list_ids_bead.copy()
                        bead.N+=1
                        bead.state=state
                        ids_aa.append(bead_id)
                        ids_chain.append(bead_id)
                        id_backbone=bead_id
                        name_backbone=residue.name

                    elif n_bead ==1:
                        
                        
                        rand_vec=generate_trialvectors(1)
                        sidechain_vec=np.cross(rand_vec,backbone_vec)
                        r = mt.sqrt(sum(x*x for x in sidechain_vec))
                        sidechain=np.array([x/r for x in sidechain_vec])
                        sidechain_vec=sidechain_vec*residue.bondl.to('sigma').magnitude
                        pos_sidechain=pos_backbone+sidechain_vec
                        system.part.add(id=[bead_id], pos=pos_sidechain, type=[bead.type[state]], q=[bead.q[state]])
                        sidechain_bond = interactions.HarmonicBond(k=residue.k.to('sigma**-2').magnitude, r_0=residue.bondl.to('sigma').magnitude)
                        system.bonded_inter.add(sidechain_bond)
                        system.part[bead_id].add_bond((sidechain_bond, id_backbone))
                        list_ids_bead.append([bead_id])
                        bead.state=state
                        bead.N+=1
                        bead.ids=list_ids_bead.copy()
                        ids_aa.append(bead_id)
                        ids_chain.append(bead_id)
                
                n_bead+=1
            
            list_ids_aa.append(ids_aa)
            residue.ids=list_ids_aa.copy()
        protein.ids.append(ids_chain)

    return 



def count_titrable_groups(protein):
    """
    Counts the number of titrable groups in the protein object

    Input:
    protein: class object as defined in this library, attributes:

    Output:
    N_ti: (int) number of titrable groups in the sequence
    """

    N_ti=0

    for residue in protein.sequence:

        if residue.pKa is not None:

            N_ti+=1

    N_ti*=protein.N

    return N_ti

def setup_protein_acidbase_reactions(RE, protein, cation):
    """
    Set up the Acid/Base reactions for acidic/basidic residues in protein. The reaction steps are done following the constant pH ensamble procedure. 

    Inputs:
    RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
    protein: class object as defined in this library, attributes:
    cation: class objects of the cations,  with the following attributes:
        q: (int) charge of the ion
        type: (int) type of the ion
    """

    reaction_absent={}

    for group in param.general.pka_hass.keys():
        
        reaction_absent[group]=True

    for residue in protein.sequence:

        if residue.pKa is not None:
        
            for bead in residue.part:

                if (reaction_absent[residue.name] and "charged" in bead.q.keys()): 

                    if (bead.q["charged"] == 1) : # Basic residue
                        
                        RE.add_reaction(gamma=10**-residue.pKa,
                            reactant_types=[bead.type["charged"]],
                            reactant_coefficients=[1],
                            product_types=[bead.type["neutral"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={bead.type["neutral"]: bead.q["neutral"],
                                             bead.type["charged"]: bead.q["charged"],
                                             cation.type: cation.q})
                        reaction_absent[residue.name] = False

                    elif (bead.q["charged"] == -1) : # Acid residue

                        RE.add_reaction(gamma=10**-residue.pKa,
                            reactant_types=[bead.type["neutral"]],
                            reactant_coefficients=[1],
                            product_types=[bead.type["charged"], cation.type],
                            product_coefficients=[1, 1],
                            default_charges={bead.type["neutral"]: bead.q["neutral"],
                                             bead.type["charged"]: bead.q["charged"],
                                             cation.type: cation.q})
                        reaction_absent[residue.name] = False

                    else:
    
                        raise ValueError("This subrutine is concived for the acid/base equilibria of monovalent ions. Charge of residue ", residue.name, " = ", bead.q["charged"])

    return

def calculate_protein_charge(system, protein):
    """ 
    Calculates the charge of the protein and its square

    Inputs:
    system: espresso class object with all system variables
    protein: class object as defined in this library

    Outputs:
    Z_prot: (float) total charge of the protein
    Z2: (float) square of Z_prot
    """
    
    Z_prot=0

    for chain in protein.ids:
        
        for id in chain:

            Z_prot+=system.part[id].q

    Z_prot=Z_prot/protein.N
    Z2=Z_prot**2

    return Z_prot, Z2

def track_ionization(system, protein):
    """
    Sets up espresso to track the average number of particles of each residue type
    
    Inputs:
    system: espresso class object with all system variables
    protein: class object as defined in this library

    """

    types=[]

    for residue in protein.sequence:
        
        for bead in residue.part:

            for type in bead.type.values():

                types.append(type)

    system.setup_type_map(types)

    return

def calculate_HH(pH, protein):
    """
    Calculates the ideal Henderson-Hassebach titration curve in the given pH range

    Inputs:
    pH: (list of floats) list of pH values
    protein: class object as defined in this library

    Outputs:
    Z_HH: (list of floats) Henderson-Hasselnach prediction of the protein charge in the given pH
    """

    Z_HH=[]

    for pH_value in pH:
        
        Z=0

        for residue in protein.sequence:

            if residue.pKa is not None:

                for bead in residue.part:

                    if "charged" in bead.q.keys():

                        Z+=bead.q["charged"]/(1+10**(bead.q["charged"]*(pH_value-residue.pKa)))

        Z_HH.append(Z)

    return Z_HH


def create_counterions(system,protein):
    """
    Adds one monovalent counter-ion (cation or anion) of opposite charge per each charge in the peptide.

    Inputs:

    system: espresso class object with all system variables.
    protein: class object as defined in this library

    Outputs:
    cation/anion: particle class objects of the counter-ions, as defined in this library
        q: (int) charge of the ion
        type: (int) type of the ion

    Assumptions:

    -The charges in system are monovalent
    -The counter-ions are monovalent

    """

    # Load the parameters

    param=parameters()

    # Create and set-up the cation particle object

    cation=particle()
    cation.name="cation"
    cation.q=param.q[cation.name]
    cation.type=param.type[cation.name]
    cation.radi=param.radi[cation.name]
    cation.state="charged"

    # Create and set-up the cation particle object

    anion=particle()
    anion.name="anion"
    anion.q=param.q[anion.name]
    anion.type=param.type[anion.name]
    anion.radi=param.radi[anion.name]
    anion.state="charged"

    # The ids of the counter ions are created to follow the ones of the peptide

    if  len(system.part[:].id) != 0:

        I_id=max(system.part[:].id)

    total_q=0

    # Create a counter-ion for each charged group in the peptide

    cation_id=[]
    anion_id=[]

    for residue in protein.sequence:

        for bead in residue.part:
           
            if bead.state is None:

                raise ValueError("The peptide chains must be created first in the system before creating their counter-ions")

            else:

                for npart in range(bead.N):

                    if bead.q[bead.state] == 1: # Create an anion

                        I_id+=1
                        I_pos=np.random.random((1, 3)) * system.box_l
                        system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[anion.q])
                        anion.N+=1
                        anion_id.append(I_id)
                        total_q+=bead.q[bead.state]+anion.q

                    elif bead.q[bead.state] == -1: # Create an anion

                        I_id+=1
                        I_pos=np.random.random((1, 3)) * system.box_l
                        system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[cation.q])
                        cation.N+=1
                        cation_id.append(I_id)
                        total_q+=bead.q[bead.state]+cation.q
                    
                    elif bead.q[bead.state] != 0:

                        raise ValueError("This subroutine only considers the case of peptide chains with monovalent charged residues")
    
    cation.ids=cation_id
    anion.ids=anion_id

    if total_q != 0:

        raise ValueError("Subrutine attempted to create a not electroneutral system, please review your peptide set-ups")


    return cation, anion

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

