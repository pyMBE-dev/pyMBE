import numpy as np
import pint
import math as mt

class particle:
    
    radi=None
    type=None
    q=None
    id=None
    name=None

class aminoacid:

    def __init__(self, name):

        self.name=name

    bondl=None
    k=None
    pKa=None
    part=[]
    ids=[]


class protein():

    N = None
    Nm = None
    ids = []
    bondl = None
    k = None

    def __init__(self, sequence, beads_per_monomer=1, pKa_set=None, pKa_custom=None):

        clean_sequence=[]
        self.sequence=[]
        param=parameters()
        self.bondl=0.388*param.ureg.nm
        self.k=100/param.ureg.sigma**2
    
        if isinstance(sequence, str):
                
            if (sequence.find("-") != -1):

                splited_sequence=sequence.split("-")

                for residue in splited_sequence:

                    if len(residue) == 1:

                        if residue == "c" or residue == "n":

                            residue_up=residue

                        else:

                            residue_up=residue.upper()

                        if (residue_up in param.one_letter_key.values()):

                            clean_sequence.append(residue_up)

                        else:

                            raise ValueError("Unknown one letter code for an aminoacid given: ", residue, " please review the input sequence")

                    else:

                        if (residue.upper() in one_letter_key.keys()):

                            clean_sequence.append(param.one_letter_key[residue.upper()])

                        else:

                            raise ValueError("Unknown three letter code for an aminoacid given: ", residue, " please review the input sequence")

            else:

                for letter in sequence:

                    if letter == "c" or letter == "n":

                        letter_up=letter
                    else:
                        letter_up=letter.upper()
                    
                    if (letter_up in param.one_letter_key.values()):

                        clean_sequence.append(letter_up)
                        
                    else:

                        raise ValueError("Unknown one letter code for an aminoacid given: ", letter, " please review the input sequence")

        if isinstance(sequence, list):

            for item in sequence:
               
                if item == "c" or item == "n":

                    item_up=item

                else:

                    item_up=item.upper()

                if item_up in param.one_letter_key.values():

                    clean_sequence.append(item_up)

                elif item_up in param.one_letter_key.keys():

                    clean_sequence.append(param.one_letter_key[item_up])
                
                else:

                    raise ValueError("Unknown code for an aminoacid given: ", item, " please review the input sequence")

        if (beads_per_monomer == 1 or beads_per_monomer == 2):

            self.beads_per_monomer=beads_per_monomer

        else:

            raise ValueError("The library is only ready for peptide models with 1 or 2 beads. Value provided for beads_per_monomer = ", beads_per_monomer , " not allowed.")

        if pKa_set is None or pKa_set.upper() == "HASS": # Default values 

            pKa=param.pka_hass

        elif pKa_set.upper() == "PLATZER": 

            pKa=param.pka_platzer

        elif pKa_set.upper() == "CRCHANDBOOK": 

            pKa=param.pka_crc

        elif ( pKa_set.upper() == "NOZAKI"): 

            pKa=param.pka_nozaki


        elif ( pKa_set.upper() == "CUSTOM"):

            if isinstance(pKa_custom, dict):
                
                pKa=param.pka_hass

                for custom_aa in pKa_custom.keys():

                    if custom_aa == "n" or custom_aa == "c":

                        custom_aa_up=custom_aa

                    else:

                        custom_aa_up=custom_aa.upper()

                    if custom_aa_up in param.one_letter_key.keys():

                        custom_aa_up=param.one_letter_key[custom_aa_up]

                    if custom_aa_up in  pKa.keys():

                        pKa[custom_aa_up]=pKa_custom[custom_aa_up]

                    else:

                        raise ValueError("Unknown aminoacid type given for custom pKa-value, please use the one letter aminoacide code as key. Key given: ",  custom_aa)

            else:

                raise ValueError("The custom pKa-values must be given as a dictionary such that pKa['one_letter_code'] = value. Given pKa = ",  pKa_custom)

        else:

            raise ValueError("Unknown option for the desired pKa set: ", pKa_set, "Valid options are 'Hass', 'Platzer', 'CRCHandbook','Nozaki' or 'custom'")

        for residue in clean_sequence:

            monomer=aminoacid(name=residue)
            monomer.bondl=param.bondl[residue]
            monomer.k=param.k[residue]
            
            bead_list=[]

            if residue in pKa.keys():

                monomer.pKa=pKa[residue]
            
            # Create an alpha carbon bead


            if beads_per_monomer == 2 and residue != "c" and residue != "n":

                alpha_carbon=particle()
                alpha_carbon.name="C_alpha"
                alpha_carbon.q={"neutral": 0}
                alpha_carbon.type={"neutral": 20}
                alpha_carbon.radi= 0.5*param.ureg.sigma
                bead_list.append(alpha_carbon)

            # Create the lateral chain bead

            lateral_chain=particle()
            lateral_chain.name=residue
            lateral_chain.q=param.q[residue]
            lateral_chain.type=param.type[residue]
            lateral_chain.radi=param.radi[residue]
            

            bead_list.append(lateral_chain)
            monomer.part=bead_list

            self.sequence.append(monomer)

class parameters:

    ureg = pint.UnitRegistry()
    ureg.define('sigma = 0.35 * nm = sig')


    one_letter_key={"ALA": "A",
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

# Values for the phosphorilated aminoacids J U and Z are always taken from Bienkiewicz & K.J. Lumb, J Biomol NMR 15: 203-206 (1999).
# Values from Hass MA, Mulder FAA. Contemporary NMR Studies of Protein Electrostatics. Annu Rev Biophys. 2015;44:53-75.

    pka_hass = { "D" : 4.0,
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

# Platzer G, Okon M, McIntosh LP. 2014. pH-dependent random coil 1 H, 13 C, and 15 N chemical shifts of the ionizable amino acids: a guide for protein pK a measurements. J. Biomol. NMR 60:109–29

    pka_platzer = { "D" : 3.86,
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

# Values from Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
            
    pka_crc = { "D" : 3.65,
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

# Y. Nozaki and C. Tanford, Methods Enzymol., 1967, 11, 715–734.

    pka_nozaki = { "D" : 4.00,
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
    
    q =    {"A": {"neutral": 0},
            "R" : {"charged": 1, "neutral": 0},
            "H" : {"charged": 1, "neutral": 0},
            "K" : {"charged": 1, "neutral": 0},
            "D" : {"neutral": 0, "charged": -1},
            "E" : {"neutral": 0, "charged": -1},
            "S" : {"neutral": 0},
            "T" : {"neutral": 0},
            "N" : {"neutral": 0},
            "Q" : {"neutral": 0},
            "C" : {"charged": 1, "neutral": 0},
            "G" : {"neutral": 0},
            "P" : {"neutral": 0},
            "A" : {"neutral": 0},
            "V" : {"neutral": 0},
            "I" : {"neutral": 0},
            "L" : {"neutral": 0},
            "M" : {"neutral": 0},
            "F" : {"neutral": 0},
            "Y" : {"neutral": 0},
            "W" : {"neutral": 0},
            "n" : {"charged": 1, "neutral": 0},
            "c" : {"neutral": 0, "charged": -1},
            "J" : {"neutral": 0, "charged": -1},
            "U" : {"neutral": 0, "charged": -1},
            "Z" : {"neutral": 0, "charged": -1},
            }

    type =    {"A": {"neutral": 21},
            "R" : {"neutral": 22, "charged": 23},
            "H" : {"neutral": 24, "charged": 25},
            "K" : {"neutral": 26, "charged": 27},
            "D" : {"neutral": 28, "charged": 29},
            "E" : {"neutral": 30, "charged": 31},
            "S" : {"neutral": 32},
            "T" : {"neutral": 33},
            "N" : {"neutral": 34},
            "Q" : {"neutral": 35},
            "C" : {"neutral": 36, "charged": 37},
            "G" : {"neutral": 38},
            "P" : {"neutral": 39},
            "A" : {"neutral": 40},
            "V" : {"neutral": 41},
            "I" : {"neutral": 42},
            "L" : {"neutral": 43},
            "M" : {"neutral": 44},
            "F" : {"neutral": 45},
            "Y" : {"neutral": 46},
            "W" : {"neutral": 47},
            "n" : {"neutral": 48, "charged": 49},
            "c" : {"neutral": 50, "charged": 51},
            "J" : {"neutral": 52, "charged": 53},
            "U" : {"neutral": 54, "charged": 55},
            "Z" : {"neutral": 56, "charged": 57},
            }

    radi =    {"A":  0.5*ureg.sigma,
            "R" :  0.5*ureg.sigma,
            "H" :  0.5*ureg.sigma, 
            "K" :  0.5*ureg.sigma, 
            "D" :  0.5*ureg.sigma, 
            "E" :  0.5*ureg.sigma, 
            "S" :  0.5*ureg.sigma,
            "T" :  0.5*ureg.sigma,
            "N" :  0.5*ureg.sigma,
            "Q" :  0.5*ureg.sigma,
            "C" :  0.5*ureg.sigma, 
            "G" :  0.5*ureg.sigma,
            "P" :  0.5*ureg.sigma,
            "A" :  0.5*ureg.sigma,
            "V" :  0.5*ureg.sigma,
            "I" :  0.5*ureg.sigma,
            "L" :  0.5*ureg.sigma,
            "M" :  0.5*ureg.sigma,
            "F" :  0.5*ureg.sigma,
            "Y" :  0.5*ureg.sigma,
            "W" :  0.5*ureg.sigma,
            "n" :  0.5*ureg.sigma, 
            "c" :  0.5*ureg.sigma, 
            "J" :  0.5*ureg.sigma, 
            "U" :  0.5*ureg.sigma, 
            "Z" :  0.5*ureg.sigma
            }

    bondl =    {"A": 1*ureg.sigma,
            "R" :  1*ureg.sigma,
            "H" :  1*ureg.sigma, 
            "K" :  1*ureg.sigma, 
            "D" :  1*ureg.sigma, 
            "E" :  1*ureg.sigma, 
            "S" :  1*ureg.sigma,
            "T" :  1*ureg.sigma,
            "N" :  1*ureg.sigma,
            "Q" :  1*ureg.sigma,
            "C" :  1*ureg.sigma,
            "G" :  1*ureg.sigma,
            "P" :  1*ureg.sigma,
            "A" :  1*ureg.sigma,
            "V" :  1*ureg.sigma,
            "I" :  1*ureg.sigma,
            "L" :  1*ureg.sigma,
            "M" :  1*ureg.sigma,
            "F" :  1*ureg.sigma,
            "Y" :  1*ureg.sigma,
            "W" :  1*ureg.sigma,
            "n" :  1*ureg.sigma,
            "c" :  1*ureg.sigma,
            "J" :  1*ureg.sigma,
            "U" :  1*ureg.sigma,
            "Z" :  1*ureg.sigma
            }
    
    k  =    {"A": 100/ureg.sigma**2,
            "R" :  100/ureg.sigma**2,
            "H" :  100/ureg.sigma**2,
            "K" :  100/ureg.sigma**2, 
            "D" :  100/ureg.sigma**2, 
            "E" :  100/ureg.sigma**2, 
            "S" :  100/ureg.sigma**2,
            "T" :  100/ureg.sigma**2,
            "N" :  100/ureg.sigma**2,
            "Q" :  100/ureg.sigma**2,
            "C" :  100/ureg.sigma**2,
            "G" :  100/ureg.sigma**2,
            "P" :  100/ureg.sigma**2,
            "A" :  100/ureg.sigma**2,
            "V" :  100/ureg.sigma**2,
            "I" :  100/ureg.sigma**2,
            "L" :  100/ureg.sigma**2,
            "M" :  100/ureg.sigma**2,
            "F" :  100/ureg.sigma**2,
            "Y" :  100/ureg.sigma**2,
            "W" :  100/ureg.sigma**2,
            "n" :  100/ureg.sigma**2,
            "c" :  100/ureg.sigma**2,
            "J" :  100/ureg.sigma**2,
            "U" :  100/ureg.sigma**2,
            "Z" :  100/ureg.sigma**2
            }

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

        raise ValueError("Unvalid key for the initial state of the aminoacid, valid options are 'neutral' or 'charged'. Key given: ",  initial_state)

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

        for aminoacid in protein.sequence:

            n_bead=0

            for bead in aminoacid.part:

                bead_id+=1

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
                    name_backbone=aminoacid.name
                    system.part.add(id=[bead_id], pos=pos_backbone, type=[bead.type[state]], q=[bead.q[state]])
                    bead.id=bead_id
                    aminoacid.ids.append(bead_id)
                    ids_chain.append(bead_id)

                else:

                    if n_bead ==0:

                        pos_backbone+=backbone_vec
                        system.part.add(id=[bead_id], pos=pos_backbone, type=[bead.type[state]], q=[bead.q[state]])
                        
                        if (protein.beads_per_monomer == 1):

                            sidechain_bond = interactions.HarmonicBond(k=aminoacid.k.to('sigma**-2').magnitude, r_0=aminoacid.bondl.to('sigma').magnitude)
                            system.bonded_inter.add(sidechain_bond)
                            system.part[bead_id].add_bond((sidechain_bond, id_backbone))

                        else:


                            if (name_backbone =="n") or ( aminoacid.name == "n"):

                                sidechain_bond = interactions.HarmonicBond(k=param.k["n"].to('sigma**-2').magnitude, r_0=param.bondl["n"].to('sigma').magnitude)
                                system.bonded_inter.add(sidechain_bond)
                                system.part[bead_id].add_bond((sidechain_bond, id_backbone))

                            elif (name_backbone =="c") or ( aminoacid.name == "c"):
                                sidechain_bond = interactions.HarmonicBond(k=param.k["c"].to('sigma**-2').magnitude, r_0=param.bondl["c"].to('sigma').magnitude)
                                system.bonded_inter.add(sidechain_bond)
                                system.part[bead_id].add_bond((sidechain_bond, id_backbone))

                            else:

                                system.part[bead_id].add_bond((backbone_bond, id_backbone))
                        
                        bead.id=bead_id
                        aminoacid.ids.append(bead_id)
                        ids_chain.append(bead_id)
                        id_backbone=bead_id
                        name_backbone=aminoacid.name

                    elif n_bead ==1:
                        
                        
                        rand_vec=generate_trialvectors(1)
                        sidechain_vec=np.cross(rand_vec,backbone_vec)
                        r = mt.sqrt(sum(x*x for x in sidechain_vec))
                        sidechain=np.array([x/r for x in sidechain_vec])
                        sidechain_vec=sidechain_vec*aminoacid.bondl.to('sigma').magnitude
                        pos_sidechain=pos_backbone+sidechain_vec
                        system.part.add(id=[bead_id], pos=pos_sidechain, type=[bead.type[state]], q=[bead.q[state]])
                        sidechain_bond = interactions.HarmonicBond(k=aminoacid.k.to('sigma**-2').magnitude, r_0=aminoacid.bondl.to('sigma').magnitude)
                        system.bonded_inter.add(sidechain_bond)
                        system.part[bead_id].add_bond((sidechain_bond, id_backbone))
                        bead.id=bead_id
                        aminoacid.ids.append(bead_id)
                        ids_chain.append(bead_id)
                
                n_bead+=1

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

    for aminoacid in protein.sequence:

        if aminoacid.pKa is not None:

            N_ti+=1

    return N_ti

def setup_protein_acidbase_reactions(RE, protein, cation):
    """
    Set up the Acid/Base reactions for acidic/basidic aminoacids in protein. The reaction steps are done following the constant pH ensamble procedure. 

    Inputs:
    RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
    protein: class object as defined in this library, attributes:
    cation: class objects of the cations,  with the following attributes:
        q: (int) charge of the ion
        type: (int) type of the ion
    """

    param=parameters()
    reaction_absent={}

    for group in param.pka_hass.keys();
        
        reaction_absent[group]=True

    for aminoacid in protein.sequence:

        if aminoacid.pKa is not None:

            if (reaction_absent[aminoacid.name]):

                if (aminoacid.q["charged"] == 1) : # Basic aminoacid

                    RE.add_reaction(gamma=10**-aminoacid.pKa,
                            reactant_types=aminoacid.type["charged"],
                            reactant_coefficients=[1],
                            product_types=[aminoacid.type["neutral"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={aminoacid.type["neutral"]: aminoacid.q["neutral"],
                                             aminoacid.type["charged"]: aminoacid.q["charged"],
                                             cation.type: cation.q})
                    reaction_absent[aminoacid.name] = False

                elif (aminoacid.q["charged"] == -1) : # Acid aminoacid


                    RE.add_reaction(gamma=10**-aminoacid.pka,
                            reactant_types=[aminoacid.type["neutral"]],
                            reactant_coefficients=[1],
                            product_types=[aminoacid.type["charged"], cation.type],
                            product_coefficients=[1, 1],
                            default_charges={aminoacid.type["neutral"]: aminoacid.q["neutral"],
                                             aminoacid.type["charged"]: aminoacid.q["charged"],
                                             cation.type: cation.q})
                    reaction_absent[aminoacid.name] = False

                else:
    
                    raise ValueError("This subrutine is concived for the acid/base equilibria of monovalent ions. Charge of aminoacid ", aminoacid.name, " = ", aminoacid.q["charged"])

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

    for id in protein.ids:

        Z_prot+=system.part[id].q

    Z_prot=Z_prot/protein.N
    Z2=Z_prot**2

    return Z_prot, Z2

def track_ionization(system, protein):
    """
    Sets up espresso to track the average number of particles of each aminoacid type
    
    Inputs:
    system: espresso class object with all system variables
    protein: class object as defined in this library

    """

    types=[]

    for aminoacid in protein.sequence:
        
        for type in aminoacid.type.values():

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

        for aminoacid in protein.sequence:

            if aminoacid.pKa is not None:

                Z+=aminoacid.q["charged"]/(1+10**(aminoacid.q["charged"]*(pH_value-protein.pka[group])))

        Z_HH.append(Z)

    return Z_HH


def create_counterions(system,cation,anion):
    """
    Adds one counter-ion (cation or anion) of opposite charge per each charge in the system.

    Inputs:

    system: espresso class object with all system variables.
    cation/anion: class objects of the ions,  with the following attributes:
        q: (int) charge of the ion
        type: (int) type of the ion

    Outputs:
    cation_ids: (list) contains the ids of the cations created
    anion_ids: (list) contains the ids of the anions created

    Assumptions:

    -The charges in system are monovalent
    -The counter-ions are monovalent

    """

    cation_ids=[]
    anion_ids=[]

    if any( abs(q) > 1 for q in system.part[:].q):

        print("ERROR: the subrutine for creating counter-ions is not prepared for multi-valent beads")
        exit()

    if (anion.q != -1) or (cation.q != 1):
        print("ERROR: the charge of the cation/anion has to be monovalent and equal to +1/-1")
        print("       Charge anion: ", anion.q ," Charge cation: ", cation.q )
        exit()

    N_cation=0
    N_anion=0

    for id in system.part[:].id:

        if (system.part[id].q == 1): # Add an anion to the system (counter-ion)

            I_id=max(system.part[:].id)+1
            I_pos=np.random.random((1, 3)) * system.box_l
            system.part.add(id=[I_id], pos=I_pos, type=[anion.type], q=[anion.q])
            N_anion+=1
            anion_ids.append(I_id)

        if (system.part[id].q == -1): # Add an anion to the system (counter-ion)

            I_id=max(system.part[:].id)+1
            I_pos=np.random.random((1, 3)) * system.box_l
            system.part.add(id=[I_id], pos=I_pos, type=[cation.type], q=[cation.q])
            N_cation+=1
            cation_ids.append(I_id)

    if (abs(sum(system.part[:].q)) > 1e-10):

        print("ERROR: System not electroneutral, consider revising the set-up unless working on ideal conditions ")
        print("Global charge: ", sum(system.part[:].q))
        exit()

    cation.N=N_cation
    anion.N=N_anion

    return cation_ids, anion_ids

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

