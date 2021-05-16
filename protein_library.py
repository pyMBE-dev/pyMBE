import numpy as np
import pint

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
    pKa=None
    part=[]
    ids=[]


class protein():

    N = None
    Nm = None
    ids = []
    
    def __init__(self, sequence, beads_per_monomer=1, pKa_set=None, pKa_custom=None):

        clean_sequence=[]
        self.sequence=[]
        param=parameters()
        self.bondl=0.388*param.ureg.nm
    
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
           
            if residue in pKa.keys():

                monomer.pKa=pKa[residue]
            
            if beads_per_monomer == 2: # Create an alpha carbon bead

                alpha_carbon=particle()
                alpha_carbon.name="C_alpha"
                alpha_carbon.q={"neutral": 0}
                alpha_carbon.type={"neutral": 20}
                alpha_carbon.radi={"neutral": 0.5*param.ureg.sigma}
                monomer.part.append(alpha_carbon)

            # Create the lateral chain bead

            lateral_chain=particle()
            lateral_chain.name=residue
            lateral_chain.q=param.q[residue]
            lateral_chain.type=param.type[residue]
            lateral_chain.radi=param.radi[residue]
            
            if beads_per_monomer == 2: # distance to the alpha carbon

                lateral_chain.bondl=param.bondl[residue]

            monomer.part.append(lateral_chain)
            
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
            "R" : {"protonated": 1, "unprotonated": 0},
            "H" : {"protonated": 1, "unprotonated": 0},
            "K" : {"protonated": 1, "unprotonated": 0},
            "D" : {"protonated": 0, "unprotonated": -1},
            "E" : {"protonated": 0, "unprotonated": -1},
            "S" : {"neutral": 0},
            "T" : {"neutral": 0},
            "N" : {"neutral": 0},
            "Q" : {"neutral": 0},
            "C" : {"protonated": 1, "unprotonated": 0},
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
            "n" : {"protonated": 1, "unprotonated": 0},
            "c" : {"protonated": 0, "unprotonated": -1},
            "J" : {"protonated": 0, "unprotonated": -1},
            "U" : {"protonated": 0, "unprotonated": -1},
            "Z" : {"protonated": 0, "unprotonated": -1},
            }

    type =    {"A": {"neutral": 21},
            "R" : {"protonated": 22, "unprotonated": 23},
            "H" : {"protonated": 24, "unprotonated": 25},
            "K" : {"protonated": 26, "unprotonated": 27},
            "D" : {"protonated": 28, "unprotonated": 29},
            "E" : {"protonated": 30, "unprotonated": 31},
            "S" : {"neutral": 32},
            "T" : {"neutral": 33},
            "N" : {"neutral": 34},
            "Q" : {"neutral": 35},
            "C" : {"protonated": 36, "unprotonated": 37},
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
            "n" : {"protonated": 48, "unprotonated": 49},
            "c" : {"protonated": 50, "unprotonated": 51},
            "J" : {"protonated": 52, "unprotonated": 53},
            "U" : {"protonated": 54, "unprotonated": 55},
            "Z" : {"protonated": 56, "unprotonated": 57},
            }

    radi =    {"A": {"neutral": 0.5*ureg.sigma},
            "R" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "H" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "K" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "D" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "E" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "S" : {"neutral": 0.5*ureg.sigma},
            "T" : {"neutral": 0.5*ureg.sigma},
            "N" : {"neutral": 0.5*ureg.sigma},
            "Q" : {"neutral": 0.5*ureg.sigma},
            "C" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "G" : {"neutral": 0.5*ureg.sigma},
            "P" : {"neutral": 0.5*ureg.sigma},
            "A" : {"neutral": 0.5*ureg.sigma},
            "V" : {"neutral": 0.5*ureg.sigma},
            "I" : {"neutral": 0.5*ureg.sigma},
            "L" : {"neutral": 0.5*ureg.sigma},
            "M" : {"neutral": 0.5*ureg.sigma},
            "F" : {"neutral": 0.5*ureg.sigma},
            "Y" : {"neutral": 0.5*ureg.sigma},
            "W" : {"neutral": 0.5*ureg.sigma},
            "n" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "c" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "J" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "U" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            "Z" : {"protonated": 0.5*ureg.sigma, "unprotonated": 0.5*ureg.sigma},
            }

    bondl =    {"A": {"neutral": 1*ureg.sigma},
            "R" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "H" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "K" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "D" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "E" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "S" : {"neutral": 1*ureg.sigma},
            "T" : {"neutral": 1*ureg.sigma},
            "N" : {"neutral": 1*ureg.sigma},
            "Q" : {"neutral": 1*ureg.sigma},
            "C" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "G" : {"neutral": 1*ureg.sigma},
            "P" : {"neutral": 1*ureg.sigma},
            "A" : {"neutral": 1*ureg.sigma},
            "V" : {"neutral": 1*ureg.sigma},
            "I" : {"neutral": 1*ureg.sigma},
            "L" : {"neutral": 1*ureg.sigma},
            "M" : {"neutral": 1*ureg.sigma},
            "F" : {"neutral": 1*ureg.sigma},
            "Y" : {"neutral": 1*ureg.sigma},
            "W" : {"neutral": 1*ureg.sigma},
            "n" : {"protonated":1*ureg.sigma, "unprotonated":  1*ureg.sigma},
            "c" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "J" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "U" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            "Z" : {"protonated": 1*ureg.sigma, "unprotonated": 1*ureg.sigma},
            }

def create_protein(system, protein, harmonic_bond):
    '''
    Creates a protein chain in the system with the sequence stored in protein, joined by harmonic_bond)
    Inputs:
    system: espresso class object with all system variables
    protein: class object as defined in this library, attributes:
        .sequence: labels of the whole sequence of the protein
        .type: dictionary containing the default types of each group
        .q: dictionary containing the charge of each group
        .Nm: number of monomers of the protein
        .N:number of protein chain
    harmonic_bond: espresso class containing the harmonic bond object

    Outputs:
    protein_ids: (list) contains all the ids of the protein particles

    '''


    # Compute the number of particles needed for each protein chain and store it

    protein.Nm=len(protein.sequence)

    # Create polymer chains in the system

    polymer=molecule(radi=protein.radi['GLY'], q=protein.q['GLY'], type=protein.type['GLY'])

    polymer.Nm=protein.Nm
    polymer.N=protein.N

    protein_ids=create_linear_polymer(system,polymer,harmonic_bond)

    # Transform the polymer chains into peptide chains by asigning the proper parameters

    for protein_chain in protein_ids:

        for index in range(protein.Nm):

            label=protein.sequence[index]
            id_group=protein_chain[index]
            q_group=protein.q[label]
            type_group=protein.type[label]

            system.part[id_group].q=q_group
            system.part[id_group].type=type_group

    return protein_ids


def check_protein_sequence(protein):
    """
    Checks that the protein sequence contained in 'protein' contains labels understandable for the library. Corrects capital letter missprints

    Input:
    protein: object by the protein class defined above
    """
    
    one_letter_key={"A": "ALA",
                    "R": "ARG",
                    "N": "ASN",
                    "D": "ASP",
                    "C": "CYS",
                    "E": "GLU",
                    "Q": "GLN",
                    "G": "GLY",
                    "H": "HIS",
                    "I": "ILE",
                    "L": "LEU",
                    "K": "LYS",
                    "M": "MET",
                    "F": "PHE",
                    "P": "PRO",
                    "S": "SER",
                    "T": "THR",
                    "W": "TRP",
                    "Y": "TYR",
                    "V": "VAL"}

    for ind_group in range(len(protein.sequence)):

        group=protein.sequence[ind_group]

        if not isinstance(group, str):

            raise ValueError("All labels in protein sequence must be strings, ", group, " is not a string")

        if not group.isupper():

            upper_label=group.upper()

        else:

            upper_label=group

        if upper_label not in protein.type.keys():

            if upper_label in one_letter_key.keys():

                correct_label= one_letter_key[upper_label]

            else:

                raise ValueError("Unknown label ", group)

        else:

            correct_label=upper_label

        protein.sequence[ind_group]=correct_label

    return

def count_titrable_groups(protein):
    """
    Counts the number of titrable groups in the protein object

    Input:
    protein: class object as defined in this library, attributes:
        .sequence: labels of the whole sequence of the protein

    Output:
    N_ti: (int) number of titrable groups in the sequence
    """

    N_ti=0

    for group in protein.sequence:

        if group in protein.pka.keys():

            N_ti+=1

    return N_ti

def setup_protein_acidbase_reactions(RE, protein, cation):
    """
    Set up the Acid/Base reactions for acidic/basidic aminoacids in protein. The reaction steps are done following the constant pH ensamble procedure. 

    Inputs:
    RE: instance of the espresso class reaction_ensemble.ConstantpHEnsemble
    protein: class object as defined in this library, attributes:
        .sequence: labels of the whole sequence of the protein
        .type: dictionary containing the default types of each group
        .q: dictionary containing the charge of each group
    cation: class objects of the counter-ions,  with the following attributes:
        q: (int) charge of the ion
        type: (int) type of the ion
    """

    reaction_present={ "ASP" : False,
            "GLU" : False,
            "HIS" : False,
            "TYR" : False,
            "LYS" : False,
            "ARG" : False,
            "CYS" : False,
            "NH2" : False,
            "COOH": False}

    for group in protein.sequence:

        if (group == "ASP" or group == "NASP") and (reaction_present["ASP"] == False):

            RE.add_reaction(gamma=10**-protein.pka["ASP"],
                            reactant_types=[protein.type["NASP"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["ASP"], cation.type],
                            product_coefficients=[1, 1],
                            default_charges={protein.type["NASP"]: protein.q["NASP"],
                                             protein.type["ASP"]: protein.q["ASP"],
                                             cation.type: cation.q})
            reaction_present["ASP"] = True

        elif (group == "GLU" or group == "NGLU") and (reaction_present["GLU"] == False):

            RE.add_reaction(gamma=10**-protein.pka["GLU"],
                            reactant_types=[protein.type["NGLU"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["GLU"], cation.type],
                            product_coefficients=[1, 1],
                            default_charges={protein.type["NGLU"]: protein.q["NGLU"],
                                             protein.type["GLU"]: protein.q["GLU"],
                                             cation.type: cation.q})
            reaction_present["GLU"] = True

        elif (group == "COOH" or group == "NCOOH") and (reaction_present["COOH"] == False):

            RE.add_reaction(gamma=10**-protein.pka["COOH"],
                            reactant_types=[protein.type["NCOOH"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["COOH"], cation.type],
                            product_coefficients=[1, 1],
                            default_charges={protein.type["NCOOH"]: protein.q["NCOOH"],
                                             protein.type["COOH"]: protein.q["COOH"],
                                             cation.type: cation.q})
            reaction_present["GLU"] = True

        elif (group == "HIS" or group == "NHIS") and (reaction_present["HIS"] == False):

            RE.add_reaction(gamma=10**-protein.pka["HIS"],
                            reactant_types=[protein.type["HIS"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["NHIS"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={protein.type["NHIS"]: protein.q["NHIS"],
                                             protein.type["HIS"]: protein.q["HIS"],
                                             cation.type: cation.q})
            reaction_present["HIS"] = True

        elif (group == "TYR" or group == "NTYR") and (reaction_present["TYR"] == False):

            RE.add_reaction(gamma=10**-protein.pka["TYR"],
                            reactant_types=[protein.type["TYR"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["NTYR"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={protein.type["NTYR"]: protein.q["NTYR"],
                                             protein.type["TYR"]: protein.q["TYR"],
                                             cation.type: cation.q})
            reaction_present["TYR"] = True

        elif (group == "LYS" or group == "NLYS") and (reaction_present["LYS"] == False):

            RE.add_reaction(gamma=10**-protein.pka["LYS"],
                            reactant_types=[protein.type["LYS"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["NLYS"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={protein.type["NLYS"]: protein.q["NLYS"],
                                             protein.type["LYS"]: protein.q["LYS"],
                                             cation.type: cation.q})
            reaction_present["LYS"] = True

        elif (group == "ARG" or group == "NARG") and (reaction_present["ARG"] == False):

            RE.add_reaction(gamma=10**-protein.pka["ARG"],
                            reactant_types=[protein.type["ARG"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["NARG"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={protein.type["NARG"]: protein.q["NARG"],
                                             protein.type["ARG"]: protein.q["ARG"],
                                             cation.type: cation.q})
            reaction_present["ARG"] = True

        elif (group == "CYS" or group == "NCYS") and (reaction_present["CYS"] == False):

            RE.add_reaction(gamma=10**-protein.pka["CYS"],
                            reactant_types=[protein.type["CYS"], anion.type],
                            reactant_coefficients=[1,1],
                            product_types=[protein.type["NCYS"]],
                            product_coefficients=[1],
                            default_charges={protein.type["NCYS"]: protein.q["NCYS"],
                                             protein.type["CYS"]: protein.q["CYS"],
                                             cation.type: cation.q})
            reaction_present["CYS"] = True

        elif (group == "NH2" or group == "NH2") and (reaction_present["NH2"] == False):

            RE.add_reaction(gamma=10**-protein.pka["NH2"],
                            reactant_types=[protein.type["NH2"]],
                            reactant_coefficients=[1],
                            product_types=[protein.type["NNH2"], cation.type],
                            product_coefficients=[1,1],
                            default_charges={protein.type["NNH2"]: protein.q["NNH2"],
                                             protein.type["NH2"]: protein.q["NH2"],
                                             cation.type: cation.q})
            reaction_present["NH2"] = True

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
    checked_groups=[]

    for group in protein.sequence:

        if group not in checked_groups:

            N_group=system.number_of_particles(type=protein.type[group])
            Z_prot+=N_group*protein.q[group]
            checked_groups.append(group)

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

    for group in protein.sequence:

        types.append(protein.type[group])

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

        for group in protein.sequence:

            if group in protein.pka.keys():

                if protein.q[group] == 1: # Basic group

                    psi=1

                elif protein.q[group] == -1: # Acidic group
                    
                    psi=-1

                else:

                    raise ValueError("Unknown label ", group) 

                Z+=psi/(1+10**(psi*(pH_value-protein.pka[group])))
        
        Z_HH.append(Z)

    return Z_HH

def create_linear_polymer(system,linear_polymer,harmonic_bond, inpt_pos=None, inpt_vec=None):
    """
    Creates a linear polymer in the system. 
    
    Inputs:
    system: espresso class object with all system variables
    linear_polymer: class variable defined above as "molecule"
    harmonic_bond: espresso class containing the harmonic bond object
    inpt_pos: (float*3) [x y z] position of the first monomer in the simulation box, by default is placed at random
    inpt_vec: determines the growing direction of the chain, by default random

    Output:
    polymer_ids:(int list) list of all the ids of the polymer particles

    Assumptions:
    -All monomers  are equal
    -All bonds are equal
    """

    # Parameters consistency check

    bondl=harmonic_bond.params['r_0']

    if ( bondl <= 0 ):

        print("ERROR: the bond length must be a positive number bondl = ", bondl)
        exit()

    elif(linear_polymer.Nm <= 1) :

        print("WARNING: The number of monomers per chain must be at least 2 to generate a linear polymer!")
        print("Number of monomers: ", linear_polymer.Nm)

        if(linear_polymer.Nm <= 0):

            print("ERROR: the number of monomers must be larger than 0!")
            exit()

    elif(linear_polymer.N < 1):

            print("ERROR: the number of chains must be 1 or more!")
            exit()

    all_ids=[]

    for chain in range(linear_polymer.N):

        if inpt_pos is None: # By default, the first monomer of the linear polymer is placed in a random position

            ini_pos=np.random.random((1, 3)) * system.box_l

        else:

            ini_pos=[inpt_pos]

        if ( (system.part[:].id).size == 0):

            M_id=-1

        else:

            M_id=max(system.part[:].id)

    # Generate the growing vector for the polymer (by default, random)

        if inpt_vec is None:

            grow_vec=generate_trialvectors(bondl)
            grow_vec=[grow_vec]

        else:

            grow_vec=[inpt_vec]

        polymer_ids=[]

        for mon in range(linear_polymer.Nm):

            M_id=M_id+1

            if not polymer_ids: # First monomer

                system.part.add(id=[M_id], pos=ini_pos, type=[linear_polymer.type], q=[linear_polymer.q])
            else: # Rest of the monomers 

                prev_id=polymer_ids[-1]
                M_pos=system.part[prev_id].pos
                M_pos=M_pos+grow_vec
                system.part.add(id=[M_id], pos=M_pos, type=[linear_polymer.type], q=[linear_polymer.q])
                system.part[prev_id].add_bond((harmonic_bond, M_id))

            polymer_ids.append(M_id)

        all_ids.append(polymer_ids)

    return all_ids

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

