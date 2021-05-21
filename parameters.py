import pint
ureg = pint.UnitRegistry()

class parameter_set:
    name=None
    type=None
    bondl=None
    radius=None
    k=None
    eps=None
    pKa=None
    beads_per_monomer=None

class general:
    name=None
    aminoacid={"ALA": "A",
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

class small_ions:
    name='small_ions'
    
    class cation:
        name="cation"
        type=18
        radius=0.15 * ureg.nm
        q=1

    class anion:
        name="anion"
        type=19
        radius=0.15 * ureg.nm
        q=-1


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


class one_bead_peptide:

    name='1beadprotein'
    beads_per_monomer=1
    
    class Alanine:
        name="A"
        radius=0.18*ureg.nm,
        type=21
        q=0

    class Aspargine:
        name="N"
        radius=0.18*ureg.nm,
        type=22
        q=0
    
    class Glutamine:
        name="Q"
        radius=0.18*ureg.nm,
        type=23
        q=0
    
    class Glycine:
        name="G"
        radius=0.18*ureg.nm,
        type=24
        q=0
    
    class Isoleucine:
        name="I"
        radius=0.18*ureg.nm,
        type=25
        q=0
    
    class Leucine:
        name="L"
        radius=0.18*ureg.nm,
        type=26
        q=0
    
    class Methionine:
        name="M"
        radius=0.18*ureg.nm,
        type=27
        q=0
    
    class Phenylalanine:
        name="F"
        radius=0.18*ureg.nm,
        type=28
        q=0
    
    class Proline:
        name="P"
        radius=0.18*ureg.nm,
        type=29
        q=0
    
    class Serine:
        name="S"
        radius=0.18*ureg.nm,
        type=30
        q=0
    
    class Threonine:
        name="T"
        radius=0.18*ureg.nm,
        type=31
        q=0
    
    class Tryptophan:
        name="W"
        radius=0.18*ureg.nm,
        type=32
        q=0
    
    class Tyrosine:
        name="Y"
        radius=0.18*ureg.nm,
        type=33
        q=0
    
    class Valine:
        name="V"
        radius=0.18*ureg.nm,
        type=34
        q=0
    
    class Lysine:
        name="K"
        radius=0.18*ureg.nm,
        type={"neutral": 35, "charged": 36}
        q={"neutral": 0, "charged": 1}
    
    class Histidine:
        name="H"
        radius=0.18*ureg.nm,
        type={"neutral": 37, "charged": 38}
        q={"neutral": 0, "charged": 1}
    
    class Cysteine:
        name="C"
        radius=0.18*ureg.nm,
        type={"neutral": 39, "charged": 40}
        q={"neutral": 0, "charged": 1}
    
    class Arginine:
        name="R"
        radius=0.18*ureg.nm,
        type={"neutral": 41, "charged": 42}
        q={"neutral": 0, "charged": 1}
    
    class Amino_end:
        name="n"
        radius=0.18*ureg.nm,
        type={"neutral": 43, "charged": 44}
        q={"neutral": 0, "charged": 1}

    class Carboxyl_end:
        name="c"
        radius=0.18*ureg.nm,
        type={"neutral": 45, "charged": 46}
        q={"neutral": 0, "charged": -1}

    class Aspartic:
        name="D"
        radius=0.18*ureg.nm,
        type={"neutral": 47, "charged": 48}
        q={"neutral": 0, "charged": -1}

    class Glutamic:
        name="E"
        radius=0.18*ureg.nm,
        type={"neutral": 49, "charged": 50}
        q={"neutral": 0, "charged": -1}
    
    class Phosphorilated_Serine:
        name="J"
        radius=0.18*ureg.nm,
        type={"neutral": 51, "charged": 52}
        q={"neutral": 0, "charged": -1}
    
    class Phosphorilated_Threonine:
        name="U"
        radius=0.18*ureg.nm,
        type={"neutral": 53, "charged": 54}
        q={"neutral": 0, "charged": -1}
    
    class Phosphorilated_Tryptophan:
        name="Z"
        radius=0.18*ureg.nm,
        type={"neutral": 55, "charged": 56}
        q={"neutral": 0, "charged": -1}

class two_bead_peptide:

    name='2beadpeptide'
    beads_per_monomer=2
    
    class Alpha_Carbon:
        name="C_alpha"
        radius=0.18*ureg.nm,
        type=20
        q=0
    
    class Alanine:
        name="A"
        radius=0.18*ureg.nm,
        type=21
        q=0

    class Aspargine:
        name="N"
        radius=0.18*ureg.nm,
        type=22
        q=0
    
    class Glutamine:
        name="Q"
        radius=0.18*ureg.nm,
        type=23
        q=0
    
    class Glycine:
        name="G"
        radius=0.18*ureg.nm,
        type=24
        q=0
    
    class Isoleucine:
        name="I"
        radius=0.18*ureg.nm,
        type=25
        q=0
    
    class Leucine:
        name="L"
        radius=0.18*ureg.nm,
        type=26
        q=0
    
    class Methionine:
        name="M"
        radius=0.18*ureg.nm,
        type=27
        q=0
    
    class Phenylalanine:
        name="F"
        radius=0.18*ureg.nm,
        type=28
        q=0
    
    class Proline:
        name="P"
        radius=0.18*ureg.nm,
        type=29
        q=0
    
    class Serine:
        name="S"
        radius=0.18*ureg.nm,
        type=30
        q=0
    
    class Threonine:
        name="T"
        radius=0.18*ureg.nm,
        type=31
        q=0
    
    class Tryptophan:
        name="W"
        radius=0.18*ureg.nm,
        type=32
        q=0
    
    class Tyrosine:
        name="Y"
        radius=0.18*ureg.nm,
        type=33
        q=0
    
    class Valine:
        name="V"
        radius=0.18*ureg.nm,
        type=34
        q=0
    
    class Lysine:
        name="K"
        radius=0.18*ureg.nm,
        type={"neutral": 35, "charged": 36}
        q={"neutral": 0, "charged": 1}
    
    class Histidine:
        name="H"
        radius=0.18*ureg.nm,
        type={"neutral": 37, "charged": 38}
        q={"neutral": 0, "charged": 1}
    
    class Cysteine:
        name="C"
        radius=0.18*ureg.nm,
        type={"neutral": 39, "charged": 40}
        q={"neutral": 0, "charged": 1}
    
    class Arginine:
        name="R"
        radius=0.18*ureg.nm,
        type={"neutral": 41, "charged": 42}
        q={"neutral": 0, "charged": 1}
    
    class Amino_end:
        name="n"
        radius=0.18*ureg.nm,
        type={"neutral": 43, "charged": 44}
        q={"neutral": 0, "charged": 1}

    class Carboxyl_end:
        name="c"
        radius=0.18*ureg.nm,
        type={"neutral": 45, "charged": 46}
        q={"neutral": 0, "charged": -1}

    class Aspartic:
        name="D"
        radius=0.18*ureg.nm,
        type={"neutral": 47, "charged": 48}
        q={"neutral": 0, "charged": -1}

    class Glutamic:
        name="E"
        radius=0.18*ureg.nm,
        type={"neutral": 49, "charged": 50}
        q={"neutral": 0, "charged": -1}
    
    class Phosphorilated_Serine:
        name="J"
        radius=0.18*ureg.nm,
        type={"neutral": 51, "charged": 52}
        q={"neutral": 0, "charged": -1}
    
    class Phosphorilated_Threonine:
        name="U"
        radius=0.18*ureg.nm,
        type={"neutral": 53, "charged": 54}
        q={"neutral": 0, "charged": -1}
    
    class Phosphorilated_Tryptophan:
        name="Z"
        radius=0.18*ureg.nm,
        type={"neutral": 55, "charged": 56}
        q={"neutral": 0, "charged": -1}



def get_subclasses(cls):
    
    import inspect    
    model_list=[]

    for attribute in cls.__dict__.items():
       
        name_attribute=attribute[0]

        if name_attribute[0:2] != "__" and inspect.isclass(attribute[1]):
            
            model_list.append(attribute[1])

    return model_list

