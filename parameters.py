import pint
ureg = pint.UnitRegistry()


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



class custom_parameters:
    name="custom"
    beads_per_residue=None
    principal_chain=None
    side_chain=None

    def newAttr(self, attr):

        setattr(self, attr.name, attr)


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
    beads_per_residue=1
    principal_chain='sequence'
    side_chain=None
    

    cation=particle()
    cation.name="cation"
    cation.type=18
    cation.radius=0.15 * ureg.nm
    cation.q=1

    anion=particle()
    anion.name="anion"
    anion.type=19
    anion.radius=0.15 * ureg.nm
    anion.q=-1


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

    name='1beadpeptide'
    beads_per_residue=1
    principal_chain='sequence'
    side_chain=None

    Alanine=particle()
    Alanine.name="A"
    Alanine.acidity="inert"
    Alanine.radius=0.18*ureg.nm,
    Alanine.type=21
    Alanine.q=0

    Aspargine=particle()
    Aspargine.name="N"
    Aspargine.radius=0.18*ureg.nm,
    Aspargine.type=22
    Aspargine.q=0
    Aspargine.acidity="inert"
    
    Glutamine=particle()
    Glutamine.name="Q"
    Glutamine.radius=0.18*ureg.nm,
    Glutamine.type=23
    Glutamine.q=0
    Glutamine.acidity="inert"
    
    Glycine=particle()
    Glycine.name="G"
    Glycine.radius=0.18*ureg.nm,
    Glycine.type=24
    Glycine.q=0
    Glycine.acidity="inert"
    
    Isoleucine=particle()
    Isoleucine.name="I"
    Isoleucine.radius=0.18*ureg.nm,
    Isoleucine.type=25
    Isoleucine.q=0
    Isoleucine.acidity="inert"
    
    Leucine=particle()
    Leucine.name="L"
    Leucine.radius=0.18*ureg.nm,
    Leucine.type=26
    Leucine.q=0
    Leucine.acidity="inert"
    
    Methionine=particle()
    Methionine.name="M"
    Methionine.radius=0.18*ureg.nm,
    Methionine.type=27
    Methionine.q=0
    Methionine.acidity="inert"
    
    Phenylalanine=particle()
    Phenylalanine.name="F"
    Phenylalanine.radius=0.18*ureg.nm,
    Phenylalanine.type=28
    Phenylalanine.q=0
    Phenylalanine.acidity="inert"
    
    Proline=particle()
    Proline.name="P"
    Proline.radius=0.18*ureg.nm,
    Proline.type=29
    Proline.q=0
    Proline.acidity="inert"
    
    Serine=particle()
    Serine.name="S"
    Serine.radius=0.18*ureg.nm,
    Serine.type=30
    Serine.q=0
    Serine.acidity="inert"
    
    Threonine=particle()
    Threonine.name="T"
    Threonine.radius=0.18*ureg.nm,
    Threonine.type=31
    Threonine.q=0
    Threonine.acidity="inert"
    
    Tryptophan=particle()
    Tryptophan.name="W"
    Tryptophan.radius=0.18*ureg.nm,
    Tryptophan.type=32
    Tryptophan.q=0
    Tryptophan.acidity="inert"
    
    Tyrosine=particle()
    Tyrosine.name="Y"
    Tyrosine.radius=0.18*ureg.nm,
    Tyrosine.type=33
    Tyrosine.q=0
    Tyrosine.acidity="inert"
    
    Valine=particle()
    Valine.name="V"
    Valine.radius=0.18*ureg.nm,
    Valine.type=34
    Valine.q=0
    Valine.acidity="inert"
    
    Lysine=particle()
    Lysine.name="K"
    Lysine.radius=0.18*ureg.nm,
    Lysine.type={"unprotonated": 35, "protonated": 36}
    Lysine.q={"unprotonated": 0, "protonated": 1}
    Lysine.acidity="basic"
    
    Histidine=particle()
    Histidine.name="H"
    Histidine.radius=0.18*ureg.nm,
    Histidine.type={"unprotonated": 37, "protonated": 38}
    Histidine.q={"unprotonated": 0, "protonated": 1}
    Histidine.acidity="basic"
    
    Cysteine=particle()
    Cysteine.name="C"
    Cysteine.radius=0.18*ureg.nm,
    Cysteine.type={"unprotonated": 39, "protonated": 40}
    Cysteine.q={"unprotonated": 0, "protonated": 1}
    Cysteine.acidity="basic"
    
    Arginine=particle()
    Arginine.name="R"
    Arginine.radius=0.18*ureg.nm,
    Arginine.type={"unprotonated": 41, "protonated": 42}
    Arginine.q={"unprotonated": 0, "protonated": 1}
    Arginine.acidity="basic"
    
    Amino_end=particle()
    Amino_end.name="n"
    Amino_end.radius=0.18*ureg.nm,
    Amino_end.type={"unprotonated": 43, "protonated": 44}
    Amino_end.q={"unprotonated": 0, "protonated": 1}
    Amino_end.acidity="basic"

    Carboxyl_end=particle()
    Carboxyl_end.name="c"
    Carboxyl_end.radius=0.18*ureg.nm,
    Carboxyl_end.type={"unprotonated": 45, "protonated": 46}
    Carboxyl_end.q={"protonated": 0, "unprotonated": -1}
    Carboxyl_end.acidity="acid"

    Aspartic=particle()
    Aspartic.name="D"
    Aspartic.radius=0.18*ureg.nm,
    Aspartic.type={"unprotonated": 47, "protonated": 48}
    Aspartic.q={"protonated": 0, "unprotonated": -1}
    Aspartic.acidity="acid"

    Glutamic=particle()
    Glutamic.name="E"
    Glutamic.radius=0.18*ureg.nm,
    Glutamic.type={"unprotonated": 49, "protonated": 50}
    Glutamic.q={"protonated": 0, "unprotonated": -1}
    Glutamic.acidity="acid"
    
    Phosphorilated_Serine=particle()
    Phosphorilated_Serine.name="J"
    Phosphorilated_Serine.radius=0.18*ureg.nm,
    Phosphorilated_Serine.type={"unprotonated": 51, "protonated": 52}
    Phosphorilated_Serine.q={"protonated": -1, "unprotonated": -2}
    Phosphorilated_Serine.acidity="acid"
    
    Phosphorilated_Threonine=particle()
    Phosphorilated_Threonine.name="U"
    Phosphorilated_Threonine.radius=0.18*ureg.nm,
    Phosphorilated_Threonine.type={"unprotonated": 53, "protonated": 54}
    Phosphorilated_Threonine.q={"protonated": -1, "unprotonated": -2}
    Phosphorilated_Threonine.acidity="acid"
    
    Phosphorilated_Tryptophan=particle()
    Phosphorilated_Tryptophan.name="Z"
    Phosphorilated_Tryptophan.radius=0.18*ureg.nm,
    Phosphorilated_Tryptophan.type={"unprotonated": 55, "protonated": 56}
    Phosphorilated_Tryptophan.q={"protonated": -1, "unprotonated": -2}
    Phosphorilated_Tryptophan.acidity="acid"

class two_bead_peptide:

    name='2beadpeptide'
    beads_per_residue=2
    principal_chain='C_alpha'
    side_chain=['sequence']
    
    Alpha_Carbon=particle()
    Alpha_Carbon.name="C_alpha"
    Alpha_Carbon.radius=0.18*ureg.nm,
    Alpha_Carbon.type=20
    Alpha_Carbon.q=0
    
    Alanine=particle()
    Alanine.name="A"
    Alanine.acidity="inert"
    Alanine.radius=0.18*ureg.nm,
    Alanine.type=21
    Alanine.q=0

    Aspargine=particle()
    Aspargine.name="N"
    Aspargine.radius=0.18*ureg.nm,
    Aspargine.type=22
    Aspargine.q=0
    Aspargine.acidity="inert"
    
    Glutamine=particle()
    Glutamine.name="Q"
    Glutamine.radius=0.18*ureg.nm,
    Glutamine.type=23
    Glutamine.q=0
    Glutamine.acidity="inert"
    
    Glycine=particle()
    Glycine.name="G"
    Glycine.radius=0.18*ureg.nm,
    Glycine.type=24
    Glycine.q=0
    Glycine.acidity="inert"
    
    Isoleucine=particle()
    Isoleucine.name="I"
    Isoleucine.radius=0.18*ureg.nm,
    Isoleucine.type=25
    Isoleucine.q=0
    Isoleucine.acidity="inert"
    
    Leucine=particle()
    Leucine.name="L"
    Leucine.radius=0.18*ureg.nm,
    Leucine.type=26
    Leucine.q=0
    Leucine.acidity="inert"
    
    Methionine=particle()
    Methionine.name="M"
    Methionine.radius=0.18*ureg.nm,
    Methionine.type=27
    Methionine.q=0
    Methionine.acidity="inert"
    
    Phenylalanine=particle()
    Phenylalanine.name="F"
    Phenylalanine.radius=0.18*ureg.nm,
    Phenylalanine.type=28
    Phenylalanine.q=0
    Phenylalanine.acidity="inert"
    
    Proline=particle()
    Proline.name="P"
    Proline.radius=0.18*ureg.nm,
    Proline.type=29
    Proline.q=0
    Proline.acidity="inert"
    
    Serine=particle()
    Serine.name="S"
    Serine.radius=0.18*ureg.nm,
    Serine.type=30
    Serine.q=0
    Serine.acidity="inert"
    
    Threonine=particle()
    Threonine.name="T"
    Threonine.radius=0.18*ureg.nm,
    Threonine.type=31
    Threonine.q=0
    Threonine.acidity="inert"
    
    Tryptophan=particle()
    Tryptophan.name="W"
    Tryptophan.radius=0.18*ureg.nm,
    Tryptophan.type=32
    Tryptophan.q=0
    Tryptophan.acidity="inert"
    
    Tyrosine=particle()
    Tyrosine.name="Y"
    Tyrosine.radius=0.18*ureg.nm,
    Tyrosine.type=33
    Tyrosine.q=0
    Tyrosine.acidity="inert"
    
    Valine=particle()
    Valine.name="V"
    Valine.radius=0.18*ureg.nm,
    Valine.type=34
    Valine.q=0
    Valine.acidity="inert"
    
    Lysine=particle()
    Lysine.name="K"
    Lysine.radius=0.18*ureg.nm,
    Lysine.type={"unprotonated": 35, "protonated": 36}
    Lysine.q={"unprotonated": 0, "protonated": 1}
    Lysine.acidity="basic"
    
    Histidine=particle()
    Histidine.name="H"
    Histidine.radius=0.18*ureg.nm,
    Histidine.type={"unprotonated": 37, "protonated": 38}
    Histidine.q={"unprotonated": 0, "protonated": 1}
    Histidine.acidity="basic"
    
    Cysteine=particle()
    Cysteine.name="C"
    Cysteine.radius=0.18*ureg.nm,
    Cysteine.type={"unprotonated": 39, "protonated": 40}
    Cysteine.q={"unprotonated": 0, "protonated": 1}
    Cysteine.acidity="basic"
    
    Arginine=particle()
    Arginine.name="R"
    Arginine.radius=0.18*ureg.nm,
    Arginine.type={"unprotonated": 41, "protonated": 42}
    Arginine.q={"unprotonated": 0, "protonated": 1}
    Arginine.acidity="basic"
    
    Amino_end=particle()
    Amino_end.name="n"
    Amino_end.radius=0.18*ureg.nm,
    Amino_end.type={"unprotonated": 43, "protonated": 44}
    Amino_end.q={"unprotonated": 0, "protonated": 1}
    Amino_end.acidity="basic"

    Carboxyl_end=particle()
    Carboxyl_end.name="c"
    Carboxyl_end.radius=0.18*ureg.nm,
    Carboxyl_end.type={"unprotonated": 45, "protonated": 46}
    Carboxyl_end.q={"protonated": 0, "unprotonated": -1}
    Carboxyl_end.acidity="acid"

    Aspartic=particle()
    Aspartic.name="D"
    Aspartic.radius=0.18*ureg.nm,
    Aspartic.type={"unprotonated": 47, "protonated": 48}
    Aspartic.q={"protonated": 0, "unprotonated": -1}
    Aspartic.acidity="acid"

    Glutamic=particle()
    Glutamic.name="E"
    Glutamic.radius=0.18*ureg.nm,
    Glutamic.type={"unprotonated": 49, "protonated": 50}
    Glutamic.q={"protonated": 0, "unprotonated": -1}
    Glutamic.acidity="acid"
    
    Phosphorilated_Serine=particle()
    Phosphorilated_Serine.name="J"
    Phosphorilated_Serine.radius=0.18*ureg.nm,
    Phosphorilated_Serine.type={"unprotonated": 51, "protonated": 52}
    Phosphorilated_Serine.q={"protonated": -1, "unprotonated": -2}
    Phosphorilated_Serine.acidity="acid"
    
    Phosphorilated_Threonine=particle()
    Phosphorilated_Threonine.name="U"
    Phosphorilated_Threonine.radius=0.18*ureg.nm,
    Phosphorilated_Threonine.type={"unprotonated": 53, "protonated": 54}
    Phosphorilated_Threonine.q={"protonated": -1, "unprotonated": -2}
    Phosphorilated_Threonine.acidity="acid"
    
    Phosphorilated_Tryptophan=particle()
    Phosphorilated_Tryptophan.name="Z"
    Phosphorilated_Tryptophan.radius=0.18*ureg.nm,
    Phosphorilated_Tryptophan.type={"unprotonated": 55, "protonated": 56}
    Phosphorilated_Tryptophan.q={"protonated": -1, "unprotonated": -2}
    Phosphorilated_Tryptophan.acidity="acid"

