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

class general:
    
    name=None
    lj_WCA = custom_lj()
    lj_WCA.name='WCA'
    lj_WCA.actors=['WCA']
    lj_WCA.epsilon= 2.5 * ureg.kJ / ureg.mol
    lj_WCA.sigma=0.35 * ureg.nm # nm
    lj_WCA.cutoff=1.12246*0.35 * ureg.nm # nm
    lj_WCA.shift='auto'
    
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
    principal_chain={"default": 'sequence'}
    side_chain=None
    
    default_bond=bond()
    default_bond.actors=['default']
    default_bond.type='harmonic'
    default_bond.bondl=0.33*ureg.nm
    default_bond.k=2000*ureg.kJ / ureg.nm**2 / ureg.mol

    A=particle()
    A.name="A"
    A.acidity="inert"
    A.radius=0.18*ureg.nm,
    A.type=21
    A.q=0

    N=particle()
    N.name="N"
    N.radius=0.18*ureg.nm,
    N.type=22
    N.q=0
    N.acidity="inert"
    
    Q=particle()
    Q.name="Q"
    Q.radius=0.18*ureg.nm,
    Q.type=23
    Q.q=0
    Q.acidity="inert"
    
    G=particle()
    G.name="G"
    G.radius=0.18*ureg.nm,
    G.type=24
    G.q=0
    G.acidity="inert"
    
    I=particle()
    I.name="I"
    I.radius=0.18*ureg.nm,
    I.type=25
    I.q=0
    I.acidity="inert"
    
    L=particle()
    L.name="L"
    L.radius=0.18*ureg.nm,
    L.type=26
    L.q=0
    L.acidity="inert"
    
    M=particle()
    M.name="M"
    M.radius=0.18*ureg.nm,
    M.type=27
    M.q=0
    M.acidity="inert"
    
    F=particle()
    F.name="F"
    F.radius=0.18*ureg.nm,
    F.type=28
    F.q=0
    F.acidity="inert"
    
    P=particle()
    P.name="P"
    P.radius=0.18*ureg.nm,
    P.type=29
    P.q=0
    P.acidity="inert"
    
    S=particle()
    S.name="S"
    S.radius=0.18*ureg.nm,
    S.type=30
    S.q=0
    S.acidity="inert"
    
    T=particle()
    T.name="T"
    T.radius=0.18*ureg.nm,
    T.type=31
    T.q=0
    T.acidity="inert"
    
    W=particle()
    W.name="W"
    W.radius=0.18*ureg.nm,
    W.type=32
    W.q=0
    W.acidity="inert"
    
    Y=particle()
    Y.name="Y"
    Y.radius=0.18*ureg.nm,
    Y.type=33
    Y.q=0
    Y.acidity="inert"
    
    V=particle()
    V.name="V"
    V.radius=0.18*ureg.nm,
    V.type=34
    V.q=0
    V.acidity="inert"
    
    K=particle()
    K.name="K"
    K.radius=0.18*ureg.nm,
    K.type={"unprotonated": 35, "protonated": 36}
    K.q={"unprotonated": 0, "protonated": 1}
    K.acidity="basic"
    
    H=particle()
    H.name="H"
    H.radius=0.18*ureg.nm,
    H.type={"unprotonated": 37, "protonated": 38}
    H.q={"unprotonated": 0, "protonated": 1}
    H.acidity="basic"
    
    C=particle()
    C.name="C"
    C.radius=0.18*ureg.nm,
    C.type={"unprotonated": 39, "protonated": 40}
    C.q={"unprotonated": 0, "protonated": 1}
    C.acidity="basic"
    
    R=particle()
    R.name="R"
    R.radius=0.18*ureg.nm,
    R.type={"unprotonated": 41, "protonated": 42}
    R.q={"unprotonated": 0, "protonated": 1}
    R.acidity="basic"
    
    n=particle()
    n.name="n"
    n.radius=0.18*ureg.nm,
    n.type={"unprotonated": 43, "protonated": 44}
    n.q={"unprotonated": 0, "protonated": 1}
    n.acidity="basic"

    c=particle()
    c.name="c"
    c.radius=0.18*ureg.nm,
    c.type={"unprotonated": 45, "protonated": 46}
    c.q={"protonated": 0, "unprotonated": -1}
    c.acidity="acid"

    D=particle()
    D.name="D"
    D.radius=0.18*ureg.nm,
    D.type={"unprotonated": 47, "protonated": 48}
    D.q={"protonated": 0, "unprotonated": -1}
    D.acidity="acid"

    E=particle()
    E.name="E"
    E.radius=0.18*ureg.nm,
    E.type={"unprotonated": 49, "protonated": 50}
    E.q={"protonated": 0, "unprotonated": -1}
    E.acidity="acid"
    
    J=particle()
    J.name="J"
    J.radius=0.18*ureg.nm,
    J.type={"unprotonated": 51, "protonated": 52}
    J.q={"protonated": -1, "unprotonated": -2}
    J.acidity="acid"
    
    U=particle()
    U.name="U"
    U.radius=0.18*ureg.nm,
    U.type={"unprotonated": 53, "protonated": 54}
    U.q={"protonated": -1, "unprotonated": -2}
    U.acidity="acid"
    
    Z=particle()
    Z.name="Z"
    Z.radius=0.18*ureg.nm,
    Z.type={"unprotonated": 55, "protonated": 56}
    Z.q={"protonated": -1, "unprotonated": -2}
    Z.acidity="acid"

class two_bead_peptide:

    name='2beadpeptide'
    principal_chain={"c": "sequence", "n": "sequence", "default": "C_alpha"}
    side_chain={"c": [], "n": [], "default": ["sequence"]}
    
    default_bond=bond()
    default_bond.actors=['default']
    default_bond.type='harmonic'
    default_bond.bondl=0.33*ureg.nm
    default_bond.k=2000*ureg.kJ / ureg.nm**2 / ureg.mol
    
    principal_bond=bond()
    principal_bond.actors=['C_alpha','C_alpha']
    principal_bond.type='harmonic'
    principal_bond.bondl=0.388*ureg.nm
    principal_bond.k=2000*ureg.kJ / ureg.nm**2 / ureg.mol

    C_alpha=particle()
    C_alpha.name="C_alpha"
    C_alpha.acidity="inert"
    C_alpha.radius=0.18*ureg.nm,
    C_alpha.type=20
    C_alpha.q=0
    
    A=particle()
    A.name="A"
    A.acidity="inert"
    A.radius=0.18*ureg.nm,
    A.type=21
    A.q=0

    N=particle()
    N.name="N"
    N.radius=0.18*ureg.nm,
    N.type=22
    N.q=0
    N.acidity="inert"
    
    Q=particle()
    Q.name="Q"
    Q.radius=0.18*ureg.nm,
    Q.type=23
    Q.q=0
    Q.acidity="inert"
    
    G=particle()
    G.name="G"
    G.radius=0.18*ureg.nm,
    G.type=24
    G.q=0
    G.acidity="inert"
    
    I=particle()
    I.name="I"
    I.radius=0.18*ureg.nm,
    I.type=25
    I.q=0
    I.acidity="inert"
    
    L=particle()
    L.name="L"
    L.radius=0.18*ureg.nm,
    L.type=26
    L.q=0
    L.acidity="inert"
    
    M=particle()
    M.name="M"
    M.radius=0.18*ureg.nm,
    M.type=27
    M.q=0
    M.acidity="inert"
    
    F=particle()
    F.name="F"
    F.radius=0.18*ureg.nm,
    F.type=28
    F.q=0
    F.acidity="inert"
    
    P=particle()
    P.name="P"
    P.radius=0.18*ureg.nm,
    P.type=29
    P.q=0
    P.acidity="inert"
    
    S=particle()
    S.name="S"
    S.radius=0.18*ureg.nm,
    S.type=30
    S.q=0
    S.acidity="inert"
    
    T=particle()
    T.name="T"
    T.radius=0.18*ureg.nm,
    T.type=31
    T.q=0
    T.acidity="inert"
    
    W=particle()
    W.name="W"
    W.radius=0.18*ureg.nm,
    W.type=32
    W.q=0
    W.acidity="inert"
    
    Y=particle()
    Y.name="Y"
    Y.radius=0.18*ureg.nm,
    Y.type=33
    Y.q=0
    Y.acidity="inert"
    
    V=particle()
    V.name="V"
    V.radius=0.18*ureg.nm,
    V.type=34
    V.q=0
    V.acidity="inert"
    
    K=particle()
    K.name="K"
    K.radius=0.18*ureg.nm,
    K.type={"unprotonated": 35, "protonated": 36}
    K.q={"unprotonated": 0, "protonated": 1}
    K.acidity="basic"
    
    H=particle()
    H.name="H"
    H.radius=0.18*ureg.nm,
    H.type={"unprotonated": 37, "protonated": 38}
    H.q={"unprotonated": 0, "protonated": 1}
    H.acidity="basic"
    
    C=particle()
    C.name="C"
    C.radius=0.18*ureg.nm,
    C.type={"unprotonated": 39, "protonated": 40}
    C.q={"unprotonated": 0, "protonated": 1}
    C.acidity="basic"
    
    R=particle()
    R.name="R"
    R.radius=0.18*ureg.nm,
    R.type={"unprotonated": 41, "protonated": 42}
    R.q={"unprotonated": 0, "protonated": 1}
    R.acidity="basic"
    
    n=particle()
    n.name="n"
    n.radius=0.18*ureg.nm,
    n.type={"unprotonated": 43, "protonated": 44}
    n.q={"unprotonated": 0, "protonated": 1}
    n.acidity="basic"

    c=particle()
    c.name="c"
    c.radius=0.18*ureg.nm,
    c.type={"unprotonated": 45, "protonated": 46}
    c.q={"protonated": 0, "unprotonated": -1}
    c.acidity="acid"

    D=particle()
    D.name="D"
    D.radius=0.18*ureg.nm,
    D.type={"unprotonated": 47, "protonated": 48}
    D.q={"protonated": 0, "unprotonated": -1}
    D.acidity="acid"

    E=particle()
    E.name="E"
    E.radius=0.18*ureg.nm,
    E.type={"unprotonated": 49, "protonated": 50}
    E.q={"protonated": 0, "unprotonated": -1}
    E.acidity="acid"
    
    J=particle()
    J.name="J"
    J.radius=0.18*ureg.nm,
    J.type={"unprotonated": 51, "protonated": 52}
    J.q={"protonated": -1, "unprotonated": -2}
    J.acidity="acid"
    
    U=particle()
    U.name="U"
    U.radius=0.18*ureg.nm,
    U.type={"unprotonated": 53, "protonated": 54}
    U.q={"protonated": -1, "unprotonated": -2}
    U.acidity="acid"
    
    Z=particle()
    Z.name="Z"
    Z.radius=0.18*ureg.nm,
    Z.type={"unprotonated": 55, "protonated": 56}
    Z.q={"protonated": -1, "unprotonated": -2}
    Z.acidity="acid"

