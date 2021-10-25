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

        class general:
            
            name=None
            lj_WCA = parameters_self.custom_lj()
            lj_WCA.name='WCA'
            lj_WCA.actors=['WCA']
            lj_WCA.epsilon= 1 * parameters_self.units('reduced_energy')
            lj_WCA.sigma=1 * parameters_self.units('reduced_length')
            lj_WCA.cutoff=2**(1.0 / 6)*parameters_self.units('reduced_length')
            lj_WCA.shift='auto'

            default_bond=parameters_self.bond()
            default_bond.actors=['default']
            default_bond.type='harmonic'
            default_bond.bondl=None
            default_bond.k=100 * parameters_self.units('reduced_energy / reduced_length**2')
        
        parameters_self.general=general 
        
        class small_ions:
            name='small_ions'
            beads_per_residue=1
            principal_chain='sequence'
            side_chain=None
            

            cation=parameters_self.particle()
            cation.name="cation"
            cation.type=18
            cation.radius=0.15 * parameters_self.units.nm
            cation.q=1

            anion=parameters_self.particle()
            anion.name="anion"
            anion.type=19
            anion.radius=0.15 * parameters_self.units.nm
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
            A.radius=0.15*parameters_self.units.nm
            A.type=21
            A.q=0

            N=parameters_self.particle()
            N.name="N"
            N.radius=0.15*parameters_self.units.nm
            N.type=22
            N.q=0
            N.acidity="inert"
            
            Q=parameters_self.particle()
            Q.name="Q"
            Q.radius=0.15*parameters_self.units.nm
            Q.type=23
            Q.q=0
            Q.acidity="inert"
            
            G=parameters_self.particle()
            G.name="G"
            G.radius=0.15*parameters_self.units.nm
            G.type=24
            G.q=0
            G.acidity="inert"
            
            I=parameters_self.particle()
            I.name="I"
            I.radius=0.15*parameters_self.units.nm
            I.type=25
            I.q=0
            I.acidity="inert"
            
            L=parameters_self.particle()
            L.name="L"
            L.radius=0.15*parameters_self.units.nm
            L.type=26
            L.q=0
            L.acidity="inert"
            
            M=parameters_self.particle()
            M.name="M"
            M.radius=0.15*parameters_self.units.nm
            M.type=27
            M.q=0
            M.acidity="inert"
            
            F=parameters_self.particle()
            F.name="F"
            F.radius=0.15*parameters_self.units.nm
            F.type=28
            F.q=0
            F.acidity="inert"
            
            P=parameters_self.particle()
            P.name="P"
            P.radius=0.15*parameters_self.units.nm
            P.type=29
            P.q=0
            P.acidity="inert"
            
            S=parameters_self.particle()
            S.name="S"
            S.radius=0.15*parameters_self.units.nm
            S.type=30
            S.q=0
            S.acidity="inert"
            
            T=parameters_self.particle()
            T.name="T"
            T.radius=0.15*parameters_self.units.nm
            T.type=31
            T.q=0
            T.acidity="inert"
            
            W=parameters_self.particle()
            W.name="W"
            W.radius=0.15*parameters_self.units.nm
            W.type=32
            W.q=0
            W.acidity="inert"
            
            Y=parameters_self.particle()
            Y.name="Y"
            Y.radius=0.15*parameters_self.units.nm
            Y.type=33
            Y.q=0
            Y.acidity="inert"
            
            V=parameters_self.particle()
            V.name="V"
            V.radius=0.15*parameters_self.units.nm
            V.type=34
            V.q=0
            V.acidity="inert"
            
            K=parameters_self.particle()
            K.name="K"
            K.radius=0.15*parameters_self.units.nm
            K.type={"unprotonated": 35, "protonated": 36}
            K.q={"unprotonated": 0, "protonated": 1}
            K.acidity="basic"
            
            H=parameters_self.particle()
            H.name="H"
            H.radius=0.15*parameters_self.units.nm
            H.type={"unprotonated": 37, "protonated": 38}
            H.q={"unprotonated": 0, "protonated": 1}
            H.acidity="basic"
            
            C=parameters_self.particle()
            C.name="C"
            C.radius=0.15*parameters_self.units.nm
            C.type={"unprotonated": 39, "protonated": 40}
            C.q={"unprotonated": 0, "protonated": 1}
            C.acidity="basic"
            
            R=parameters_self.particle()
            R.name="R"
            R.radius=0.15*parameters_self.units.nm
            R.type={"unprotonated": 41, "protonated": 42}
            R.q={"unprotonated": 0, "protonated": 1}
            R.acidity="basic"
            
            n=parameters_self.particle()
            n.name="n"
            n.radius=0.15*parameters_self.units.nm
            n.type={"unprotonated": 43, "protonated": 44}
            n.q={"unprotonated": 0, "protonated": 1}
            n.acidity="basic"

            c=parameters_self.particle()
            c.name="c"
            c.radius=0.15*parameters_self.units.nm
            c.type={"unprotonated": 45, "protonated": 46}
            c.q={"protonated": 0, "unprotonated": -1}
            c.acidity="acid"

            D=parameters_self.particle()
            D.name="D"
            D.radius=0.15*parameters_self.units.nm
            D.type={"unprotonated": 47, "protonated": 48}
            D.q={"protonated": 0, "unprotonated": -1}
            D.acidity="acid"

            E=parameters_self.particle()
            E.name="E"
            E.radius=0.15*parameters_self.units.nm
            E.type={"unprotonated": 49, "protonated": 50}
            E.q={"protonated": 0, "unprotonated": -1}
            E.acidity="acid"
            
            J=parameters_self.particle()
            J.name="J"
            J.radius=0.15*parameters_self.units.nm
            J.type={"unprotonated": 51, "protonated": 52}
            J.q={"protonated": -1, "unprotonated": -2}
            J.acidity="acid"
            
            U=parameters_self.particle()
            U.name="U"
            U.radius=0.15*parameters_self.units.nm
            U.type={"unprotonated": 53, "protonated": 54}
            U.q={"protonated": -1, "unprotonated": -2}
            U.acidity="acid"
            
            Z=parameters_self.particle()
            Z.name="Z"
            Z.radius=0.15*parameters_self.units.nm
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
            C_alpha.radius=0.15*parameters_self.units.nm
            C_alpha.type=20
            C_alpha.q=0
            
            A=parameters_self.particle()
            A.name="A"
            A.acidity="inert"
            A.radius=0.15*parameters_self.units.nm
            A.type=21
            A.q=0

            N=parameters_self.particle()
            N.name="N"
            N.radius=0.15*parameters_self.units.nm
            N.type=22
            N.q=0
            N.acidity="inert"
            
            Q=parameters_self.particle()
            Q.name="Q"
            Q.radius=0.15*parameters_self.units.nm
            Q.type=23
            Q.q=0
            Q.acidity="inert"
            
            G=parameters_self.particle()
            G.name="G"
            G.radius=0.15*parameters_self.units.nm
            G.type=24
            G.q=0
            G.acidity="inert"
            
            I=parameters_self.particle()
            I.name="I"
            I.radius=0.15*parameters_self.units.nm
            I.type=25
            I.q=0
            I.acidity="inert"
            
            L=parameters_self.particle()
            L.name="L"
            L.radius=0.15*parameters_self.units.nm
            L.type=26
            L.q=0
            L.acidity="inert"
            
            M=parameters_self.particle()
            M.name="M"
            M.radius=0.15*parameters_self.units.nm
            M.type=27
            M.q=0
            M.acidity="inert"
            
            F=parameters_self.particle()
            F.name="F"
            F.radius=0.15*parameters_self.units.nm
            F.type=28
            F.q=0
            F.acidity="inert"
            
            P=parameters_self.particle()
            P.name="P"
            P.radius=0.15*parameters_self.units.nm
            P.type=29
            P.q=0
            P.acidity="inert"
            
            S=parameters_self.particle()
            S.name="S"
            S.radius=0.15*parameters_self.units.nm
            S.type=30
            S.q=0
            S.acidity="inert"
            
            T=parameters_self.particle()
            T.name="T"
            T.radius=0.15*parameters_self.units.nm
            T.type=31
            T.q=0
            T.acidity="inert"
            
            W=parameters_self.particle()
            W.name="W"
            W.radius=0.15*parameters_self.units.nm
            W.type=32
            W.q=0
            W.acidity="inert"
            
            Y=parameters_self.particle()
            Y.name="Y"
            Y.radius=0.15*parameters_self.units.nm
            Y.type=33
            Y.q=0
            Y.acidity="inert"
            
            V=parameters_self.particle()
            V.name="V"
            V.radius=0.15*parameters_self.units.nm
            V.type=34
            V.q=0
            V.acidity="inert"
            
            K=parameters_self.particle()
            K.name="K"
            K.radius=0.15*parameters_self.units.nm
            K.type={"unprotonated": 35, "protonated": 36}
            K.q={"unprotonated": 0, "protonated": 1}
            K.acidity="basic"
            
            H=parameters_self.particle()
            H.name="H"
            H.radius=0.15*parameters_self.units.nm
            H.type={"unprotonated": 37, "protonated": 38}
            H.q={"unprotonated": 0, "protonated": 1}
            H.acidity="basic"
            
            C=parameters_self.particle()
            C.name="C"
            C.radius=0.15*parameters_self.units.nm
            C.type={"unprotonated": 39, "protonated": 40}
            C.q={"unprotonated": 0, "protonated": 1}
            C.acidity="basic"
            
            R=parameters_self.particle()
            R.name="R"
            R.radius=0.15*parameters_self.units.nm
            R.type={"unprotonated": 41, "protonated": 42}
            R.q={"unprotonated": 0, "protonated": 1}
            R.acidity="basic"
            
            n=parameters_self.particle()
            n.name="n"
            n.radius=0.15*parameters_self.units.nm
            n.type={"unprotonated": 43, "protonated": 44}
            n.q={"unprotonated": 0, "protonated": 1}
            n.acidity="basic"

            c=parameters_self.particle()
            c.name="c"
            c.radius=0.15*parameters_self.units.nm
            c.type={"unprotonated": 45, "protonated": 46}
            c.q={"protonated": 0, "unprotonated": -1}
            c.acidity="acid"

            D=parameters_self.particle()
            D.name="D"
            D.radius=0.15*parameters_self.units.nm
            D.type={"unprotonated": 47, "protonated": 48}
            D.q={"protonated": 0, "unprotonated": -1}
            D.acidity="acid"

            E=parameters_self.particle()
            E.name="E"
            E.radius=0.15*parameters_self.units.nm
            E.type={"unprotonated": 49, "protonated": 50}
            E.q={"protonated": 0, "unprotonated": -1}
            E.acidity="acid"
            
            J=parameters_self.particle()
            J.name="J"
            J.radius=0.15*parameters_self.units.nm
            J.type={"unprotonated": 51, "protonated": 52}
            J.q={"protonated": -1, "unprotonated": -2}
            J.acidity="acid"
            
            U=parameters_self.particle()
            U.name="U"
            U.radius=0.15*parameters_self.units.nm
            U.type={"unprotonated": 53, "protonated": 54}
            U.q={"protonated": -1, "unprotonated": -2}
            U.acidity="acid"
            
            Z=parameters_self.particle()
            Z.name="Z"
            Z.radius=0.15*parameters_self.units.nm
            Z.type={"unprotonated": 55, "protonated": 56}
            Z.q={"protonated": -1, "unprotonated": -2}
            Z.acidity="acid"

        parameters_self.two_bead_peptide=two_bead_peptide
