#
# Copyright (C) 2024-2026 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import tempfile
import espressomd
import pandas as pd
import numpy as np
import unittest as ut
import pyMBE
from pyMBE.lib.lattice import DiamondLattice
import pyMBE.lib.handy_functions as hf


espresso_system=espressomd.System (box_l = [100]*3)

class Test(ut.TestCase):

    def test_io_particles_and_particle_states_templates(self):
        """
        Checks that information in the pyMBE database about 
        particle and particle_state templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        units = pmb.units
        pmb.define_particle(name="Z",
                    sigma=3.5 * units.reduced_length,
                    cutoff=4 * units.reduced_length,
                    offset=0 * units.reduced_length,
                    epsilon=0.2 * units.reduced_energy,
                    acidity="acidic",
                    pka=4.25)
        pmb.define_particle(name="X",
                        sigma=3.5 * units.reduced_length,
                        cutoff=4 * units.reduced_length,
                        offset=0 * units.reduced_length,
                        epsilon=0.2 * units.reduced_energy,
                        z=1)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded particle_state templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="particle_state"),
                                      new_pmb.get_templates_df(pmb_type="particle_state"))
        # Test that the loaded particle templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="particle"),
                                      new_pmb.get_templates_df(pmb_type="particle"))
        
    def test_io_lj_templates(self):
        """
        Checks that information in the pyMBE database about 
        LennardJOnes templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        units = pmb.units
        pmb.define_particle(name="Z",
                            sigma=3.5 * units.reduced_length,
                            cutoff=4 * units.reduced_length,
                            offset=0 * units.reduced_length,
                            epsilon=0.2 * units.reduced_energy,
                            acidity="acidic",
                            pka=4.25)
        pmb.define_particle(name="X",
                            sigma=3.5 * units.reduced_length,
                            cutoff=4 * units.reduced_length,
                            offset=0 * units.reduced_length,
                            epsilon=0.2 * units.reduced_energy,
                            z=1)
        pmb.setup_lj_interactions(espresso_system=espresso_system)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded LJ templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="lj"),
                                      new_pmb.get_templates_df(pmb_type="lj"))
    
    def test_io_bond_templates(self):
        """
        Checks that information in the pyMBE database about 
        residue templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        parameters1 = {"k":  100.0 * pmb.units.reduced_energy / (pmb.units.reduced_length**2),
                       "r_0": 1.0  * pmb.units.reduced_length}
        parameters2 = {'r_0'    : 0.4 * pmb.units.nm,
                        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                        'd_r_max': 0.8 * pmb.units.nm}
        pmb.define_bond(bond_type="harmonic",
                    bond_parameters=parameters1,
                    particle_pairs=[["Z","Z"], 
                                    ["Z","X"],
                                    ["X","X"]])
        pmb.define_default_bond(bond_type="FENE",
                                bond_parameters=parameters2)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded bond templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="bond"),
                                      new_pmb.get_templates_df(pmb_type="bond"))

    def test_io_residues_templates(self):
        """
        Checks that information in the pyMBE database about 
        residue templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_residue(name="R1", 
                           central_bead="Z", 
                           side_chains=["X","Z"])
        pmb.define_residue(name="R2", 
                           central_bead="X", 
                           side_chains=[])
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded residue templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="residue"),
                                      new_pmb.get_templates_df(pmb_type="residue"))
        
    def test_io_molecule_templates(self):
        """
        Checks that information in the pyMBE database about 
        molecule templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_molecule(name="M1", 
                            residue_list=["R1","R2"]*2)
        pmb.define_molecule(name="M2", 
                            residue_list=["R2","R1"]*20)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded molecule templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="molecule"),
                                      new_pmb.get_templates_df(pmb_type="molecule"))
        
    def test_io_peptide_templates(self):
        """
        Checks that information in the pyMBE database about 
        peptide templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_peptide(name="Peptide1",
                           model="1beadAA",
                           sequence="KKKKDDDD")
        pmb.define_peptide(name="Peptide2",
                           model="2beadAA",
                           sequence="RRRREEEE")
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded peptide templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="peptide"),
                                      new_pmb.get_templates_df(pmb_type="peptide"))
        
    def test_io_protein_templates(self):
        """
        Checks that information in the pyMBE database about 
        protein templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_protein(name="1beb",
                            model="2beadAA",
                            sequence="KKKKKK")
        pmb.define_protein(name="Avidin",
                            model="1beadAA",
                            sequence="KKKKKK")
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded protein templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="protein"),
                                      new_pmb.get_templates_df(pmb_type="protein"))
        
    def test_io_hydrogel_templates(self):
        """
        Checks that information in the pyMBE database about 
        hydrogel templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        diamond_lattice = DiamondLattice(4, 3.5 * pmb.units.reduced_length)
        lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)
        # Setting up node topology
        indices = diamond_lattice.indices
        node_topology = []
        for index in range(len(indices)):
            node_topology.append({"particle_name": "X",
                                "lattice_index": indices[index]})
        # Setting up chain topology
        node_labels = lattice_builder.node_labels
        chain_labels = lattice_builder.chain_labels
        reverse_node_labels = {v: k for k, v in node_labels.items()}
        chain_topology = []

        for chain_data in chain_labels.items():
            node_label_pair = chain_data[0]
            node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
            chain_topology.append({'node_start':reverse_node_labels[node_label_s],
                                    'node_end': reverse_node_labels[node_label_e],
                                    'molecule_name':"M1"})

        pmb.define_hydrogel(name="my_hydrogel", 
                            node_map=node_topology, 
                            chain_map=chain_topology)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded hydrogel templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="hydrogel"),
                                      new_pmb.get_templates_df(pmb_type="hydrogel"))

    def test_io_reaction_templates(self):
        """
        Checks that information in the pyMBE database about 
        reaction templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        path_to_pka=pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json"
        pmb.load_pka_set(filename=path_to_pka)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded protein templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_reactions_df(),
                                      new_pmb.get_reactions_df())
        
    def test_io_instances(self):
        """
        Checks that information in the pyMBE database about 
        instances created into espresso is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        # Test instances of a hydrogel (tests hydrogel, molecule, residue, bond and particle instances)
        pmb.define_particle(name="Z",
                    sigma=3.5 * pmb.units.reduced_length,
                    cutoff=4  * pmb.units.reduced_length,
                    offset=0  * pmb.units.reduced_length,
                    epsilon=0.2 * pmb.units.reduced_energy,
                    acidity="acidic",
                    pka=4.25)
        pmb.define_particle(name="X",
                            sigma=3.5 * pmb.units.reduced_length,
                            cutoff=4 * pmb.units.reduced_length,
                            offset=0 * pmb.units.reduced_length,
                            epsilon=0.2 * pmb.units.reduced_energy,
                            z=1)
        pmb.define_residue(name="R1", 
                   central_bead="Z", 
                   side_chains=["X","Z"])
        parameters = {"k":  100.0 * pmb.units.reduced_energy / (pmb.units.reduced_length**2),
                      "r_0": 1.0  * pmb.units.reduced_length}
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=parameters,
                        particle_pairs=[["Z","Z"], 
                                        ["Z","X"],
                                        ["X","X"]])
        pmb.define_molecule(name="M1", 
                            residue_list=["R1"]*1)
        diamond_lattice = DiamondLattice(4, 3.5 * pmb.units.reduced_length)
        lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)
        # Setting up node topology --> Nodes are particles of type "X"
        indices = diamond_lattice.indices
        node_topology = []
        for index in range(len(indices)):
            node_topology.append({"particle_name": "X",
                                "lattice_index": indices[index]})
        # Setting up chain topology --> Chains are molecules of type "M1"
        node_labels = lattice_builder.node_labels
        chain_labels = lattice_builder.chain_labels
        reverse_node_labels = {v: k for k, v in node_labels.items()}
        chain_topology = []
        for chain_data in chain_labels.items():
            node_label_pair = chain_data[0]
            node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
            chain_topology.append({'node_start': reverse_node_labels[node_label_s],
                                    'node_end': reverse_node_labels[node_label_e],
                                    'molecule_name': "M1"})
        pmb.define_hydrogel("my_hydrogel", 
                            node_topology, 
                            chain_topology)
        assembly_id = pmb.create_hydrogel(name="my_hydrogel",
                                          espresso_system=espresso_system)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded instances are equal to the original
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="hydrogel"),
                                      new_pmb.get_instances_df(pmb_type="hydrogel"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="molecule"),
                                      new_pmb.get_instances_df(pmb_type="molecule"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="residues"),
                                      new_pmb.get_instances_df(pmb_type="residues"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="bond"),
                                      new_pmb.get_instances_df(pmb_type="bond"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="particle"),
                                      new_pmb.get_instances_df(pmb_type="particle"))
        # Clean up before the next test
        pmb.delete_instances_in_system(espresso_system=espresso_system,
                                       instance_id=assembly_id,
                                       pmb_type="hydrogel")
        # Test instances of a peptide (tests peptide, residue, bond and particle instances)
        path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021"
        path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"
        pmb.load_database (folder=path_to_interactions) # Defines particles
        pmb.load_pka_set(filename=path_to_pka)
        pka_set = pmb.get_pka_set()
        for particle_name in pka_set.keys():
            pmb.define_monoprototic_particle_states(particle_name=particle_name,
                                                    acidity=pka_set[particle_name]["acidity"])
        hf.define_peptide_AA_residues(sequence="KKKDDD",
                                      model="1beadAA",
                                      pmb=pmb)
        parameters2 = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}
        pmb.define_default_bond(bond_type="FENE",
                                bond_parameters=parameters2)
        pmb.define_peptide(name="Peptide1",
                            model="1beadAA",
                            sequence="KKKKDDDD")
        pep_ids = pmb.create_molecule(name="Peptide1",
                            number_of_molecules=2,
                            espresso_system=espresso_system,
                            use_default_bond=True)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded instances are equal to the original
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="peptide"),
                                      new_pmb.get_instances_df(pmb_type="peptide"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="residues"),
                                      new_pmb.get_instances_df(pmb_type="residues"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="bond"),
                                      new_pmb.get_instances_df(pmb_type="bond"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="particle"),
                                      new_pmb.get_instances_df(pmb_type="particle"))
        # Clean up before the next test
        for pepid in pep_ids:
            pmb.delete_instances_in_system(espresso_system=espresso_system,
                                        instance_id=pepid,
                                        pmb_type="peptide")
        pmb.db.delete_templates(pmb_type="particle")
        pmb.db.delete_templates(pmb_type="particle_state")
        pmb.db.delete_templates(pmb_type="residue")
        pmb.db.delete_reactions()
        # Test instances of a protein (tests protein, residue and particle instances)
        path_to_protein_structure = pmb.root / "parameters" / "globular_proteins" / f"1beb.vtf",    
        topology_dict, sequence = pmb.read_protein_vtf (filename=path_to_protein_structure[0])
        pmb.load_pka_set(filename=path_to_pka)
        # Define AA particles and residues
        hf.define_protein_AA_particles(topology_dict=topology_dict,
                                       pmb=pmb,
                                       pka_set=pka_set)
        hf.define_protein_AA_residues(sequence=sequence,
                                      model="2beadAA",
                                      pmb=pmb)
        pmb.define_protein(name="1beb",
                           model="2beadAA",
                           sequence="KKKKKK")
        prot_ids = pmb.create_protein(name="1beb",
                                    number_of_proteins=1,
                                    espresso_system=espresso_system,
                                    topology_dict=topology_dict)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded instances are equal to the original
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="protein"),
                                      new_pmb.get_instances_df(pmb_type="protein"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="residues"),
                                      new_pmb.get_instances_df(pmb_type="residues"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="bond"),
                                      new_pmb.get_instances_df(pmb_type="bond"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="particle"),
                                      new_pmb.get_instances_df(pmb_type="particle"))
        # Clean up 
        for protid in prot_ids:
            pmb.delete_instances_in_system(espresso_system=espresso_system,
                                        instance_id=protid,
                                        pmb_type="protein")
        

if __name__ == '__main__':
    ut.main()

