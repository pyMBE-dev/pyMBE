#
# Copyright (C) 2026 pyMBE-dev team
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

# Import pyMBE and other libraries
import pyMBE
from pyMBE.lib.lattice import DiamondLattice
import numpy as np
import unittest as ut
import espressomd
import itertools

# Create an instance of the ESPResSo system
box_l=[10]*3
espresso_system = espressomd.System(box_l=box_l)
espresso_system.time_step = 0.002
espresso_system.cell_system.skin = 0.4


class Test(ut.TestCase):
    def setUp(self):
        espresso_system.part.clear()
        espresso_system.bonded_inter.clear()

    def define_templates(self, pmb):
        for name in ('A', 'B', 'C', 'D'):
            pmb.define_particle(name=name,
                z=0,
                sigma=0.4*pmb.units.nm,
                epsilon=1*pmb.units('reduced_energy'))

        self.periodic_dihedral_params = {
            'k': 5 * pmb.units('reduced_energy'),
            'phi_s': np.pi * pmb.units(''),
            'mult': 1,
        }

        self.periodic_dihedral_params_alt = {
            'k': 3 * pmb.units('reduced_energy'),
            'phi_s': np.pi / 2 * pmb.units(''),
            'mult': 2,
        }

        self.multi_term_dihedral_params = [
            {'k': 4 * pmb.units('reduced_energy'), 'phi_s': 0.0 * pmb.units(''), 'mult': 1},
            {'k': 2 * pmb.units('reduced_energy'), 'phi_s': np.pi * pmb.units(''), 'mult': 3},
        ]

        self.harmonic_bond_params = {
            'r_0': 0.4 * pmb.units.nm,
            'k': 400 * pmb.units('reduced_energy / reduced_length**2'),
        }

    def get_dihedral_objects(self, particle_id2):
        """
        Returns the list of dihedral interaction objects stored in ESPResSo for a
        given second particle. Dihedral bonds are added to the second particle of
        the quadruplet in ESPResSo and carry three partner particles (tuple arity 4).
        """
        bonds = espresso_system.part.by_id(particle_id2).bonds
        dihedral_objects = []
        for bond_tuple in bonds:
            if len(bond_tuple) == 4:
                dihedral_objects.append(bond_tuple[0])
        return dihedral_objects

    def check_dihedral_setup(self, dihedral_objects, input_terms):
        """
        Checks that pyMBE sets up the dihedral interaction objects correctly.
        pyMBE takes GROMACS type 1 parameters; the ESPResSo phase is -(phi_s + pi).

        Args:
            dihedral_objects: list of ESPResSo dihedral interaction objects.
            input_terms('list[dict]'): list of the term dictionaries passed to pyMBE.
        """
        self.assertEqual(len(dihedral_objects), len(input_terms),
                         msg="pyMBE did not register one ESPResSo Dihedral object per term")
        for dihedral_object in dihedral_objects:
            self.assertIn("dihedral", str(type(dihedral_object)).lower(),
                          msg="pyMBE does not store the correct type of dihedral interaction")
        actual = sorted((d.get_params()['bend'],
                         d.get_params()['mult'],
                         round(d.get_params()['phase'], 10)) for d in dihedral_objects)
        expected = sorted((term['k'].m_as('reduced_energy'),
                           int(term['mult']),
                           round(-(float(term["phi_s"].magnitude) + np.pi), 10)) for term in input_terms)
        for (b_a, m_a, p_a), (b_e, m_e, p_e) in zip(actual, expected):
            np.testing.assert_almost_equal(actual=b_a, desired=b_e)
            self.assertEqual(m_a, m_e)
            np.testing.assert_almost_equal(actual=p_a, desired=p_e)

    def build_linear_chain(self, pmb):
        """Create particles A-B-C-D and the three backbone bonds."""
        pids = []
        for name in ('A', 'B', 'C', 'D'):
            pids.append(pmb.create_particle(name=name,
                                            box_l=box_l,
                                            number_of_particles=1)[0])
        pmb.create_bond(particle_id1=pids[0], particle_id2=pids[1])
        pmb.create_bond(particle_id1=pids[1], particle_id2=pids[2])
        pmb.create_bond(particle_id1=pids[2], particle_id2=pids[3])
        return pids

    def test_dihedral_setup_periodic(self):
        """
        Unit test to check the setup of a periodic dihedral potential in pyMBE.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])

        pids = self.build_linear_chain(pmb)
        pmb.create_dihedral_potential(particle_id1=pids[0],
                                      particle_id2=pids[1],
                                      particle_id3=pids[2],
                                      particle_id4=pids[3])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        dihedral_objects = self.get_dihedral_objects(particle_id2=pids[1])
        self.assertEqual(len(dihedral_objects), 1, "No dihedral object found on the first particle")
        self.check_dihedral_setup(dihedral_objects=dihedral_objects,
                                  input_terms=[self.periodic_dihedral_params])

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_multi_term_dihedral(self):
        """
        Unit test to check that a multi-term dihedral registers one ESPResSo
        Dihedral object per term, all attached to the same quadruplet.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.multi_term_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])

        pids = self.build_linear_chain(pmb)
        pmb.create_dihedral_potential(particle_id1=pids[0],
                                      particle_id2=pids[1],
                                      particle_id3=pids[2],
                                      particle_id4=pids[3])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        dihedral_objects = self.get_dihedral_objects(particle_id2=pids[1])
        self.assertEqual(len(dihedral_objects), 2, "Expected two Dihedral objects for a two-term dihedral")
        self.check_dihedral_setup(dihedral_objects=dihedral_objects,
                                  input_terms=self.multi_term_dihedral_params)
        self.assertEqual(len(pmb.get_instances_df(pmb_type="dihedral")), 1)

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_dihedral_instance_reuses_cached_object(self):
        """
        Test that cached ESPResSo dihedral interactions are reused for the same template.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])
        dihedral_template = pmb.get_dihedral_template(name1="A", name2="B", name3="C", name4="D")

        pmb.set_simulation_engine(espresso_system)
        first = pmb.simulation_engine._get_dihedral_instance(dihedral_template=dihedral_template)
        second = pmb.simulation_engine._get_dihedral_instance(dihedral_template=dihedral_template)

        self.assertIs(first, second)
        self.assertEqual(len(pmb.db.espresso_dihedral_instances), 1)

        pmb.db.delete_templates(pmb_type="dihedral")

    def test_default_dihedral(self):
        """
        Unit test to check the setup of a default dihedral potential.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_default_dihedral_potential(dihedral_type="periodic",
                                             dihedral_parameters=self.periodic_dihedral_params)

        pids = self.build_linear_chain(pmb)
        pmb.create_dihedral_potential(particle_id1=pids[0],
                                      particle_id2=pids[1],
                                      particle_id3=pids[2],
                                      particle_id4=pids[3],
                                      use_default_dihedral=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        dihedral_objects = self.get_dihedral_objects(particle_id2=pids[1])
        self.assertEqual(len(dihedral_objects), 1)
        self.check_dihedral_setup(dihedral_objects=dihedral_objects,
                                  input_terms=[self.periodic_dihedral_params])

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_specific_dihedral_over_default(self):
        """
        Test that a specific dihedral template is preferred over the default.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_default_dihedral_potential(dihedral_type="periodic",
                                             dihedral_parameters=self.periodic_dihedral_params_alt)
        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])

        pids = self.build_linear_chain(pmb)
        pmb.create_dihedral_potential(particle_id1=pids[0],
                                      particle_id2=pids[1],
                                      particle_id3=pids[2],
                                      particle_id4=pids[3],
                                      use_default_dihedral=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        dihedral_objects = self.get_dihedral_objects(particle_id2=pids[1])
        self.check_dihedral_setup(dihedral_objects=dihedral_objects,
                                  input_terms=[self.periodic_dihedral_params])

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_dihedral_raised_exceptions(self):
        """
        Unit test to check that dihedral-related methods raise the correct exceptions.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        for callback in [pmb.define_dihedral_potential, pmb.define_default_dihedral_potential]:
            with self.subTest(msg=f'using method {callback.__qualname__}()'):
                self.check_dihedral_exceptions(callback, pmb)

    def check_dihedral_exceptions(self, callback, pmb):
        good_term = {'k': 5 * pmb.units('reduced_energy'),
                     'phi_s': np.pi * pmb.units(''),
                     'mult': 1}

        def call(dihedral_type, dihedral_parameters):
            kwargs = {"dihedral_type": dihedral_type,
                      "dihedral_parameters": dihedral_parameters}
            if callback == pmb.define_dihedral_potential:
                kwargs["particle_quadruplets"] = [('A', 'B', 'C', 'D')]
            return kwargs

        # unknown dihedral type
        np.testing.assert_raises(NotImplementedError, callback,
                                 **call('quartic', good_term))
        # tabulated is reserved but not implemented yet
        np.testing.assert_raises(NotImplementedError, callback,
                                 **call('tabulated', good_term))
        # missing k
        np.testing.assert_raises(ValueError, callback,
                                 **call('periodic', {'phi_s': np.pi * pmb.units(''), 'mult': 1}))
        # missing phi_s
        np.testing.assert_raises(ValueError, callback,
                                 **call('periodic', {'k': 5 * pmb.units('reduced_energy'), 'mult': 1}))
        # missing mult
        np.testing.assert_raises(ValueError, callback,
                                 **call('periodic', {'k': 5 * pmb.units('reduced_energy'),
                                                     'phi_s': np.pi * pmb.units('')}))
        # duplicate quadruplet in define_dihedral_potential
        if callback == pmb.define_dihedral_potential:
            np.testing.assert_raises(RuntimeError,
                                     pmb.define_dihedral_potential,
                                     dihedral_type="periodic",
                                     dihedral_parameters=good_term,
                                     particle_quadruplets=[('A', 'B', 'C', 'D'),
                                                           ('D', 'C', 'B', 'A')])

    def test_sanity_get_dihedral_template(self):
        """
        Tests that get_dihedral_template raises ValueError when no template is found.
        """
        pmb = pyMBE.pymbe_library(seed=51)
        np.testing.assert_raises(ValueError,
                                 pmb.get_dihedral_template,
                                 name1="W", name2="X", name3="Y", name4="Z")

    def test_canonical_dihedral_name(self):
        """
        Test that dihedral names are canonical (invariant under full reversal).
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('D', 'C', 'B', 'A')])

        tpl = pmb.get_dihedral_template(name1="A", name2="B", name3="C", name4="D")
        self.assertEqual(tpl.name, "A-B-C-D")
        tpl2 = pmb.get_dihedral_template(name1="D", name2="C", name3="B", name4="A")
        self.assertEqual(tpl2.name, "A-B-C-D")

        pmb.db.delete_templates(pmb_type="dihedral")

    def test_dihedral_without_bonds_raises_error(self):
        """
        Test that creating a dihedral without bonds between the particles raises a ValueError.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])

        pids = []
        for name in ('A', 'B', 'C', 'D'):
            pids.append(pmb.create_particle(name=name,
                                            box_l=box_l,
                                            number_of_particles=1)[0])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        with self.assertRaises(ValueError):
            pmb.create_dihedral_potential(particle_id1=pids[0],
                                          particle_id2=pids[1],
                                          particle_id3=pids[2],
                                          particle_id4=pids[3])

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")

    def test_dihedral_with_partial_bonds_raises_error(self):
        """
        Test that creating a dihedral with a missing bond in the i-j-k-l chain raises a ValueError.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])
        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])

        pids = []
        for name in ('A', 'B', 'C', 'D'):
            pids.append(pmb.create_particle(name=name,
                                            box_l=box_l,
                                            number_of_particles=1)[0])
        pmb.create_bond(particle_id1=pids[0], particle_id2=pids[1])
        pmb.create_bond(particle_id1=pids[1], particle_id2=pids[2])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        with self.assertRaises(ValueError):
            pmb.create_dihedral_potential(particle_id1=pids[0],
                                          particle_id2=pids[1],
                                          particle_id3=pids[2],
                                          particle_id4=pids[3])

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_generate_dihedrals_for_entity_without_particles_returns(self):
        """
        Test that auto-generating dihedrals for an empty entity returns without changes.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb._generate_dihedrals_for_entity(entity_id=999,
                                           entity_id_col="residue_id")
        self.assertEqual(len(pmb.get_instances_df(pmb_type="dihedral")), 0)

    def test_generate_dihedrals_for_entity_creates_dihedral_from_bonds(self):
        """
        Test that auto-generating dihedrals creates one dihedral from a bonded quadruplet.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])

        pids = self.build_linear_chain(pmb)
        for particle_id in pids:
            pmb.db._update_instance(instance_id=particle_id,
                                    pmb_type="particle",
                                    attribute="residue_id",
                                    value=0)

        pmb._generate_dihedrals_for_entity(entity_id=0,
                                           entity_id_col="residue_id")
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        dihedral_objects = self.get_dihedral_objects(particle_id2=pids[1])
        self.assertEqual(len(dihedral_objects), 1)
        self.assertEqual(len(pmb.get_instances_df(pmb_type="dihedral")), 1)

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_generate_dihedrals_for_entity_skips_missing_templates(self):
        """
        Test that auto-generating dihedrals skips bonded quadruplets without a matching template.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])

        pids = self.build_linear_chain(pmb)
        for particle_id in pids:
            pmb.db._update_instance(instance_id=particle_id,
                                    pmb_type="particle",
                                    attribute="residue_id",
                                    value=0)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        pmb._generate_dihedrals_for_entity(entity_id=0,
                                           entity_id_col="residue_id")

        self.assertEqual(len(self.get_dihedral_objects(particle_id2=pids[1])), 0)
        self.assertEqual(len(pmb.get_instances_df(pmb_type="dihedral")), 0)

        for inst_id in pids:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_create_molecule_with_gen_dihedral_generates_dihedrals(self):
        """
        Test that create_molecule generates backbone dihedrals automatically when gen_dihedral is enabled.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_dihedral_potential(dihedral_type="periodic",
                                     dihedral_parameters=self.periodic_dihedral_params,
                                     particle_quadruplets=[('A', 'B', 'C', 'D')])
        pmb.define_residue(name="R_A", central_bead="A", side_chains=[])
        pmb.define_residue(name="R_B", central_bead="B", side_chains=[])
        pmb.define_residue(name="R_C", central_bead="C", side_chains=[])
        pmb.define_residue(name="R_D", central_bead="D", side_chains=[])
        pmb.define_molecule(name="M_dihedral",
                            residue_list=["R_A", "R_B", "R_C", "R_D"])

        molecule_ids = pmb.create_molecule(name="M_dihedral",
                                           number_of_molecules=1,
                                           box_l=box_l,
                                           backbone_vector=[1.0, 0.0, 0.0],
                                           gen_dihedral=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        self.assertEqual(len(pmb.get_instances_df(pmb_type="dihedral")), 1)

        pmb.delete_instances_in_system(instance_id=molecule_ids[0],
                                       pmb_type="molecule")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_create_nested_residue_with_gen_dihedral_generates_dihedrals(self):
        """
        Test that create_residue generates the dihedral of a chain A-B-C-D formed by nested residues.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_default_dihedral_potential(dihedral_type="periodic",
                                              dihedral_parameters=self.periodic_dihedral_params)
        pmb.define_residue(name="R3", central_bead="C", side_chains=["D"])
        pmb.define_residue(name="R2", central_bead="B", side_chains=["R3"])
        pmb.define_residue(name="R1", central_bead="A", side_chains=["R2"])

        residue_ids = pmb.create_residue(name="R1",
                                         box_l=box_l,
                                         gen_dihedral=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        self.assertEqual(len(pmb.get_instances_df(pmb_type="dihedral")), 1)
        n_engine_dihedrals = sum(1 for particle in espresso_system.part.all()
                                 for bond in particle.bonds if len(bond) == 4)
        self.assertEqual(n_engine_dihedrals, 1)

        pmb.delete_instances_in_system(instance_id=residue_ids,
                                       pmb_type="residue")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def test_create_hydrogel_with_gen_dihedral_generates_chain_dihedrals(self):
        """
        Test that create_hydrogel generates the dihedrals internal to every chain, i.e. mpc-3 per chain,
        and none involving the crosslinker nodes.
        """
        mpc = 5
        pmb = pyMBE.pymbe_library(seed=42)
        bond_length = 0.355 * pmb.units.nm
        pmb.define_particle(name="N", sigma=bond_length, epsilon=1 * pmb.units("reduced_energy"))
        pmb.define_particle(name="C", sigma=bond_length, epsilon=1 * pmb.units("reduced_energy"))
        pmb.define_residue(name="Res", central_bead="C", side_chains=[])
        pmb.define_molecule(name="chain", residue_list=["Res"] * mpc)
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters={"r_0": bond_length,
                                         "k": 400 * pmb.units("reduced_energy / reduced_length**2")},
                        particle_pairs=[["C", "C"], ["C", "N"]])
        pmb.define_default_dihedral_potential(dihedral_type="periodic",
                                              dihedral_parameters={'k': 5 * pmb.units('reduced_energy'),
                                                                   'phi_s': np.pi * pmb.units(''),
                                                                   'mult': 1})
        lattice = DiamondLattice(mpc, bond_length)
        hydrogel_box_l = [lattice.box_l] * 3
        try:
            espresso_system.box_l = hydrogel_box_l
            pmb.initialize_lattice_builder(lattice)
            node_topology = [{"particle_name": "N", "lattice_index": index} for index in lattice.indices]
            chain_topology = [{"node_start": str(lattice.indices[connectivity[0]]),
                               "node_end": str(lattice.indices[connectivity[1]]),
                               "molecule_name": "chain"} for connectivity in lattice.connectivity]
            pmb.define_hydrogel("gel", node_topology, chain_topology)
            pmb.create_hydrogel("gel", hydrogel_box_l, gen_dihedral=True)
            pmb.set_simulation_engine(espresso_system)
            pmb.add_instances_to_engine()

            dihedrals = pmb.db.get_instances(pmb_type="dihedral")
            particle_names = {pid: inst.name for pid, inst in pmb.db.get_instances(pmb_type="particle").items()}
            self.assertEqual(len(dihedrals), len(lattice.connectivity) * (mpc - 3))
            for dihedral in dihedrals.values():
                names = [particle_names[getattr(dihedral, f"particle_id{i}")] for i in (1, 2, 3, 4)]
                self.assertNotIn("N", names, msg="crosslinker-adjacent dihedrals are not generated")
            n_engine_dihedrals = sum(1 for particle in espresso_system.part.all()
                                     for bond in particle.bonds if len(bond) == 4)
            self.assertEqual(n_engine_dihedrals, len(dihedrals))
        finally:
            espresso_system.part.clear()
            espresso_system.bonded_inter.clear()
            espresso_system.box_l = box_l


def iupac_dihedral(p1, p2, p3, p4):
    """
    Independent reference for the dihedral angle of the quadruplet p1-p2-p3-p4
    in the IUPAC/GROMACS convention (cis = 0, trans = pi, right-handed positive).
    """
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))


def periodic_energy(positions, terms):
    """Analytic GROMACS type 1 energy sum_t k_t [1 + cos(n_t phi - phi_s_t)] of a quadruplet, terms = [(k, n, phi_s), ...]."""
    phi = iupac_dihedral(*positions)
    return sum(k * (1 + np.cos(n * phi - phi_s)) for k, n, phi_s in terms)


class DihedralPhysicsTest(ut.TestCase):
    """
    Checks that the dihedral potential set up by pyMBE in ESPResSo is the right
    physical potential, using references computed independently with numpy:
    energies, forces, momentum conservation, energy conservation in MD,
    relaxation to the minimum, symmetry, GROMACS mapping and periodic images.
    """
    K = 5.0
    TERMS = [(1, 0.0), (1, 1.0), (2, 0.3), (3, 2.0), (3, np.pi)]  # (mult, phi_s)

    def setUp(self):
        espresso_system.part.clear()
        espresso_system.bonded_inter.clear()
        espresso_system.integrator.set_vv()
        espresso_system.thermostat.turn_off()
        self.rng = np.random.default_rng(42)

    def build(self, terms):
        """
        Builds A-B-C-D with harmonic backbone bonds and a dihedral with the given
        terms [(k, n, phi_s), ...] and returns (pmb, pids, r_0 in reduced units).
        """
        pmb = pyMBE.pymbe_library(seed=42)
        for name in ('A', 'B', 'C', 'D'):
            pmb.define_particle(name=name, z=0, sigma=0.4 * pmb.units.nm,
                                epsilon=1 * pmb.units('reduced_energy'))
        r_0 = 0.4 * pmb.units.nm
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters={'r_0': r_0,
                                         'k': 400 * pmb.units('reduced_energy / reduced_length**2')},
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_dihedral_potential(
            dihedral_type="periodic",
            dihedral_parameters=[{'k': k * pmb.units('reduced_energy'),
                                  'phi_s': phi_s * pmb.units(''),
                                  'mult': n} for k, n, phi_s in terms],
            particle_quadruplets=[('A', 'B', 'C', 'D')])
        pids = [pmb.create_particle(name=name, box_l=box_l, number_of_particles=1)[0]
                for name in ('A', 'B', 'C', 'D')]
        for i in range(3):
            pmb.create_bond(particle_id1=pids[i], particle_id2=pids[i + 1])
        pmb.create_dihedral_potential(particle_id1=pids[0], particle_id2=pids[1],
                                      particle_id3=pids[2], particle_id4=pids[3])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        return pmb, pids, r_0.m_as('reduced_length')

    @staticmethod
    def geometry(r_0, theta, alpha=1.7, beta=1.9, center=(5., 5., 5.)):
        """
        Positions of A, B, C, D with all three bond lengths exactly r_0 (so the harmonic
        bonds contribute no energy or force). alpha/beta are the polar angles of the
        A and D bonds relative to the B->C axis, theta is the rotation about it.
        """
        B = np.array(center, dtype=float)
        C = B + r_0 * np.array([0., 0., 1.])
        A = B + r_0 * np.array([np.sin(alpha), 0., np.cos(alpha)])
        D = C + r_0 * np.array([np.sin(beta) * np.cos(theta), np.sin(beta) * np.sin(theta), np.cos(beta)])
        return np.array([A, B, C, D])

    @staticmethod
    def place(pids, positions):
        for pid, pos in zip(pids, positions):
            espresso_system.part.by_id(pid).pos = pos

    @staticmethod
    def cleanup(pmb, pids):
        for pid in pids:
            pmb.delete_instances_in_system(instance_id=pid, pmb_type="particle")
        pmb.db.delete_templates(pmb_type="dihedral")
        pmb.db.delete_templates(pmb_type="bond")

    def energy_and_forces(self, pids):
        espresso_system.integrator.run(0, recalc_forces=True)
        energy = espresso_system.analysis.energy()['bonded']
        forces = np.array([espresso_system.part.by_id(pid).f for pid in pids])
        return energy, forces

    def test_energy_matches_analytic_formula(self):
        """ESPResSo energy equals the GROMACS expression k[1 + cos(n*phi - phi_s)] for several (n, phi_s) and geometries."""
        for n, phi_s in self.TERMS:
            pmb, pids, r_0 = self.build([(self.K, n, phi_s)])
            for theta in np.linspace(-2.9, 2.9, 11):
                pos = self.geometry(r_0, theta)
                self.place(pids, pos)
                energy, _ = self.energy_and_forces(pids)
                np.testing.assert_allclose(energy, periodic_energy(pos, [(self.K, n, phi_s)]),
                                           atol=1e-8, err_msg=f"n={n}, phi_s={phi_s}, theta={theta}")
            self.cleanup(pmb, pids)
            espresso_system.part.clear()
            espresso_system.bonded_inter.clear()

    def test_energy_is_not_symmetric_in_sign_of_phi(self):
        """For phi_s not in {0, pi}, phi and -phi differ in energy: guards the sign convention."""
        phi_s = 1.0
        pmb, pids, r_0 = self.build([(self.K, 1, phi_s)])
        energies = []
        for theta in (0.8, -0.8):
            pos = self.geometry(r_0, theta)
            self.place(pids, pos)
            energies.append(self.energy_and_forces(pids)[0])
            np.testing.assert_allclose(energies[-1], periodic_energy(pos, [(self.K, 1, phi_s)]), atol=1e-8)
        self.assertGreater(abs(energies[0] - energies[1]), 0.1)
        self.cleanup(pmb, pids)

    def test_multi_term_energy_is_the_sum_of_terms(self):
        """A stacked multi-term dihedral adds the energies of its terms."""
        terms = [(4.0, 1, 0.0), (2.0, 3, np.pi), (1.5, 2, 0.7)]
        pmb, pids, r_0 = self.build(terms)
        for theta in np.linspace(-2.8, 2.8, 9):
            pos = self.geometry(r_0, theta)
            self.place(pids, pos)
            np.testing.assert_allclose(self.energy_and_forces(pids)[0], periodic_energy(pos, terms), atol=1e-8)
        self.cleanup(pmb, pids)

    def test_forces_match_numerical_gradient(self):
        """ESPResSo forces equal -grad of the analytic energy (central finite differences)."""
        terms = [(self.K, 2, 0.3), (1.5, 1, 1.0)]
        pmb, pids, r_0 = self.build(terms)
        h = 1e-6
        for theta in (-2.3, -1.1, 0.7, 1.9, 2.6):
            pos = self.geometry(r_0, theta, alpha=1.5, beta=2.0)
            self.place(pids, pos)
            _, forces = self.energy_and_forces(pids)
            reference = np.zeros_like(pos)
            for i in range(4):
                for d in range(3):
                    plus, minus = pos.copy(), pos.copy()
                    plus[i, d] += h
                    minus[i, d] -= h
                    reference[i, d] = -(periodic_energy(plus, terms) - periodic_energy(minus, terms)) / (2 * h)
            np.testing.assert_allclose(forces, reference, atol=1e-5, err_msg=f"theta={theta}")
        self.cleanup(pmb, pids)

    def test_forces_conserve_momentum_and_angular_momentum(self):
        """Net force and net torque of an internal torsion vanish."""
        pmb, pids, r_0 = self.build([(self.K, 3, 0.4)])
        for theta in (-2.0, 0.6, 1.4):
            pos = self.geometry(r_0, theta, alpha=1.3, beta=2.2)
            self.place(pids, pos)
            _, forces = self.energy_and_forces(pids)
            self.assertGreater(np.abs(forces).max(), 1e-3, "test needs non-zero forces to be meaningful")
            np.testing.assert_allclose(forces.sum(axis=0), 0., atol=1e-9)
            np.testing.assert_allclose(np.cross(pos, forces).sum(axis=0), 0., atol=1e-9)
        self.cleanup(pmb, pids)

    def test_nve_energy_conservation(self):
        """Velocity-Verlet NVE dynamics conserves total energy (forces consistent with the energy)."""
        pmb, pids, r_0 = self.build([(self.K, 2, 0.3), (2.0, 1, 1.0)])
        self.place(pids, self.geometry(r_0, 1.2, alpha=1.7, beta=1.9))
        for pid in pids:
            espresso_system.part.by_id(pid).v = 0.3 * self.rng.normal(size=3)
        espresso_system.integrator.run(0, recalc_forces=True)
        e_0 = espresso_system.analysis.energy()['total']
        drift = 0.
        for _ in range(40):
            espresso_system.integrator.run(100)
            drift = max(drift, abs(espresso_system.analysis.energy()['total'] - e_0))
        self.assertLess(drift, 5e-3 * max(1., abs(e_0)), msg=f"energy drift {drift} from E_0={e_0}")
        self.cleanup(pmb, pids)

    def test_relaxation_reaches_the_minimum_with_the_right_sign(self):
        """
        Steepest descent from a twisted conformation ends at the minimum of the GROMACS
        type 1 potential, phi_min = phi_s + pi for n = 1 (IUPAC angle), not at -phi_min.
        This depends on the sign convention and on the forces, not on the energy formula.
        """
        for phi_min in (2.0, -2.0):
            phi_s = phi_min - np.pi if phi_min > 0 else phi_min + np.pi
            pmb, pids, r_0 = self.build([(self.K, 1, phi_s)])
            start = 1.0 if phi_min > 0 else -1.0
            self.place(pids, self.geometry(r_0, start))
            espresso_system.integrator.set_steepest_descent(f_max=1e-4, gamma=1e-3, max_displacement=0.01)  # gamma*k_bond << 2 for stability
            espresso_system.integrator.run(5000)
            espresso_system.integrator.set_vv()
            pos = np.array([espresso_system.part.by_id(pid).pos for pid in pids])
            # unfold the chain around B (bond lengths are far below half the box)
            pos = pos[1] + (pos - pos[1] + np.array(box_l) / 2) % np.array(box_l) - np.array(box_l) / 2
            phi = iupac_dihedral(*pos)
            self.assertAlmostEqual(phi, phi_min, places=3)
            self.assertLess(espresso_system.analysis.energy()['bonded'], 1e-6)
            self.cleanup(pmb, pids)
            espresso_system.part.clear()
            espresso_system.bonded_inter.clear()

    def test_reversed_quadruplet_gives_the_same_physics(self):
        """Creating the dihedral on (D, C, B, A) hits the same template and gives the same energy and forces."""
        terms = [(self.K, 2, 0.3), (1.5, 1, 1.0)]
        pmb, pids, r_0 = self.build(terms)
        pos = self.geometry(r_0, 0.9, alpha=1.4, beta=2.1)
        self.place(pids, pos)
        e_forward, f_forward = self.energy_and_forces(pids)
        self.cleanup(pmb, pids)
        espresso_system.part.clear()
        espresso_system.bonded_inter.clear()

        pmb = pyMBE.pymbe_library(seed=42)
        for name in ('A', 'B', 'C', 'D'):
            pmb.define_particle(name=name, z=0, sigma=0.4 * pmb.units.nm,
                                epsilon=1 * pmb.units('reduced_energy'))
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters={'r_0': 0.4 * pmb.units.nm,
                                         'k': 400 * pmb.units('reduced_energy / reduced_length**2')},
                        particle_pairs=[['A', 'B'], ['B', 'C'], ['C', 'D']])
        pmb.define_dihedral_potential(
            dihedral_type="periodic",
            dihedral_parameters=[{'k': k * pmb.units('reduced_energy'), 'phi_s': phi_s * pmb.units(''), 'mult': n}
                                 for k, n, phi_s in terms],
            particle_quadruplets=[('A', 'B', 'C', 'D')])
        pids = [pmb.create_particle(name=name, box_l=box_l, number_of_particles=1)[0]
                for name in ('A', 'B', 'C', 'D')]
        for i in range(3):
            pmb.create_bond(particle_id1=pids[i], particle_id2=pids[i + 1])
        pmb.create_dihedral_potential(particle_id1=pids[3], particle_id2=pids[2],
                                      particle_id3=pids[1], particle_id4=pids[0])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        self.place(pids, pos)
        e_reverse, f_reverse = self.energy_and_forces(pids)
        np.testing.assert_allclose(e_reverse, e_forward, atol=1e-10)
        np.testing.assert_allclose(f_reverse, f_forward, atol=1e-9)
        np.testing.assert_allclose(e_reverse, periodic_energy(pos, terms), atol=1e-8)
        self.cleanup(pmb, pids)

    def test_gromacs_proper_dihedral_mapping(self):
        """
        A GROMACS proper dihedral (type 1), V = k[1 + cos(n*phi - phi_s)], given to pyMBE exactly
        as in a topology (k, mult = n, phi_s in radians) is reproduced by the ESPResSo engine.
        """
        for n in (1, 2, 3):
            for phi_s_deg in (0., 180., 60., -90.):
                phi_s = np.radians(phi_s_deg)
                pmb, pids, r_0 = self.build([(self.K, n, phi_s)])
                for theta in (-2.4, -0.9, 0.5, 1.6, 2.7):
                    pos = self.geometry(r_0, theta)
                    self.place(pids, pos)
                    gromacs = self.K * (1 + np.cos(n * iupac_dihedral(*pos) - phi_s))
                    np.testing.assert_allclose(self.energy_and_forces(pids)[0], gromacs, atol=1e-8,
                                               err_msg=f"n={n}, phi_s={phi_s_deg}, theta={theta}")
                self.cleanup(pmb, pids)
                espresso_system.part.clear()
                espresso_system.bonded_inter.clear()

    def test_dihedral_across_periodic_boundary(self):
        """A torsion whose particles are on opposite sides of the box has the same energy and forces."""
        terms = [(self.K, 2, 0.3)]
        pmb, pids, r_0 = self.build(terms)
        results = []
        for center in ((5., 5., 5.), (box_l[0] - 0.4 * r_0, 0.3 * r_0, box_l[2] - 0.2 * r_0)):
            pos = self.geometry(r_0, 0.9, alpha=1.4, beta=2.1, center=center)
            self.place(pids, pos)
            if center[0] != 5.:
                self.assertTrue(np.any(pos > box_l[0]) or np.any(pos < 0.),
                                "test geometry should straddle the periodic boundary")
            results.append(self.energy_and_forces(pids))
            np.testing.assert_allclose(results[-1][0], periodic_energy(pos, terms), atol=1e-8)
        np.testing.assert_allclose(results[1][0], results[0][0], atol=1e-8)
        np.testing.assert_allclose(results[1][1], results[0][1], atol=1e-8)
        self.cleanup(pmb, pids)

    def count_generated_dihedrals(self, n_particles, bonds):
        """Builds a bonded graph with a default dihedral and returns the generated dihedral quadruplets."""
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_particle(name='A', z=0, sigma=0.4 * pmb.units.nm, epsilon=1 * pmb.units('reduced_energy'))
        pmb.define_default_bond(bond_type="harmonic",
                                bond_parameters={'r_0': 0.4 * pmb.units.nm,
                                                 'k': 400 * pmb.units('reduced_energy / reduced_length**2')})
        pmb.define_default_dihedral_potential(dihedral_type="periodic",
                                              dihedral_parameters={'k': 5 * pmb.units('reduced_energy'),
                                                                   'phi_s': np.pi * pmb.units(''),
                                                                   'mult': 1})
        pids = pmb.create_particle(name='A', box_l=box_l, number_of_particles=n_particles)
        for i, j in bonds:
            pmb.create_bond(particle_id1=pids[i], particle_id2=pids[j], use_default_bond=True)
            pmb.db._update_instance(instance_id=pids[i], pmb_type="particle", attribute="residue_id", value=0)
            pmb.db._update_instance(instance_id=pids[j], pmb_type="particle", attribute="residue_id", value=0)
        pmb._generate_dihedrals_for_entity(entity_id=0, entity_id_col="residue_id")
        index = {pid: i for i, pid in enumerate(pids)}
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        found = [tuple(index[getattr(d, f"particle_id{k}")] for k in (1, 2, 3, 4))
                 for d in pmb.db.get_instances(pmb_type="dihedral").values()]
        self.cleanup(pmb, pids)
        return found

    @staticmethod
    def reference_torsions(n_particles, bonds):
        """Independent count: simple 4-vertex paths over the bond graph, a path and its reverse being one torsion."""
        edges = {frozenset(b) for b in bonds}
        paths = {p for p in itertools.permutations(range(n_particles), 4)
                 if all(frozenset(p[i:i + 2]) in edges for i in range(3))}
        return {min(p, p[::-1]) for p in paths}

    def test_gen_dihedral_matches_graph_enumeration(self):
        """Auto-generation finds exactly the torsions of the bond graph: chain, branch, rings, K4."""
        topologies = {"linear chain of 6": (6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]),
                      "branched": (5, [(0, 1), (1, 2), (2, 3), (2, 4)]),
                      "4-ring": (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
                      "triangle (no torsion)": (3, [(0, 1), (1, 2), (2, 0)]),
                      "K4": (4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])}
        for label, (n_particles, bonds) in topologies.items():
            found = self.count_generated_dihedrals(n_particles, bonds)
            reference = self.reference_torsions(n_particles, bonds)
            self.assertEqual(len(found), len(reference), msg=f"{label}: number of dihedrals")
            self.assertEqual({min(q, q[::-1]) for q in found}, reference, msg=f"{label}: quadruplets")
            espresso_system.part.clear()


if __name__ == '__main__':
    ut.main()
