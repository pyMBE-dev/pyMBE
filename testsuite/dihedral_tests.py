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
import numpy as np
import unittest as ut
import espressomd

# Create an instance of the ESPResSo system
box_l=[10]*3
espresso_system = espressomd.System(box_l=box_l)


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
            'phi_0': np.pi * pmb.units(''),
            'mult': 1,
        }

        self.periodic_dihedral_params_alt = {
            'k': 3 * pmb.units('reduced_energy'),
            'phi_0': np.pi / 2 * pmb.units(''),
            'mult': 2,
        }

        self.multi_term_dihedral_params = [
            {'k': 4 * pmb.units('reduced_energy'), 'phi_0': 0.0 * pmb.units(''), 'mult': 1},
            {'k': 2 * pmb.units('reduced_energy'), 'phi_0': np.pi * pmb.units(''), 'mult': 3},
        ]

        self.harmonic_bond_params = {
            'r_0': 0.4 * pmb.units.nm,
            'k': 400 * pmb.units('reduced_energy / reduced_length**2'),
        }

    def get_dihedral_objects(self, particle_id1):
        """
        Returns the list of dihedral interaction objects stored in ESPResSo for a
        given first particle. Dihedral bonds are added to the first particle in
        ESPResSo and carry three partner particles (tuple arity 4).
        """
        bonds = espresso_system.part.by_id(particle_id1).bonds
        dihedral_objects = []
        for bond_tuple in bonds:
            if len(bond_tuple) == 4:
                dihedral_objects.append(bond_tuple[0])
        return dihedral_objects

    def check_dihedral_setup(self, dihedral_objects, input_terms):
        """
        Checks that pyMBE sets up the dihedral interaction objects correctly.

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
                           round(float(term['phi_0'].magnitude), 10)) for term in input_terms)
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

        dihedral_objects = self.get_dihedral_objects(particle_id1=pids[0])
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

        dihedral_objects = self.get_dihedral_objects(particle_id1=pids[0])
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

        dihedral_objects = self.get_dihedral_objects(particle_id1=pids[0])
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

        dihedral_objects = self.get_dihedral_objects(particle_id1=pids[0])
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
                     'phi_0': np.pi * pmb.units(''),
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
                                 **call('periodic', {'phi_0': np.pi * pmb.units(''), 'mult': 1}))
        # missing phi_0
        np.testing.assert_raises(ValueError, callback,
                                 **call('periodic', {'k': 5 * pmb.units('reduced_energy'), 'mult': 1}))
        # missing mult
        np.testing.assert_raises(ValueError, callback,
                                 **call('periodic', {'k': 5 * pmb.units('reduced_energy'),
                                                     'phi_0': np.pi * pmb.units('')}))
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

        dihedral_objects = self.get_dihedral_objects(particle_id1=pids[0])
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

        self.assertEqual(len(self.get_dihedral_objects(particle_id1=pids[0])), 0)
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


if __name__ == '__main__':
    ut.main()
