# test.py
from pyMBE.storage.manager import Manager
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleState
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
from pyMBE.storage.pint_quantity import PintQuantity
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.templates.hydrogel import HydrogelTemplate, HydrogelNode, HydrogelChain
from pyMBE.storage.instances.hydrogel import HydrogelInstance

import pyMBE.storage.io as io

import pint 
import scipy.constants

def main():

    units = pint.UnitRegistry()
    unit_length= 0.355*units.nm
    temperature = 298.15 * units.K
    kB=scipy.constants.k * units.J / units.K
    kT=temperature*kB
    units.define(f'reduced_energy = {kT} ')
    units.define(f'reduced_length = {unit_length}')   

    db = Manager(units=units)

    # ============================================================
    # 1. CREATE PARTICLE TEMPLATES + STATES
    # ============================================================
    
    # A particle (acid)
    tpl_A = ParticleTemplate(name="A", 
                             sigma=PintQuantity.from_quantity(q=3.5 * units.reduced_length, expected_dimension="length", ureg=units),
                             cutoff=PintQuantity.from_quantity(q=4 * units.reduced_length, expected_dimension="length", ureg=units),
                             offset=PintQuantity.from_quantity(q=0 * units.reduced_length, expected_dimension="length", ureg=units),
                             epsilon=PintQuantity.from_quantity(q=0.2 * units.reduced_energy, expected_dimension="energy", ureg=units))
                             
    tpl_A.add_state(ParticleState(name="HA", z=0, es_type=0))
    tpl_A.add_state(ParticleState(name="A-", z=-1, es_type=1))

    # H+ particle (single-state)
    tpl_H = ParticleTemplate(name="H", sigma=PintQuantity.from_quantity(q=3.5 * units.reduced_length, expected_dimension="length", ureg=units),
                             cutoff=PintQuantity.from_quantity(q=4 * units.reduced_length, expected_dimension="length", ureg=units),
                                offset=PintQuantity.from_quantity(q=0 * units.reduced_length, expected_dimension="length", ureg=units),
                                epsilon=PintQuantity(magnitude=0.2, units="J", dimension="energy"))
    tpl_H.add_state(ParticleState(name="H+", z=+1, es_type=2))

    # Register templates
    db._register_template(tpl_A)
    db._register_template(tpl_H)
    print("\n=== Particle Templates DataFrame ===")
    print(db._get_templates_df(pmb_type="particle"))


    tpl_R1 = ResidueTemplate(name="R1", central_bead="A", side_chains=["H","A"])
    tpl_R2 = ResidueTemplate(name="R2", central_bead="HA", side_chains=["H","HA"])
    db._register_template(tpl_R1)
    db._register_template(tpl_R2)
    print("\n=== Residue Templates DataFrame ===")
    print(db._get_templates_df(pmb_type="residue"))


    tpl_M1 = MoleculeTemplate(name="M1", residue_list=["R1","R2"])
    db._register_template(tpl_M1)
    print("\n=== Molecule Templates DataFrame ===")
    print(db._get_templates_df(pmb_type="molecule"))

    parameters = {"k":  PintQuantity.from_quantity(q=100.0 * units.reduced_energy / (units.reduced_length**2), expected_dimension="energy/length**2", ureg=units),
                  "r0": PintQuantity.from_quantity(q=1.0  * units.reduced_length, expected_dimension="length", ureg=units),}
    
    tpl_bond = BondTemplate(name="A1-A2",
                            bond_type="harmonic", 
                            parameters=parameters,
                            l0=PintQuantity.from_quantity(q=1.0 * units.reduced_length, 
                                                          expected_dimension="length", 
                                                          ureg=units))
    db._register_template(tpl_bond)
    print("\n=== Bond Templates DataFrame ===")
    print(db._get_templates_df(pmb_type="bond"))

    print("\n=== Peptide Templates DataFrame ===")
    tpl_P1 = PeptideTemplate(name="Peptide1",
                             model="Model1",
                             residue_list=["R1","R2"],
                             sequence=["R1","R2"])
    db._register_template(tpl_P1)
    print(db._get_templates_df(pmb_type="peptide"))

    print("\n=== Protein Templates DataFrame ===")
    tpl_PR1 = ProteinTemplate(name="Protein1",
                             model="ModelP1",
                             residue_list=["R1","R2"],
                             sequence=["R1","R2"])
    db._register_template(tpl_PR1)
    print(db._get_templates_df(pmb_type="protein"))

    
    print("\n=== Hydrogel Templates DataFrame ===")
    node1 = HydrogelNode(particle_name="A", lattice_index=[0,0,0])
    node2 = HydrogelNode(particle_name="HA", lattice_index=[1,0,0])
    chain1 = HydrogelChain(node_start="A", node_end="HA", residue_list=["R1","R2"])
    tpl_HG1 = HydrogelTemplate(name="Hydrogel1",
                                 node_map=[node1, node2],
                                 chain_map=[chain1])
    db._register_template(tpl_HG1)
    print(db._get_templates_df(pmb_type="hydrogel"))

    # ============================================================
    # 2. CREATE INSTANCES (optional for testing)
    # ============================================================

    inst1 = ParticleInstance(name="A", particle_id=1, initial_state="HA")
    inst2 = ParticleInstance(name="A", particle_id=2, initial_state="A-",residue_id=0)
    inst3 = ParticleInstance(name="H", particle_id=3, initial_state="H+")

    db._register_instance(inst1)
    db._register_instance(inst2)
    db._register_instance(inst3)
    
    print("\n=== Particle Instances DataFrame ===")
    print(db._get_instances_df(pmb_type="particle"))


    db._update_instance(pmb_type="particle", instance_id=1, attribute="residue_id", value=int(0))
    print("\n=== Particle Instances DataFrame (after update) ===")
    print(db._get_instances_df(pmb_type="particle"))

    inst1 = ResidueInstance(name="R1", 
                            residue_id=1)
    inst2 = ResidueInstance(name="R2", 
                            residue_id=2)
    inst3 = ResidueInstance(name="R1", 
                            residue_id=3, 
                            molecule_id=0)

    db._register_instance(inst1)
    db._register_instance(inst2)
    db._register_instance(inst3)
    
    print("\n=== Residue Instances DataFrame ===")
    print(db._get_instances_df(pmb_type="residue"))


    db._update_instance(pmb_type="residue",instance_id=1, attribute="molecule_id", value=int(0))
    print("\n=== Residue Instances DataFrame (after update)===")
    print(db._get_instances_df(pmb_type="residue"))


    inst1 = MoleculeInstance(name="M1", molecule_id=1)
    inst2 = MoleculeInstance(name="M1", molecule_id=2)
    db._register_instance(inst1)
    db._register_instance(inst2)
    print("\n=== Molecule Instances DataFrame ===")
    print(db._get_instances_df(pmb_type="molecule"))

    inst_bond = BondInstance(name="A1-A2", bond_id=1, particle_id1=1, particle_id2=2)
    db._register_instance(inst_bond)
    print("\n=== Bond Instances DataFrame ===")
    print(db._get_instances_df(pmb_type="bond"))     

    print("\n=== Peptide Instances DataFrame ===")
    inst_peptide1 = PeptideInstance(name="Peptide1", molecule_id=3)
    db._register_instance(inst_peptide1)
    print(db._get_instances_df(pmb_type="peptide"))

    print("\n=== Protein Instances DataFrame ===")
    inst_protein1 = ProteinInstance(name="Protein1", molecule_id=4)
    db._register_instance(inst_protein1)
    print(db._get_instances_df(pmb_type="protein"))

    print("\n=== Hydrogel Instances DataFrame ===")
    inst_hydrogel1 = HydrogelInstance(name="Hydrogel1", hydrogel_id=1, molecule_ids=["1","2","3"])
    db._register_instance(inst_hydrogel1)
    print(db._get_instances_df(pmb_type="hydrogel"))


    # ============================================================
    # 3. DEFINE A REACTION:  HA <-> A- + H+
    # ============================================================

    rx = Reaction(
        name="acid_dissociation",
        pK=4.75,
        reaction_type="acid/base",
        participants=[
            ReactionParticipant(particle_name="A", state_name="HA", coefficient=-1),
            ReactionParticipant(particle_name="A", state_name="A-", coefficient=+1),
            ReactionParticipant(particle_name="H", state_name="H+", coefficient=+1),
        ],
    )

    db._register_reaction(rx)

    # ============================================================
    # 4. PRINT DATAFRAMES
    # ============================================================

    
    print("\n=== Instances DataFrame ===")
    print(db._get_instances_df(pmb_type="particle"))

    print("\n=== Reactions DataFrame ===")
    print(db._get_reactions_df())

    # -------------------------
    # Now create a different registry with different reduced unit definitions
    # and re-create a DFManager with that registry. The DB still stores SI values,
    # so conversions are consistent.
    # -------------------------
    ureg2 = pint.UnitRegistry()
    # define different reduced units (different numeric size)
    unit_length2 = 0.2 * ureg2.nanometer
    temperature2 = 310.0 * ureg2.kelvin
    kB2 = scipy.constants.k * ureg2.joule / ureg2.kelvin
    kT2 = temperature2 * kB2
    ureg2.define(f"reduced_length = {unit_length2}")
    ureg2.define(f"reduced_energy = {kT2}")

    # create a new DFManager that uses the same stored templates but different ureg
    # For this demo we will copy the stored templates (in real use you would re-load from serialized storage)
    db2 = Manager(units=ureg2)
    # re-insert templates by transferring stored representation (simulate loading)
    for ptype, tdict in db.templates.items():
        for tname, t in tdict.items():
            db2._register_template(t)

    print("\nTemplates shown with registry 2 (different reduced units):")
    print(db2._get_templates_df("particle"))

    io._save_database_csv(db, folder="test_db_csv")

    db3 = Manager(units=ureg2)

    io._load_database_csv(db3, folder="test_db_csv")
    print("\nLoaded DB3 Templates DataFrame:")
    print(db3._get_templates_df("particle"))
    print(db3._get_templates_df("residue"))
    print(db3._get_templates_df("molecule"))
    print(db3._get_templates_df("bond"))
    print(db3._get_templates_df("peptide"))
    print(db3._get_templates_df("protein"))
    print(db3._get_templates_df("hydrogel"))
    print("\nLoaded DB3 Instances DataFrame:")
    print(db3._get_instances_df("particle"))
    print(db3._get_instances_df("residue"))
    print(db3._get_instances_df("molecule"))
    print(db3._get_instances_df("bond"))
    print(db3._get_instances_df("peptide"))
    print(db3._get_instances_df("protein"))
    print(db3._get_instances_df("hydrogel"))
    print("\nLoaded DB3 Reactions DataFrame:")
    print(db3._get_reactions_df())


if __name__ == "__main__":
    main()

