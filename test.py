
import pyMBE
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
from pyMBE.storage.templates.lj import LJInteractionTemplate

from pyMBE.lib.lattice import DiamondLattice
import importlib.resources

import pyMBE.storage.io as io

import pint 
import scipy
import espressomd

def main():

    units = pint.UnitRegistry()
    unit_length= 0.355*units.nm
    temperature = 298.15 * units.K
    kB=scipy.constants.k * units.J / units.K
    kT=temperature*kB
    units.define(f'reduced_energy = {kT} ')
    units.define(f'reduced_length = {unit_length}')   
    espresso_system=espressomd.System (box_l = [10]*3)
    db = Manager(units=units)

    # ============================================================
    # 1. CREATE PARTICLE TEMPLATES + STATES
    # ============================================================
    
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
                        epsilon=0.2 * units.reduced_energy)
    
    pmb.define_residue(name="R1", central_bead="Z", side_chains=["X","Z"])
    pmb.define_residue(name="R2", central_bead="Z", side_chains=["X","R1"])
    
    print("\n=== Residue Templates DataFrame ===")
    print(pmb.db._get_templates_df(pmb_type="residue"))

    pmb.define_molecule(name="M1", residue_list=["R1","R2"])
    print("\n=== Molecule Templates DataFrame ===")
    print(pmb.db._get_templates_df(pmb_type="molecule"))


    print("\n=== Hydrogel Templates DataFrame ===")
    diamond_lattice = DiamondLattice(30, 3.5 * units.reduced_length)
    lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

    # Setting up node topology
    indices = diamond_lattice.indices
    node_topology = []

    for index in range(len(indices)):
        node_topology.append({"particle_name": "A",
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

    pmb.define_hydrogel("my_hydrogel", node_topology, chain_topology)
    print(pmb.db._get_templates_df(pmb_type="hydrogel"))

    print("\n=== Molecule Templates DataFrame ===")
    print(pmb.db._get_templates_df(pmb_type="molecule"))

    tpl = LJInteractionTemplate(state1  = "A",
                                state2  = "AH",
                                sigma   = PintQuantity.from_quantity(q=1.0  * units.reduced_length, expected_dimension="length", ureg=units),
                                cutoff  = PintQuantity.from_quantity(q=1.2  * units.reduced_length, expected_dimension="length", ureg=units),
                                offset  = PintQuantity.from_quantity(q=0  * units.reduced_length, expected_dimension="length", ureg=units),
                                epsilon = PintQuantity.from_quantity(q=1.0  * units.reduced_energy, expected_dimension="energy", ureg=units),
                                shift   = "auto"
                                )
    db._register_template(tpl)
    print(db._get_templates_df(pmb_type="lj"))                     
    
    print("\n=== Particle Templates DataFrame ===")
    print(pmb.db._get_templates_df(pmb_type="particle"))
    print(pmb.db._get_reactions_df())

    # Update reaction
    """
    pmb.db._update_reaction_participant(reaction_name="AH <-> A",
                                        particle_name="H",
                                        state_name="H",
                                        coefficient=1)
    print(pmb.db._get_reactions_df())
    """


    
    parameters = {"k":  100.0 * units.reduced_energy / (units.reduced_length**2),
                  "r_0": 1.0  * units.reduced_length}
    
    
    pmb.define_bond(bond_type="harmonic",
                    bond_parameters=parameters,
                    particle_pairs=[["Z","Z"], 
                                    ["Z","X"],
                                    ["X","X"]])
    
    pmb.define_default_bond(bond_type="harmonic",
                            bond_parameters=parameters)
    
    print("\n=== Bond Templates DataFrame ===")
    print(pmb.db._get_templates_df(pmb_type="bond"))

    print("\n=== Peptide Templates DataFrame ===")
    pmb.define_peptide(name="Peptide1",
                       model="1beadAA",
                       sequence="KKKKDDDD")
    
    print(pmb.db._get_templates_df(pmb_type="peptide"))

    print("\n=== Protein Templates DataFrame ===")
    path = importlib.resources.files(pyMBE) / "parameters" / "globular_proteins" / f"1beb.vtf",
    
    topology_dict = pmb.read_protein_vtf_in_df (filename=path[0])

    pmb.define_protein(name="blabla",
                       model="2beadAA",
                       sequence="KKKKKK")
    
    print(db._get_templates_df(pmb_type="protein"))

    
    

    # ============================================================
    # 2. CREATE INSTANCES (optional for testing)
    # ============================================================
    pmb.create_particle(name="Z",
                        espresso_system=espresso_system,
                        number_of_particles=3)
    pmb.create_particle(name="X",
                        espresso_system=espresso_system,
                        number_of_particles=1)
    
    print("\n=== Particle Instances DataFrame ===")
    print(pmb.db._get_instances_df(pmb_type="particle"))


    pmb.db._update_instance(pmb_type="particle", instance_id=1, attribute="residue_id", value=int(0))
    print("\n=== Particle Instances DataFrame (after update) ===")
    print(pmb.db._get_instances_df(pmb_type="particle"))

    pmb.create_residue(name="R1",
                       espresso_system=espresso_system)
    pmb.create_residue(name="R2",
                       espresso_system=espresso_system)

    print("\n=== Residue Instances DataFrame ===")
    print(pmb.db._get_instances_df(pmb_type="residue"))


    pmb.db._update_instance(pmb_type="residue",instance_id=0, attribute="molecule_id", value=int(0))
    print("\n=== Residue Instances DataFrame (after update)===")
    print(pmb.db._get_instances_df(pmb_type="residue"))


    inst1 = MoleculeInstance(name="M1", molecule_id=1)
    inst2 = MoleculeInstance(name="M1", molecule_id=2)
#    db._register_instance(inst1)
#    db._register_instance(inst2)
    print("\n=== Molecule Instances DataFrame ===")
    print(db._get_instances_df(pmb_type="molecule"))

    
    print("\n=== Bond Instances DataFrame ===")
    print(pmb.db._get_instances_df(pmb_type="bond"))     

    print("\n=== Peptide Instances DataFrame ===")
 #   inst_peptide1 = PeptideInstance(name="Peptide1", molecule_id=3)
 #   db._register_instance(inst_peptide1)
    print(db._get_instances_df(pmb_type="peptide"))

    print("\n=== Protein Instances DataFrame ===")
 #   inst_protein1 = ProteinInstance(name="Protein1", molecule_id=4)
 #   db._register_instance(inst_protein1)
    print(db._get_instances_df(pmb_type="protein"))

    print("\n=== Hydrogel Instances DataFrame ===")
 #   inst_hydrogel1 = HydrogelInstance(name="Hydrogel1", hydrogel_id=1, molecule_ids=["1","2","3"])
 #   db._register_instance(inst_hydrogel1)
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
    for ptype, tdict in db._templates.items():
        for tname, t in tdict.items():
            db2._register_template(t)

    print("\nTemplates shown with registry 2 (different reduced units):")
    print(db2._get_templates_df("particle"))

    io._save_database_csv(pmb.db, folder="test_db_csv")

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
    print(db3._get_templates_df("lj"))
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

