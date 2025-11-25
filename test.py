# test.py
from pyMBE.storage.df_management import _DFManagement
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleState
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
from pyMBE.storage.pint_quantity import PintQuantity

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

    db = _DFManagement(units=units)


    # ============================================================
    # 1. CREATE PARTICLE TEMPLATES + STATES
    # ============================================================
    
    
    # A particle (acid)
    tpl_A = ParticleTemplate(name="A", sigma=PintQuantity.from_quantity(q=3.5 * units.reduced_length, expected_dimension="length", ureg=units),
                                cutoff=PintQuantity.from_quantity(q=4 * units.reduced_length, expected_dimension="length", ureg=units),
                                offset=PintQuantity.from_quantity(q=0 * units.reduced_length, expected_dimension="length", ureg=units),
                                epsilon=PintQuantity.from_quantity(q=0.2 * units.reduced_energy, expected_dimension="energy", ureg=units))
                             
    tpl_A.add_state(ParticleState(name="HA", z=0, es_type=0))
    tpl_A.add_state(ParticleState(name="A-", z=-1, es_type=1))

    # H+ particle (single-state)
    tpl_H = ParticleTemplate(name="H", sigma=PintQuantity(magnitude=3.5, units="nm", dimension="length"),
                             cutoff=PintQuantity.from_quantity(q=4 * units.reduced_length, expected_dimension="length", ureg=units),
                                offset=PintQuantity.from_quantity(q=0 * units.reduced_length, expected_dimension="length", ureg=units),
                                epsilon=PintQuantity(magnitude=0.2, units="J", dimension="energy"))
    tpl_H.add_state(ParticleState(name="H+", z=+1, es_type=2))

    # Register templates
    db.register_template(tpl_A)
    db.register_template(tpl_H)

    # ============================================================
    # 2. CREATE INSTANCES (optional for testing)
    # ============================================================

    inst1 = ParticleInstance(name="A", particle_id=1, initial_state="HA")
    inst2 = ParticleInstance(name="A", particle_id=2, initial_state="A-",residue_id=0)
    inst3 = ParticleInstance(name="H", particle_id=3, initial_state="H+")

    db.register_instance(inst1)
    db.register_instance(inst2)
    db.register_instance(inst3)
    
    print("\n=== Instances DataFrame ===")
    print(db.get_instances_df())


    db.update_particle_instance(particle_id=1, attribute="residue_id", value=int(0))

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

    db.register_reaction(rx)

    # ============================================================
    # 4. PRINT DATAFRAMES
    # ============================================================

    print("\n=== Templates DataFrame ===")
    print(db.get_templates_df())

    print("\n=== Instances DataFrame ===")
    print(db.get_instances_df())

    print("\n=== Reactions DataFrame ===")
    print(db.get_reactions_df())

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
    db2 = _DFManagement(units=ureg2)
    # re-insert templates by transferring stored representation (simulate loading)
    for name, tpl_obj in db.templates.items():
        db2.register_template(tpl_obj)  # tpl_obj stores SI/base units internally

    print("\nUsing registry 2 (different reduced units):")
    print(db2.get_templates_df())

if __name__ == "__main__":
    main()

