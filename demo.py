import pyMBE
import espressomd
from pyMBE.lib.lattice import DiamondLattice

# Setup
pmb = pyMBE.pymbe_library(seed=42)
units = pmb.units
espresso_system = espressomd.System(box_l=[10, 10, 10])
# Define some particle templates

pmb.define_particle(name="Z",
                    sigma=3.5 * units.reduced_length,
                    cutoff=4  * units.reduced_length,
                    offset=0  * units.reduced_length,
                    epsilon=0.2 * units.reduced_energy,
                    acidity="acidic",
                    pka=4.25)

pmb.define_particle(name="X",
                    sigma=3.5 * units.reduced_length,
                    cutoff=4 * units.reduced_length,
                    offset=0 * units.reduced_length,
                    epsilon=0.2 * units.reduced_energy,
                    z=1)

print("\n=== Particle Templates DataFrame ===")
print(pmb.get_templates_df(pmb_type="particle"))

# Access some data in the database
tpl_particle_Z = pmb.db.get_template(name="Z", pmb_type="particle")
tpl_particle_X = pmb.db.get_template(name="X", pmb_type="particle")

# PintQuantity usage example
print("\n=== PintQuantity Usage Example ===")
print(f"PintQuantity class stored in the pyMBE database: {tpl_particle_Z.sigma}")
# Convert to Pint Quantity
sigma_Z = tpl_particle_Z.sigma.to_quantity(units)
print(f"Converted sigma_Z: {sigma_Z} ({sigma_Z.to('reduced_length')})")
# Operate with Pint Quantity
sigma_X = tpl_particle_X.sigma.to_quantity(units)
print(sigma_Z+sigma_X)

# Setup LJ interactions
pmb.setup_lj_interactions(espresso_system=espresso_system)
print("\n=== LJ Templates DataFrame ===")
print(pmb.get_templates_df(pmb_type="lj"))

# Create instances of particles
pmb.create_particle(name="Z",
                        espresso_system=espresso_system,
                        number_of_particles=3)
pmb.create_particle(name="X",
                    espresso_system=espresso_system,
                    number_of_particles=1)
print("\n=== Particle Instances DataFrame ===")
print(pmb.get_instances_df(pmb_type="particle"))

# Delete instances of particles 0-2
for i in range(3):
    pmb.delete_instances_in_system(espresso_system=espresso_system,
                                pmb_type="particle",
                                instance_id=i)

print("\n=== Particle Instances DataFrame After Deletion ===")
print(pmb.get_instances_df(pmb_type="particle"))
pmb.delete_instances_in_system(espresso_system=espresso_system,
                                pmb_type="particle",
                                instance_id=3)

# Create residue
##  Define residues and bonds
pmb.define_residue(name="R1", central_bead="Z", side_chains=["X","Z"])
parameters = {"k":  100.0 * units.reduced_energy / (units.reduced_length**2),
                  "r_0": 1.0  * units.reduced_length}
pmb.define_bond(bond_type="harmonic",
                bond_parameters=parameters,
                particle_pairs=[["Z","Z"], 
                                ["Z","X"],
                                ["X","X"]])

print("\n=== Residue Templates DataFrame ===")
print(pmb.get_templates_df(pmb_type="residue"))
print("\n=== Bond Templates DataFrame ===")
print(pmb.get_templates_df(pmb_type="bond"))

# Create residue instance
pmb.create_residue(name="R1",
                   espresso_system=espresso_system)


print("\n=== Particle Instances DataFrame ===")
print(pmb.get_instances_df(pmb_type="particle"))
print("\n=== Residue Instances DataFrame ===")
print(pmb.get_instances_df(pmb_type="residue"))
print("\n=== Bond Instances DataFrame ===")
print(pmb.get_instances_df(pmb_type="bond"))

# Save database
pmb.save_database("demo_csv")

### Now create a new pyMBE instance with another set of reduced units
pmb2 = pyMBE.pymbe_library(seed=24)
pmb2.set_reduced_units(unit_length=0.6*pmb2.units.nanometer,)

pmb2.load_database("demo_csv")
print("\n=== Original Particle Templates DataFrame  ===")
print(pmb.get_templates_df(pmb_type="particle"))
print("\n=== Loaded Particle Templates DataFrame  ===")
print(pmb2.get_templates_df(pmb_type="particle"))

# Access some data in the database
tpl_particle_Z = pmb.db.get_template(name="Z", pmb_type="particle")
tpl_particle_Z_loaded = pmb2.db.get_template(name="Z", pmb_type="particle")

print("\n=== PintQuantity Usage Example After Loading Database ===")
original_sigma_Z = tpl_particle_Z.sigma.to_quantity(pmb.units)
loaded_sigma_Z   = tpl_particle_Z_loaded.sigma.to_quantity(pmb2.units)
print(f"Original sigma_Z: {original_sigma_Z.to('nanometer')} {original_sigma_Z.to('reduced_length')}")
print(f"Loaded sigma_Z: {loaded_sigma_Z.to('nanometer')} {loaded_sigma_Z.to('reduced_length')}")

# Delete the residue before proceding to the last example
pmb.delete_instances_in_system(espresso_system=espresso_system,
                              pmb_type="residue",
                              instance_id=0)
print("\n=== Particle Instances DataFrame After Deletion ===")
print(pmb.get_instances_df(pmb_type="particle"))
print("\n=== Residue Instances DataFrame After Deletion ===")
print(pmb.get_instances_df(pmb_type="residue"))
print("\n=== Bond Instances DataFrame After Deletion ===")
print(pmb.get_instances_df(pmb_type="bond"))

# Final example: let's create a hydrogel
## First define a molecule for the chains of the hydrogel
pmb.define_molecule(name="M1", 
                    residue_list=["R1"]*1)
diamond_lattice = DiamondLattice(4, 3.5 * units.reduced_length)
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
    chain_topology.append({'node_start':reverse_node_labels[node_label_s],
                            'node_end': reverse_node_labels[node_label_e],
                            'molecule_name':"M1"})

pmb.define_hydrogel("my_hydrogel", node_topology, chain_topology)

print("\n=== Molecule Templates DataFrame ===")
print(pmb.get_templates_df(pmb_type="molecule"))
print("\n=== Hydrogel Templates DataFrame ===")
print(pmb.get_templates_df(pmb_type="hydrogel"))

pmb.create_hydrogel(name="my_hydrogel",
                    espresso_system=espresso_system)
print("\n=== Particle Instances DataFrame After Hydrogel Creation ===")
print(pmb.get_instances_df(pmb_type="particle"))
print("\n=== Residue Instances DataFrame After Hydrogel Creation ===")
print(pmb.get_instances_df(pmb_type="residue"))
print("\n=== Bond Instances DataFrame After Hydrogel Creation ===")
print(pmb.get_instances_df(pmb_type="bond"))
print("\n=== Molecule Instances DataFrame After Hydrogel Creation ===")
print(pmb.get_instances_df(pmb_type="molecule"))
print("\n=== Hydrogel Instances DataFrame After Hydrogel Creation ===")
print(pmb.get_instances_df(pmb_type="hydrogel"))
pmb.save_database("demo_csv")
