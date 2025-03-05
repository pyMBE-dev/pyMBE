import espressomd
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import interpolate
import argparse
from lib.lattice import DiamondLattice
import pickle
import pyMBE
from lib import analysis

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# Load some functions from the handy_scripts library for convinience
from lib.handy_functions import setup_electrostatic_interactions
from lib.handy_functions import minimize_espresso_system_energy
from lib.handy_functions import setup_langevin_dynamics
#######################################################
# Setting parameters for the simulation
#######################################################

##### Read command-line arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument("--MPC",
                    type=int,
                    help="Number of monomers per chain ")
parser.add_argument("--L_target",
                    type=float,
                    help="Target box length to achieve starting from L_max")
parser.add_argument("--c_salt_res",
                    type=float,
                    help="Concentration of NaCl in the reservoir in mol/l.")
parser.add_argument("--pH_res",
                    type=float,
                    help="pH-value in the reservoir.")
parser.add_argument("--pKa_value",
                    type=float,
                    help="pKa-value of the hydrogel monomers")
parser.add_argument('--mode',
                    type=str,
                    default= "short-run",
                    help='sets for how long the simulation runs, valid modes are {valid_modes}')
parser.add_argument('--output',
                    type=str,
                    required= False,
                    help='output directory')
parser.add_argument('--no_verbose',
                    action='store_false',
                    help="Switch to deactivate verbose",
                    default=True)

args = parser.parse_args()
mode=args.mode
c_salt_res = args.c_salt_res * pmb.units.mol/ pmb.units.L
solvent_permittivity = 78.9

NodeType = "N"
pmb.define_particle(name=NodeType, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))

BeadType = "C"
pmb.define_particle(name=BeadType, 
                    acidity="acidic", 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'), 
                    pka=args.pKa_value)

pmb.define_residue(name='Res', 
                   central_bead=BeadType, 
                   side_chains=[])

bond_type = 'FENE'
bond_length = 0.966 * pmb.units("reduced_length")
fene_spring_constant = 30 * pmb.units('reduced_energy / reduced_length**2')
fene_r_max = 1.5 * pmb.units('reduced_length')

fene_bond = {'k'      : fene_spring_constant,
             'd_r_max': fene_r_max,
            }

pmb.define_bond(bond_type = bond_type, 
                bond_parameters = fene_bond, 
                particle_pairs = [[BeadType,NodeType],[BeadType, BeadType]])

# Parameters of the small ions
proton_name = 'Hplus'
hydroxide_name = 'OHminus'
sodium_name = 'Na'
chloride_name = 'Cl'

pmb.define_particle(name=proton_name, 
                    z=1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=hydroxide_name, 
                    z=-1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=sodium_name, 
                    z=1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=chloride_name, 
                    z=-1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))

diamond_lattice = DiamondLattice(args.MPC, bond_length)
espresso_system = espressomd.System(box_l = [diamond_lattice.BOXL]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

# Setting up node topology
indices = diamond_lattice.indices
node_topology = []

for index in range(len(indices)):
    node_topology.append({"particle_name": NodeType,
                          "lattice_index": indices[index]})

# Setting up chain topology
node_labels = lattice_builder.node_labels
chain_labels = lattice_builder.chain_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
chain_topology = []
residue_list = ["Res"]*args.MPC

for chain_data in chain_labels.items():
    chain_label = chain_data[1]
    node_label_pair = chain_data[0]
    node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
    chain_topology.append({'node_start':reverse_node_labels[node_label_s],
                              'node_end': reverse_node_labels[node_label_e],
                              'residue_list':residue_list})

pmb.define_hydrogel("my_hydrogel", node_topology, chain_topology)
hydrogel_info = pmb.create_hydrogel("my_hydrogel", espresso_system)

pmb.create_counterions(object_name="my_hydrogel", 
                       cation_name=proton_name, 
                       anion_name=hydroxide_name, 
                       espresso_system=espresso_system)

print("*** Finished adding counterions ***")

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system, 
                                          cation_name=sodium_name, 
                                          anion_name=chloride_name, 
                                          c_salt=c_salt_res)

print(f"Salt concentration {c_salt_calculated} added")

print("*** Setting LJ interactions ***")

DT = 0.01  # Timestep 
espresso_system.time_step = DT
pmb.setup_lj_interactions(espresso_system=espresso_system)

setup_langevin_dynamics(espresso_system=espresso_system,
                        kT = pmb.kT,
                        SEED = 87643,
                        time_step=DT,
                        tune_skin=False)

print("*** Force capping.. ***")
lj_cap = 50
espresso_system.force_cap = lj_cap
i=0
while i < 100:
    espresso_system.integrator.run(steps=500+args.MPC*2)
    i += 1
    lj_cap = lj_cap + 10
    espresso_system.force_cap = lj_cap

print("min dist",espresso_system.analysis.min_dist())

espresso_system.force_cap = 0
L_max = diamond_lattice.BOXL
L_target = args.L_target * L_max
steps_size = 0.5
steps_needed = int(abs((L_max - L_target)/steps_size))

print("*** BOX Dimension changing... ***")

for j in tqdm(np.arange(0,steps_needed)):
    espresso_system.change_volume_and_rescale_particles(d_new = espresso_system.box_l[0]-steps_size, dir = "xyz")
    espresso_system.integrator.run(1000)  #1000 here
espresso_system.change_volume_and_rescale_particles(d_new = L_target, dir = "xyz")
espresso_system.thermostat.turn_off()
minimize_espresso_system_energy(espresso_system=espresso_system, Nsteps=1e4, max_displacement=0.01, skin=0.4)

print("*** Setting up the reactions... ***")

# Set up the reactions
path_to_ex_pot=pmb.get_resource("parameters/salt/")
ionic_strength, excess_chemical_potential_monovalent_pairs_in_bulk_data, bjerrums, excess_chemical_potential_monovalent_pairs_in_bulk_data_error =np.loadtxt(f"{path_to_ex_pot}/monovalent_salt_excess_chemical_potential.dat", unpack=True)
excess_chemical_potential_monovalent_pair_interpolated = interpolate.interp1d(ionic_strength, excess_chemical_potential_monovalent_pairs_in_bulk_data)
activity_coefficient_monovalent_pair = lambda x: np.exp(excess_chemical_potential_monovalent_pair_interpolated(x.to('1/(reduced_length**3 * N_A)').magnitude))
pka_set = {BeadType: {"pka_value": args.pKa_value, "acidity": "acidic"}}

grxmc, labels, ionic_strength_res = pmb.setup_grxmc_reactions(pH_res=args.pH_res, c_salt_res=c_salt_res, proton_name=proton_name, hydroxide_name=hydroxide_name, salt_cation_name=sodium_name, salt_anion_name=chloride_name, activity_coefficient=activity_coefficient_monovalent_pair, pka_set=pka_set)

# Setup espresso to track the ionization of the acid groups
type_map = pmb.get_type_map()
types = list(type_map.values())
espresso_system.setup_type_map(type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
grxmc.set_non_interacting_type (type=non_interacting_type)

for i in tqdm(range(100)):
    espresso_system.integrator.run(steps=1000)
    grxmc.reaction(1000)

setup_electrostatic_interactions(units=pmb.units,
                                espresso_system=espresso_system,
                                kT=pmb.kT,
                                solvent_permittivity=solvent_permittivity)
espresso_system.thermostat.turn_off()
minimize_espresso_system_energy(espresso_system=espresso_system, Nsteps=1e4, max_displacement=0.01, skin=0.4)

setup_langevin_dynamics(espresso_system=espresso_system,
                                    kT = pmb.kT,
                                    SEED = 7653,
                                    time_step=DT,
                                    tune_skin=False)

print("*** before electrostatics ***")

if mode == "long-run":
    N_warmup_loops = 1000
else:
    N_warmup_loops = 100

print("*** Running warmup with electrostatics ***")
for i in tqdm(range(N_warmup_loops)):
    espresso_system.integrator.run(steps=1000)
    grxmc.reaction(100)


# Main loop
print("Started production run.")

labels_obs=["time", "alpha", "pressure"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

if mode == "long-run":
    N_production_loops = 5000
else:
    N_production_loops = 100

for i in tqdm(range(N_production_loops)):
    espresso_system.integrator.run(steps=1000)
    grxmc.reaction(200)

    # Measure time
    time_series["time"].append(espresso_system.time)

    # Measure degree of ionization
    time_series["alpha"].append(len(espresso_system.part.select(type=2,q=-1))/(16*args.MPC))
    time_series["pressure"].append(espresso_system.analysis.pressure()["total"])

with open("time_series.pkl", "wb") as f:
    pickle.dump(time_series, f)

inputs={"csalt": c_salt_res,
        "L_target": args.L_target,
        "pH": args.pH_res,
        "pKa": args.pKa_value}

data_path = args.output
if data_path is None:
    data_path=pmb.get_resource(path="samples/Landsgessel22a")+"/time_series/"

Path(data_path).mkdir(parents=True,
                       exist_ok=True)

time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(f"{data_path}/{filename}_time_series.csv", index=False)

