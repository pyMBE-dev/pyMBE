#
# Copyright (C) 2024 pyMBE-dev team
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

# Load python modules
import os 
import espressomd
import numpy as np
import pandas as pd
from scipy import interpolate
import argparse
from tqdm import tqdm

# Import pyMBE
import pyMBE
from lib import analysis

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)



# Command line arguments
parser = argparse.ArgumentParser(description='Script that runs a simulation of a salt solution in the grand-canonical ensemble using pyMBE and ESPResSo.')
parser.add_argument("--c_salt_res", 
                    type=float, 
                    help="Concentration of salt in the reservoir in mol/l.")
parser.add_argument('--mode',
                    type=str,
                    default= "ideal",
                    help='Set if an ideal or interacting system is simulated.')
parser.add_argument('--output',
                    type=str,
                    required= False,
                    help='output directory')
parser.add_argument('--no_verbose', 
                    action='store_false', 
                    help="Switch to deactivate verbose")
args = parser.parse_args()

# Inputs
inputs={"csalt": args.c_salt_res,
        "mode": args.mode}

mode=args.mode
valid_modes=["ideal", "interacting"]
verbose=args.no_verbose

if args.mode not in valid_modes:
    raise ValueError(f"Mode {mode} is not currently supported, valid modes are {valid_modes}")

#Import functions from handy_functions script 
from lib.handy_functions import minimize_espresso_system_energy
from lib.handy_functions import setup_electrostatic_interactions
from lib.handy_functions import setup_langevin_dynamics



# Units and general parameters
pmb.set_reduced_units(unit_length=0.355*pmb.units.nm)
solvent_permittivity = 78.9

# Integration parameters
DT = 0.01
LANGEVIN_SEED = 42

# Define salt
cation_name = 'Na'
anion_name = 'Cl'
pmb.define_particle(name=cation_name, q=1, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, sigma=0.355*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
c_salt_res = args.c_salt_res * pmb.units.mol/ pmb.units.L
N_SALT_ION_PAIRS = 200
volume = N_SALT_ION_PAIRS/(pmb.N_A*c_salt_res)
L = volume ** (1./3.) # Side of the simulation box

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)
if verbose:
    print("Created espresso object")

# Add salt
c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=0.5*c_salt_res)
if verbose:
    print("Added salt")

# Set up reactions
if args.mode == "interacting":
    path_to_ex_pot=pmb.get_resource("testsuite/data/src/")
    ionic_strength, excess_chemical_potential_monovalent_pairs_in_bulk_data, bjerrums, excess_chemical_potential_monovalent_pairs_in_bulk_data_error =np.loadtxt(f"{path_to_ex_pot}/excess_chemical_potential.dat", unpack=True)
    excess_chemical_potential_monovalent_pair_interpolated = interpolate.interp1d(ionic_strength, excess_chemical_potential_monovalent_pairs_in_bulk_data)
    activity_coefficient_monovalent_pair = lambda x: np.exp(excess_chemical_potential_monovalent_pair_interpolated(x.to('1/(reduced_length**3 * N_A)').magnitude))
    RE = pmb.setup_gcmc(c_salt_res=c_salt_res, salt_cation_name=cation_name, salt_anion_name=anion_name, activity_coefficient=activity_coefficient_monovalent_pair)
elif args.mode == "ideal":
    RE = pmb.setup_gcmc(c_salt_res=c_salt_res, salt_cation_name=cation_name, salt_anion_name=anion_name)
if verbose:
    print("Set up GCMC...")

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map(type_list = types)
print(type_map)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type(type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

espresso_system.time_step = DT

if args.mode == "interacting":
    #Set up the short-range interactions
    pmb.setup_lj_interactions(espresso_system=espresso_system)

# Minimzation
minimize_espresso_system_energy(espresso_system=espresso_system, Nsteps=1e4, max_displacement=0.01, skin=0.4)
setup_langevin_dynamics(espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = LANGEVIN_SEED,
                                    time_step=DT,
                                    tune_skin=False)

if verbose:
    print("Running warmup without electrostatics")
for i in tqdm(range(100),disable=not verbose):
    espresso_system.integrator.run(steps=100)
    RE.reaction(reaction_steps=100)

if args.mode == "interacting":
    setup_electrostatic_interactions(units=pmb.units,
                                    espresso_system=espresso_system,
                                    kT=pmb.kT,
                                    solvent_permittivity=solvent_permittivity)

espresso_system.thermostat.turn_off()
minimize_espresso_system_energy(espresso_system=espresso_system, Nsteps=1e4, max_displacement=0.01, skin=0.4)
setup_langevin_dynamics(espresso_system=espresso_system, 
                                    kT = pmb.kT, 
                                    SEED = LANGEVIN_SEED,
                                    time_step=DT,
                                    tune_skin=False)

if verbose:
    print("Running warmup with electrostatics")

N_warmup_loops = 100
for i in tqdm(range(N_warmup_loops),disable=not verbose):
    espresso_system.integrator.run(steps=100)
    RE.reaction(reaction_steps=100)

# Main loop
print("Started production run.")

labels_obs=["time", "c_salt"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

N_production_loops = 100
for i in tqdm(range(N_production_loops),disable=not verbose):
    espresso_system.integrator.run(steps=100)
    RE.reaction(reaction_steps=100)

    # Measure time
    time_series["time"].append(espresso_system.time)

    # Measure degree of ionization
    number_of_ion_pairs = espresso_system.number_of_particles(type_map[cation_name])
    time_series["c_salt"].append((number_of_ion_pairs/(volume * pmb.N_A)).magnitude)

data_path = args.output
if data_path is None:
    data_path=pmb.get_resource(path="samples/Beyer2024")+"/time_series/gcmc"

if not os.path.exists(data_path):
    os.makedirs(data_path)

time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(f"{data_path}/{filename}_time_series.csv", index=False)
particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df(filename=f'{data_path}/df.csv')
