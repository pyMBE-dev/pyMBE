import numpy as np
import pyMBE

pmb = pyMBE.pymbe_library(seed=42)
reduced_unit_set = pmb.get_reduced_units()
particle_parameters={"S1":{"name": "S1",
                        "sigma":0.355*pmb.units.nm, 
                        "epsilon":1*pmb.units('reduced_energy'),
                            "z":0},
                    "S2":{"name": "S2",
                        "sigma":0.355*pmb.units.nm, 
                        "epsilon":1*pmb.units('reduced_energy'),
                            "z":1},
                    "S3":{"name": "S3",
                        "sigma":0.355*pmb.units.nm, 
                        "epsilon":1*pmb.units('reduced_energy'),
                            "z":2}}

pmb.define_particles(parameters=particle_parameters)

generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_length = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_length,
                       'k'      : generic_harmonic_constant}
pmb.define_bond(bond_type = 'harmonic',
                        bond_parameters = HARMONIC_parameters, particle_pairs = [["S1", "S2"]])

pmb.delete_entries_in_df(entry_name="S1-S2")
assert pmb.df[pmb.df["name"]=="S1-S2"].empty
pmb.delete_entries_in_df(entry_name="S1")
assert pmb.df[pmb.df["name"]=="S1"].empty

residue_parameters={"R1":{"name": "R1",
                        "central_bead": "S2",
                        "side_chains": []},
                    "R2":{"name": "R2",
                        "central_bead": "S2",
                        "side_chains": ["S2","S3"]}}

for parameter_set in residue_parameters.values():
    pmb.define_residue(**parameter_set)

pmb.delete_entries_in_df(entry_name="R1")
assert pmb.df[pmb.df["name"]=="R1"].empty

molecule_parameters={"M1":{"name": "M1",
                    "residue_list": ["R2","R2","R2"]}}

for parameter_set in molecule_parameters.values():
    pmb.define_molecule(**parameter_set)

pmb.delete_entries_in_df(entry_name="M1")
assert pmb.df[pmb.df["name"]=="M1"].empty
