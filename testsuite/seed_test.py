import numpy as np 
import espressomd
import warnings
from espressomd import interactions
import pyMBE

espresso_system = espressomd.System(box_l = [100]*3)

def build_peptide_in_espresso(SEED):
    pmb = pyMBE.pymbe_library(SEED=SEED)

    # Simulation parameters
    pmb.set_reduced_units(unit_length=0.4*pmb.units.nm,
                          verbose=False)

    # Peptide parameters
    sequence = 'EEEEEEE'
    model = '2beadAA'  # Model with 2 beads per each aminoacid

    # Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
    path_to_interactions=pmb.get_resource("parameters/peptides/Lunkad2021.json")
    path_to_pka=pmb.get_resource("parameters/pka_sets/CRC1991.json")
    pmb.load_interaction_parameters(filename=path_to_interactions, verbose=False) 
    pmb.load_pka_set(path_to_pka, verbose=False)

    # Defines the peptide in the pyMBE data frame
    peptide_name = 'generic_peptide'
    pmb.define_peptide(name=peptide_name, sequence=sequence, model=model)

    # Bond parameters
    generic_bond_lenght=0.4 * pmb.units.nm
    generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

    HARMONIC_parameters = {'r_0'    : generic_bond_lenght,
                           'k'      : generic_harmonic_constant}

    pmb.define_default_bond(bond_type = 'harmonic',
                            bond_parameters = HARMONIC_parameters)

    # Add all bonds to espresso system
    pmb.add_bonds_to_espresso(espresso_system=espresso_system)

    # Create molecule in the espresso system
    pmb.create_pmb_object(name=peptide_name, number_of_objects=1, espresso_system=espresso_system, use_default_bond=True)

    # Extract positions of particles in the peptide
    positions = []
    molecule_id = pmb.df.loc[pmb.df['name']==peptide_name].molecule_id.values[0]
    particle_id_list = pmb.df.loc[pmb.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
    for pid in particle_id_list:
        positions.append(espresso_system.part.by_id(pid).pos)

    return np.asarray(positions)


print(f"*** Check that the using the same SEED results in the same initial particle positions***")
positions1 = build_peptide_in_espresso(42)
positions2 = build_peptide_in_espresso(42)

np.testing.assert_almost_equal(positions1, positions2)

print(f"*** Test passed ***")
