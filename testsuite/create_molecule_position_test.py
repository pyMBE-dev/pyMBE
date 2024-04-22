import numpy as np 
import espressomd
from espressomd import interactions
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

print(f"***create_molecule with input position list unit test ***")
print(f"*** Unit test: Check that the positions of the central bead of the first residue in the generated molecules are equal to the input positions***")
# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
SEED = 100
solvent_permitivity = 78.3
N_molecules = 3
molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

pos_list = [[10,10,10], [20,20,20], [30,30,30]]
pmb.define_particle(name='central_mon',
                        acidity='inert',
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name='side_mon',
                        acidity='inert',
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))

pmb.define_residue(
    name = 'res1',
    central_bead = 'central_mon',
    side_chains = ['side_mon', 'side_mon']
    )

bond_type = 'harmonic'
generic_bond_lenght=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_lenght,
                 'k'      : generic_harmonic_constant,
                 }

pmb.define_default_bond(bond_type = bond_type, bond_parameters = harmonic_bond)

# Defines the peptine in the pyMBE data frame
molecule_name = 'generic_molecule'
pmb.define_molecule (name=molecule_name, residue_list = ['res1'])

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, q=1, sigma=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  q=-1, sigma=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
volume = N_molecules/(pmb.N_A*molecule_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_molecules/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
molecule = pmb.create_molecule(name=molecule_name,
                        number_of_molecules= N_molecules,
                        espresso_system=espresso_system,
                        use_default_bond=True,
                        list_of_first_residue_positions = pos_list)

###Running unit test here. Use np.testing.assert_almost_equal of the input position list and the central_bead_pos list under here.
central_bead_pos = []
for molecule_map in molecule:
            for molecule_id, info in molecule_map.items():
                central_bead_id = info['central_bead_id']
                side_chain_ids = info['side_chain_ids']
                central_bead_pos.append(espresso_system.part.by_id(central_bead_id).pos.tolist())


np.testing.assert_almost_equal(pos_list, central_bead_pos)

print(f"*** Unit test passed ***")
