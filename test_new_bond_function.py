import old_pyMBE
import pyMBE
import espressomd
from espressomd import interactions
import pandas as pd

## Script to compare

scripts = [old_pyMBE,pyMBE]

## Create espresso_system

box_l = 5 # in reduced units
espresso_system = espressomd.System(box_l = [box_l]*3)


## Generating a polymer with two side chains, storing the outputs using both scripts 

output = []

for script in scripts:

    ## Creating an instance of pyMBE

    pmb = script.pymbe_library()
    
    ## Particle list

    part_list = ['A','B','C']

    ## Generic parameters
    
    generic_q                 = 0
    generic_diameter          = 0.4*pmb.units.nm
    generic_epsilon           = 1*pmb.units('reduced_energy')
    generic_bond_lenght       = 0.5*pmb.units.nm
    generic_harmonic_constant = 400*pmb.units('reduced_energy / nm**2')


    ## Defining particles
    
    for part in part_list:

        pmb.define_particle(name     = part, 
                            q        = generic_q, 
                            diameter = generic_diameter,
                            epsilon  = generic_epsilon)
    
    ## Defining residue

    pmb.define_residue(name         = 'R',
                       central_bead = 'A',
                       side_chains  = ['B', 'C'])
    
    ## Defining polymer

    pmb.define_molecule(name = 'P',
                        residue_list = ['R']*5)

    ## Defining bonds

    if script == scripts[0]:
        
        generic_bond = interactions.HarmonicBond(k = generic_harmonic_constant.to('reduced_energy / reduced_length**2').magnitude, r_0 = generic_bond_lenght.to('reduced_length').magnitude)
        
        for part in part_list:

            pmb.define_bond(bond_object = generic_bond,
                            bond_type = 'harmonic',
                            particle_name1 = 'A',
                            particle_name2 = part)

    elif script == scripts[1]:
        
        pmb.define_bond(bond_type      = 'harmonic',
                        k              = generic_harmonic_constant, 
                        r_0            = generic_bond_lenght,
                        particle_pairs = [['A', 'A'],
                                          ['A', 'B'],
                                          ['A', 'C']])
    else:
        0

    ## Adding bonds to espresso
    
    pmb.add_bonds_to_espresso(espresso_system = espresso_system)

    ## Creating polymer in espresso
    
    pmb.create_pmb_object(name = 'P', 
                          number_of_objects = 1,
                          espresso_system = espresso_system, 
                          position = [[box_l/2]*3]) 
    
    output.append(pmb.filter_df(pmb_type = 'bond').astype(str))
    
    ## Deleting the particles 

    pmb.destroy_pmb_object_in_system(name = 'P', 
                                     espresso_system = espresso_system)

if output[0].equals(output[1]):
    print("The DataFrames are exactly similar. TEST PASSED")
else:
    comparison= output[0].compare(output[1])
    print("The DataFrames are not exactly similar. TEST FAILED")
    print(f"The differences are: {comparison}\n")
