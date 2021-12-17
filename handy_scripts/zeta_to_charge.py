import numpy as np
import os
import sys
import inspect

# For loading sugar from parent folder

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import sugar

# It estimates the charge density of a colloidal object from the zeta potential.
# It is calculated using the Gouy−Chapman equation 
# See Eq. 3 from J. Phys. Chem. B 2017, 121, 3394−3402 DOI:  10.1021/acs.jpcb.6b08589
# Developed by:
# MsC. Albert Martinez (Royal College of Surgeons in Ireland)
# Dr. Pablo M. Blanco (Charles University) 

sg=sugar.sugar_library()

##### INPUT PARAMETERS #####

radius_NP = 11 * sg.units.nm
zeta_potential=-30 * sg.units.mV
conc=1e-3 * sg.units('mol/L')
er=78.4
e0=8.8541878176e-12 * sg.units('F/m')
T=298.15 * sg.units.K

surface_area=4*np.pi*(radius_NP)**2

root=(8*conc*sg.N_A*er*e0*sg.Kb*T)**(1./2.)
sinh=np.sinh(sg.e*zeta_potential/(2*sg.Kb*T))
charge_surface_density=root*sinh

print('the charge density of the surface is ', charge_surface_density.to('C/nm**2'),
        ' or ', (charge_surface_density/sg.e).to('nm**-2'), 'in elementary charge units')

surface_area=4*np.pi*(radius_NP)**2
N_charges=int(abs(surface_area*charge_surface_density/sg.e))

print('Punctual charges:', N_charges)
