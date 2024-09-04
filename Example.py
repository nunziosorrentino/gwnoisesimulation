from gwskysim.sources.gwpointsource import GWPointSource
from gwskysim.sources.gwbinarysource import GWBinarySource
from gwskysim.sources.gwsignal import GWSignal

import numpy as np
import time

"""
This is an example of how gwskysim spurce simulation works, aimed to 
developers who want to contribut to the project. Here we present the 
simulation of two kinf of source without noise presence: the first is a 
close encounter source (CE) contained in this package and the second
one is a NSBH system which template comes from the pycbc banck. 
"""

# You are free to change teh source parameters
m1 = 1.14 # in solar masses
m2 = 10 # in solar masses
ra = 83.63 # in deg
dec = 22.01 # in deg
dist = 1.4*1e3*((27/14)*(10**16)) # 1.4 Mpc in solar masses
pol = 0 # in deg
incl = 0 # in deg

# Parameters for NSBH
spin1z, spin2z, coa_phase = (0.05, 0.05, 0.05) # no dim
# Parameters for CE (eccentricity and semi latus rectum)
e0, p0 = (0.99, 2*(m1 + m2)) # the second one is in solar mass unit
t_per = 0. # time of the periastron passage

# Setup CE
binary_source = GWBinarySource('CEncounter', ra=ra, dec=dec, polarization=pol, 
                               inclination=incl, distance=dist, e0=e0, p0=p0,
                               m1=m1, m2=m2, tp=t_per)
print(binary_source)
# Setup NSBH
nsbh_source = GWPointSource('NSBH', ra=ra, dec=dec, polarization=pol, 
                            inclination=incl, mass1=m1, mass2=m2, 
                            spin1z=spin1z, spin2z=spin2z, coa_phase=coa_phase)
print(nsbh_source)
    
# Configuration parameters used for the simulation (can be changed)
t0 = -50
duration = 100
fs = 4096
det = 'V1'
approx = 'SEOBNRv4'
l_freq=20

# gwskysim doesn't generate noise yes, so we set the backgound to zero
zero_signal = GWSignal('No noise', t0, duration, fs, detector = det)
h_v = zero_signal

# Generate and add the first signal from CE
print('Generating CE')
init_time1 = time.time()
h_v1 = binary_source.generate_signal(t0, duration, fs, detectors = [det])[0]
h_v += h_v1 
end_time1 = time.time()
print('Done!')

# Same with NSBH
print('Generating NSBH')
init_time2 = time.time()
h_v2 = nsbh_source.generate_signal(t0, duration, fs, detectors = [det], 
                   approximant=approx, lower_freq=l_freq)[0]
h_v += h_v2
end_time2 = time.time()
print('Done!')

# Check out run time for a comparison
print('CE "generate_signal" RUN TIME {:.4} s.'.format(end_time1 - init_time1))
print('NSBH "generate_signal" RUN TIME {:.4} s.'.format(end_time2 - init_time2))

# Plot and see if it satifies you!

import pylab

h_v1.plot()
#pylab.xlim((t0, t0+duration))
#pylab.ylim((-12*1e-22, 7*1e-22))

h_v2.plot()
pylab.xlim((t0, t0+duration))

h_v.plot()
#pylab.xlim((t0, t0+duration))
    
pylab.show()

# Thank you for the attention!
