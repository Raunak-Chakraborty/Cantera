# Calculating 1-D Flamelets

import cantera as ct
import numpy as np
import csv
import math
import os

class FlameExtinguished(Exception):
    pass

air = 'O2:0.21, N2:0.79'
fuel = 'CH4:1'
reaction_mechanism = 'gri30.xml' #'mech.cti'
gas = ct.Solution(reaction_mechanism)
f = ct.CounterflowDiffusionFlame(gas, width=1.)
f.transport_model = 'Multi'
f.P = 101325.
f.fuel_inlet.mdot = 0.01;
f.fuel_inlet.X = fuel;
f.fuel_inlet.T = 298.;
f.oxidizer_inlet.mdot = 0.025;
f.oxidizer_inlet.X = air;
f.oxidizer_inlet.T = 298.;
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
temperature_limit_extinction = 700

def interrupt_extinction(t):
    if np.max(f.T) < temperature_limit_extinction:
        raise FlameExtinguished('Flame extinguished')
    return 0
f.set_interrupt(interrupt_extinction)

# strain rate loop
strain_factor = 1.3 #increase strain rate by 30\% each step
exp_mdot_a = 1/2 #coeff from Fiala and Sattelmayer 2014
data_directory = '01_flamelet_results/'
file_name = 'initial_solution.xml'
n = 0

while (np.max(f.T) > temperature_limit_extinction or n==0):
    n += 1
    print('strain rate iteration', n)
    f = ct.CounterflowDiffusionFlame(gas, width=1)
    f.transport_model = 'Multi'
    f.P = 101325
    f.fuel_inlet.X = fuel;
    f.fuel_inlet.T = 298;
    f.oxidizer_inlet.X = air;
    f.oxidizer_inlet.T = 298;
    f.fuel_inlet.mdot = 0.01 * strain_factor ** (n*exp_mdot_a)
    f.oxidizer_inlet.mdot = 0.025 * strain_factor ** (n*exp_mdot_a)
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    
    try:
        f.solve(loglevel=0, auto=True)
        file_name = 'strain_loop_' + format(n, '02d') + '.xml'
        f.save(data_directory + file_name, name='diff1D', loglevel=0, description='Cantera version ' + ct.__version__ +', reaction mechanism ' + reaction_mechanism)
    
    except FlameExtinguished:
        print('Flame extinguished')
        break
    
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        break
    
    print("Solved on {:d} points.".format(len(f.T)))
    print("T max {:.0f} K.".format(max(f.T)))