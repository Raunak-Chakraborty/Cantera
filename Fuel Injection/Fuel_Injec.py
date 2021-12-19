# The following Program simulates Fuel Injection of n-dodecane into a mixture
# of vitiated air (air in which Oxygen-content has been reduced) :

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

##################################################################################

# We use a reduced n-dodecane mechanism with PAH pathways
gas = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_IG')
gas.case_sensitive_species_names = True

# Setting up Inlet Conditions
gas.TPX = 300, 20*ct.one_atm, 'c12h26 : 1.0'
inlet = ct.Reservoir(gas)

# Setting up a reactor & specifying it's initial conditions
# to be the product of a lean reaction
gas.TP = 1000, 20*ct.one_atm
gas.set_equivalence_ratio(0.30, 'c12h26 : 1.0', 'n2:3.76, o2:1.0')
gas.equilibrate('TP', 'auto')

r = ct.IdealGasReactor(gas)
r.volume = 0.001    # in m^3

# Defining a function for the inlet mass flow rate of the fuel,
# supplied as a Gaussian Pulse
def m_dot_fuel(t) :

    total = 3.0e-3    # Mass of Fuel [kg]
    width = 0.5       # Width of Pulse [s]
    t0 = 2.0          # Time of Pulse Peak [s]
    amplitude = total / (width * np.sqrt(2 * np.pi))
    
    return amplitude * np.exp(-(t-t0)**2 / (2*width**2))

# Inlet Mass Flow Controller
inlet_mfc = ct.MassFlowController(inlet, r, mdot = m_dot_fuel)

sim = ct.ReactorNet([r])

states = ct.SolutionArray(gas, extra = ['t'])

# Time Integration
for t in np.linspace(0, 10, 1000) :
    sim.advance(t)
    states.append(r.thermo.state, t=t)

##################################################################################

species_aliases = { 
    'o2' : 'O$_2$',
    'h2o' : '$H_2$O',
    'o2' : 'CO',
    'co2' : 'CO$_2$',
    'h2' : 'H$_2$',
    'ch4' : 'CH$_4$'}

for name, alias in species_aliases.items() :
    gas.add_species_alias(name, alias)

# PAH (Poly Aromatic HydroCarbons) are considered to be precursors to soot formation
PAH_aliases = { 
    'A1c2h' : 'phenylacetylene [C$_8$H$_6$]',
    'A1c2h3' : 'styrene [C$_8$H$_8$]',
    'A1' : 'benzene [C$_6$H$_6$]',
    'A2' : 'napthalene [C$_10$H$_8$]',
    'A2r5' : 'acenapthylene [C$_12$H$_8$]',
    'A3' : 'phenanthrene [C$_14$H$_10$]',
    'A4' : 'pyrene [C$_16$H$_10$]'}

for name, alias in PAH_aliases.items() :
    gas.add_species_alias(name, alias)

# Plotting Mole Fractions of different Species with Time
f, ax = plt.subplots(1, 2)
plt.suptitle('Variation of Mole Fraction of different Species with Time', fontweight='bold')

for s in species_aliases.values() :
    ax[0].plot(states.t, states(s).X, label=s)

for s in PAH_aliases.values() :
    ax[1].plot(states.t, states(s).X, label=s)

for a in ax :

    a.legend(loc='best')
    a.set_xlabel('Time [s]', fontweight='bold')
    a.set_ylabel('Mole Fraction', fontweight='bold')
    a.set_xlim([0, 10])

plt.show()



