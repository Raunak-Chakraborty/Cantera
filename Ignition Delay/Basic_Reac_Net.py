import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.xml')

gas.TPX = 1250, 5*101325, {'CH4' : 1, 'O2' : 2, 'N2' : 7.52}

# In the Ideal Gas Reactor Model, the Total Temparature T is used instead
# of Total Internal Energy U, as a State Variable 
reac = ct.IdealGasReactor(gas)

# Reactors primarily define Basic Governing Equations. The Time Integration
# is handled by a Reactor Network
simulation = ct.ReactorNet([reac])

states = ct.SolutionArray(gas)

time = 0

for n in range(1000) :

    time += 1e-3        # Final Time, not dt
    simulation.advance(time)
    states.append(reac.thermo.state, t = time)

print(states.P)
print(states.T)

plt.plot(states.T)
plt.show()


