# This Program compares the 1-D Reaction Scenarios for Methane
# & Hydrogen Combustion using the Free Flame Object in Cantera

# CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2
# 2H2O + O2 -> 2H2O

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

###################################################################

# Instantating Object of Solutions Class
gas_CH4 = ct.Solution('gri30.cti', 'gri30_mix')
gas_H2 = ct.Solution('H2_mech.cti')

# Specifying Initial Conditions
gas_CH4.TPX = 300, ct.one_atm, {'CH4':1, 'O2':2, 'N2':7.52}
gas_H2.TPX = 300, ct.one_atm, {'H2':2, 'O2':1}

# Creating the Free Flame Object and passing the respective  
# gas objects & width of the domain [m] as parameters 
flame_CH4 = ct.FreeFlame(gas_CH4, width=0.03)
flame_H2 = ct.FreeFlame(gas_H2, width=0.03)

# Cantera will keep adding Grid Points till the values of these
# 3 parameters drop below the stated upper limits
flame_CH4.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
flame_H2.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)

# The following command actually solves the problem by repeating the
# following steps : Finite Difference Discretization, Newton-Raphson
# Method, Grid Refinement, repeat.
flame_CH4.solve(loglevel=1, refine_grid=True, auto=True)
flame_H2.solve(loglevel=1, refine_grid=True, auto=True)

# The observer sits on the Flame frame of referrence, & thus 
# Flame Speed = Speed of Unburnt Gases entering the Flame
flame_speed_CH4 = flame_CH4.velocity
flame_speed_H2 = flame_H2.velocity

# This calculates the Volumetric Heat Released for both Reactions
# standard_enthalpies_RT returns an array containing [Ideal Gas Enthalpies / (R.T)]
# net_production_rates method returns the net rate of production of each species (kmol/m^3/s)
Q_CH4 = abs(sum(flame_CH4.standard_enthalpies_RT * flame_CH4.net_production_rates, 0) * ct.gas_constant * flame_CH4.T)
Q_H2 = abs(sum(flame_H2.standard_enthalpies_RT * flame_H2.net_production_rates, 0) * ct.gas_constant * flame_H2.T)

# Extracting Concentration of CO2 for both reactions
CO2_conc_CH4 = flame_CH4.concentrations[gas_CH4.species_index('CO2'), :]
CO2_conc_H2 = flame_H2.concentrations[gas_H2.species_index('CO2'), :]

###################################################################

# Plotting Temperature VS Axial Distance
plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(flame_CH4.grid, flame_CH4.T)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Temperature [K]', fontweight = 'bold')
plt.title('Variation Of Temperature at each Grid Point for CH4 Combustion', fontweight = 'bold')

plt.subplot(2, 1, 2)
plt.plot(flame_H2.grid, flame_H2.T)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Temperature [K]', fontweight = 'bold')
plt.title('Variation Of Temperature at each Grid Point for H2 Combustion', fontweight = 'bold')

# Plotting Flame Speed VS Axial Distance
plt.figure(2)

plt.subplot(2, 1, 1)
plt.plot(flame_CH4.grid, flame_speed_CH4)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Flame Speed [m/s]', fontweight = 'bold')
plt.title('Variation Of Flame Speed at each Grid Point for CH4 Combustion', fontweight = 'bold')

plt.subplot(2, 1, 2)
plt.plot(flame_H2.grid, flame_speed_H2)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Flame Speed [m/s]', fontweight = 'bold')
plt.title('Variation Of Flame Speed at each Grid Point for H2 Combustion', fontweight = 'bold')

# Plotting Volumetric Heat Release VS Axial Distance
plt.figure(3)

plt.subplot(2, 1, 1)
plt.plot(flame_CH4.grid, Q_CH4)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Volumetric Heat Released [J/m^3]', fontweight = 'bold')
plt.title('Variation Of Volumetric Heat Release at each Grid Point for CH4 Combustion', fontweight = 'bold')

plt.subplot(2, 1, 2)
plt.plot(flame_H2.grid, Q_H2)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Volumetric Heat Released [J/m^3]', fontweight = 'bold')
plt.title('Variation Of Volumetric Heat Release at each Grid Point for H2 Combustion', fontweight = 'bold')

# Plotting Concentration of CO2 VS Axial Distance
plt.figure(4)

plt.subplot(2, 1, 1)
plt.plot(flame_CH4.grid, CO2_conc_CH4)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Concentration of CO2 [Kmol/m^3]', fontweight = 'bold')
plt.title('Variation Of Concentration of CO2 at each Grid Point for CH4 Combustion', fontweight = 'bold')

plt.subplot(2, 1, 2)
plt.plot(flame_H2.grid, CO2_conc_H2)
plt.tight_layout()
plt.grid('on')
plt.xlabel('grid [m]', fontweight = 'bold')
plt.ylabel('Concentration of CO2 [Kmol/m^3]', fontweight = 'bold')
plt.title('Variation Of Concentration of CO2 at each Grid Point for H2 Combustion', fontweight = 'bold')

plt.show()