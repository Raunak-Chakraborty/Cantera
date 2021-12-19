# Program to demonestrate the variation of Adiabatic Flame Temparature with a range
# of Inlet Temparatures (essentially embodying the Idea of Exhaust Gas Recycling)

# We consider Stoichiometric Combustion of Methane at Constant Pressure inside a Furnace
# CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2  

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

gas = ct.Solution('gri30.cti')

# Defining Fuel

Fuel = ct.Quantity(gas)
Fuel.TPX = 298.15, ct.one_atm, 'CH4 : 1'
Fuel.moles = 1

# Defining Air, the Inlet Temparature of which, is to be varied

Air = ct.Quantity(gas)
Inlet_Temp = np.linspace(298, 600, 200)     # Desired reange of Inlet Temparatures
AFT = []                                    # Empty Array to store each Iteration Value

for i in Inlet_Temp :

    Air.TPX = i, ct.one_atm, {'O2' : 0.21, 'N2' : 0.79}     # ct.one_atm sets the P as 1 atm = 101325 pa
    Air.moles = 2  + (2 * 3.67)

    Mixture = Air + Fuel        

    Mixture.equilibrate('HP', 'auto')

    AFT.append(Mixture.T)       # Inserting each value in the AFT array 

print("\nAFT when Inlet Air Temparature  = 298.15 K : ", AFT[0], 'K')
print("\nAFT when Inlet Air Temparature  = 600 K : ", AFT[-1], 'K', "\n")

# Plotting the Results

plt.plot(Inlet_Temp, AFT, '.', color = 'black')
plt.grid()
plt.xlabel('Inlet Temparature of Air [K]', fontsize = 12, fontweight = 'bold')
plt.ylabel('Adiabatic Flame Temparature [K]', fontsize = 12, fontweight = 'bold')
plt.title('Variation of Adiabatic Flame Temparature wrt. Inlet Temparature of Pre-Heated Air', fontsize = 13, fontweight = 'bold')
plt.show()