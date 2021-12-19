# Program to study the effect of Pre-Heating on Efficiency & Fuel Savings

# We consider Stoichiometric Combustion of Methane at Constant Pressure inside a Furnace
# CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2  

"""
(m_air * h_air) + (m_fuel * h_fuel) = Q_out + Q_loss + (m_exhaust * h_exhaust)
=> Q_out = (m_air * h_air) + (m_fuel * h_fuel) - (m_exhaust * h_exhaust)
=> Q_out = (m_fuel * h_fuel) + (m_air * h_air) - (m_air + m_fuel) * h_exhaust
=> Q_out = m_fuel * {h_fuel + [(m_air/m_fuel) * h_air] - [(1 + (m_air/m_fuel)) * h_exhaust]}
=> Q_specific = h_fuel + [AF_ratio * h_air] - [(1 + AF_ratio) * h_exhaust] 

   Combustion Efficiency = (Q_specific / LHV)

   Fuel Savings = (Fuel used with pre-heating - Fuel used without pre-heating) / (Fuel used with preheating)
=> Fuel Savings = 1 - (Fuel used without pre-heating / Fuel used with preheating)
=> Fuel Savings = 1 - (Q_specific without pre-heating / Q_specific with pre-heating)
=> Fuel Savings = 1 - (Combustion Efficiency without pre-heating / Combustion Efficiency with pre-heating)

"""
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

gas = ct.Solution('gri30.cti')

############################################################################################################

# Defining Fuel
Fuel = ct.Quantity(gas)
Fuel.TPX = 298.15, ct.one_atm, 'CH4 : 1' 
Fuel.moles = 1

# Defining Prouct
Product = ct.Quantity(gas)
Product.TPX = 1700, ct.one_atm, {'CO2' : 1/10.52, 'H2O' : 2/10.52, 'N2' : 7.52/10.52}
Product.moles = (1 + 2 + 7.52)

# Defining Air
Air = ct.Quantity(gas)
Inlet_Temp = np.linspace(298.15, 600, 150)
LHV = 50e6 # J/kg         # Lower Heating Value is the Maximum Amount of Heat that can be released (Ideal)
Efficiency = []           # during the reaction. We shall use it as referrence to calculate Efficiency
Fuel_Saving = []

for i in Inlet_Temp :

    Air.TPX = i, ct.one_atm, {'O2' : 0.21, 'N2' : 0.79}
    Air.moles = 2 + (2 * 3.76)

    # <Speicies_Name>.enthalpy_mass returns the Specific Enthalpy of that Species
    AF_ratio = Air.mass/Fuel.mass
    h_Fuel = Fuel.enthalpy_mass
    h_Air = Air.enthalpy_mass
    h_Product = Product.enthalpy_mass

    # Heat Released in our case (Actual)
    Q_specific = h_Fuel + (AF_ratio * h_Air) - (1 + AF_ratio) * h_Product

    ita = (Q_specific/LHV) * 100        # Efficiency = Actual/Ideal
    Efficiency.append(ita)

    #Fuel Saved = 1 - (Efficiency without Pre-Heating/Efficiency with Pre-Heating) 
    
    fuel_saved = (1 - (Efficiency[0]/Efficiency[-1])) * 100
    Fuel_Saving.append(fuel_saved)

print("\nEfficiency without Pre-Heating (T_inlet = 298.15 k) : ", Efficiency[0], '%')
print("\nEfficiency with Pre-Heating (T_inlet = 600 k) : ", Efficiency[-1], '%')
print("\nMaximum Fuel Saving with Pre-Heating (T_inlet = 600 k) : ", Fuel_Saving[-1], '%', "\n")

############################################################################################################

# Plotting the Results

plt.plot(Inlet_Temp, Efficiency, '.', color = 'blue')
plt.grid()
plt.xlabel('Inlet Temparature of Air', fontsize = 12, fontweight = 'bold')
plt.ylabel('Combustion Efficiency [%]', fontsize = 12, fontweight = 'bold')
plt.title('Variation of Combustion Efficiency wrt. Inlet Temparature of Pre-Heated Air', fontsize = 13, fontweight = 'bold')
plt.show()

plt.plot(Inlet_Temp, Fuel_Saving, '.', color = 'red')
plt.grid()
plt.xlabel('Inlet Temparature of Air [K]', fontsize = 12, fontweight = 'bold')
plt.ylabel('Fuel Saving [%]', fontsize = 12, fontweight = 'bold')
plt.title('Variation of Fuel Savings wrt. Inlet Temparature of Pre-Heated Air', fontsize = 13, fontweight = 'bold')
plt.show()

