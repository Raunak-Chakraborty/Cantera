# Program to observe dependence of Adiabatic Flame Temparature on number of bonds between 2 Carbon atoms in fuel :

# General Hydro Carbon Combustion Equation :

# CxHyOz + (1/phi)*(x + (y/4) - (z/2))(O2 + 3.76N2)        -------> 
# xCO2 + (y/2)H2O + (3.76/phi)*(x + (y/4) - (z/2))N2 + (x + (y/4) - (z/2))*((1/phi)-1)O2

# for z = 0 & phi = 1

# Ethane (single bond) : C2H6 + 3.5(O2 + 3.76N2) --> 2CO2 + 3H2O + 13.16N2       for x = 2, y = 6
# Ethene (double bond) : C2H4 + 3(O2 + 3.76N2) --> 2CO2 + 2H2O + 11.28N2         for x = 2, y = 4
# Ethyne (tripe bond) : C2H2 + 2.5(O2 + 3.76N2) --> 2CO2 + H2O + 9.4 N2          for x = 2, y = 2

##########################################################################################################################

# Solution using Python :

import matplotlib.pyplot as plt
import math

R = 8.314 # J/mol-K

# Function to evaluate Enthalpy using NASA Polynomials

def h(T, coeffs) :
        
    a1 = coeffs[0]
    a2 = coeffs[1]
    a3 = coeffs[2]
    a4 = coeffs[3]
    a5 = coeffs[4]
    a6 = coeffs[5]

    return (a1 + a2*T/2 + a3*pow(T, 2)/3 + a4*pow(T, 3)/4 + a5*pow(T, 4)/5 + a6/T)*R*T

# list of coefficients

O2_coeffs_l = [ 3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00 ]

N2_coeffs_l = [ 0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04, 0.03950372E+02 ]

C2H6_coeffs_l = [4.29142492E+00, -5.50154270E-03, 5.99438288E-05, -7.08466285E-08, 2.68685771E-11, -1.15222055E+04, 2.66682316E+00 ]

C2H4_coeffs_l = [3.95920148E+00, -7.57052247E-03, 5.70990292E-05, -6.91588753E-08, 2.69884373E-11, 5.08977593E+03, 4.09733096E+00 ]

C2H2_coeffs_l = [8.08681094E-01, 2.33615629E-02, -3.55171815E-05, 2.80152437E-08, -8.50072974E-12, 2.64289807E+04, 1.39397051E+01 ]

N2_coeffs_h = [ 0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02 ]

CO2_coeffs_h = [ 3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00 ]

H2O_coeffs_h =  [ 3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00 ]

# Switch between C2H6, C2H4 & C2H2

T_std = 298.15  # K

def f(T, n_bonds) :

    h_O2_r = h(T_std, O2_coeffs_l)
    h_N2_r = h(T_std, N2_coeffs_l)
    h_C2H6_r = h(T_std, C2H6_coeffs_l)
    h_C2H4_r = h(T_std, C2H4_coeffs_l)
    h_C2H2_r = h(T_std, C2H2_coeffs_l)
    h_N2_p = h(T, N2_coeffs_h)
    h_CO2_p = h(T, CO2_coeffs_h)
    h_H2O_p = h(T, H2O_coeffs_h)

    if n_bonds==1 :
        x = 2
        y = 6
        h_fuel = h_C2H6_r

    elif n_bonds==2 :
        x = 2
        y = 4
        h_fuel = h_C2H4_r

    elif n_bonds==3 :
        x = 2
        y = 2
        h_fuel = h_C2H2_r

    H_r = h_fuel + (x + (y/4)) * (h_O2_r + (3.76 * h_N2_r))
    H_p = x * h_CO2_p + (y/2) * h_H2O_p + (3.76 * (x + (y/4))) * h_N2_p

    return (H_p - H_r)

# Numerical Derivative of the above function using Forward Differencing

def f_dash(T, n_bonds) :
    return (f(T + 1e-06, n_bonds) - f(T, n_bonds))/1e-06


# Iterating using Newton-Raphson

tol = 1e-08       # Desired Tolerance
iter = 0
alpha = 0.5       # Relaxation Factor
Temp_python = []
n_bonds = [1, 2, 3]

for i in n_bonds :
    
    T_guess = 1500      
    while(abs(f(T_guess, i)) > tol) :
        T_guess = T_guess - alpha*(f(T_guess, i)/f_dash(T_guess, i))      
        iter = iter + 1
    Temp_python.append(T_guess)

print(Temp_python)

# Plotting Bar Graph for Python Data

plt.bar(['Ethane', 'Ethene', 'Ethyne'], Temp_python, color='red')
plt.ylabel('Adiabatic Flame Temparature [K]', fontsize = 12)
plt.title('Variation of AFT with Number of Bonds between Carbon\Atoms for 2 Carbon-Chain Fuels in Python', fontsize = 13, fontweight = 'bold')
plt.grid('on')
plt.show()

#####################################################################################################################

# Solution using Cantera :

import cantera as ct

gas = ct.Solution('gri30.xml')

Temp_cantera = []

for i in range(0, len(n_bonds)) :

    if n_bonds[i]==1 :
        x = 2
        y = 6
        fuel = 'C2H6'

    elif n_bonds[i]==2 :
        x = 2
        y = 4
        fuel = 'C2H4'

    elif n_bonds[i]==3 :
        x = 2
        y = 2
        fuel = 'C2H2'

    # Calculating Mole Fractions of Reactants

    n_total =  1 + (x + (y/4)) + (3.76 * (x + (y/4)))
    n_fuel = 1/n_total
    n_O2 = (x + (y/4))/n_total
    n_N2 = (3.76 * (x + (y/4)))/n_total

    # Defining the Reactant Mixture

    gas.TPX = 298.15, 101325, {fuel:n_fuel, 'O2':n_O2, 'N2':n_N2}

    gas.equilibrate('HP', 'auto')

    Temp_cantera.append(gas.T)

print(Temp_cantera)

# Plotting Bar Graph for Cantera Data

plt.bar(['Ethane', 'Ethene', 'Ethyne'], Temp_cantera)
plt.ylabel('Adiabatic Flame Temparature [K]', fontsize = 12)
plt.title('Variation of AFT with Number of Bonds between Carbon\Atoms for 2 Carbon-Chain Fuels in Cantera', fontsize = 13, fontweight = 'bold')
plt.grid('on')
plt.show()

#########################################################################################################################

# Tabulating Results from both Python & Cantera

Reactant = ['Ethane', 'Ethene', 'Ethyne']                                
print("---------------------------------------------------")
print("Reactant", 't', "Python", 'tt', "Cantera")
print("---------------------------------------------------")

for x in range(3) :
    print(Reactant[x], 't', Temp_python[x], 'tt', Temp_cantera[x])
print("---------------------------------------------------")









