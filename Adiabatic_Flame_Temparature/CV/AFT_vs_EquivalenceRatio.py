# Program to calculate Adiabatic Flame Temparature of Methane at different Equivalence Ratios :

# CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2       in case of Stoichiometric Combustion

# CH4 + (2/phi)(O2 + 3.76N2) -> CO2 + 2H2O + (7.52/phi)N2 + ((2/phi)-2)O2     in case of Lean Combustion

# CH4 + (2/phi)(O2 + 3.76N2) -> ((4/phi)-3)CO2 + 2H2O + (7.52/phi)N2 + (4-(4/phi))CO     in case of Rich Combustion


import matplotlib.pyplot as plt
import math
import numpy as np

R = 8.314 # J/mol-K

# Solution using Python : 

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

CH4_coeffs_l = [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04, -4.64130376E+00 ] 

O2_coeffs_l = [ 3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00 ]

O2_coeffs_h = [ 3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00 ]

N2_coeffs_l = [ 0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04, 0.03950372E+02 ]

N2_coeffs_h = [ 0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02 ]

CO2_coeffs_h = [ 3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00 ]

H2O_coeffs_h =  [ 3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00 ]

CO_coeffs_h = [ 2.71518561E+00, 2.06252743E-03, -9.98825771E-07, 2.30053008E-10, -2.03647716E-14, -1.41518724E+04, 7.81868772E+00 ]


# In a Constant Volume Process, the Internal Energy remains constant, which is embodied by the function f
# U_r = U_p => H_r - (PV)_R = H_p - (PV)_p => H_r - (nRT)_R = H_p - (nRT)_p => (H_p - H_r) - ((nRT)_p - (nRT)_r) = 0

def f(T, phi) :

    # Calculating the Enthalpy of the Products
    h_CO2_p = h(T, CO2_coeffs_h)
    h_N2_p = h(T, N2_coeffs_h)                      
    h_H2O_p = h(T, H2O_coeffs_h)
    h_CO_p = h(T, CO_coeffs_h)
    h_O2_p = h(T, O2_coeffs_h)

    # Switch between Lean & Rich Mixtures

    if phi<1 :                                                                 # in case of a Lean Mixture
        H_p = h_CO2_p + 2*h_H2O_p + (7.52/phi)*h_N2_p + ((2/phi)-2)*h_O2_p     # Enthalpy_Products_Lean 
        N_p = (1 + 2 + (7.52/phi) + (2/phi)-2)*R*T                             # No. of Moles_Products_Lean     

    elif phi>=1 :                                                                       # in case of a Rich Mixture
        H_p = ((4/phi)-3)*h_CO2_p + 2*h_H2O_p + (7.52/phi)*h_N2_p + (4-(4/phi))*h_CO_p  # Enthalpy_Products_Rich
        N_p = ((4/phi) - 3 + 2 + (7.52/phi) + 4 - (4/phi))*R*T                          # No. of Moles_Products_Rich


    T_std = 298.15 # K         # Reactants exist at STP


    # Calculating the Enthalpy of Reactants
    h_CH4_r = h(T_std, CH4_coeffs_l)
    h_O2_r = h(T_std, O2_coeffs_l)
    h_N2_r = h(T_std, N2_coeffs_l)

    H_r = h_CH4_r + (2/phi)*h_O2_r + (7.52/phi)*h_N2_r            # Total Enthalpy of Reactants
    N_r = (1 + (2/phi) + (7.52/phi))*R*T_std                      # Total Number of Moles of Reactants
    
    return (H_p - H_r) - (N_p - N_r)                              # This is the function that must be minimized


# Numerical Derivative of the above function using Forward Differencing

def f_dash(T, phi) :
    return (f(T + 1e-06, phi) - f(T, phi))/1e-06


# Iterating using Newton-Raphson

tol = 1e-08               # Desired Tolerance
iter = 0
alpha = 0.5               # Relaxation Factor
AFT_Python = []           # This initially-empty array will store all T_guess values

phi = np.linspace(0.5, 2, 16)      # Creating an array of length 16 with equally-spaced elements, from 0.5 to 2

for i in phi :
    T_guess = 1500      

    while(abs(f(T_guess, i)) > tol) :
        T_guess = T_guess - alpha*(f(T_guess, i)/f_dash(T_guess, i))
        iter = iter + 1

    AFT_Python.append(T_guess)             # Insering current T_guess at the end of Temp array

print(max(AFT_Python))
print(AFT_Python)

#############################################################################################################################

# Solution using Cantera :

import cantera as ct

gas = ct.Solution ('gri30.xml')

gas.name = 'Methane + Oxygen + Nitrogen'

phi = np.linspace(0.5, 2, 16)

AFT_Cantera = []

for i in range (len(phi)) :

    gas.TP = 298.15, 101325
    gas.set_equivalence_ratio (phi[i], "CH4:1", "O2:2, N2:7.52" )

    gas.equilibrate('UV', 'auto')

    AFT_Cantera.append(gas.T)

print(AFT_Cantera)
print(max(AFT_Cantera))

#########################################################################################################

# Plotting the Results

plt.plot(phi, AFT_Python, '-o', color = 'red', label='Python')
plt.legend()
plt.plot(phi, AFT_Cantera, '-o', color='blue', label='Cantera')
plt.legend()

plt.xlabel('Equivalence Ratio($phi$)', fontsize=12)
plt.ylabel('Adiabatic Flame Temparature [K]', fontsize=12)

plt.grid('on')
plt.title('AFT vs Equivalence Ratio Plot for Methane-Air Combustion', fontsize=13, fontweight='bold')

############################################################################################################

# Tabulating Results from both Python & Cantera

print("---------------------------------------------------")
print("phi", '      ', "Python", '      ', "Cantera")
print("---------------------------------------------------")

for x in range(len(phi)) :
    print(phi[x], '    ', AFT_Python[x], '    ', AFT_Cantera[x])
print("---------------------------------------------------")

plt.show()