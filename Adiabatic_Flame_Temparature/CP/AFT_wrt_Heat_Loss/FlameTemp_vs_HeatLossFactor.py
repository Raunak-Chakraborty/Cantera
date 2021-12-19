# Program to observe dependence of Flame Temparature on Heat Loss Factor :

# CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2       in case of Stoichiometric Combustion

###################################################################################################################################

import matplotlib.pyplot as plt
import math
import numpy as np

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

CH4_coeffs_l = [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04, -4.64130376E+00 ] 

O2_coeffs_l = [ 3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00 ]

N2_coeffs_l = [ 0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04, 0.03950372E+02 ]

N2_coeffs_h = [ 0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02 ]

CO2_coeffs_h = [ 3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00 ]

H2O_coeffs_h =  [ 3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00 ]

# In a Constant Pressure Process, the Enthalpy remains constant, which is embodied by the function f
# H_r + (Heat Loss Factor * Lower Heating Value) = H_p => H_p - H_r + (H_loss * LHV) = 0

def f(T, H_loss) :

    # Calculating the Enthalpy of the Products
    h_CO2_p = h(T, CO2_coeffs_h)
    h_N2_p = h(T, N2_coeffs_h)                      
    h_H2O_p = h(T, H2O_coeffs_h)
   
    H_p = h_CO2_p + 2*h_H2O_p + 7.52*h_N2_p

    T_std = 298.15 # K         # Reactants exist at STP

    # Calculating the Enthalpy of Reactants
    h_CH4_r = h(T_std, CH4_coeffs_l)
    h_O2_r = h(T_std, O2_coeffs_l)
    h_N2_r = h(T_std, N2_coeffs_l)

    H_r = h_CH4_r + 2*h_O2_r + 7.52*h_N2_r
    
    LHV = h_CH4_r + 2*h_O2_r - h_CO2_p - 2*h_H2O_p   # Maximun possible Heat Loss
    
    return (H_p - H_r + (H_loss*LHV))                # This is the function that must be minimized

# Numerical Derivative of the above function using Forward Differencing

def f_dash(T, H_loss) :
    return (f(T + 1e-06, H_loss) - f(T, H_loss))/1e-06


# Iterating using Newton-Raphson

tol = 1e-08       # Desired Tolerance
iter = 0
alpha = 0.5       # Relaxation Factor
FT = []           # This initially-empty array will store all T_guess values

H_loss = np.linspace(0, 1, 11)      # This is the Heat Loss Factor (a coefficient to be varied from 0 to 1)

for i in H_loss :
    T_guess = 1500      

    while(abs(f(T_guess, i)) > tol) :
        T_guess = T_guess - alpha*(f(T_guess, i)/f_dash(T_guess, i))      
        iter = iter + 1

    FT.append(T_guess)             # Insering current T_guess at the end of Temp array

print(FT)

#####################################################################################################################

# Plotting the Results

plt.plot(H_loss, FT, '-o', color = 'orange')


plt.xlabel('Heat Loss Factor', fontsize=12, fontweight = 'bold')
plt.ylabel('Flame Temparature [K]', fontsize=12, fontweight = 'bold')

plt.grid('on')
plt.title('FT vs Heat Loss Factor Plot for Methane-Air Combustion', fontsize=13, fontweight='bold')

############################################################################################################

# Tabulating Results from both Python & Cantera

print("---------------------------------------------------")
print("Heat Loss Factor", 'tt', "Python")
print("---------------------------------------------------")

for x in range(len(H_loss)) :
    print(H_loss[x], 'tt', FT[x])
print("---------------------------------------------------")


plt.show()