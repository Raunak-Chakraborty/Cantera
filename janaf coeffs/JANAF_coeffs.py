import cantera as ct
import numpy as np
import math
import matplotlib.pyplot as plt

def cp(cant,t,p,sol):
    cant.TPX=t,p,sol
    return(cant.cp_mass)

def janafPolynomials(p,sol,defn,Tlow,Tcommon,Thigh,deltaT):
    '''
    p - pressure
    sol - solution description
    defn - phase definition, species and reactions definition file
    t* - low, common, high temperatures
    deltat - temperature step between values
    '''
    #define cantera model
    gas = ct.Solution(defn);
    
    #calculate polynomials
    T_low=np.arange(Tlow,Tcommon,deltaT)
    T_high=np.arange(Tcommon,Thigh,deltaT)
    
    cp_low=[cp(gas,x,p,sol) for x in T_low]
    cp_high=[cp(gas,x,p,sol) for x in T_high]
    
    cp_low_polynomials = np.polyfit(T_low, cp_low, 4)
    cp_high_polynomials = np.polyfit(T_high, cp_high, 4)
    
    #calculate integration constants
    
    #get enthalpy and entropy at 298.15K
    T_ref=298.15
    gas.TPX=T_ref,p,sol
    h=gas.enthalpy_mass
    S=gas.entropy_mass
    
    R=ct.gas_constant/gas.mean_molecular_weight
    
    #reverse order and divide by R (np.ployfit has calculated coeffs*R so we must divide by R to get just the coeffs)
    cp_low_poly_rev=np.flip(cp_low_polynomials)/R
    cp_high_poly_rev=np.flip(cp_high_polynomials)/R
    
    h_off_low=h/R-np.sum(np.array([cp_low_poly_rev[i]*T_ref**(i+1)/(i+1) for i in range(len(cp_low_poly_rev))]))
    h_off_high=h/R-np.sum(np.array([cp_high_poly_rev[i]*T_ref**(i+1)/(i+1) for i in range(len(cp_high_poly_rev))]))
    
    s_off_low=S/R-np.sum(np.array([cp_low_poly_rev[i]*T_ref**(i)/(i) for i in range(1,len(cp_low_poly_rev))]))-cp_low_poly_rev[0]*math.log(T_ref)
    s_off_high=S/R-np.sum(np.array([cp_high_poly_rev[i]*T_ref**(i)/(i) for i in range(1,len(cp_high_poly_rev))]))-cp_high_poly_rev[0]*math.log(T_ref)
        
    return (cp_low_poly_rev, cp_high_poly_rev, h_off_low, h_off_high, s_off_low, s_off_high)

cpl,cph,h_l,h_h,s_l,s_h=janafPolynomials(101325,'OH:1','gri30.cti',200,1000,3500,100)

print('highCpCoeffs: ',cph, h_h, s_h)
print('lowCpCoeffs: ',cpl, h_l, s_l)