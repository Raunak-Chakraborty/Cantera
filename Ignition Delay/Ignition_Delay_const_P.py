# Program to depict the variation of Auto Ignition characteristics of Methane with
# varying Initial Temparature (950 - 1450 K), while Initial Pressure is held constant (5 atm)

"""
CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2

We shall employ Cantera's Ideal Gas Reactor Model & run the Simulation for 10 seconds

We assume the following :
Auto-Ignition occurs when the Reactor Temparature reaches (T_initial + 400) Kelvin

"""
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.xml')

temparature_range = np.linspace(950, 1450, 5)

# The following Array will store the values Ignition Delay for each Temparature
ignition_delay = []

##########################################################################################

# Loop 1 : Looping through Range of Temparatures
for i in temparature_range :

    # Specifying the Initial Conditions
    gas.TPX = i, 5*101325, {'CH4' : 1, 'O2' : 2, 'N2' : 7.52}
    
    # In the Ideal Gas Reactor Model, T is used instead of the Total Internal Energy U
    # as a State Variable
    reac = ct.IdealGasReactor(gas)

    # The following creates a Network of Reactors for the purpose of Time Integration
    reac_net = ct.ReactorNet([reac])

    # The SolutionArray Class helps organize output data. Additionally, we want time in ms
    states = ct.SolutionArray(gas, extra = ['t_ms', 't'])
    time = 0
    
    # The counter has been used to ensure that, only the first value of Ignition Delay
    # is inserted into the Array
    counter = 0

    # Loop 2 : Time Integration Loop 
    for j in range(0, 10000) :
        
        # This is the final time, rather than time-step. For each value of j, our simulaton
        # will advance to the updated value of final time
        time += 1e-3

        # The Advance Method in Cantera calculates the State of our System at the specified
        # final time. Cantera uses an appropriate time-step on it's own 
        reac_net.advance(time)

        # Here, we append the results into the states object
        states.append(reac.thermo.state, t_ms = time*1e3, t = time)
        
        # For each loop, we assign the calculated Temparature to a new variable
        T_new = states.T
        
        # This if statement is used to determine when Auto-Ignition first occurs
        # We are only interested in the first instance when Ignition occurs & 
        # hence the use of the counter
        if (T_new[j] - 1250 > 400 and counter==0) :
            ignition_delay.append(time * 1e3)
            counter += 1
        # End of Loop 2

    # Plotting Temparature VS Time(in ms)
    plt.figure(1)       
    plt.plot(states.t_ms, states.T)
    plt.grid('on')
    plt.legend([str(T) + '[K]' for T in temparature_range])
    plt.xlabel('Time [ms]', fontsize = 13, fontweight = 'bold')
    plt.ylabel('Temparature [K]', fontsize = 13, fontweight = 'bold')
    plt.title('Auto-Ignition Characteristics of Methane with different \n Initial Temparatures (Initial Pressure = 5 atm)', fontsize = 15, fontweight = 'bold')
# End of Loop 1

##########################################################################################

# Plotting Ignition Delay(in ms) VS Initial Temparature
plt.figure(2)
plt.plot(temparature_range, ignition_delay, '-o', color = 'black')
plt.grid('on')
plt.xlabel('Temparature [K]', fontsize = 13, fontweight = 'bold')
plt.ylabel('Ignition Delay [ms]', fontsize = 13, fontweight = 'bold')
plt.title('Variation of Ignition Delay with Initial \n Temparature (Initial Pressure = 5 atm)', fontsize = 15, fontweight = 'bold')

plt.show()
        



