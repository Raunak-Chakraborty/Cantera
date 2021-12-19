# Program to study the Rate of Change of Molar Concentrations of H2O, O2 & OH for the
# Initial Temparatures of 500 & 1000 K, while Initial Pressure is held constant (5 atm)

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

Temparature_range = [500, 1000]

##########################################################################################

# Loop 1 : Looping through Range of Temparatures
for i in Temparature_range :

    # Specifying the Initial Conditions
    gas.TPX = i, 5*ct.one_atm, {'CH4' : 1, 'O2' : 2, 'N2' : 7.52}

    # In the Ideal Gas Reactor Model, T is used instead of the Total Internal Energy U
    # as a State Variable
    reactor = ct.IdealGasReactor(gas)

    # The following creates a Network of Reactors for the purpose of Time Integration
    reactor_network = ct.ReactorNet([reactor])

    # The SolutionArray Class helps organize output data. Additionally, we want time in ms
    states = ct.SolutionArray(gas, extra = ['t_ms', 't'])
    time = 0

    # Loop 2 : Time Integration Loop 
    for j in range(0, 10000) :

        time += 1e-3
        reactor_network.advance(time)
        states.append(reactor.thermo.state, t_ms = time*1e3, t = time)
    # End of Loop 2 
    
    plt.subplot(1, 3, 1)    # plt.subplot(n_rows, n_columns, index)
    plt.plot(states.t_ms, states.X[:, gas.species_index('H2O')])
    plt.grid('on')
    plt.xlabel('Time [ms]', fontsize = 10, fontweight = 'bold')
    plt.ylabel('Mole Fraction', fontsize = 10, fontweight = 'bold')
    plt.title('Mole Fraction of H2O VS Time', fontsize = 12, fontweight = 'bold')

    plt.subplot(1, 3, 2)      
    plt.plot(states.t_ms, states.X[:, gas.species_index('O2')])
    plt.grid('on')
    plt.xlabel('Time [ms]', fontsize = 10, fontweight = 'bold')
    plt.title('Mole Fraction of O2 VS Time', fontsize = 12, fontweight = 'bold')

    plt.subplot(1, 3, 3)       
    plt.plot(states.t_ms, states.X[:, gas.species_index('OH')])
    plt.grid('on')
    plt.xlabel('Time [ms]', fontsize = 10, fontweight = 'bold')
    plt.title('Mole Fraction of OH VS Time', fontsize = 12, fontweight = 'bold')

    plt.show()
# End of Loop 1
