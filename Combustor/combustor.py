# Program to demonestrate the Variation of Temperature & Volumetric
# Heat Release Rate with Residence Time for a Constant Pressure
# Combustor using Methane as Fuel

# Residence Time may be defined as the amount of time, a particle
# spends in the concerned Control Volume

# Residence Time = Mass / Mass Flow Rate

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

#############################################################################

gas = ct.Solution('gri30.yaml')

# Specifying Inlet Conditions as a mixture of Methane & Air
gas.TP = 300.0, ct.one_atm
gas.set_equivalence_ratio(0.5, 'CH4:1.0', 'O2:1.0, N2:3.76')
inlet = ct.Reservoir(gas)

gas.equilibrate('HP', 'auto')

combustor = ct.IdealGasReactor(gas)
combustor.volume = 1.0    # in m^3

# Storing Exhaust Conditions after Equillibrium is reached
exhaust = ct.Reservoir(gas)

# Here, we shall use a Variable Mass Flow Rate to mantain the
# Residence Time inside the Reactor as a constant
def m_dot(t) :
    return (combustor.mass/residence_time)

# Defining Inlet & Outlet Mass Flow Controllers
inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=m_dot)
outlet_mfc = ct.PressureController(combustor, exhaust, master = inlet_mfc, K=0.01)

sim = ct.ReactorNet([combustor])
states = ct.SolutionArray(gas, extra=['t_res'])

# Initial Residence Time
residence_time = 0.1

# Considering 500 K to be the Auto-Ignition Temperature of the Mixture
while combustor.T > 500 :

    sim.set_initial_time(0.0)
    sim.advance_to_steady_state()
    
    print('t_res = {:.2e}, T = {:.1f}'.format(residence_time, combustor.T))
    states.append(combustor.thermo.state, t_res = residence_time)

    residence_time *= 0.9

#############################################################################

# Plotting the Results
plt.figure(1)
plt.suptitle('Variation of Q_dot & Temperature with Residence Time', fontweight='bold')

plt.subplot(2, 1, 1)
plt.plot(states.t_res, states.heat_release_rate, '-o', color='black')
plt.xlabel('Residence Time [s]', fontweight='bold')
plt.ylabel('Vol. Heat Release Rate [W/m$^3$]', fontweight='bold')
plt.grid('on')

plt.subplot(2, 1, 2)
plt.plot(states.t_res, states.T, '-o', color='red')
plt.xlabel('Residence Time [s]', fontweight='bold')
plt.ylabel('Temperature [K]', fontweight='bold')
plt.grid('on')

plt.show()

