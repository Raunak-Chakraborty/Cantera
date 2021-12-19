import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.xml')

gas.TPX = 1500, ct.one_atm, {'CH4' : 1, 'O2' : 2, 'N2' : 7.52}

reac = ct.IdealGasConstPressureReactor(gas, name = 'R1')

reac_net = ct.ReactorNet([reac])

for i in range(10) :
    reac.add_sensitivity_reaction(i)

reac_net.rtol = 1e-6
reac_net.atol = 1e-15

reac_net.rtol_sensitivity = 1e-6
reac_net.atol_sensitivity = 1e-6

states = ct.SolutionArray(gas, extra = ['t_ms', 's2', 's3'])

for t in np.arange(0, 2e-3, 5e-6) :

    reac_net.advance(t)
    s2 = reac_net.sensitivity('OH', 2)   
    s3 = reac_net.sensitivity('OH', 3)
    states.append(reac.thermo.state, t_ms = 1000*t, s2 = s2, s3 = s3)

    print('%10.3e %10.3f %10.3f %14.6e %10.3f %10.3f' % 
    (reac_net.time, reac.T, reac.thermo.P, reac.thermo.u, s2, s3))

plt.figure(1)

plt.subplot(2, 2, 1)
plt.plot(states.t_ms, states.T)

plt.subplot(2, 2, 2)
plt.plot(states.t_ms, states('OH').X)

plt.subplot(2, 2, 3)
plt.plot(states.t_ms, states('CH4').X)

plt.subplot(2, 2, 1)
plt.plot(states.t_ms, states('H').X)

plt.figure(2)
plt.plot(states.t_ms, states.s2)
plt.plot(states.t_ms, states.s3)

plt.show()




