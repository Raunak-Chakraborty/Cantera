# This Program compares the results of various Explicit & Implicit Methods (for solving an ODE)
# with the Analytical Solution

"""
We consider the following ODE :

    f_dash(y , t) = dy/dt = -1000y - (e)^(-t)

Explicit :

Euler's Forward Difference Method : 
   (y_(n+1) - y_n) / dt = f_dash(y_n , t_n)        Here, n is the old value
=> y_(n+1) = y_n + f_dash(y_n , t_n)*dt
Stability Value = (1 - 100*dt)
 
Implicit :

Euler's Backward Difference Method : 
   (y_n - y_(n-1))/dt = f_dash(y_n , t_n)          Here, n is the new value
=> y_n = y_(n-1) + f_dash(y_n , t_n)*dt
Stability Value = 1/(1 + 1000*dt)

"""

import numpy as np
import matplotlib.pyplot as plt

t_final = 0.1       # Maximum Time
dt = 0.0001         # Time-Step

# Defining the ODE
def f_dash(t, y) :      
    return (-1000*y) - (np.exp(-t))

# Defining function for Explicit Euler (Forward Difference)  
def explicit_euler(t_final, dt) :
    
    t = np.arange(0, t_final, dt)   # linspace takes in the number of steps as the 3rd parameter,
    y = np.zeros(len(t))            # whereas, arange requires the step-size.
    
    for i in range(1, len(t)) :
        y[i] = y[i-1] + (dt * f_dash(t[i-1], y[i-1]))
    return [t, y]       # a 2-D Array containing values of t & y

# Defining function for Implicit Euler (Backward Difference)
def implicit_euler(t_final, dt) :
    
    t = np.arange(0, t_final, dt)
    y = np.zeros(len(t))
    
    for j in range(1, len(t)) :
        y[j] = (y[j-1] - (dt * np.exp(-t[j]))) / (1 + 1000*dt)
    return [t, y]

# Defining the Analytical Solution
def analytical(t_final, dt) :

    t = np.arange(0, t_final, dt)
    y = np.zeros(len(t))

    for k in range(1, len(t)) :
        y[k] =  (np.exp(-1000*t[k]) * (1 - np.exp(999*t[k])))/999
    return [t, y]

t, y_explicit = explicit_euler(t_final, dt)     
t, y_implicit = implicit_euler(t_final, dt)
t, y_analytical = analytical(t_final, dt)

# Plotting the Results
plt.plot(t, y_implicit, color = 'blue', label = 'Implicit')
plt.plot(t, y_explicit, color = 'red', label = 'Explicit')
plt.plot(t, y_analytical, color = 'black', label = 'Analytical')
plt.legend()
plt.xlabel('t', fontweight = 'bold')
plt.ylabel('y', fontweight = 'bold')
plt.title('dt = 0.0001', fontweight = 'bold')
plt.grid('on')

plt.show()