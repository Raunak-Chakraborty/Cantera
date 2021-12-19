# Program to solve a Coupled System of ODEs using Multivariate Newton-Raphson Solver

"""
We consider the following Equations :
    
    d(y1)/dt = (-0.04 * y1) + (1e4 * y2 * y3)
    d(y2)/dt = (0.04 * y1) - (1e4 * y2 * y3) - (3 * 1e7 * (y2)^2)
    d(y3)/dt = (3 * 1e7 * (y2)^2)

After applying Backward Differencing, we arrive at the following :

    f1 = y1_old - y1 + dt * ((-0.04 * y1) + (1e4 * y2 * y3)) = 0
    f2 = y2_old - y2 + dt * ((0.04 * y1) - (1e4 * y2 * y3) - (3 * 1e7 * pow(y2, 2))) = 0
    f3 = y3_old - y3 + dt * (3e7 * pow(y2, 2)) = 0

We shall try to solve for y1, y2, y3 using the Multivariate Newton-Raphson Solver, where the 
update occurs as follows :

    X_new = X_old - (Relaxation Factor) * (Inverse of the Jacobian Matrix) * (Matrix containing old values)
 => x_(n+1) = x_n - (apha) * [J]^(-1) * f(x_n)

"""

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

####################################################################################################

# Defining the 3 Non-Linear Equations & their Jacobian Matrix

def f1(y1, y2, y3, y1_old, dt) :
    return (y1_old - y1 + dt * ((-0.04 * y1) + (1e4 * y2 * y3)))

def f2(y1, y2, y3, y2_old, dt) :
    return (y2_old - y2 + dt * ((0.04 * y1) - (1e4 * y2 * y3) - (3 * 1e7 * pow(y2, 2))))

def f3(y1, y2, y3, y3_old, dt) :
    return (y3_old - y3 + dt * (3e7 * pow(y2, 2)))

# Jacobain Matrix calculated using Numerical Differentiation
def Jacobian(y1, y2, y3, y1_old, y2_old, y3_old, dt) :
    h = 1e-8        # Step-Size for Numerical Differentiation   
    J = np.ones((3, 3))
    # Row 1 :
    J[0, 0] = (f1(y1 + h, y2, y3, y1_old, dt) - f1(y1, y2, y3, y1_old, dt)) / h
    J[0, 1] = (f1(y1, y2 + h, y3, y1_old, dt) - f1(y1, y2, y3, y1_old, dt)) / h
    J[0, 2] = (f1(y1, y2, y3 + h, y1_old, dt) - f1(y1, y2, y3, y1_old, dt)) / h
    # Row 2 :
    J[1, 0] = (f2(y1 + h, y2, y3, y2_old, dt) - f2(y1, y2, y3, y2_old, dt)) / h
    J[1, 1] = (f2(y1 , y2 + h, y3, y2_old, dt) - f2(y1, y2, y3, y2_old, dt)) / h
    J[1, 2] = (f2(y1 , y2, y3 + h, y2_old, dt) - f2(y1, y2, y3, y2_old, dt)) / h
    # Row 3 :
    J[2, 0] = (f3(y1 + h, y2, y3, y3_old, dt) - f3(y1, y2, y3, y3_old, dt)) / h
    J[2, 1] = (f3(y1, y2 + h, y3, y3_old, dt) - f3(y1, y2, y3, y3_old, dt)) / h
    J[2, 2] = (f3(y1, y2, y3 + h, y3_old, dt) - f3(y1, y2, y3, y3_old, dt)) / h

    return J

####################################################################################################

# Following are the Initial Guess Values for y1, y2 & y3
y1_old = 1
y2_old = 0
y3_old = 0

# Creating a Matrix for the above values
Y_old = np.ones((3, 1))
Y_old[0] = y1_old
Y_old[1] = y2_old
Y_old[2] = y3_old

F = np.copy(Y_old)

dt = 0.1       # Time-Step Size
t = np.arange(0, 600, dt)   # creating the Time Array

alpha = 0.8     # Relaxation Factor
iter = 1        # Initial Iteration Value

y1_new = []     # Array for storing newly calculated y1 value for each time-step
y2_new = []     # Array for storing newly calculated y2 value for each time-step
y3_new = []     # Array for storing newly calculated y3 value for each time-step
error = []      # Array for storing newly calculated error value for each time-step

####################################################################################################

# Outer Time Loop
for i in range(0, len(t)) :
    
    err = 1e-09
    tol = 1e-12

    # Inner Newton-Raphson Loop
    while err > tol :

        # Constructing array F by sending old values to functions f1, f2, f3
        F[0] = f1(Y_old[0], Y_old[1], Y_old[2], y1_old, dt)
        F[1] = f2(Y_old[0], Y_old[1], Y_old[2], y2_old, dt)
        F[2] = f3(Y_old[0], Y_old[1], Y_old[2], y3_old, dt)

        # Constructing array J by sending old values to function Jacobian
        J = Jacobian(Y_old[0], Y_old[1], Y_old[2], y1_old, y2_old, y3_old, dt)

        # Matrix Y_new is calculated using Matrices Y_old, J & F
        Y_new = Y_old - alpha * np.matmul(inv(J), F)

        # Defining Error as the Maximum(Absolute(Y_new - Y_old))
        err = np.max(np.abs(Y_new - Y_old))

        # After each Iteration of the Inner Loop, the resulting Y_new Matrix
        # is set to Y_old, to be used by the next Iteration
        Y_old = Y_new

        iter = iter + 1
    # Inner Loop End
    
    # Only after the Inner Loop converges, this message will be printed
    # Therefore, the converged results of each time-step will be printed
    log_message = 'Time = {0}, y1 = {1}, y2 = {2}, y3 = {3}, error = {4}'.format(t[i], Y_new[0], Y_new[1], Y_new[2], err)
    print(log_message)

    # Value Update fotr the Outer Loop
    Y_old = Y_new           

    # Storing newly-calculated values in Arrays 
    y1_new.append(Y_new[0])
    y2_new.append(Y_new[1])
    y3_new.append(Y_new[2])
    error.append(err)

    # Updating the values which make up the Y_old Matrix
    y1_old = Y_new[0]
    y2_old = Y_new[1]
    y3_old = Y_new[2]
# Outer Loop End

####################################################################################################

# Plotting the values of y1, y2 & y3 with Time
plt.plot(t, y1_new, '-', color = 'blue', label = 'y1')
plt.plot(t, y2_new, '-', color = 'red', label = 'y2')
plt.plot(t, y3_new, '-', color = 'green', label = 'y3')
plt.legend()
plt.grid('on')
plt.xlabel('Time [seconds]', fontsize = '13', fontweight = 'bold')
plt.ylabel('Value of y', fontsize = '13', fontweight = 'bold')
plt.title('Variation of y1, y2 & y3 with Time (dt = 0.1 s)', fontsize = '14', fontweight = 'bold')
plt.show()

# Plotting values of y2 with Time
plt.plot(t, y2_new, '-', color = 'red', label = 'y2')
plt.legend()
plt.grid('on')
plt.xlabel('Time [seconds]', fontsize = '13', fontweight = 'bold')
plt.ylabel('Value of y2', fontsize = '13', fontweight = 'bold')
plt.title('Variation of y2 with Time (dt = 0.1 s)', fontsize = '14', fontweight = 'bold')
plt.show()

# Plotting Error with Time
plt.plot(t, error, '-', color = 'black', label = 'error')
plt.legend()
plt.grid('on')
plt.xlabel('Time [seconds]', fontsize = '13', fontweight = 'bold')
plt.ylabel('Maximum Error', fontsize = '13', fontweight = 'bold')
plt.title('Variation of Error with Time (dt = 0.1 s)', fontsize = '14', fontweight = 'bold')
plt.show()







         