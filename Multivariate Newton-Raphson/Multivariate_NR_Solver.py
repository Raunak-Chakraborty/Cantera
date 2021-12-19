# Program to solve a Coupled System of Non-Linear Equations using Multivariate Newton-Raphson Solver

# (3 * x1) - cos(x2 . x3) - 1.5 = 0
# {4 * (x1)^2} - {625 * (x2)^2} + (2 * x3) - 1 = 0
# (20 * x3) + {(e)^(-x1 * x2)} + 9 = 0

# We shall solve for x1, x2, x3
# In a Multivariate Newton-Raphson Solver, the update occurs as follows :

# X_new = X_old - (Relaxation Factor) * (Inverse of the Jacobian Matrix) * (Matrix containing old values)
# => x_(n+1) = x_n - (apha) * [J]^(-1) * f(x_n)

# Here, the Jacobian Matrix is calculated analytically

import numpy as np
from numpy.linalg import inv

####################################################################################################

# Defining the 3 Non-Linear Equations & their Jacobian Matrix

def f1(x1, x2, x3) :
    return ((3 * x1) - (np.cos(x2 * x3))  -1.5)

def f2(x1, x2, x3) :
    return ((4 * pow(x1, 2)) - (625 * pow(x2 , 2)) + (2 * x3) - 1)

def f3(x1, x2, x3) :
    return ((20 * x3) + (np.exp(-x1 * x2)) + 9)

def Jacobian(x1, x2, x3) :
    J = np.ones((3, 3))
    # Row 1 :
    J[0, 0] = 3
    J[0, 1] = np.sin(x2 * x3) * x3
    J[0, 2] = np.sin(x2 * x3) * x2
    # Row 2 :
    J[1, 0] = 8 * x1
    J[1, 1] = 625 * 2 * x2
    J[1, 2] = 2
    # Row 3 :
    J[2, 0] = np.exp(-x1 * x2) * (-x2)
    J[2, 1] = np.exp(-x1 * x2) * (-x1)
    J[2, 2] = 20

    return J

####################################################################################################

# Specifying Guess Values for x1, x2, x3 & Tolerance for the Newton-Raphson Solver

x1 = 1
x2 = 0
x3 = 0

error = 9e-09

X_old = np.ones((3, 1))
X_old[0] = x1
X_old[1] = x2
X_old[2] = x3

F = np.copy(X_old)

alpha = 0.6
tolerance = 1e-9
iter = 1

# Setting up the Newton-Raphson Solver

while(error > tolerance) :

    F[0] = f1(X_old[0], X_old[1], X_old[2])
    F[1] = f2(X_old[0], X_old[1], X_old[2])
    F[2] = f3(X_old[0], X_old[1], X_old[2])

    J = Jacobian(X_old[0], X_old[1], X_old[2])

    X_new = X_old - alpha * np.matmul(inv(J), F)

    error = np.max(np.abs(X_new - X_old))

    X_old = X_new

    log_message = 'iteration # = {0}, x1 = {1}, x2 = {2}, x3 = {3}'.format(iter, X_new[0], X_new[1], X_new[2])
    
    print(log_message)

    iter = iter + 1

# Final Results

x1_final = X_new[0]
x2_final = X_new[1]
x3_final = X_new[2]

print("The Final Value of x1 is : ", x1_final)
print("The Final Value of x2 is : ", x2_final)
print("The Final Value of x3 is : ", x3_final)








