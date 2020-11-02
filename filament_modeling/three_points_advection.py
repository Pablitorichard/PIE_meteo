# -*- coding: utf-8 -*-
"""
Spyder Editor

Three time points semi-Lagrangian scheme. 
Corresponds to Section 2a of Staniforth.
"""

import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt


# GEOMETRY ------------------------------------------------------
dt = 0.1 # Time increment
Nt = 30 # number of time steps

dx = 1 # Spatial increment
Nx = 100 # Number of horizontal steps
x_grid = np.arange(0,Nx,dx) # Abscissas


# INPUTS --------------------------------------------------------
# All types of data will be stored with arrays, with space as the first 
# dimension for space and the second for time. Note that in a 
# computationally effective implementation you can avoid storing all
# all the previous time steps.

F_m1 = np.zeros(Nx)
F_m1[Nx//2] = 1  # Advected field at time t0-dt is a Dirac

F_0 = np.zeros(Nx)
F_0[Nx//2] = 1  # Advected field at time t0 is still a Dirac

F = np.zeros((Nx,Nt+1)) # We need a +1 in dimension as we need an initial
# condition at t0-dt. The indices are weird for this array but either 
# hime or U has to be weird.
F [:,0] = F_m1
F [:,1] = F_0


alpha_0 = np.zeros(Nx) # Initial guess for the displacement

alpha = np.zeros((Nx,Nt))
alpha [:,0] = alpha_0


U = 10*np.ones((Nx,Nt)) # Array of velocities. The velocity is uniform
# in time and space.


# LOOP ----------------------------------------------------------


for k in range (0, Nt-1):
    # U_cubic is a function which when called on a grid realizes
    # the actual interpolation. The coarguments implies Dirichlet boundary
    # conditions (The filed is equal to zero outside).
    U_cubic = interpolate.interp1d(x_grid, U[:,k], kind='cubic',
                                   bounds_error = False, fill_value = 0) 
    
    
    # The approximate displacement is updated
    alpha [:,k+1] = dt * U_cubic(x_grid - alpha[:,k])
    
    # F_cubic is a function which when called on a grid realizes
    # the actual interpolation
    F_cubic = interpolate.interp1d(x_grid,F[:,k], kind='cubic',
                                   bounds_error = False, fill_value = 0)  
    
    F [:,k+2] = F_cubic(x_grid - 2*alpha[:,k+1])
    
    
# DISPLAY -------------------------------------------------------

plt.plot(x_grid, F[:,1], label= 't=0')
plt.plot(x_grid, F[:,11], label= 't=10')
plt.legend()
plt.show()