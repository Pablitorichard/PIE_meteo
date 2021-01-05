# -*- coding: utf-8 -*-
"""
MWE for 2D interpolation function

@author: Olivier Goux
"""

import numpy as np
from upstream_interp import upstream_interp
import matplotlib.pyplot as plt
import time

Nx = 256
Ny = 128
ratio_x = Nx//10 # 1 arrow every <ratio> point
ratio_y = Ny//10
[x_grid, y_grid] = np.mgrid[0:Nx, 0:Ny]

F = np.zeros((Nx, Ny))
bubble_indices = np.where ( (x_grid - Nx//2)**2 +
                                (y_grid - Ny//2)**2 < 20**2 )
[X,Y] = np.mgrid[0:Nx,0:Ny]
F = 300*np.ones((Nx, Ny))
depth = -15
# Half width of the stripe
dY = Ny//10
# Distance between the side of the X axis and the first V-profile
dX = Nx//8


F[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] = \
    F[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] +\
    (np.abs(Y[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] - Ny//2)/dY -1)*depth

left_ind= np.where(np.logical_and(\
                (X - dX)**2 + (Y - Ny//2)**2 < (dY)**2,
                X < dX ))
F[left_ind] = F[left_ind] +\
         (np.sqrt((Y[left_ind]-Ny//2)**2 + (X[left_ind]-dX)**2 ) \
         /dY -1)*depth
             
right_ind= np.where(np.logical_and(\
                (X - (Nx - dX - 1))**2 + (Y - Ny//2)**2 < dY**2,
                X > (Nx - dX - 1) ))
F[right_ind] = F[right_ind] +\
         (np.sqrt((Y[right_ind]-Ny//2)**2 + \
                  (X[right_ind]-(Nx-dX-1))**2 )/dY -1)*depth
#F[60:65,60:65]=1

# alpha_x = 0.5*np.ones((Nx,Ny)) + 10*np.random.rand(Nx,Ny)
# alpha_y = 0.5*np.ones((Nx,Ny)) + 10*np.random.rand(Nx,Ny)

alpha_x = 30.5*np.ones((Nx,Ny)) 
alpha_y = 30.5*np.ones((Nx,Ny))

"""
Problem: Finding the field F_upstream such that its value in (x,y) 
         represents the value of the field F in (x-alpha_x, y-alpha_y) 
"""

t0=time.time()
F_cubic = upstream_interp(alpha_x, alpha_y, F, method='bicubic')
t1=time.time() - t0

t0=time.time()
F_linear = upstream_interp(alpha_x, alpha_y, F, method='linear')
# F_failure = interpolate.griddata(
#             (x_grid.flatten(), y_grid.flatten()), F.flatten(),
#             (x_grid - alpha_x, y_grid - alpha_y),
#             method='cubic',fill_value=0) 
t2=time.time() - t0


plt.subplot(2,2,1)
plt.pcolormesh(F.T)
plt.colorbar()
plt.xlabel('Before interpolation')
plt.subplot(2,2,2)
plt.pcolormesh(F_linear.T)
plt.colorbar()
plt.quiver(x_grid[0:Nx:ratio_x,0],y_grid[0,0:Ny:ratio_y],
                  alpha_x[0:Nx:ratio_x, 0:Ny:ratio_y],
                  alpha_y[0:Nx:ratio_x, 0:Ny:ratio_y],
                  color='b', width=0.003)
plt.xlabel('Linear interpolation, dt = ' + str(round(t2,4)))
plt.subplot(2,2,3)
plt.pcolormesh(F_cubic.T)
plt.xlabel('Cubic interpolation, dt = ' + str(round(t1,4)))
plt.colorbar()
plt.subplot(2,2,4)
plt.pcolormesh((F_cubic-F_linear).T)
plt.colorbar()
plt.xlabel('Cubic - Linear')

print("Computation time for upstream_interp: ", t1, " ; for griddata: ",t2)