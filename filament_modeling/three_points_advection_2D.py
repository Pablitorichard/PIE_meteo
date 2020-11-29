# -*- coding: utf-8 -*-
"""
Spyder Editor

Three time points semi-Lagrangian scheme in 2D. 
Corresponds to Section 2b of Staniforth.
"""

import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time

# GEOMETRY ------------------------------------------------------
dt = 0.1 # Time increment
Nt = 31 # number of time steps


Lx = 100 # Zonal length of the domain 
Ly = 100 # Meridional length of the domain

Nx = 300 # Number of zonal points.
Ny = 300 # Number of meridional points.

dx = Lx/Nx # Zonal resolution
dy = Ly/Ny # Meridional resolution

grid = np.mgrid[0:Lx:dx, 0:Ly:dy]# Array of dimensions (2, Nx, Ny)
# grid[0,:,:] is the X-coordinate of each point, grid[1,:,:] is the 
# Y-coordinate for each point
X_grid = grid[0,:,:]
Y_grid = grid[1,:,:]
X_coord = grid[0,:,0]
Y_coord = grid[1,0,:]

cx = Lx//2 #Zonal coordinates of the center of the bubble
cy = Ly//2 #Meridional coordinates of the center of the bubble
radius = 10 # radius of the bubble

ratio = Nx//10 # 1 arrow every <ratio> point

# INPUTS --------------------------------------------------------
"""
All types of data will be stored with arrays, with the zonal spatial
coordinate as the first dimension, the meridional spatial coordinate
as the second dimension, and the time coordinate as the third dimension.
Note that in a computationally effective implementation you should avoid 
storing all the previous time steps.
"""
bubble_indices = np.where ( (grid[0,:,:] - cx)**2 + (grid[1,:,:] - cy)**2
                           < radius**2 )

F_m1 = np.zeros((Nx,Ny))
F_m1[bubble_indices] = 1  # Advected field at time t0-dt is a Dirac

F_0 = np.zeros((Nx,Ny))
F_0[bubble_indices] = 1  # Advected field at time t0 is still a Dirac

F = np.zeros((Nx,Ny,Nt+1)) # F is the field advected by U. (theta_tp in the 
# Wirth article). We need a +1 in dimension as we need an initial
# condition at t0-dt. The indices are weird for this array but either 
# hime or U has to be weird.
F [:,:,0] = F_m1
F [:,:,1] = F_0


alpha_0 = np.zeros((2,Nx,Ny)) # Initial guess for the displacement
alpha = np.zeros((2,Nx,Ny,Nt)) 
alpha [:,:,:,0] = alpha_0
# [0,:,:,:] -> zonal displacement ; U[1,:,:,:] -> meridional displacement

U = 10*np.ones((2,Nx,Ny,Nt)) # Array of velocities. The velocity is uniform
# in time and space. 
# U[0,:,:,:] -> zonal wind ; U[1,:,:,:] -> meridional wind

# distort_x = np.cos(2*np.pi*X_grid/Lx)
# distort_x = np.repeat(distort_x[:, :, np.newaxis], Nt, axis=2)
# distort_y = np.sin(2*np.pi*Y_grid/Ly)
# distort_y = np.repeat(distort_y[:, :, np.newaxis], Nt, axis=2)
# U[0,:,:,:] =U[0,:,:,:] * distort_x 
# U[1,:,:,:] =U[1,:,:,:] * distort_y


# LOOP ----------------------------------------------------------

t0 = time.time()
for k in range (0, Nt-1):
    
    #   INSERT PSEUDO SPECTRAL FUNCTION WHICH TAKES F AND SEND U[k]
    
    # The displacement is updated by interpolating the velocity field at time 
    # tk at the points x - alpha^k. alpha[0,:,:,k+1] is the zonal component, 
    # alpha[1,:,:,] is the meridional component.
    alpha [0,:,:,k+1] = dt * interpolate.griddata( (X_grid.flatten(),
                                    Y_grid.flatten()), U[0,:,:,k].flatten(), 
                                    (X_grid - alpha[0,:,:,k], 
                                    Y_grid - alpha[1,:,:,k]), method='cubic',
                                    fill_value = 0) 
    
    alpha [1,:,:,k+1] = dt * interpolate.griddata((X_grid.flatten(),
                                    Y_grid.flatten()), U[1,:,:,k].flatten(),
                                    (X_grid - alpha[0,:,:,k], 
                                    Y_grid - alpha[1,:,:,k]), method='cubic',
                                    fill_value = 0) 
    
    #F at time tk+ 2*dt is udpated by interpolating F at time tk at the 
    # locations x - 2* alpha
    F [:,:,k+2] = interpolate.griddata( (X_grid.flatten(), Y_grid.flatten()) ,
                                    F [:,:,k].flatten(), 
                                    (X_grid - 2*alpha[0,:,:,k+1], 
                                    Y_grid - 2*alpha[1,:,:,k+1]), method='cubic',
                                    fill_value = 0) 
    
    
cpu_time = time.time() - t0 
print(cpu_time)
# # DISPLAY -------------------------------------------------------
fig, (( ax_t0, ax_t10,), (ax_t20, ax_t30)) =  plt.subplots(nrows=2, ncols =2,
                                                        figsize=[6.8/0.8,6.8])
im_t0 = ax_t0.pcolormesh(grid[0,:,:], grid[1,:,:], F[:,:,1], cmap=cm.binary)
ax_t0.quiver(X_coord[0:Nx:ratio], Y_coord[0:Ny:ratio],
             U[0, 0:Nx:ratio, 0:Ny:ratio, 0], U[1, 0:Nx:ratio, 0:Ny:ratio, 0],
             color='b', width=0.003)
ax_t0.set_xlabel('t=0')

ax_t10.pcolormesh(grid[0,:,:], grid[1,:,:], F[:,:,11], cmap=cm.binary)
ax_t10.quiver(X_coord[0:Nx:ratio], Y_coord[0:Ny:ratio],
             U[0, 0:Nx:ratio, 0:Ny:ratio, 10], U[1, 0:Nx:ratio, 0:Ny:ratio, 10],
             color='b', width=0.003)
ax_t10.set_xlabel('t=10')

ax_t20.pcolormesh(grid[0,:,:], grid[1,:,:], F[:,:,21], cmap=cm.binary)
ax_t20.quiver(X_coord[0:Nx:ratio], Y_coord[0:Ny:ratio],
             U[0, 0:Nx:ratio, 0:Ny:ratio, 20], U[1, 0:Nx:ratio, 0:Ny:ratio, 20],
             color='b', width=0.003)
ax_t20.set_xlabel('t=20')

ax_t30.pcolormesh(grid[0,:,:], grid[1,:,:], F[:,:,31], cmap=cm.binary)
ax_t30.quiver(X_coord[0:Nx:ratio], Y_coord[0:Ny:ratio],
             U[0, 0:Nx:ratio, 0:Ny:ratio, 30], U[1, 0:Nx:ratio, 0:Ny:ratio, 30],
             color='b', width=0.003)
ax_t30.set_xlabel('t=30')

plt.tight_layout()
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.10, 0.05, 0.7])
fig.colorbar(im_t0, cax=cbar_ax)

# plt.plot(x_grid, F[:,11], label= 't=10')
# plt.plot(x_grid, F[:,21], label= 't=20')
# plt.plot(x_grid, F[:,31], label= 't=30')
# plt.legend()
# plt.show()