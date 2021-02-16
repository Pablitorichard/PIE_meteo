# -*- coding: utf-8 -*-
"""
MWE for 2D interpolation function

@author: Olivier Goux
"""

import numpy as np
from upstream_interp import upstream_interp
import matplotlib.pyplot as plt
import time

Nx = 20
Ny = 10

ratio_x = Nx//10 # 1 arrow every <ratio> point
ratio_y = Ny//10
[x_grid, y_grid] = np.mgrid[0:Nx, 0:Ny]

F = np.zeros((Nx, Ny))
bubble_indices = np.where ( (x_grid - Nx//2)**2 +
                                (y_grid - Ny//2)**2 < 20**2 )
[X,Y] = np.mgrid[0:Nx,0:Ny]

depth = -15
# Half width of the stripe
dY = Ny//10
# Distance between the side of the X axis and the first V-profile
dX = Nx//6


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


#F[Nx//2: Nx//2+1, Ny//2: Ny//2 +1] = 1

# F = np.tile(F, (2,1,1))
# F[0,...] *= 2

[X, Y] = np.mgrid[0:Nx, 0:Ny]

alpha_x = 0.2*np.random.randint(5, size=(Nx,Ny))
alpha_y = 0.65*np.random.randint(5, size=(Nx,Ny))
# alpha_x = 0.9 * np.ones((Nx,Ny))
# alpha_x[X>Nx//2] *= 0
# alpha_y = 0*np.ones((Nx,Ny))
"""
Problem: Finding the field F_upstream such that its value in (x,y) 
         represents the value of the field F in (x-alpha_x, y-alpha_y) 
"""

t0=time.time()
F_cubic = upstream_interp(alpha_x, alpha_y, F, method='bicubic')
t1=time.time() - t0

t0=time.time()
F_cubic_V1 = upstream_interp(alpha_x, alpha_y, F, method='bicubic_V1')
t2=time.time() - t0

t0=time.time()
F_cubic_V2 = upstream_interp(alpha_x, alpha_y, F, method='bicubic_V2')
t3=time.time() - t0



plt.subplot(2,2,1)
plt.pcolormesh(F.T)


plt.subplot(2,2,2)
plt.pcolormesh(F_cubic.T)
plt.colorbar()
plt.xlabel('cubic naive, dt = ' + str(round(t1,6)))


plt.subplot(2,2,3)
plt.pcolormesh(F_cubic_V1.T)
plt.colorbar()
plt.xlabel('cubic optimized 1, dt = ' + str(round(t2,6)))


plt.subplot(2,2,4)
plt.pcolormesh(F_cubic_V2.T)
plt.colorbar()
plt.xlabel('cubic optimized 2 = ' + str(round(t3,6)))



print (np.max(np.abs(F_cubic_V2 - F_cubic)))