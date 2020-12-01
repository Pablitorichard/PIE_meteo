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
from bubble_test import bubble_test
from netCDF4 import Dataset

# My cocktail to make path stuff work in Python
from pathlib import Path
import os
import inspect
# Path to 
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))



path = self_path / "out.nc"

Lx = 2048
Ly = 1024
Nx = 256
Ny = 128
dt = 30

T = 1*900
#T =48*3600
Nt = T//dt

cx = 500
cy = 500
radius = 100
wind_norm = 15
bubble_test("out.nc", Lx, Ly, Nx, Ny, T, Nt, cx, cy, radius, wind_norm)
print("init done")

ratio_x = Nx//10 # 1 arrow every <ratio> point
ratio_y = Ny//10

def three_points_advection(path, pseudo_spectral_wind):
    
    #UNPACK DATA ------------------------------------------------
    handle = Dataset(path, 'r+',format='NETCDF4')  
    Nt = handle.Nt
    dt = handle.dt
    ut_minus = handle['ut'][:,:,0]
    vt_minus = handle['vt'][:,:,0]
    x_grid = handle['x_grid'][:,:]
    y_grid = handle['y_grid'][:,:]
    alpha_u_minus = handle['alpha_u'][:,:,0]
    alpha_v_minus = handle['alpha_v'][:,:,0]
    theta_minus = handle['theta_t'][:,:,0]
    theta = handle['theta_t'][:,:,1]
    t0 = time.time()
    for k in range (1, Nt):
        
        #Update of the tropopause wind
        if (pseudo_spectral_wind == True):
            #function taking theta as input and returning ut, vt
            print('no pseudo spectral function yet, ut,vt kept contants')
        else:
            ut = ut_minus
            vt = vt_minus
        
        # The displacement is updated by interpolating the velocity field 
        alpha_u = dt * interpolate.griddata( 
            (x_grid.flatten(), y_grid.flatten()),ut.flatten(),
            (x_grid - alpha_u_minus, y_grid - alpha_v_minus),
             method='cubic', fill_value = 0) 
        
        alpha_v = dt * interpolate.griddata( 
            (x_grid.flatten(), y_grid.flatten()), vt.flatten(),
            (x_grid - alpha_u_minus, y_grid - alpha_v_minus),
             method='cubic', fill_value = 0)
        
        # theta_t at time tk+ 2*dt is udpated by interpolating F at time tk at the 
        # locations x - 2* alpha
        theta_plus = interpolate.griddata(
            (x_grid.flatten(), y_grid.flatten()), theta_minus.flatten(),
            (x_grid - 2*alpha_u, y_grid - 2*alpha_v),
            method='cubic', fill_value = 0) 
        
        # Write new quantitiesto netCDF
        handle['alpha_u'][:,:,k] = alpha_u
        handle['alpha_v'][:,:,k] = alpha_v
        handle['theta_t'][:,:,k+1] = theta_plus
        handle['ut'][:,:,k] = ut
        handle['vt'][:,:,k] = vt
        handle['t'][k+1] = handle ['t'][k] + dt 
        
        # Update temporary arays
        alpha_u_minus = alpha_u
        alpha_v_minus = alpha_v
        theta_minus = theta
        theta = theta_plus
        ut_minus = ut
        vt_minus = vt
        
    cpu_time = time.time() - t0 
    print("CPU time = ",cpu_time," seconds")
    handle.close()

three_points_advection(path, 0)

# # DISPLAY -------------------------------------------------------

def ugly_plot(path, ratio_x, ratio_y):
    # UNPACK
    #UNPACK DATA ------------------------------------------------
    handle = Dataset(path, 'r+',format='NETCDF4')  
    Nt = handle.Nt
    dt = handle.dt
    Nx = handle.Nx
    Ny = handle.Ny
    x_grid = handle['x_grid'][:,:]
    y_grid = handle['y_grid'][:,:]
    X_coord = x_grid[:,0]
    Y_coord = y_grid[0,:]
  
    
    # DISPLAY
    fig, (( ax_t0, ax_t1,), (ax_t2, ax_t3)) =  plt.subplots(nrows=2, ncols =2,
                                                            figsize=[6.8/0.8,6.8])
    im_t0 = ax_t0.pcolormesh(x_grid, y_grid, handle['theta_t'][:,:,1],
                             cmap=cm.binary)
    ax_t0.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, 0],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, 0],
                  color='b', width=0.003)
    t0_label = 't = '+ str(handle['t'][1])
    ax_t0.set_xlabel(t0_label)
    
    
    ax_t1.pcolormesh(x_grid, y_grid, handle['theta_t'][:,:,Nt//3],
                     cmap=cm.binary)
    ax_t1.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt//3-1],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt//3-1],
                  color='b', width=0.003)
    t1_label = 't = '+ str(handle['t'][Nt//3])
    ax_t1.set_xlabel(t1_label)
    
    ax_t2.pcolormesh(x_grid, y_grid, handle['theta_t'][:,:,2*Nt//3],
                     cmap=cm.binary)
    ax_t2.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, 2*Nt//3-1],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, 2*Nt//3-1],
                  color='b', width=0.003)
    t2_label = 't = '+ str(handle['t'][2*Nt//3])
    ax_t2.set_xlabel(t2_label)
    
    ax_t3.pcolormesh(x_grid, y_grid, handle['theta_t'][:,:,Nt],
                     cmap=cm.binary)
    ax_t3.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt-1],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt-1],
                  color='b', width=0.003)
    t3_label = 't = '+ str(handle['t'][Nt-1 ])
    ax_t3.set_xlabel(t3_label)
    
    plt.tight_layout()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.10, 0.05, 0.7])
    fig.colorbar(im_t0, cax=cbar_ax)
    handle.close()
    plt.show()
    
ugly_plot(path, ratio_x,ratio_y)