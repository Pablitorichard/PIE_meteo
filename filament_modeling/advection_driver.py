# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 11:31:54 2020

@author: 33676
"""

import time
from netCDF4 import Dataset
# My cocktail to make path stuff work in Python


from test_cases import bubble_test, gaussian_test, v_stripe_test
from advection_step_3P import advection_step_3P 
from spectral import streamfunc, geostwind




def advection_driver(path, pseudo_spectral_wind=1, 
                     alpha_method = 'mix', F_method='bicubic'):
    
    #UNPACK DATA ------------------------------------------------
    handle = Dataset(path, 'r+',format='NETCDF4')  
    Nt = handle.Nt
    dt = handle.dt
    Lx = handle.Lx
    Ly = handle.Ly
    Nx = handle.Nx
    Ny = handle.Ny
    dx=handle.dx
    dy=handle.dx
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
        print("Progress: ",round(100*k/Nt,2)," %")
        
        # UPDATE OF THE TROPOPAUSE WIND---------------------------------
        if (pseudo_spectral_wind == True):
            psi = streamfunc(Lx, Ly, Nx, Ny, theta)
            vt,ut = geostwind(Lx, Ly, Nx, Ny, psi)
           
        else:
            ut = ut_minus
            vt = vt_minus 
        #---------------------------------------------------------------
                
        
        #ADVECTION AT THE TROPOPAUSE------------------------------------
        alpha_u, alpha_v, theta_plus = \
        advection_step_3P(alpha_u_minus, alpha_v_minus,theta_minus,
                          dt, ut, vt, dx, dy, alpha_method = alpha_method,
                          F_method=F_method)
        #---------------------------------------------------------------
        
        
        
        
        # Write new quantities to netCDF
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
