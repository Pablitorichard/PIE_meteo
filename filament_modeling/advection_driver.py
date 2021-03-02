# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 11:31:54 2020

@author: 33676
"""

import time
from netCDF4 import Dataset
# My cocktail to make path stuff work in Python

import numpy as np
from test_cases import bubble_test, gaussian_test, v_stripe_test
from advection_step_3P import advection_step_3P 
from spectral import geostwind, vertwind




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
    
    
    #k_hour = int(3600/dt)
    k_hour = 2
    
    ut_minus = handle['ut'][:,:,0]  
    vt_minus = handle['vt'][:,:,0]
    us_minus = handle['us'][:,:,0]  
    vs_minus = handle['vs'][:,:,0]
    
    x_grid = handle['x_grid'][:,:]
    y_grid = handle['y_grid'][:,:]
    
    alpha_ut_minus = handle['alpha_ut'][:,:,0]
    alpha_vt_minus = handle['alpha_vt'][:,:,0]
    alpha_us_minus = handle['alpha_us'][:,:,0]
    alpha_vs_minus = handle['alpha_vs'][:,:,0]
    
    Delta_z_minus = handle['Delta_z'][:,:,0]
    Delta_z = handle['Delta_z'][:,:,1]
    Delta_T_hist_minus = handle['Delta_T_hist'][:,:,0]
    Delta_T_hist = handle['Delta_T_hist'][:,:,1]
    theta_minus = handle['theta_t'][:,:,0]
    theta = handle['theta_t'][:,:,1]
    
    t0 = time.time()
    for k in range (1, Nt):
        print("Progress: ",round(100*k/Nt,2)," %")
        
        # UPDATE OF THE TROPOPAUSE WIND---------------------------------
        if (pseudo_spectral_wind == True):
            ut, vt = geostwind(Lx,Ly,theta, z=0)
            
            us, vs = geostwind(Lx,Ly,theta, z=handle.z_star)
           
        else:
            ut = ut_minus
            vt = vt_minus 
            us = us_minus
            vs = vs_minus 
        #---------------------------------------------------------------
                

        
        
        #ADVECTION AT THE TROPOPAUSE AND MORE------------------------------------
        alpha_ut, alpha_vt, theta_plus = \
        advection_step_3P(alpha_ut_minus, alpha_vt_minus,theta_minus,
                          dt, ut, vt, dx, dy, alpha_method = alpha_method,
                          F_method=F_method)
        
        alpha_us, alpha_vs, Delta_z_plus = \
        advection_step_3P(alpha_us_minus, alpha_vs_minus, Delta_z_minus,
                          dt, us, vs, dx, dy, alpha_method = alpha_method,
                          F_method=F_method)
        
        alpha_us, alpha_vs, Delta_T_hist_plus = \
        advection_step_3P(alpha_us_minus, alpha_vs_minus, Delta_T_hist_minus,
                          dt, us, vs, dx, dy, alpha_method = alpha_method,
                          F_method=F_method)
        #---------------------------------------------------------------
        
                #UPDATE OF W ---------------------------------------------------
        if ((k-1)%k_hour==0 ):# k-1 -> k ?  
            w = vertwind(Lx, Ly,theta, theta_minus, dt, z=handle.z_star)
            w_plus = vertwind(Lx, Ly,theta_plus, theta, dt, z=handle.z_star)
            w_mean = (w + w_plus)/2
            print("w: ",np.max(w)," , ", np.min(w))
            Delta_z += k_hour * dt * w_mean
            Delta_z_plus += k_hour * dt * w_mean
            
        #---------------------------------------------------------------
        print(k, ": ",np.max(Delta_z)," , ", np.min(Delta_z))
        
        
        Delta_T_disp = handle.gamma_2 * Delta_z
 
        Delta_T_cloud = handle.Delta_Tc * ( Delta_z > handle.Delta_zc ) 

        
        Delta_T_bb = Delta_T_hist + Delta_T_disp + Delta_T_cloud
        
        
        # Write new quantities to netCDF
        handle['alpha_ut'][:,:,k] = alpha_ut
        handle['alpha_vt'][:,:,k] = alpha_vt
        handle['alpha_us'][:,:,k] = alpha_us
        handle['alpha_vs'][:,:,k] = alpha_vs
        
        handle['theta_t'][:,:,k+1] = theta_plus
        
        handle['Delta_z'][:,:,k+1] = Delta_z_plus
        handle['Delta_T_hist'][:,:,k+1] = Delta_T_hist_plus
        
        handle['Delta_T_bb'][:,:,k] = Delta_T_bb
        
        handle['ut'][:,:,k] = ut
        handle['vt'][:,:,k] = vt
        handle['us'][:,:,k] = us
        handle['vs'][:,:,k] = vs
        
        handle['t'][k+1] = handle ['t'][k] + dt 
        
        # Update temporary arays
        alpha_ut_minus = alpha_ut
        alpha_vt_minus = alpha_vt
        alpha_us_minus = alpha_us
        alpha_vs_minus = alpha_vs
        
        theta_minus = theta
        theta = theta_plus
        Delta_z_minus = Delta_z
        Delta_z = Delta_z_plus
        Delta_T_hist_minus = Delta_T_hist
        Delta_T_hist = Delta_T_hist_plus
        
        ut_minus = ut
        vt_minus = vt
        us_minus = us
        vs_minus = vs
        
    cpu_time = time.time() - t0 
    print("CPU time = ",cpu_time," seconds")
    handle.close()
