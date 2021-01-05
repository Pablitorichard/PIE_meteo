# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 15:07:29 2020

@author: 33676
"""

#LIBRARY ----------------------------------------------------------------------

# Path handling
from pathlib import Path
import os
import inspect
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))

#Other libraries
import numpy as np

#My own functions
from create_netcdf import create_netcdf

#------------------------------------------------------------------------------


def v_stripe_test(path, Lx, Ly, Nx, Ny, T, Nt):
    
    handle = create_netcdf(path, Lx, Ly,T,  Nx, Ny, Nt)
#ATTRIBUTES -------------------------------------------------------------------

    handle.T = T
    handle.dt = T/Nt
    handle.Nt = Nt
#------------------------------------------------------------------------------

#INITIALIZATION ---------------------------------------------------------------
        
    #Bubble creation
    [X,Y] = np.mgrid[0:Nx,0:Ny]
    F = 300*np.ones((Nx, Ny))
    depth = -15
    # Half width of the stripe
    dY = Ny//10
    # Distance between the side of the X axis and the first V-profile
    dX = Nx//8
    Y_ind =  np.arange(Ny//2 - dY,Ny//2 + dY)
    X_ind =  np.arange(dX, Nx - dX)

    F[X_ind, Y_ind] = F[X_ind, Y_ind] +\
            (np.abs(Y[X_ind, Y_ind] - Ny//2)/dY -1)*depth
    
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
                 

    
    
    #Potential temperature
    handle['theta_t'][:,:,0] = F
    handle['theta_t'][:,:,1] = F
    
    #Initial displacement guess for advection
    handle['alpha_u'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_v'][:,:,0] = np.zeros((Nx,Ny))
    
    #Uniform wind
    handle['ut'][:,:,0] = np.zeros((Nx,Ny))
    handle['vt'][:,:,0] = np.zeros((Nx,Ny))
    
    #time storage
    handle['t'][0] = 0
    handle['t'][1] = handle.dt
#------------------------------------------------------------------------------
    
    handle.close()
    
#------------------------------------------------------------------------------
   