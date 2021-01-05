# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 18:34:12 2020

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


def bubble_test(path, Lx, Ly, Nx, Ny, T, Nt, cx, cy, radius):
    
    handle = create_netcdf(path, Lx, Ly,T,  Nx, Ny, Nt)
#ATTRIBUTES -------------------------------------------------------------------

    handle.T = T
    handle.dt = T/Nt
    handle.Nt = Nt
#------------------------------------------------------------------------------

#INITIALIZATION ---------------------------------------------------------------
        
    #Bubble creation
    bubble_indices = np.where ( np.abs(handle['x_grid'][:,:] - cx) + \
                               np.abs(handle['y_grid'][:,:] - cy) < radius )
    F = np.zeros((Nx, Ny))
    F[bubble_indices] = 1
    
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
   