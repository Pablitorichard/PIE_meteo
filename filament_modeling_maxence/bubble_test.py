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
from netCDF4 import Dataset

#My own functions
from create_netcdf import create_netcdf

#------------------------------------------------------------------------------


def bubble_test(path, Lx, Ly, Nx, Ny, T, Nt, cx, cy, radius, wind_norm):
    
    handle = create_netcdf(path, Lx, Ly,T,  Nx, Ny, Nt)
#ATTRIBUTES -------------------------------------------------------------------

    handle.T = T
    handle.dt = T/Nt
    handle.Nt = Nt
#------------------------------------------------------------------------------

#INITIALIZATION ---------------------------------------------------------------
        
    #Bubble creation
    bubble_indices = np.where ( (handle['x_grid'][:,:] - cx)**2 +
                               (handle['y_grid'][:,:] - cy)**2 < radius**2 )
    F = np.zeros((Nx, Ny))
    F[bubble_indices] = 1
    
    #Potential temperature
    handle['theta_t'][:,:,0] = F
    handle['theta_t'][:,:,1] = F
    
    #Initial displacement guess for advection
    handle['alpha_u'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_v'][:,:,0] = np.zeros((Nx,Ny))
    
    #Uniform wind
    u = wind_norm/np.sqrt(handle.dx**2+handle.dy**2)
    handle['ut'][:,:,0] = u * np.ones((Nx,Ny))
    handle['vt'][:,:,0] = u * np.ones((Nx,Ny))
    
    #time storage
    handle['t'][0] = 0
    handle['t'][1] = handle.dt
#------------------------------------------------------------------------------
    
    handle.close()
    
#------------------------------------------------------------------------------
   