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


def filament_test(path, Lx, Ly, Nx, Ny, T, Nt):
    
    handle = create_netcdf(path, Lx, Ly,T,  Nx, Ny, Nt)
#ATTRIBUTES -------------------------------------------------------------------

    handle.T = T
    handle.dt = T/Nt
    handle.Nt = Nt
#------------------------------------------------------------------------------

#INITIALIZATION ---------------------------------------------------------------
        
    #Filament creation
    thetatp = np.zeros((Nx,Ny))
    import scipy.ndimage as spnd 
    Px = 8 
    Py = 64
    thetaanom = 15
    for i in range(Nx): 
        for j in range(Ny):
            if abs(i-Nx/2) < Px and abs(j-Ny/2) < Py:
                thetatp[i,j] = thetaanom
    spnd.gaussian_filter(thetatp, 10, output=thetatp)
    
    #Potential temperature
    handle['theta_t'][:,:,0] = thetatp
    handle['theta_t'][:,:,1] = thetatp
    
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
   