import numpy as np
import time
from netCDF4 import Dataset

from pathlib import Path
import os
import inspect
# Path to 
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))

def create_netcdf(path, Lx, Ly, T, Nx, Ny, Nt):

#CREATION OF THE NETCDF FILE --------------------------------------------------

    file_path = self_path / path
    handle = Dataset(file_path, 'w',format='NETCDF4')

#------------------------------------------------------------------------------

#DIMENSIONS -------------------------------------------------------------------

    handle.createDimension("Nx", Nx)
    handle.createDimension("Ny", Ny)
    handle.createDimension("Nt", None) 

#------------------------------------------------------------------------------

#ATTRIBUTE --------------------------------------------------------------------
    handle.T = T    
    handle.Lx = Lx
    handle.Ly = Ly
    handle.Nx = Nx
    handle.Ny = Ny
    handle.Nt = Nt
    handle.dx = Lx / Nx
    handle.dy = Ly/ Ny  
    handle.dt = T/Nt
#------------------------------------------------------------------------------

#VARIABLES --------------------------------------------------------------------
    # "f8" is a data type: 64-bit floating point variable
    handle.createVariable("ut","f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("theta_t", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_u", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_v", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("t", "f8", ("Nt"))
    handle.createVariable("x_grid", "f8", ("Nx", "Ny"))
    handle.createVariable("y_grid", "f8", ("Nx", "Ny"))
#------------------------------------------------------------------------------

#GEOMETRY INITIALIZATION ------------------------------------------------------
    
    grid = np.mgrid[0:handle.Lx:handle.dx, 0:handle.Ly:handle.dy]
    handle['x_grid'][:,:] = grid[0,:,:]
    handle['y_grid'][:,:] = grid[1,:,:]

#------------------------------------------------------------------------------
    return(handle)



def bubble_test(path, Lx, Ly, Nx, Ny, T, Nt, cx, cy, radius, wind_norm):
    
    handle = create_netcdf(path, Lx, Ly,T,  Nx, Ny, Nt)
    print("netcdf created")
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