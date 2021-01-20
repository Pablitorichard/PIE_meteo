# -*- coding: utf-8 -*-
"""
dirver_netcdf(path, Nx, Ny)

init_cond takes a relative path as an input, and creates a netcdf file at
 this spot. The initial condition is a 'bubble'
"""


#LIBRARY ----------------------------------------------------------------------

# My cocktail to make path stuff work in Python
from pathlib import Path
import os
import inspect
# Path to 
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))

#Other libraries
import numpy as np
from netCDF4 import Dataset

#------------------------------------------------------------------------------

def create_netcdf(path, Lx, Ly, T, Nx, Ny, Nt):

#CREATION OF THE NETCDF FILE --------------------------------------------------

    file_path= self_path / path
    handle = Dataset(file_path, 'w',format='NETCDF4')

#------------------------------------------------------------------------------

#DIMENSIONS -------------------------------------------------------------------

    handle.createDimension("Nx", Nx)
    handle.createDimension("Ny", Ny)
    handle.createDimension("Nt", None) 
    handle.createDimension("one",1)

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
    handle.z_star = -500
    handle.gamma_1 = -4
    handle.gamma_2 = -8.5
    handle.Delta_zc = 500
    handle.Delta_Tc = -5
    handle.g = 9.81
    handle.Nt = 0.01
    handle.Ns = 2E-2
    handle.theta_00 = 300
    
#------------------------------------------------------------------------------

#VARIABLES --------------------------------------------------------------------
    # "f8" is a data type: 64-bit floating point variable
    handle.createVariable("ut","f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("us","f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vs", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("w", "f8", ("Nx", "Ny", "Nt"))
    

    handle.createVariable("theta_t", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("Delta_T_bb", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("Delta_z", "f8",("Nx", "Ny", "Nt"))
    
    handle.createVariable("alpha_ut", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_us", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_vs", "f8", ("Nx", "Ny", "Nt"))
    
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



