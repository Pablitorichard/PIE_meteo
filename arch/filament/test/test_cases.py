import numpy as np
from netCDF4 import Dataset

#------------------------------------------------------------------------------

def create_results_netcdf(path, grid, T, Nt):

    handle = Dataset(path, 'w', format='NETCDF4', parallel=False)

    handle.createDimension("Nx", grid.Nx)
    handle.createDimension("Ny", grid.Ny)
    handle.createDimension("Nt", None) 

    handle.T = T    
    handle.Lx = grid.Lx
    handle.Ly = grid.Ly
    handle.Nx = grid.Nx
    handle.Ny = grid.Ny
    handle.Nt = Nt
    handle.dx = grid.dx
    handle.dy = grid.dy
    handle.dt = T/Nt

    # "f8" is a data type: 64-bit floating point variable
    handle.createVariable("ut","f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("us","f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vs", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("theta_t", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("theta_s", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("delta_z", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("Delta_z", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("w", "f8",("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_ut", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_us", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_vs", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("t", "f8", ("Nt"))
    
    handle.createVariable("x_grid", "f8", ("Nx", "Ny"))
    handle.createVariable("y_grid", "f8", ("Nx", "Ny"))
    handle['x_grid'][:,:] = grid.x_grid
    handle['y_grid'][:,:] = grid.y_grid

    return handle

#------------------------------------------------------------------------------

def bubble_test(path, grid, T, Nt, cx, cy, radius, wind_norm):
    
    handle = create_results_netcdf(path, grid, T, Nt)
        
    #Bubble creation
    bubble_indices = np.where ((handle['x_grid'][:,:] - cx)**2 +
                               (handle['y_grid'][:,:] - cy)**2 < radius**2 )
    
    F = np.zeros((grid.Nx, grid.Ny))
    F[bubble_indices] = 1
    
    #Potential temperature
    handle['theta_t'][:,:,0] = F
    handle['theta_t'][:,:,1] = F
    
    #Initial displacement guess for advection
    handle['alpha_ut'][:,:,0] = np.zeros((grid.Nx,grid.Ny))
    handle['alpha_vt'][:,:,0] = np.zeros((grid.Nx,grid.Ny))
    
    #Uniform wind
    u = wind_norm/np.sqrt(handle.dx**2+handle.dy**2)
    handle['ut'][:,:,0] = u * np.ones((grid.Nx,grid.Ny))
    handle['vt'][:,:,0] = u * np.ones((grid.Nx,grid.Ny))
    
    #time storage
    handle['t'][0] = 0
    handle['t'][1] = handle.dt

    return handle

#------------------------------------------------------------------------------

def gaussian_test(path, grid, T, Nt):
    
    handle = create_results_netcdf(path, grid, T, Nt)
     
    #Gaussian creation
    thetatp = np.zeros((grid.Nx, grid.Ny))
    import scipy.ndimage as spnd 
    Px = 8 
    Py = 64
    thetaanom = 15
    for i in range(grid.Nx): 
        for j in range(grid.Ny):
            if abs(i-grid.Nx/2) < Py and abs(j-grid.Ny/2) < Px:
                thetatp[i,j] = thetaanom
    spnd.gaussian_filter(thetatp, 10, output=thetatp)
    
    #Potential temperature
    handle['theta_t'][:,:,0] = thetatp
    handle['theta_t'][:,:,1] = thetatp
    
    #Initial displacement guess for advection
    handle['alpha_ut'][:,:,0] = np.zeros((grid.Nx, grid.Ny))
    handle['alpha_vt'][:,:,0] = np.zeros((grid.Nx, grid.Ny))
    
    #Uniform wind
    handle['ut'][:,:,0] = np.zeros((grid.Nx, grid.Ny))
    handle['vt'][:,:,0] = np.zeros((grid.Nx, grid.Ny))
    
    #time storage
    handle['t'][0] = 0
    handle['t'][1] = handle.dt

    return handle
    
#------------------------------------------------------------------------------