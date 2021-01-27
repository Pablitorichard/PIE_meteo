import numpy as np
from netCDF4 import Dataset

#------------------------------------------------------------------------------

def create_results_netcdf(path, grid, T, params, **kwargs):

    handle = Dataset(path, 'w', format='NETCDF4', parallel=False)

    handle.createDimension("Nx", grid.Nx)
    handle.createDimension("Ny", grid.Ny)
    handle.createDimension("Nt", None) 
	
	handle.T = T
    list(map(lambda item: handle.setncattr(*item), params.items())) 

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
    
    
    
    handle['x_grid'][:,:] = grid.x_grid
    handle['y_grid'][:,:] = grid.y_grid
    
    handle.close()
    
#------------------------------------------------------------------------------

def create_initial_netcdf(path, Lx, Ly, Nx, Ny, dt, nb_state):

#CREATION OF THE NETCDF FILE --------------------------------------------------
    handle = Dataset(path, 'w',format='NETCDF4')

#DIMENSIONS -------------------------------------------------------------------
    handle.createDimension("Nx", Nx)
    handle.createDimension("Ny", Ny)
    handle.createDimension("Nt", nb_state) 
    handle.createDimension("one",1)

#ATTRIBUTE --------------------------------------------------------------------   
    # Geometry
    handle.Lx = Lx
    handle.Ly = Ly
    handle.Nx = Nx
    handle.Ny = Ny
    handle.dx = Lx / Nx
    handle.dy = Ly / Ny 
    
    ## Parameters
    handle.z_star = -500
    handle.gamma_1 = -4E-3 #K.m-1
    handle.gamma_2 = -8.5E-3 #K.m-1
    handle.Delta_zc = 500 
    handle.Delta_Tc = -5
    handle.g = 9.81
    handle.N_t = 0.01 # Brunt-Vaisala frequency of the troposphere (s^{-1})
    handle.N_s = 2E-2 # Brunt-Vaisala frequency of the stratosphere (s^{-1})
    handle.theta_00 = 300

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

#GEOMETRY INITIALIZATION ------------------------------------------------------
    grid = np.mgrid[0:handle.Lx:handle.dx, 0:handle.Ly:handle.dy]
    handle['x_grid'][:,:] = grid[0,:,:]
    handle['y_grid'][:,:] = grid[1,:,:]
    
#TIME INITIALIZATION-----------------------------------------------------------
	handle['t'][:] = dt * np.arange(nb_state)

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
