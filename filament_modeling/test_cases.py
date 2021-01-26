# Path handling
from pathlib import Path
import os
import inspect
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))

#Other libraries
import numpy as np

#My own functions
from create_netcdf import create_netcdf


def v_stripe_test(path, Lx, Ly, Nx, Ny, T, Nt, dX, dY):
    
    handle = create_netcdf(path, Lx, Ly,T,  Nx, Ny, Nt)
#ATTRIBUTES -------------------------------------------------------------------

    handle.T = T
    handle.dt = T/Nt
    handle.Nt = Nt
#------------------------------------------------------------------------------

#INITIALIZATION ---------------------------------------------------------------
        
    #Bubble creation
    [X,Y] = np.mgrid[0:Nx,0:Ny]
    F = np.zeros((Nx, Ny)) #handle.theta_00*np.ones((Nx, Ny))
    depth = 15#- 2800 * (handle.N_t*handle.N_s*handle.theta_00)/ handle.g

    F[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] = \
        F[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] +\
        (np.abs(Y[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] - Ny//2)/dY -1)*depth
    
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
    
    handle['Delta_T_hist'][:,:,0] = handle.gamma_1 \
        * F * handle.g/(handle.N_t*handle.N_s*handle.theta_00)
    handle['Delta_T_hist'][:,:,1] =  handle['Delta_T_hist'][:,:,0] 
    handle['Delta_T_bb'][:,:,0] = handle['Delta_T_hist'][:,:,0]
    handle['Delta_z'][:,:,0] = np.zeros((Nx,Ny))
    handle['Delta_z'][:,:,1] = np.zeros((Nx,Ny))
    
    #Initial displacement guess for advection
    handle['alpha_ut'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_vt'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_us'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_vs'][:,:,0] = np.zeros((Nx,Ny))
    
    #Uniform wind
    handle['ut'][:,:,0] = np.zeros((Nx,Ny))
    handle['vt'][:,:,0] = np.zeros((Nx,Ny))
    handle['us'][:,:,0] = np.zeros((Nx,Ny))
    handle['vs'][:,:,0] = np.zeros((Nx,Ny))
    handle['w'][:,:,0] = np.zeros((Nx,Ny))
    
    #time storage
    handle['t'][0] = 0
    handle['t'][1] = handle.dt
#------------------------------------------------------------------------------
    
    handle.close()
    
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
    
#-------------------------------------------------------------------------


def gaussian_test(path, Lx, Ly, Nx, Ny, T, Nt):
    
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
            if abs(i-Nx/2) < Py and abs(j-Ny/2) < Px:
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