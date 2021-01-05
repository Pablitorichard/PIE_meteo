#Path management
from pathlib import Path
import os
import inspect
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))



import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from matplotlib.ticker import MaxNLocator
from celluloid import Camera

def make_video(path, lvl_num=10, colormap='Greys', ratio=25):
    
    #Unpack NetCDF
    handle = Dataset(path, 'r+',format='NETCDF4')  
    Nt = handle.Nt
    dt = handle.dt
    Nx = handle.Nx
    Ny = handle.Ny
    x_grid = handle['x_grid'][:,:]
    y_grid = handle['y_grid'][:,:]
    X = x_grid[:,0]
    Y = y_grid[0,:]
    
    # Create figure
    fig, (( ax_F), (ax_V)) =  \
        plt.subplots(nrows=2, ncols =1, figsize=[4.8, 6.4])
    
    camera = Camera(fig)
    
    for t in range (0,Nt,20):
        # Potential temperature plot
        F=handle['theta_t'][:,:,t]
        lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
        im_F = ax_F.contourf(F.T,cmap=colormap,levels=lvl)
        F_label = 't = '+ str(round(handle['t'][t]/3600,3)) + " h"
        ax_F.set_xlabel(F_label)
        
        
        #Velocity plots
        Xr = X[0:Nx:Nx//ratio]
        Yr = Y[0:Ny:Ny//ratio]
        
        #Fetch velocity components and their norm
        U = handle['ut'][:,:,t]
        V = handle['vt'][:,:,t]
        N = np.sqrt(U**2 + V**2)
        #Take a subset and transpose for quiver.
        Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
        Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
        Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
        # Normalize the arrows to just get a readable direction
        Ur = np.divide(Ur,Nr)
        Vr = np.divide(Vr,Nr)
        #Plot
        im_V = ax_V.pcolormesh(x_grid, y_grid, N)
        ax_V.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
        
        
        # Layout settings and color bar
        # plt.tight_layout()
        # fig.subplots_adjust(right=0.85)
        # cbar1_ax = fig.add_axes([0.9, 0.55, 0.04, 0.4])
        # fig.colorbar(im_F, cax=cbar1_ax)
        # cbar2_ax = fig.add_axes([0.9, 0.05, 0.04, 0.4])
        # fig.colorbar(im_V, cax=cbar2_ax)
        plt.pause(0.001)
        camera.snap()
   
    handle.close()
    animation = camera.animate()
    animation.save(self_path / 'animation.gif',writer='PillowWriter', fps=Nt//10)
    
path = self_path / "out.nc" 
make_video(path)