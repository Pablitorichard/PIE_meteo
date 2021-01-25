import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from matplotlib.ticker import MaxNLocator

def ugly_plot(path, lvl_num=5, colormap='Greys', ratio=8):
    
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
    fig, (( ax_t0, ax_t1, ax_t2), (ax_t3, ax_t4,ax_t5)) =  \
        plt.subplots(nrows=2, ncols =3)
    
    # Potential temperature plot
    F=handle['theta_t'][:,:,1]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    im_t0=ax_t0.contourf(F.T,cmap=colormap,levels=lvl)
    t0_label = 't = '+ str(round(handle['t'][0]/3600,3)) + " h"
    ax_t0.set_xlabel(t0_label)
    
    F=handle['theta_t'][:,:,Nt//2]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    ax_t1.contourf(F.T,cmap=colormap,levels=lvl)
    t1_label = 't = '+ str(round(handle['t'][Nt//2]/3600,3)) + " h"
    ax_t1.set_xlabel(t1_label)
   
    F=handle['theta_t'][:,:,Nt-1]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    ax_t2.contourf(F.T,cmap=colormap,levels=lvl)
    t2_label = 't = '+ str(round(handle['t'][Nt]/3600,3)) + " h"
    ax_t2.set_xlabel(t2_label)
    
    #Velocity plots
    Xr = X[0:Nx:Nx//ratio]
    Yr = Y[0:Ny:Ny//ratio]
    
    #Fetch velocity components and their norm
    U = handle['ut'][:,:,1]
    V = handle['vt'][:,:,1]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t3 = ax_t3.pcolormesh(x_grid, y_grid, N)
    ax_t3.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    
    #Fetch velocity components and their norm
    U = handle['ut'][:,:,Nt//2]
    V = handle['vt'][:,:,Nt//2]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t4 = ax_t4.pcolormesh(x_grid, y_grid, N)
    ax_t4.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    
    #Fetch velocity components and their norm
    U = handle['ut'][:,:,Nt-1]
    V = handle['vt'][:,:,Nt-1]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t5 = ax_t5.pcolormesh(x_grid, y_grid, N)
    ax_t5.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    # Layout settings and color bar
    plt.tight_layout()
    fig.subplots_adjust(right=0.85)
    cbar1_ax = fig.add_axes([0.9, 0.55, 0.04, 0.4])
    fig.colorbar(im_t0, cax=cbar1_ax)
    cbar2_ax = fig.add_axes([0.9, 0.05, 0.04, 0.4])
    fig.colorbar(im_t5, cax=cbar2_ax)
    handle.close()
    plt.show()
    


def ugly_WV(path, lvl_num=5, colormap='Greys', ratio=8):
    
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
    fig, (( ax_t0, ax_t1, ax_t2), (ax_t3, ax_t4,ax_t5)) =  \
        plt.subplots(nrows=2, ncols =3)
    
    # Potential temperature plot
    F=handle['Delta_T_bb'][:,:,0]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    im_t0=ax_t0.contourf(F.T,cmap=colormap,levels=lvl)
    t0_label = 't = '+ str(round(handle['t'][0]/3600,3)) + " h"
    ax_t0.set_xlabel(t0_label)
    
    F=handle['Delta_T_bb'][:,:,Nt//2]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    ax_t1.contourf(F.T,cmap=colormap,levels=lvl)
    t1_label = 't = '+ str(round(handle['t'][Nt//2]/3600,3)) + " h"
    ax_t1.set_xlabel(t1_label)
   
    F=handle['Delta_T_bb'][:,:,Nt-1]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    ax_t2.contourf(F.T,cmap=colormap,levels=lvl)
    t2_label = 't = '+ str(round(handle['t'][Nt]/3600,3)) + " h"
    ax_t2.set_xlabel(t2_label)
    
    #Velocity plots
    Xr = X[0:Nx:Nx//ratio]
    Yr = Y[0:Ny:Ny//ratio]
    
    #Fetch velocity components and their norm
    U = handle['us'][:,:,1]
    V = handle['vs'][:,:,1]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t3 = ax_t3.pcolormesh(x_grid, y_grid, N)
    ax_t3.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    
    #Fetch velocity components and their norm
    U = handle['us'][:,:,Nt//2]
    V = handle['vs'][:,:,Nt//2]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t4 = ax_t4.pcolormesh(x_grid, y_grid, N)
    ax_t4.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    
    #Fetch velocity components and their norm
    U = handle['us'][:,:,Nt-1]
    V = handle['vs'][:,:,Nt-1]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t5 = ax_t5.pcolormesh(x_grid, y_grid, N)
    ax_t5.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    # Layout settings and color bar
    plt.tight_layout()
    fig.subplots_adjust(right=0.85)
    cbar1_ax = fig.add_axes([0.9, 0.55, 0.04, 0.4])
    fig.colorbar(im_t0, cax=cbar1_ax)
    cbar2_ax = fig.add_axes([0.9, 0.05, 0.04, 0.4])
    fig.colorbar(im_t5, cax=cbar2_ax)
    handle.close()
    plt.show()
    
def ugly_Delta_z(path, lvl_num=10, colormap='Greys', ratio=8):
    
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
    fig, (( ax_t0, ax_t1, ax_t2), (ax_t3, ax_t4,ax_t5)) =  \
        plt.subplots(nrows=2, ncols =3)
    
    # Potential temperature plot
    F=handle['Delta_z'][:,:,0]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    im_t0=ax_t0.contourf(F.T,cmap=colormap,levels=lvl)
    t0_label = 't = '+ str(round(handle['t'][0]/3600,3)) + " h"
    ax_t0.set_xlabel(t0_label)
    
    F=handle['Delta_z'][:,:,Nt//2]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    ax_t1.contourf(F.T,cmap=colormap,levels=lvl)
    t1_label = 't = '+ str(round(handle['t'][Nt//2]/3600,3)) + " h"
    ax_t1.set_xlabel(t1_label)
   
    F=handle['Delta_z'][:,:,Nt-1]
    lvl=MaxNLocator(nbins=lvl_num).tick_values(F.min(), F.max())       
    im_t2=ax_t2.contourf(F.T,cmap=colormap,levels=lvl)
    t2_label = 't = '+ str(round(handle['t'][Nt]/3600,3)) + " h"
    ax_t2.set_xlabel(t2_label)
    
    #Velocity plots
    Xr = X[0:Nx:Nx//ratio]
    Yr = Y[0:Ny:Ny//ratio]
    
    #Fetch velocity components and their norm
    U = handle['us'][:,:,1]
    V = handle['vs'][:,:,1]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t3 = ax_t3.pcolormesh(x_grid, y_grid, N)
    ax_t3.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    
    #Fetch velocity components and their norm
    U = handle['us'][:,:,Nt//2]
    V = handle['vs'][:,:,Nt//2]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t4 = ax_t4.pcolormesh(x_grid, y_grid, N)
    ax_t4.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    
    #Fetch velocity components and their norm
    U = handle['us'][:,:,Nt-1]
    V = handle['vs'][:,:,Nt-1]
    N = np.sqrt(U**2 + V**2)
    #Take a subset and transpose for quiver.
    Ur =U[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Vr = V[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    Nr = N[0:Nx:Nx//ratio,0:Ny:Ny//ratio].T
    # Normalize the arrows to just get a readable direction
    Ur = np.divide(Ur,Nr)
    Vr = np.divide(Vr,Nr)
    #Plot
    im_t5 = ax_t5.pcolormesh(x_grid, y_grid, N)
    ax_t5.quiver(Xr, Yr, Ur,Vr, width=0.002, color='w')
    
    # Layout settings and color bar
    plt.tight_layout()
    fig.subplots_adjust(right=0.85)
    cbar1_ax = fig.add_axes([0.9, 0.55, 0.04, 0.4])
    fig.colorbar(im_t2, cax=cbar1_ax)
    cbar2_ax = fig.add_axes([0.9, 0.05, 0.04, 0.4])
    fig.colorbar(im_t5, cax=cbar2_ax)
    handle.close()
    plt.show()