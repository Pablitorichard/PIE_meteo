import matplotlib.pyplot as plt
import matplotlib.cm as cm
from netCDF4 import Dataset

def ugly_plot(path, ratio_x, ratio_y):
    # UNPACK
    #UNPACK DATA ------------------------------------------------
    handle = Dataset(path, 'r+',format='NETCDF4')  
    Nt = handle.Nt
    dt = handle.dt
    Nx = handle.Nx
    Ny = handle.Ny
    x_grid = handle['x_grid'][:,:]
    y_grid = handle['y_grid'][:,:]
    X_coord = x_grid[:,0]
    Y_coord = y_grid[0,:]

    # DISPLAY
    fig, (( ax_t0, ax_t1,), (ax_t2, ax_t3)) =  plt.subplots(nrows=2, ncols =2,
                                                            figsize=[6.8/0.8,6.8])

    ax_t0.pcolormesh(x_grid, y_grid,handle['theta_t'][:,:,1],
                             cmap=cm.gist_yarg)
    ax_t0.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, 0],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, 0],
                  color='b', width=0.003)
    t0_label = 't = '+ str(handle['t'][1])
    ax_t0.set_xlabel(t0_label)
    

    ax_t1.pcolormesh(x_grid, y_grid, handle['theta_t'][:,:,Nt//3],
                     cmap=cm.gist_yarg)
    ax_t1.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt//3-1],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt//3-1],
                  color='b', width=0.003)
    t1_label = 't = '+ str(handle['t'][Nt//3])
    ax_t1.set_xlabel(t1_label)
    
    

    ax_t2.pcolormesh(x_grid, y_grid,handle['theta_t'][:,:,2*Nt//3],
                     cmap=cm.gist_yarg)
    ax_t2.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, 2*Nt//3-1],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, 2*Nt//3-1],
                  color='b', width=0.003)
    t2_label = 't = '+ str(handle['t'][2*Nt//3])
    ax_t2.set_xlabel(t2_label)
    
    

    im_t3=ax_t3.pcolormesh(x_grid, y_grid ,handle['theta_t'][:,:,Nt],
                     cmap=cm.gist_yarg)
    ax_t3.quiver(X_coord[0:Nx:ratio_x], Y_coord[0:Ny:ratio_y],
                  handle['ut'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt-1],
                  handle['vt'][0:Nx:ratio_x, 0:Ny:ratio_y, Nt-1],
                  color='b', width=0.003)
    t3_label = 't = '+ str(handle['t'][Nt-1 ])
    ax_t3.set_xlabel(t3_label)
    
    plt.tight_layout()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.10, 0.05, 0.7])
    fig.colorbar(im_t3, cax=cbar_ax)
    handle.close()
    plt.show()