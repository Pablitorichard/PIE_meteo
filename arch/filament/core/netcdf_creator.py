import numpy as np
from netCDF4 import Dataset
from copy import deepcopy

from .state import forced_variables
#------------------------------------------------------------------------------

def create_results_netcdf(path, initialCDF, params, grid, T, Nt, 
                            methods, methods_kwargs, save_rate, backup_rate, **kwargs):

    handle = Dataset(path, 'w', format='NETCDF4', parallel=False)

    handle.createDimension("Nx", grid.Nx)
    handle.createDimension("Ny", grid.Ny)
    handle.createDimension("Nt", None) 

    handle.T = deepcopy(T)
    list(map(lambda item: handle.setncattr(*item), params.items())) 
    handle.Nt = deepcopy(Nt)

    
    handle.methods = [m.__name__ for m in methods]
    # Probleme a la sauvegarde des arguments en temps qu'attributs a regler
    #meth_kwar = []
    #for i,mk in enumerate(methods_kwargs):
    #    meth_kwar.append(i)
    #    meth_kwar.append(list(mk))
    #handle.methods_kwargs = meth_kwar

    handle.save_rate = deepcopy(save_rate)
    handle.backup_rate = deepcopy(backup_rate) 

    # "f8" is a data type: 64-bit floating point variable
    for var in initialCDF.variables:
        if (var not in forced_variables):
            handle.createVariable(var, "f8", ("Nx", "Ny", "Nt"))

    handle.createVariable("t", "f8", ("Nt"))
    handle.createVariable("x_grid", "f8", ("Nx", "Ny"))
    handle.createVariable("y_grid", "f8", ("Nx", "Ny"))
    
    handle['x_grid'][:,:] = grid.x_grid
    handle['y_grid'][:,:] = grid.y_grid
    
    handle.close()

def results_netcdf_frombackup(path, backupCDF, pre_resultCDF, **kwargs):

    handle = Dataset(path, 'r+', format='NETCDF4', parallel=False)

    t_tocopy = np.where(pre_resultCDF['t'][:].data < backupCDF['t'][0])
    kmax = np.argmax(t_tocopy) if t_tocopy else 0
    
    for var in handle.variables:
        if (var not in forced_variables):
            for k in range(kmax):
                handle[var][:,:,k] = pre_resultCDF[var][:,:,k]

    for k in range(kmax):
        handle['t'][k] = pre_resultCDF['t'][k]

    handle.close()


    


    
