import numpy as np
from copy import deepcopy
        
forced_variables = ['x_grid','y_grid','t'] # General variables     

class State():
    
    def __init__(self, t, vrs={}):
        self.t = t
        self.vrs = deepcopy(vrs)

    @classmethod
    def fromCDF(cls, netCDF_file, k=None):
        t = np.float(netCDF_file['t'][k].data) if k is not None else None

        variables = [var for var in netCDF_file.variables if var not in forced_variables]
        vrs = {var: netCDF_file[var][:,:,k].data if k is not None else netCDF_file[var][:].data
               for var in variables}

        return cls(t, vrs)

    @classmethod
    def copy(cls, otherState):
        return cls(otherState.t, otherState.vrs)

    def save(self, netCDF_file, saved_vrs=None, backup=False, k=None):
        
        # current location to fill
        k = netCDF_file.dimensions['Nt'].size if not backup else k

        # filling the variables
        netCDF_file['t'][k] = self.t

        variables = [var for var in netCDF_file.variables if var not in forced_variables]
        export_variables = saved_vrs if saved_vrs is not None else variables
        for var in export_variables:
            netCDF_file[var][:,:,k] = self.vrs[var]
