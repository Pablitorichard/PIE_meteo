import numpy as np
        
variables = ['delta_z', 'theta_t', 'ut', 'vt', 'alpha_ut', 'alpha_vt', 'w', 
             'Delta_z', 'theta_s', 'us', 'vs', 'alpha_ut', 'alpha_vt']

saved_variables = ['theta_t', 'ut', 'vt', 'Delta_z', 'theta_s']
        
class State():
    
    def __init__(self, t, vrs={}):
        self.t = t
        self.vrs = vrs

    @classmethod
    def fromCDF(cls, netCDF_file, k=None):
        t = netCDF_file['t'][k].data if k is not None else None
        vrs = {variable: netCDF_file[variable][:,:,k].data if k is not None else netCDF_file[variable][:].data
                for variable in variables}
        return cls(t, vrs)

    @classmethod
    def copy(cls, otherState):
        return cls(otherState.t, otherState.vrs)

    def save(self, netCDF_file, backup=False):
        
        # current location to fill
        k = netCDF_file.dimensions['Nt'].size

        # filling the variables
        netCDF_file['t'][k] = self.t

        export_variables = variables if backup else saved_variables
        for variable in export_variables:
            netCDF_file[variable][:,:,k] = self.vrs[variable]
