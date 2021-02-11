import numpy as np
import copy
        
variables = ['ut', 'vt', 'alpha_ut', 'alpha_vt', 'theta_t', 
             'us', 'vs', 'alpha_us', 'alpha_vs', 
             'w', 'Delta_z', 'Delta_T_bb', 'Delta_T_hist']

saved_variables = ['theta_t', 'ut', 'vt', 'Delta_z', 'Delta_T_bb', 'us', 'vs', 'w'] ##### to change
        
class State():
    
    def __init__(self, t, vrs={}):
        self.t = t
        self.vrs = copy.deepcopy(vrs)

    @classmethod
    def fromCDF(cls, netCDF_file, k=None):
        t = np.float(netCDF_file['t'][k].data) if k is not None else None
        vrs = {variable: netCDF_file[variable][:,:,k].data if k is not None else netCDF_file[variable][:].data
               for variable in variables}
        return cls(t, vrs)

    @classmethod
    def copy(cls, otherState):
        return cls(otherState.t, otherState.vrs)

    def save(self, netCDF_file, backup=False, k=None):
        
        # current location to fill
        k = netCDF_file.dimensions['Nt'].size if not backup else k

        # filling the variables
        netCDF_file['t'][k] = self.t

        export_variables = variables if backup else saved_variables
        for variable in export_variables:
            netCDF_file[variable][:,:,k] = self.vrs[variable]
