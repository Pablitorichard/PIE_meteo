import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
from netCDF4 import Dataset

## Interactive plot (specific to IPython)
import ipywidgets as widgets

def plot_data(pathCDF, variable, time, cmap='magma'):
    
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)
    
    min_value = np.min(resultsCDF[variable][:])
    max_value = np.max(resultsCDF[variable][:])

    if(len(np.shape(resultsCDF[variable])) == 3):
        fig = plt.figure(figsize=(12,8))
        ax = plt.subplot(111)
        ax.set_title(variable, fontsize=25, pad=15)
        cbar = fig.colorbar(ScalarMappable(cmap=cmap, norm=Normalize(vmin=min_value, vmax=max_value)), ax=ax)
        cbar_ticks = cbar.get_ticks()
        cbar.set_ticks(cbar_ticks)
        plt.imshow(resultsCDF[variable][:,:,time].T,
                    origin='lower', 
                            cmap=cmap,
                            vmin=min_value,
                            vmax=max_value) 
    else :
        print(variable + ' is not a variable evolving on the 2D-grid through time.')
            
    resultsCDF.close()
    
def interactive_plot(pathCDF):
	
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)

    times = resultsCDF['t'][:].data
    times_ind = np.arange(len(times))

    disp_times = [ '{}'.format(int(np.floor(times[k]//3600))) + 'h {}min'.format(int(np.floor(times[k]%3600//60))) for k in times_ind]
    
    vars_opts = list(resultsCDF.variables)
    cmap_opts = ['magma','Greys','hot','viridis','plasma','inferno','cividis']
    time_opts = [(disp_times[k], k) for k in times_ind]
    
    resultsCDF.close()

    vars_widg = widgets.ToggleButtons(options=vars_opts, 
                                    value=vars_opts[0], 
                                    description='variable :',
                                    disabled=False)
    
    time_widg = widgets.SelectionSlider(options=time_opts, 
    									value=0, 
    									description='t = ',
    									continuous_update=False, 
    									disabled=False)

    cmap_widg = widgets.Select(options=cmap_opts,
                                value='magma',
                                description='colormap :',
                                disabled=False)
    
    inter = widgets.interactive(plot_data, 
                                pathCDF=widgets.fixed(pathCDF),
                                cmap=cmap_widg,
                                variable=vars_widg,
                                time=time_widg)
    return inter