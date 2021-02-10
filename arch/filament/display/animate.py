import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
from netCDF4 import Dataset
from IPython.display import HTML
from base64 import b64encode
import os

def make_video(pathCDF, save_path, variable, cmap='magma'):
    try:
        os.mkdir(save_path+'/videos')
    except FileExistsError:
        pass
    
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)
    
    min_value = np.min(resultsCDF[variable][:].data)
    max_value = np.max(resultsCDF[variable][:].data)

    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    ax.set_title(variable, fontsize=25, pad=15)
    cbar = fig.colorbar(ScalarMappable(cmap=cmap, norm=Normalize(vmin=min_value, vmax=max_value)), ax=ax)
    cbar_ticks = cbar.get_ticks()
    cbar.set_ticks(cbar_ticks)
    frames = [[plt.imshow(resultsCDF[variable][:,:,iteration_nb].T, 
                          origin='lower', 
                          cmap=cmap,
                          vmin=min_value,
                          vmax=max_value)]
              for iteration_nb in range(resultsCDF.dimensions['Nt'].size)]
    resultsCDF.close()
    anim = ArtistAnimation(fig, frames, interval=200, blit=False, repeat=False)
    anim.save(save_path+'/videos/'+variable+'.mp4')
    plt.close()
    
    mp4 = open(save_path+'/videos/'+variable+'.mp4','rb').read()
    data_url = "data:"+save_path+"/videos/"+variable+"/mp4;base64," + b64encode(mp4).decode()
    return HTML("""<video width=800 controls>
                         <source src="%s" type="video/mp4">
                   </video>""" % data_url)