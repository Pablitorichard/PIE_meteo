import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from netCDF4 import Dataset
from IPython.display import HTML
from base64 import b64encode
import os

def make_video(pathCDF, save_path, variable):
    try:
        os.mkdir(save_path+'/videos')
    except FileExistsError:
        pass
    
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)

    fig = plt.figure(figsize=(12,8))
    frames = [[plt.imshow(resultsCDF[variable][:,:,iteration_nb].T, 
                          origin='lower', 
                          cmap='magma')]
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