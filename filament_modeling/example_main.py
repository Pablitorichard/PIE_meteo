
#Path management
from pathlib import Path
import os
import inspect
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))


#Custom functions
from test_cases import bubble_test, gaussian_test, v_stripe_test
from advection_driver import advection_driver
from ugly_plot import ugly_plot, ugly_WV, ugly_Delta_z

# Pick an adress for the netCDF file (it will over write existing files)
path = self_path / "outputs/out_dt300_WV.nc" 

# Parameters
Lx = 2048E3
Ly = 1024E3
Nx = 256
Ny = 128

# Nt will drive the computation time. On my standard laptop (when plugged in), 
# I measure around 500 time cycles per hour for linear interpolation.
dt = 300
T = 48*3600
Nt = int(T//dt)



# Pick a test case

# Half width of the stripe
dY = Ny//15
# Distance between the side of the X axis and the first V-profile
dX = Nx//8
v_stripe_test(path, Lx, Ly, Nx, Ny, T, Nt, dX, dY)




# The advection driver will propagate the solution in time
advection_driver(path, pseudo_spectral_wind = 1,
                alpha_method = 'linear', F_method='bicubic')

# the name is accurate, it is a plot, it is ugly
ugly_plot(path, ratio=20, lvl_num=50)
ugly_WV(path, ratio=20, lvl_num=50)
ugly_Delta_z(path, ratio=20, lvl_num=50,colormap='bone')