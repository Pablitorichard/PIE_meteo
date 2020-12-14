
#Path management
from pathlib import Path
import os
import inspect
self_path = Path(os.path.dirname(inspect.getfile(lambda: None)))


#Custom functions
from filament_test import filament_test
from advection_driver import advection_driver
from ugly_plot import ugly_plot

# Parameters
Lx = 2048
Ly = 1024
Nx = 256
Ny = 128
dt = 30

T = 1*300
#T =48*3600
Nt = T//dt

ratio_x = Nx//10 # 1 arrow every <ratio> point
ratio_y = Ny//10

# Pick an adress for the netCDF file (it will over write existing files)
path = self_path / "out.nc" 

# Initialize the bubbble test case
filament_test("out.nc", Lx, Ly, Nx, Ny, T, Nt)

# The advection driver will propagate the solution in time
advection_driver(path, 1)

# the name is accurate, it is a plot, it is ugly
ugly_plot(path, ratio_x,ratio_y)