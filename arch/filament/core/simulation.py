import numpy as np
from netCDF4 import Dataset
import os


from .history import History
from .grid import Grid
from ..test.test_cases import create_results_netcdf

class Simulation():

    def __init__(self, initialCDF, Lx, Ly, Nx, Ny, T, Nt, methods, output_folder, save_rate, backup_rate, verbose):
       
        # store the initial data
        self.history = History.fromCDF(initialCDF)
        initialCDF.close()

        # store the physical dimensions
        self.grid = Grid(Lx, Ly, Nx, Ny)
        self.timeline = np.array([0]+(Nt-1)*[None]) #np.linspace(0, T, Nt, endpoint=False)
        self.T = T
        self.Nt = Nt
        self.dt = T/Nt ############@

        # checkout the methods
        self.methods = methods

        # handling the output netCDF files (save & backup)
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            pass
        self.output_folder = output_folder
        self.save_rate = save_rate
        self.backup_rate = backup_rate
        resultsCDF = create_results_netcdf(output_folder + '/results.nc', self.grid, T, Nt) ##### to change of course
        resultsCDF.close()
        backupCDF = create_results_netcdf(output_folder + '/backup.nc', self.grid, T, Nt) ######
        backupCDF.close()

        # other parameters
        self.verbose = verbose


    def run(self):

        for iter_nb in range(self.Nt):

            # first handle saving
            if iter_nb % self.backup_rate == 0:
                backupCDF = Dataset(self.output_folder + '/backup.nc', 'r+', format='NETCDF4', parallel=False)
                self.history.save(backupCDF, backup=True)
                backupCDF.close()
                print("backup "+str(iter_nb)+" done") if self.verbose else None
            if iter_nb % self.save_rate == 0:
                resultsCDF = Dataset(self.output_folder + '/results.nc', 'r+', format='NETCDF4', parallel=False)
                self.history.save(resultsCDF, backup=False)  
                resultsCDF.close()  
                print("save "+str(iter_nb)+" done") if self.verbose else None

            # then perform forward
            self.forward()
    

    def forward(self):
        for method in self.methods:
            method(**self.__dict__)