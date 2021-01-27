import numpy as np
from netCDF4 import Dataset
import os

from .history import History
from .grid import Grid
from ..test.test_cases import create_results_netcdf

#-------------------------------------------------------------------------------------------

class Simulation():

    def __init__(self, initialCDF, T, Nt, methods, output_folder, save_rate, backup_rate, verbose=True, methods_kwargs=None):
       
        # store the initial data
        self.history = History.fromCDF(initialCDF)
        self.params = initialCDF.__dict__
        self.grid = Grid(**self.params)
        
        tinit = initialCDF['t'][:].data
        tnext = np.array((Nt-len(tinit)) * [None])
        self.timeline = np.concatenate((tinit,tnext), axis=0)
        self.T = T
        self.Nt = Nt
        
        initialCDF.close()
        
        # checkout the methods
        self.methods = methods
        self.methods_kwargs = methods_kwargs

        # handling the output netCDF files (save & backup)
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            pass
        self.output_folder = output_folder
        self.save_rate = save_rate
        self.backup_rate = backup_rate
        create_results_netcdf(output_folder + '/results.nc', **self)
        create_results_netcdf(output_folder + '/backup.nc', **self) 

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
