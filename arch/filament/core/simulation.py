import numpy as np
from netCDF4 import Dataset
import os
import time

from .history import History
from .grid import Grid
from ..test.test_cases import create_results_netcdf

#-------------------------------------------------------------------------------------------

class Simulation():

    def __init__(self, initialCDF, T, Nt, methods, output_folder, save_rate, backup_rate, verbose=0, methods_kwargs=None):
       
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
        self.methods_kwargs = methods_kwargs if methods_kwargs is not None else len(methods)*[None]

        # handling the output netCDF files (save & backup)
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            pass
        self.output_folder = output_folder
        self.save_rate = save_rate
        self.backup_rate = backup_rate
        create_results_netcdf(output_folder + '/results.nc', **self.__dict__)
        create_results_netcdf(output_folder + '/backup.nc', **self.__dict__) 

        # other parameters
        self.verbose = verbose


    def run(self):
        cpu_tot_time = np.zeros(len(self.methods))
        if self.verbose:
            print("          ------------------------")
            print("          |  RUNNING SIMULATION  |")
            print("          ------------------------")
        for iter_nb in range(self.Nt):
            print("\n\nIteration ", iter_nb, "...") if self.verbose else None
            # first handle saving
            if iter_nb % self.backup_rate == 0:
                backupCDF = Dataset(self.output_folder + '/backup.nc', 'r+', format='NETCDF4', parallel=False)
                self.history.save(backupCDF, backup=True)
                backupCDF.close()
                print("---> backup refreshed at iteration "+str(iter_nb)) if self.verbose else None
            if iter_nb % self.save_rate == 0:
                resultsCDF = Dataset(self.output_folder + '/results.nc', 'r+', format='NETCDF4', parallel=False)
                self.history.save(resultsCDF, backup=False)  
                resultsCDF.close()  
                print("---> saved results of iteration "+str(iter_nb)) if self.verbose else None

            # then perform forward
            cpu_time = self.forward()
            cpu_tot_time += cpu_time
            
            simu_time = 0
        
        # Print Total and Mean CPU time per method
        for ind, method in enumerate(self.methods):
            print("\n\nTotal CPU time for method ", method.__name__, " = {:.2f}".format(cpu_tot_time[ind]), " seconds") if self.verbose else None
            print("Mean CPU time for method ", method.__name__, " per call = {:.2f}".format(cpu_tot_time[ind]/self.Nt), " seconds") if self.verbose else None
            simu_time += cpu_tot_time[ind] + cpu_tot_time[ind]/self.Nt
            
        print("\n**************************************************\n")
        print("TOTAL SIMULATION TIME = {:.2f}".format(simu_time), " seconds")

    def forward(self):
        cpu_time = np.zeros(len(self.methods))
        for ind, (method, kwargs) in enumerate(zip(self.methods, self.methods_kwargs)):
            t0 = time.time()
            print("      *** Proceeding to method: "+method.__name__) if self.verbose > 1 else None
            method(**self.__dict__, **kwargs)
            cpu_time[ind] = time.time() - t0 
            print("      *** CPU time = {:.2f}".format(cpu_time[ind]), " seconds") if self.verbose > 1 else None
        return cpu_time
