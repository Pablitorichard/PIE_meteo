import numpy as np

class Grid():
    
    def __init__(self, Lx, Ly, Nx, Ny):
        
        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
        
        self.dx = self.Lx / self.Nx
        self.dy = self.Ly / self.Ny     
        
        grid = np.mgrid[0:self.Lx:self.dx, 0:self.Ly:self.dy]
        self.x_grid = grid[0,:,:]
        self.y_grid = grid[1,:,:]