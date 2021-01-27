"""
 alpha_interp interpolates a multidimensionnal field F from a 2D grid to 
 to an 'upstream' unstructured mesh defined by the displacements 
 alpha_u, alpha_v: If F were a continuous field,we would have:
     F_int(x,y) = F(x-alpha, y-alpha)
 Inputs:
     - alpha_x, double(Nx,Ny): Displacement to apply to the x coordinate
     of each grid point to get the upstream mesh
     
     - alpha_x, double(Nx,Ny): Displacement to apply to the x coordinate
     of each grid point to get the upstream mesh
     
     - F, double(Nx,Ny,dim): Multidimensional field on a 2D grid of 
     dimensions Nx,Ny
     
     - method, string: Interpolation method, set by default to 'linear'. 
     This is the only available option right now.

Output:
    - F_int, double(Nx,Ny,dim) : Value of the field F interpolated on the
     upstream mesh. 
"""


import numpy as np


def upstream_interp(alpha_x, alpha_y, F, method='linear'):
    
    if len(F.shape)==2:
        F = np.array([F])
        
    [dim,Nx,Ny] = F.shape
    F_int = np.zeros((dim,Nx,Ny))
    if method=='nearest':
         # Coordinates of the points of the upstream mesh
         [x_grid, y_grid] = np.mgrid[0:Nx-1, 0:Ny-1]
         x_grid = x_grid -alpha_x
         y_grid = y_grid- alpha_y
         
         # Apply periodic boundary coniditions. Note that nothing is done
         # for x,y<0 as the indices are already periodic in Python. 
         x_grid[np.where(x_grid > Nx-1)] -= (Nx-1)
         y_grid[np.where(y_grid > Ny-1)] -= (Ny-1)
         
         # Fetch the closest neighbor
         F_int[:,:,:] = F[np.round(x_grid), np.round(y_grid),:]
         
    elif method=='linear':
        #We loop on the elements of alpha_u/v
        for x in range(Nx):
            for y in range(Ny):
                #Upstream coordinates
                xi = x - alpha_x[x,y]
                yi = y - alpha_y[x,y]
                
                # Periodic boundary conditions
                #print("xi,yi : ", xi,yi)
                if (xi > 2*(Nx-1) or x< -Nx-1 
                    or yi > 2*(Ny-1) or yi< -Ny-1): 
                    print('Displacement is too large: alpha = [',\
                          alpha_x[x,y],", ", alpha_y[x,y],"] grid points")
                
                if xi<0 :
                    xi += Nx-1 
                elif xi> Nx-1:
                    xi -= (Nx-1)
                if xi<0 :
                    yi += Ny-1
                elif yi> Ny-1:
                    yi -= (Ny-1)
                
                
                #Coordinates of the top right corner, and distance   
                x_top = int(np.ceil(xi))
                y_top = int(np.ceil(yi))
                
                
                x_coeff = x_top - xi
                y_coeff = y_top - yi
                
                #Update of F_int
                F_int[:,x,y] = x_coeff * y_coeff * F[:,x_top-1,y_top-1] \
                    + x_coeff * (1-y_coeff) * F[:,x_top-1,y_top] \
                    + (1-x_coeff) * (1-y_coeff) * F[:,x_top,y_top] \
                    + (1-x_coeff) * y_coeff * F[:,x_top,y_top-1] 
                    
    elif method=='bicubic':
        # For a given point, bi-cubic interpolation fits the value of the 
        # four surrounding points as well as the slope at each of these 
        # points. The slope is evaluated using centered finite differences.
        # We can reduce this operation to a weight for each of the sixteen
        # surrounding points. We compute these weights for each grid cell.
        # Each of the 16 is thus stored in an array. A halo is added to 
        # the grid to handle boundary conditions.
        
        #Boundary conditions ---------------------------------------------
        F_halo=np.zeros((dim,Nx + 4, Ny + 4))
        F_halo[:, 2:-2, 2:-2] = F # Center
        
        F_halo[:, :2, 2:-2]  = F[:, -2:, :]#West band
        F_halo[:, -2:, 2:-2] = F[:, :2, :]#East band
        F_halo[:, 2:-2, :2]  = F[:, :, -2:]#South band
        F_halo[:, 2:-2, -2:] = F[:, :, :2]#North band
        
        F_halo[:, :2, :2] = F[:, -2:, -2:] #South-West corner
        F_halo[:, -2:, :2] = F[:, :2, -2:] #South-East corner
        F_halo[:, -2:, -2:] = F[:, :2, :2] #North-East corner
        F_halo[:, :2, -2:] = F[:, -2:, :2] #NOrth-West corner
        
        
        # Computation of the weights -------------------------------------
        a00 = F_halo[:, 1:-2, 1:-2]
        
        a01 = -0.5 * F_halo[:, 1:-2, 0:-3] + 0.5 * F_halo[:, 1:-2, 2:-1]
        
        a02 = F_halo[:, 1:-2, 0:-3] - 2.5 * F_halo[:, 1:-2, 1:-2] + \
            2 * F_halo[:, 1:-2, 2:-1] -  0.5 * F_halo[:, 1:-2, 3:]
            
        a03 = -.5 * F_halo[:, 1:-2, 0:-3] + 1.5 * F_halo[:, 1:-2, 1:-2] - \
            1.5 * F_halo[:, 1:-2, 2:-1] + 0.5 * F_halo[:, 1:-2, 3:]
            
        a10 = -.5 * F_halo[:, 0:-3, 1:-2] + 0.5 * F_halo[:, 2:-1, 1:-2]
        
        a11 = 0.25 * F_halo[:, 0:-3, 0:-3] - 0.25 * F_halo[:, 0:-3, 2:-1] - \
            0.25 * F_halo[:, 2:-1, 0:-3] + 0.25 * F_halo[:, 2:-1, 2:-1]
            
        a12 = -.5 * F_halo[:, 0:-3, 0:-3] + 1.25 * F_halo[:, 0:-3, 1:-2] - \
            F_halo[:, 0:-3, 2:-1] +  0.25 * F_halo[:, 0:-3, 3:] + \
                0.5 * F_halo[:, 2:-1, 0:-3] - 1.25 * F_halo[:, 2:-1, 1:-2] + \
                    F_halo[:, 2:-1, 2:-1] - 0.25 * F_halo[:, 2:-1, 3:]
                    
        a13 = 0.25 * F_halo[:, 0:-3, 0:-3] - 0.75 * F_halo[:, 0:-3, 1:-2] + \
            0.75 * F_halo[:, 0:-3, 2:-1] -  \
            0.25 * F_halo[:, 0:-3, 3:] - \
            0.25 * F_halo[:, 2:-1, 0:-3] + 0.75 * F_halo[:, 2:-1, 1:-2] - \
            0.75 * F_halo[:, 2:-1, 2:-1] + 0.25 * F_halo[:, 2:-1, 3:]
                
        a20 = F_halo[:, 0:-3, 1:-2] - 2.5 * F_halo[:, 1:-2, 1:-2] + \
            2 * F_halo[:, 2:-1, 1:-2] -  0.5 * F_halo[:, 3:, 1:-2]
            
        a21 = -.5 * F_halo[:, 0:-3, 0:-3] + 0.5 * F_halo[:, 0:-3, 2:-1] + \
            1.25 * F_halo[:, 1:-2, 0:-3] -  \
            1.25 * F_halo[:, 1:-2, 2:-1] - F_halo[:, 2:-1, 0:-3] + \
                F_halo[:, 2:-1, 2:-1] + 0.25 * F_halo[:, 3:, 0:-3] -  \
                    0.25 * F_halo[:, 3:, 2:-1]
                    
        a22 = F_halo[:, 0:-3, 0:-3] - 2.5 * F_halo[:, 0:-3, 1:-2] + \
            2 * F_halo[:, 0:-3, 2:-1] -  \
            0.5 * F_halo[:, 0:-3, 3:] - 2.5 * F_halo[:, 1:-2, 0:-3] + \
                6.25 * F_halo[:, 1:-2, 1:-2] - 5 * F_halo[:, 1:-2, 2:-1] + \
                1.25 * F_halo[:, 1:-2, 3:] + 2 * F_halo[:, 2:-1, 0:-3] \
                - 5 * F_halo[:, 2:-1, 1:-2] + 4 * F_halo[:, 2:-1, 2:-1] - \
                    F_halo[:, 2:-1, 3:] - 0.5 * F_halo[:, 3:, 0:-3] + \
                        1.25 * F_halo[:, 3:, 1:-2]  - F_halo[:, 3:, 2:-1] + \
                            0.25 * F_halo[:, 3:, 3:]
                            
        a23 = -.5 * F_halo[:, 0:-3, 0:-3] + 1.5 * F_halo[:, 0:-3, 1:-2] - \
            1.5 * F_halo[:, 0:-3, 2:-1] +  0.5 * F_halo[:, 0:-3, 3:] + \
            1.25 * F_halo[:, 1:-2, 0:-3] - 3.75 * F_halo[:, 1:-2, 1:-2] + \
                3.75 * F_halo[:, 1:-2, 2:-1] -  1.25 * F_halo[:, 1:-2, 3:] - \
                F_halo[:, 2:-1, 0:-3] + 3 * F_halo[:, 2:-1, 1:-2] - \
                    3 * F_halo[:, 2:-1, 2:-1] +  F_halo[:, 2:-1, 3:] + \
                        0.25 * F_halo[:, 3:, 0:-3] - \
                            0.75 * F_halo[:, 3:, 1:-2] +  \
                                0.75 * F_halo[:, 3:, 2:-1] - \
                                    0.25 * F_halo[:, 3:, 3:]
                    
        a30 = -.5 * F_halo[:, 0:-3, 1:-2] + 1.5 * F_halo[:, 1:-2, 1:-2] - \
            1.5 * F_halo[:, 2:-1, 1:-2] + 0.5 * F_halo[:, 3:, 1:-2]
            
        a31 = 0.25 * F_halo[:, 0:-3, 0:-3] - 0.25 * F_halo[:, 0:-3, 2:-1] - \
            0.75 * F_halo[:, 1:-2, 0:-3] + 0.75 * F_halo[:, 1:-2, 2:-1] + \
                0.75 * F_halo[:, 2:-1, 0:-3] - 0.75 * F_halo[:, 2:-1, 2:-1] - \
                    0.25 * F_halo[:, 3:, 0:-3] + 0.25 * F_halo[:, 3:, 2:-1];
        
        a32 = -.5 * F_halo[:, 0:-3, 0:-3] + 1.25 * F_halo[:, 0:-3, 1:-2] - \
            F_halo[:, 0:-3, 2:-1] + 0.25 * F_halo[:, 0:-3, 3:] + \
                1.5 * F_halo[:, 1:-2, 0:-3] - 3.75 * F_halo[:, 1:-2, 1:-2] + \
                    3 * F_halo[:, 1:-2, 2:-1] - 0.75 * F_halo[:, 1:-2, 3:] - \
                        1.5 * F_halo[:, 2:-1, 0:-3] + \
                            3.75 * F_halo[:, 2:-1, 1:-2] - \
                                3 * F_halo[:, 2:-1, 2:-1] +\
                                    0.75 * F_halo[:, 2:-1, 3:] + \
                                        0.5 * F_halo[:, 3:, 0:-3] - \
                                            1.25 * F_halo[:, 3:, 1:-2] + \
                                                F_halo[:, 3:, 2:-1] - \
                                                    0.25 * F_halo[:, 3:, 3:]
		
        a33 = 0.25 * F_halo[:, 0:-3, 0:-3] - 0.75 * F_halo[:, 0:-3, 1:-2] + \
            0.75 * F_halo[:, 0:-3, 2:-1] - 0.25 * F_halo[:, 0:-3, 3:] -  \
            0.75 * F_halo[:, 1:-2, 0:-3] + 2.25 * F_halo[:, 1:-2, 1:-2] - \
                2.25 * F_halo[:, 1:-2, 2:-1] + 0.75 * F_halo[:, 1:-2, 3:] + \
                0.75 * F_halo[:, 2:-1, 0:-3] - 2.25 * F_halo[:, 2:-1, 1:-2] +\
                    2.25 * F_halo[:, 2:-1, 2:-1] - 0.75 * F_halo[:, 2:-1, 3:]\
                        - 0.25 * F_halo[:, 3:, 0:-3] + \
                            0.75 * F_halo[:, 3:, 1:-2] - \
                                0.75 * F_halo[:, 3:, 2:-1] + \
                                    0.25 * F_halo[:, 3:, 3:]
    
    # Loop through the mesh points
        for x in range(Nx):
            for y in range(Ny):
                #Upstream coordinates
                xi = x - alpha_x[x,y]
                yi = y - alpha_y[x,y]
                                
                # Periodic boundary conditions
                #print("xi,yi : ", xi,yi)
                if (xi > 2*(Nx-1) or yi > 2*(Ny-1)):
                    print('Displacement too large: xi,yi=',xi,yi)
                
                if xi<0 :
                    xi += Nx-1 
                elif xi> Nx-1:
                    xi -= (Nx-1)
                if xi<0 :
                    yi += Ny-1
                elif yi> Ny-1:
                    yi -= (Ny-1)
                
                #Normalized coordinates inside the cell
                xb = xi - np.floor(xi)
                yb = yi - np.floor(yi)
                
                #Indice of the weights corresponding to the cell
                xf = int(np.floor(xi))+1
                yf = int(np.floor(yi))+1
                
                #Update of F
                F_int[:, x, y] = \
                    (a00[:, xf, yf] + a01[:, xf, yf] * yb +\
                     a02[:, xf, yf] * yb**2 + a03[:, xf, yf] * yb**3) +\
                    (a10[:, xf, yf] + a11[:, xf, yf] * yb +\
                     a12[:, xf, yf] * yb**2 + a13[:, xf, yf] * yb**3) * xb +\
                    (a20[:, xf, yf] + a21[:, xf, yf] * yb +\
                     a22[:, xf, yf] * yb**2 + a23[:, xf, yf] * yb**3) * xb**2 +\
                    (a30[:, xf, yf] + a31[:, xf, yf] * yb +\
                     a32[:, xf, yf] * yb**2 + a33[:, xf, yf] * yb**3) * xb**3 
    else:
        print('Unkwon method for interpolation: ',method)
    if dim==1:
        F_int = F_int[0,:,:]
    
    return F_int

