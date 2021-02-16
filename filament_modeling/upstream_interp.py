"""
 alpha_interp intrpolates a multidimensionnal field F from a 2D grid to 
 to an 'upstream' unstrctured mesh defined by the displacements 
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
    elif method=='linear_V1':

        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = np.mod(X - alpha_x, Nx - 1)
        Yi = np.mod(Y - alpha_y, Ny - 1)
        
        Xt = np.ceil(Xi).astype(int)
        Yt = np.ceil(Yi).astype(int)
        
        Xc = Xt - Xi
        Yc = Yt - Yi
        
        F_int[:, X, Y] = Xc * Yc * F[:, Xt - 1, Yt -1] \
                   + Xc * (1 - Yc) * F[:, Xt - 1, Yt] \
                   + (1 - Xc) * (1 - Yc) * F[:, Xt, Yt] \
                   + (1 - Xc) * Yc * F[:, Xt, Yt - 1] 
                   
    elif method=='linear_V2':

        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = X - alpha_x
        Yi = Y - alpha_y
        
        Xt = np.ceil(Xi).astype(int)
        Yt = np.ceil(Yi).astype(int)
        
        Xc = Xt - Xi
        Yc = Yt - Yi
        
        Xt = np.mod(Xt, Nx)
        Yt = np.mod(Yt, Ny)
        
        Ft = F[:, Xt, Yt]
        
        F_int = Xc * Yc * np.roll(Ft, (1, 1), (1,2)) \
                   + Xc * (1 - Yc) * np.roll(Ft, 1, 1) \
                   + (1 - Xc) * (1 - Yc) * Ft \
                   + (1 - Xc) * Yc * np.roll(Ft, 1, 2) 
                   
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
    
    elif method=='bicubic_V1':
        # For a given point, bi-cubic interpolation fits the value of the 
        # four surrounding points as well as the slope at each of these 
        # points. The slope is evaluated using centered finite differences.
        # We can reduce this operation to a weight for each of the sixteen
        # surrounding points. We compute these weights for each grid cell.
        # Each of the 16 is thus stored in an array..
        
        #-----------------------------------------------------------------
        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = X - alpha_x
        Yi = Y - alpha_y
        
        Xf = np.mod( np.floor(Xi).astype(int), Nx)
        Yf = np.mod( np.floor(Yi).astype(int), Ny)
        
        Xb = np.abs(Xi - np.floor(Xi))
        Yb = np.abs(Yi - np.floor(Yi))
        
        Ft = F[:, Xf, Yf]
        
        # We define 'rolled' arrays that appears multiple times 
        F00 = np.roll(Ft, (1,1), (1,2))
        F01 = np.roll(Ft, (1,0), (1,2))
        F02 = np.roll(Ft, (1,-1), (1,2))
        F03 = np.roll(Ft, (1,-2), (1,2))
        F10 = np.roll(Ft, (0,1), (1,2))
        F11 = np.roll(Ft, (0,0), (1,2))
        F12 = np.roll(Ft, (0,-1), (1,2))
        F13 = np.roll(Ft, (0,-2), (1,2))
        F20 = np.roll(Ft, (-1,1), (1,2))
        F21 = np.roll(Ft, (-1,0), (1,2))
        F22 = np.roll(Ft, (-1,-1), (1,2))
        F23 = np.roll(Ft, (-1,-2), (1,2))
        F30 = np.roll(Ft, (-2,1), (1,2))
        F31 = np.roll(Ft, (-2,0), (1,2))
        F32 = np.roll(Ft, (-2,-1), (1,2))
        F33 = np.roll(Ft, (-2,-2), (1,2))
        
        # Computation of the weights -------------------------------------
        A00 = F11
        
        A01 = -0.5 * F10 + 0.5 * F12
        
        A02 = F10 - 2.5 * F11 + 2 * F12 -  0.5 * F13
            
        A03 = -.5 * F10 + 1.5 * F11 - 1.5 * F12 + 0.5 * F13
            
        A10 = -.5 * F01 + 0.5 * F21
        
        A11 = 0.25 * F00 - 0.25 * F02 - 0.25 * F20 + 0.25 * F22
            
        A12 = -.5 * F00 + 1.25 * F01 - F02 +  0.25 * F03 + 0.5 * F20 - 1.25  \
            * F21 + F22 - 0.25 * F23
                    
        A13 = 0.25 * F00 - 0.75 * F01 + 0.75 * F02 - 0.25 * F03 - 0.25 * F20 \
            + 0.75 * F21 - 0.75 * F22 + 0.25 * F23
                
        A20 = F01 - 2.5 * F11 + 2 * F21 -  0.5 * F31
            
        A21 = -.5 * F00 + 0.5 * F02 + 1.25 * F10 - 1.25 * F12 - F20 + F22 + \
            0.25 * F30 - 0.25 * F32
                    
        A22 = F00 - 2.5 * F01 + 2 * F02 - 0.5 * F03 - 2.5 * F10 +  6.25 * \
            F11 - 5 * F12 + 1.25 * F13 + 2 * F20  - 5 * F21 + 4 * F22 - \
            F23 - 0.5 * F30 + 1.25 * F31  - F32 + 0.25 * F33
                            
        A23 = -.5 * F00 + 1.5 * F01 - 1.5 * F02 +  0.5 * F03 + 1.25 * F10 - \
            3.75 * F11 + 3.75 * F12 -  1.25 * F13 - F20 + 3 * F21 - 3 * F22 \
                +  F23 + 0.25 * F30 - 0.75 * F31 + 0.75 * F32 - 0.25 * F33
                    
        A30 = -.5 * F01 + 1.5 * F11 - 1.5 * F21 + 0.5 * F31
            
        A31 = 0.25 * F00 - 0.25 * F02 - 0.75 * F10 + 0.75 * F12 + 0.75 * F20 \
            - 0.75 * F22 - 0.25 * F30 + 0.25 * F32;
        
        A32 = -.5 * F00 + 1.25 * F01 - F02 + 0.25 * F03 + 1.5 * F10 - 3.75 *\
            F11 + 3 * F12 - 0.75 * F13 - 1.5 * F20 + 3.75 * F21 - 3 * F22 +\
            0.75 * F23 + 0.5 * F30 - 1.25 * F31 + F32 - 0.25 * F33
		
        A33 = 0.25 * F00 - 0.75 * F01 + 0.75 * F02 - 0.25 * F03 - 0.75 * F10 \
            + 2.25 * F11 - 2.25 * F12 + 0.75 * F13 + 0.75 * F20 - 2.25 * F21 +\
            2.25 * F22 - 0.75 * F23 - 0.25 * F30 + 0.75 * F31 - 0.75 * F32 +\
            0.25 * F33
    
        
        #Update of F
        F_int = (A00 + A01 * Yb + A02 * Yb**2 + A03 * Yb**3) + \
                    (A10 + A11 * Yb + A12 * Yb**2 + A13 * Yb**3) * Xb + \
                    (A20 + A21 * Yb + A22 * Yb**2 + A23 * Yb**3) * Xb**2 + \
                    (A30 + A31 * Yb + A32 * Yb**2 + A33 * Yb**3) * Xb**3 
    
    elif method=='bicubic_V2':
        # For a given point, bi-cubic interpolation fits the value of the 
        # four surrounding points as well as the slope at each of these 
        # points. The slope is evaluated using centered finite differences.
        # We can reduce this operation to a weight for each of the sixteen
        # surrounding points. We compute these weights for each grid cell.
        # Each of the 16 is thus stored in an array..
        
        #-----------------------------------------------------------------
        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = X - alpha_x
        Yi = Y - alpha_y
        
        Xf = np.mod( np.floor(Xi).astype(int), Nx)
        Yf = np.mod( np.floor(Yi).astype(int), Ny)
        
        Xb = np.abs(Xi - np.floor(Xi))
        Yb = np.abs(Yi - np.floor(Yi))
        
        
        
        # We define 'rolled' arrays that appears multiple times 
        F00 = np.roll(F, (1,1), (1,2))
        F01 = np.roll(F, (1,0), (1,2))
        F02 = np.roll(F, (1,-1), (1,2))
        F03 = np.roll(F, (1,-2), (1,2))
        F10 = np.roll(F, (0,1), (1,2))
        F11 = np.roll(F, (0,0), (1,2))
        F12 = np.roll(F, (0,-1), (1,2))
        F13 = np.roll(F, (0,-2), (1,2))
        F20 = np.roll(F, (-1,1), (1,2))
        F21 = np.roll(F, (-1,0), (1,2))
        F22 = np.roll(F, (-1,-1), (1,2))
        F23 = np.roll(F, (-1,-2), (1,2))
        F30 = np.roll(F, (-2,1), (1,2))
        F31 = np.roll(F, (-2,0), (1,2))
        F32 = np.roll(F, (-2,-1), (1,2))
        F33 = np.roll(F, (-2,-2), (1,2))
        
        # Computation of the weights -------------------------------------
        A00 = F11
        
        A01 = -0.5 * F10 + 0.5 * F12
        
        A02 = F10 - 2.5 * F11 + 2 * F12 -  0.5 * F13
            
        A03 = -.5 * F10 + 1.5 * F11 - 1.5 * F12 + 0.5 * F13
            
        A10 = -.5 * F01 + 0.5 * F21
        
        A11 = 0.25 * F00 - 0.25 * F02 - 0.25 * F20 + 0.25 * F22
            
        A12 = -.5 * F00 + 1.25 * F01 - F02 +  0.25 * F03 + 0.5 * F20 - 1.25  \
            * F21 + F22 - 0.25 * F23
                    
        A13 = 0.25 * F00 - 0.75 * F01 + 0.75 * F02 - 0.25 * F03 - 0.25 * F20 \
            + 0.75 * F21 - 0.75 * F22 + 0.25 * F23
                
        A20 = F01 - 2.5 * F11 + 2 * F21 -  0.5 * F31
            
        A21 = -.5 * F00 + 0.5 * F02 + 1.25 * F10 - 1.25 * F12 - F20 + F22 + \
            0.25 * F30 - 0.25 * F32
                    
        A22 = F00 - 2.5 * F01 + 2 * F02 - 0.5 * F03 - 2.5 * F10 +  6.25 * \
            F11 - 5 * F12 + 1.25 * F13 + 2 * F20  - 5 * F21 + 4 * F22 - \
            F23 - 0.5 * F30 + 1.25 * F31  - F32 + 0.25 * F33
                            
        A23 = -.5 * F00 + 1.5 * F01 - 1.5 * F02 +  0.5 * F03 + 1.25 * F10 - \
            3.75 * F11 + 3.75 * F12 -  1.25 * F13 - F20 + 3 * F21 - 3 * F22 \
                +  F23 + 0.25 * F30 - 0.75 * F31 + 0.75 * F32 - 0.25 * F33
                    
        A30 = -.5 * F01 + 1.5 * F11 - 1.5 * F21 + 0.5 * F31
            
        A31 = 0.25 * F00 - 0.25 * F02 - 0.75 * F10 + 0.75 * F12 + 0.75 * F20 \
            - 0.75 * F22 - 0.25 * F30 + 0.25 * F32;
        
        A32 = -.5 * F00 + 1.25 * F01 - F02 + 0.25 * F03 + 1.5 * F10 - 3.75 *\
            F11 + 3 * F12 - 0.75 * F13 - 1.5 * F20 + 3.75 * F21 - 3 * F22 +\
            0.75 * F23 + 0.5 * F30 - 1.25 * F31 + F32 - 0.25 * F33
		
        A33 = 0.25 * F00 - 0.75 * F01 + 0.75 * F02 - 0.25 * F03 - 0.75 * F10 \
            + 2.25 * F11 - 2.25 * F12 + 0.75 * F13 + 0.75 * F20 - 2.25 * F21 +\
            2.25 * F22 - 0.75 * F23 - 0.25 * F30 + 0.75 * F31 - 0.75 * F32 +\
            0.25 * F33
    
        
        #Update of F
        F_int = (A00[:,Xf, Yf] + A01[:,Xf, Yf] * Yb + A02[:,Xf, Yf] * Yb**2 + \
                 A03[:,Xf, Yf] * Yb**3) + \
                (A10[:,Xf, Yf] + A11[:,Xf, Yf] * Yb + A12[:,Xf, Yf] * Yb**2 + \
                 A13[:,Xf, Yf] * Yb**3) * Xb + \
                (A20[:,Xf, Yf] + A21[:,Xf, Yf] * Yb + A22[:,Xf, Yf] * Yb**2 + \
                 A23[:,Xf, Yf] * Yb**3) * Xb**2 + \
                (A30[:,Xf, Yf] + A31[:,Xf, Yf] * Yb + A32[:,Xf, Yf] * Yb**2 + \
                 A33[:,Xf, Yf] * Yb**3) * Xb**3 
    else:
        print('Unkwon method for interpolation: ',method)
    if dim==1:
        F_int = F_int[0,:,:]
    
    return F_int

