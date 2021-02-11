import numpy as np

from .upstream_interp import upstream_interp

def advection_step_3P(alpha_u_minus, alpha_v_minus, field_minus,
                      dt, u, v, dx, dy,
                      alpha_method,
                      order_alpha,
                      F_method,
                      verbose=0):
    
    # if alpha_method == 'damped_bicubic' or F_method == 'damped_bicubic':
    #     d = 0.5 * np.sqrt(
    #         np.square( (np.roll(u,-1,0) - np.roll(u,1,0)) / (2 * dx)   - \
    #                    (np.roll(v,-1,1) - np.roll(v,1,1)) / (2 * dy) ) + \
    #         np.square( (np.roll(u,-1,1) - np.roll(u,1,1)) / (2 * dy)   + \
    #                    (np.roll(v,-1,0) - np.roll(v,1,0)) / (2 * dx) ) )
    
    # ITERATIVE ESTIMATION OF THE DISPLACEMENT-------------------------------
    # The displacement alpha is estimated iteratively by interpolating the wind 
    # in x -alpha_minus, where alpha is the previous estimate. At each time
    # step, a number of iterations order_alpha is used, and the iterative 
    # scheme is initialized with the estimate at the previous time step.
    for k in range(order_alpha):
        # Staniforth et Al. states that for the interpolation of the 
        # estimated displacement, linear interpolation is usually 
        # sufficient. Bicubic interpolation is more accurate but more 
        # expensive. Damped bicubic interpolation also combines the accuracy 
        # of bicubic interpolation with a relaxation based on the linear 
        # interpolation to limit small scale numerical noise.
        if k < order_alpha-1 and alpha_method !='linear':
            method = 'linear'
        else:
            method = alpha_method
            
        print("      advection_step_3P with alpha order "+str(k)+" and method "\
              +method) if verbose > 2 else None
        
        [alpha_u, alpha_v] = (dt/dx)*upstream_interp(alpha_u_minus,
                                                     alpha_v_minus,
                                                     np.array([u,v]),
                                                     method=method,
                                                     verbose=verbose)
        alpha_u_minus = alpha_u
        alpha_v_minus = alpha_v
        
#------------------------------------------------------------------------
    
    # The field at time tk + dt is udpated by interpolating the field at 
    # time tk -dt at the locations x - 2* alpha.
    field_plus = upstream_interp(2*alpha_u, 2*alpha_v, field_minus,
                                 method=F_method, verbose=verbose)

    return alpha_u, alpha_v, field_plus
    

