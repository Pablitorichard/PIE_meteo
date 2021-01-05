from upstream_interp import upstream_interp
import numpy as np

def advection_step_3P(alpha_u_minus, alpha_v_minus,field_minus, dt, u, v,
                      dx, dy, alpha_method='linear',
                      order_alpha=2, F_method='bicubic'):
    
# ITERATIVE ESTIMATION OF THE DISPLACEMENT-------------------------------
# The displacement alpha is estimated iteratively by interpolating the wind 
# in x -alpha_minus, where alpha is the previous estimate. At each time
# step, a number of iterations order_alpha is used, and the iterative 
# scheme is initialized with the estimate at the previous time step.
    for k in range (order_alpha):
        # Staniforth et Al. states that for the interpolation of the 
        # estimated displacement, linear interpolation is usually 
        # sufficient. Cubic interpolation is more accurate but more 
        # expensive. An interesting compromise is to use the option 'mix',
        # which uses linear inteprolation for all iterations but the last
        # one, for which cubic interpolation is used.
        if k==order_alpha-1 and alpha_method =='mix':
            method = 'bicubic'
        elif k < order_alpha-1 and alpha_method =='mix':
            method='linear'
        else:
            method=alpha_method
            
        [alpha_u, alpha_v] = (dt/dx)*upstream_interp(alpha_u_minus,
                                        alpha_v_minus,
                                        np.array([u,v]),method=method)
        alpha_u_minus = alpha_u
        alpha_v_minus = alpha_v
        
#------------------------------------------------------------------------
    
    # The field at time tk + dt is udpated by interpolating the field at 
    # time tk -dt at the locations x - 2* alpha.
    field_plus = upstream_interp(2*alpha_u, 2*alpha_v, field_minus,
                                 method = F_method)

    return alpha_u, alpha_v, field_plus
    

