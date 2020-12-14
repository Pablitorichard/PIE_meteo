import scipy.interpolate as interpolate


def advection_step_3P(alpha_u_minus, alpha_v_minus,field_minus, dt, u, v,
                      x_grid, y_grid, alpha_method='linear',
                      order_alpha=2, F_method='cubic'):
    
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
            method = 'cubic'
        elif k==order_alpha-1 and alpha_method =='mix':
            method='linear'
        else:
            method=alpha_method
  
        alpha_u = dt * interpolate.griddata( 
            (x_grid.flatten(), y_grid.flatten()),u.flatten(),
            (x_grid - alpha_u_minus, y_grid - alpha_v_minus),
             method=method, fill_value = 0)
        alpha_v = dt * interpolate.griddata( 
            (x_grid.flatten(), y_grid.flatten()), v.flatten(),
            (x_grid - alpha_u_minus, y_grid - alpha_v_minus),
             method=method, fill_value = 0)
        
        alpha_u_minus = alpha_u
        alpha_v_minus = alpha_v
        
#------------------------------------------------------------------------


    # The field at time tk + dt is udpated by interpolating the field at 
    # time tk -dt at the locations x - 2* alpha.
    field_plus = interpolate.griddata(
            (x_grid.flatten(), y_grid.flatten()), field_minus.flatten(),
            (x_grid - 2*alpha_u, y_grid - 2*alpha_v),
            method='cubic', fill_value = 0) 

    return alpha_u, alpha_v, field_plus
    

