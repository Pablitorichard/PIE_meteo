import scipy.interpolate as interpolate


def advection_step_3P(alpha_u_minus, alpha_v_minus,field_minus, dt, u, v,
                      x_grid, y_grid):
    # The displacement is updated by interpolating the velocity field 
    alpha_u = dt * interpolate.griddata( 
            (x_grid.flatten(), y_grid.flatten()),u.flatten(),
            (x_grid - alpha_u_minus, y_grid - alpha_v_minus),
             method='cubic', fill_value = 0) 
        
    alpha_v = dt * interpolate.griddata( 
            (x_grid.flatten(), y_grid.flatten()), v.flatten(),
            (x_grid - alpha_u_minus, y_grid - alpha_v_minus),
             method='cubic', fill_value = 0)
        
    # The field at time tk + dt is udpated by interpolating the field at time
    # tk -dt at the locations x - 2* alpha
    field_plus = interpolate.griddata(
            (x_grid.flatten(), y_grid.flatten()), field_minus.flatten(),
            (x_grid - 2*alpha_u, y_grid - 2*alpha_v),
            method='cubic', fill_value = 0) 

    return alpha_u, alpha_v, field_plus
    

