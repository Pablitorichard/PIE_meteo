from .advection_step_3P import advection_step_3P
# from .spectral import vertwind
from ..core.state import State #, variables

def wrap_advection_step_3P(history, grid, params, alpha_method, order_alpha, F_method, verbose=0, **kwargs):

    assert history.size > 1
    pre_state = history.state_list[-2]
    cur_state = history.state_list[-1]

    dt = cur_state.t - pre_state.t # constant step

    new_state = State.copy(cur_state)
    new_state.t += dt              # constant step
    
    a_ut, a_vt, theta_new = advection_step_3P(pre_state.vrs['alpha_ut'],
                                              pre_state.vrs['alpha_vt'],
                                              pre_state.vrs['theta_t'],
                                              dt,
                                              cur_state.vrs['ut'],
                                              cur_state.vrs['vt'],
                                              grid.dx,
                                              grid.dy,
                                              alpha_method,
                                              order_alpha,
                                              F_method,
                                              verbose)
    print("      ut vt done") if verbose > 2 else None
    a_us, a_vs, dz_new = advection_step_3P(pre_state.vrs['alpha_us'],
                                           pre_state.vrs['alpha_vs'],
                                           pre_state.vrs['Delta_z'],
                                           dt,
                                           cur_state.vrs['us'],
                                           cur_state.vrs['vs'],
                                           grid.dx,
                                           grid.dy,
                                           alpha_method,
                                           order_alpha,
                                           F_method,
                                           verbose)
    print("      us vs done") if verbose > 2 else None
    
    #UPDATE OF W ---------------------------------------------------
    #k_hour = int(3600/dt)
    #if ((cur_state.t-2)%k_hour==0 ):# k-1 -> k ?
        #w = vertwind(grid.Lx, grid.Ly, cur_state.vrs['theta'], pre_state.vrs['theta'], dt, z=params.z_star)
        #print("w: ",np.max(w)," , ", np.min(w))
        #cur_state.vrs['Delta_z'] += k_hour * dt * w
        #dz_next += k_hour * dt * w
        
    #dT_disp = params.gamma_2 * cur_state.vrs['Delta_z']
    #dT_cloud = params.Delta_Tc * ( cur_state.vrs['Delta_z'] > params.Delta_zc ) 
            
    #params.Delta_T_bb 
    
    cur_state.vrs['alpha_ut'] = a_ut
    cur_state.vrs['alpha_vt'] = a_vt
    cur_state.vrs['alpha_us'] = a_us
    cur_state.vrs['alpha_vs'] = a_vs
    
    new_state.vrs['theta_t'] = theta_new
    new_state.vrs['Delta_z'] = dz_new # Pas sûr !
    
    
    #curr_state['Delta_z'] = Delta_z
    #curr_state['Delta_T_bb'] = Delta_T_bb
    
    history.append(new_state)
    history.pop(0)
