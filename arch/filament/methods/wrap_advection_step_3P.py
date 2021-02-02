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
  
    
    cur_state.vrs['alpha_ut'] = a_ut
    cur_state.vrs['alpha_vt'] = a_vt
    
    new_state.vrs['theta_t'] = theta_new
    
    history.append(new_state)
