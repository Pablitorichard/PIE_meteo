from .spectral import geostwind

def pseudo_spectral_wind(history, grid, params, **kwargs):
	
	assert history.size > 0
	current_state = history.state_list[-1]
	
	ut, vt = geostwind(grid.Lx, grid.Ly, current_state.vrs['theta'], z=0)
	us, vs = geostwind(grid.Lx, grid.Ly, current_state.vrs['theta'], z=params['z_star'])
	
	
	current_state.vrs['ut'] = ut
	current_state.vrs['vt'] = vt
	current_state.vrs['us'] = us
	current_state.vrs['vt'] = vt	
	
