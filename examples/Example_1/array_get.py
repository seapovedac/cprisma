def array_get():

	# Tuple of residues 
	target_residues = ('D', 'E', 'R', 'H', 'K', 'Y')

	# Tuple of sequences name
	name_sequence = ('UG_chA_x-ray', 'UG_chB_x-ray', 'UG_chA_mod1', 'UG_chA_mod2', 'SE_chA_mod1', 'CAR_chA_mod1', 'MA_chA_mod1', 'TH_chA_mod1', 'YAP_chA_mod1', 'BR_chA_x-ray', 'BR_chB_x-ray', 'BR_chA_mod1', 'BR_chA_mod2', '')

	# Operations
	descriptor_ope = { 'n' }

	# Color
	descriptor_col = { 'nc' }

	# Visualization
	descriptor_vis = { 'ReY','DeN','LeY' }

	# Maximum restriction
	descriptor_mxr = { 'nr' }

	# Array for comparisons
	dict_ref = { 'main_0': [ {'0' : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] } ] ,
		         'main_13': [ {'13' : [] } ] }

	# Array for operations
	dict_ope = { 'main_0': [ { 'd' } ] ,
		         'main_13': [ { 'n' } ] }

	# Array for color
	dict_col = { 'main_0' : [ { 'pkac': True } ] ,
				 'main_13' : [ { 'ssc' } ] }

	# Array for visualization
	dict_vis = { 'main_0': [ { 'ReN','DeY','LeY' } ] ,
		         'main_13': [ { 'ReY','DeY','LeN' } ] }

	# Array for maximum restriction
	dict_mxr = { 'main_0': [ { 'rm' } ] ,                 
		         'main_13': [ { 'nr' } ] }

	return target_residues, name_sequence, descriptor_ope, descriptor_col, descriptor_vis, descriptor_mxr, dict_ref, dict_ope, dict_col, dict_vis, dict_mxr
