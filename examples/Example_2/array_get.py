def array_get():

	# Tuple of residues 
	target_residues = ('D', 'E', 'R', 'H', 'K', 'Y')

	# Tuple of sequences name
	name_sequence = ('SASA', 'pKa_NS1 (alone)','PROCEEDpKa', 'PISA', 'DiscoTope', 'ElliPro', 'SEPPA2', 'SEPPA3', 'Consensus', '')

	# Operations
	descriptor_ope = { 'n' }

	# Color
	descriptor_col = { 'nc' }

	# Visualization
	descriptor_vis = { 'ReY','DeN','LeY' }

	# Maximum restriction
	descriptor_mxr = { 'nr' }

	# Array for comparisons
	dict_ref = { 'main_0': [ {'0' : [] } ] ,
		         'main_1': [ {'1' : [2] } ] ,
		         'main_3': [ {'3' : [] } ] ,
		         'main_4': [ {'4' : [] } ] ,
		         'main_5': [ {'5' : [] } ] ,
		         'main_6': [ {'6' : [] } ] ,
		         'main_7': [ {'7' : [] } ] ,
		         'main_8': [ {'8' : [] } ] ,
		         'main_9': [ {'9' : [] } ] }

	# Array for operations
	dict_ope = { 'main_0': [ { 'n' } ] ,
		         'main_1': [ { 'da' } ] ,
		         'main_3': [ { 'n' } ] ,
		         'main_4': [ { 'n' } ] ,
		         'main_5': [ { 'n' } ] ,
		         'main_6': [ { 'n' } ] ,
		         'main_7': [ { 'n' } ] ,
		         'main_8': [ { 'n' } ] ,
		         'main_9': [ { 'n' } ] }

	# Array for color

	dict_col = { 'main_0': [ { 'tc' : [54, 0, '>'] } ] ,
		         'main_1': [ { 'tc' : [6, 0.007, '>'] } ] ,
		         'main_3': [ { 'tc' : [75, 0, '>'] } ] ,
		         'main_4': [ { 'tc' : [20, -6.85, '>='] } ] ,
		         'main_5': [ { 'tc' : [64, 0, '>'] } ] ,
		         'main_6': [ { 'tc' : [22, 0.049, '>='] } ] ,
		         'main_7': [ { 'tc' : [96, 0.093, '>='] } ] ,
		         'main_8': [ { 'tc' : [92, 0, '>'] } ] ,
		         'main_9': [ { 'ssc' } ] }

	# Array for visualization
	dict_vis = { 'main_0': [ { 'ReY','DeN','LeY' } ] ,
		         'main_1': [ { 'ReN','DeY','LeY' } ] ,
		         'main_3': [ { 'ReY','DeN','LeY' } ] ,
		         'main_4': [ { 'ReY','DeN','LeY' } ] ,
		         'main_5': [ { 'ReY','DeN','LeY' } ] ,
		         'main_6': [ { 'ReY','DeN','LeY' } ] ,
		         'main_7': [ { 'ReY','DeN','LeY' } ] ,
		         'main_8': [ { 'ReY','DeN','LeY' } ] ,
		         'main_9': [ { 'ReY','DeY','LeN' } ] }

	# Array for maximum restriction
	dict_mxr={ }

	return target_residues, name_sequence, descriptor_ope, descriptor_col, descriptor_vis, descriptor_mxr, dict_ref, dict_ope, dict_col, dict_vis, dict_mxr
