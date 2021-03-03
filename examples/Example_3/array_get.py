def array_get():

	# Tuple of residues 
	target_residues = ()

	# Tuple of sequences name
	name_sequence = ('', 'ZIKV-UG', 'ZIKV-BR')

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
		         'main_2': [ {'2' : [1] } ] ,
				 'main_1': [ {'1' : [2] } ] }

	# Array for operations
	dict_ope = {}

	# Array for color
	dict_col = { 'main_0' : [ { 'pic' : [ { '0' : [81, 'n', '5' ,'7', '12', '14', '27-30', '114-119', '121-122', '159', '161-164'] } , 
										  { '0' : [35, 'n', '207-209', '231-239', '252', '255-261', '263-264', '286', '288', '292-294', '313-316', '332', '350'] } ] } ] ,
				 'main_2' : [ { 'pimc' : [ { '1' : [20, 'n', '0-28'] } , 
										   { '1' : [6, 'n', '29-37'] } , 
										   { '1' : [94, 'n', '38-149'] } , 
										   { '1' : [6, 'n', '150-179'] } , 
										   { '1' : [64, 'n', '180-351'] } ] } ] ,
				 'main_1' : [ { 'pimc' : [ { '2' : [20, 'n', '0-28'] } , 
		                                   { '2' : [6, 'n', '29-37'] } , 
										   { '2' : [94, 'n', '38-149'] } , 
										   { '2' : [6, 'n', '150-179'] } , 
										   { '2' : [64, 'n', '180-351'] } ] } ] }

	# Array for visualization
	dict_vis = { 'main_0': [ { 'ReY','DeN','LeN' } ] ,
		         'main_2': [ { 'ReN','DeN','LeY' } ] ,
				 'main_1': [ { 'ReN','DeN','LeY' } ] }

	# Array for maximum restriction
	dict_mxr = {}

	return target_residues, name_sequence, descriptor_ope, descriptor_col, descriptor_vis, descriptor_mxr, dict_ref, dict_ope, dict_col, dict_vis, dict_mxr
