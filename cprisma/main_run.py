from cprisma.main_cprisma import *
from cprisma.cover import *
from pathlib import Path
import argparse
import sys
import os

# To get array_get.py
try:
	cwd = os.path.abspath(os.getcwd())
	sys.path.insert(0, cwd)
	from array_get import array_get
except ModuleNotFoundError:
	pass

def str2bool(v):

	if isinstance(v, bool):
	   return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def main():

	# Starting parser

	parser = argparse.ArgumentParser(prog='cprisma', description="CPRISMA - Coloring PRoteins by Inputs and Sets of Multiple Alignments")

	parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
	parser.add_argument('-ns', metavar='bool', help="rename the sequences based on a tuple", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-va', metavar='[int]', help="method to visualize the alignment", type=int, default=1, choices=[1, 2])
	parser.add_argument('-j', metavar='bool', help="join sequences (only available for va = 2)", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-hc', metavar='bool', help="hide the conservation line", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-n', metavar='[int]', help="change the first and subsequent numbers in the alignment", type=int, default=1)
	parser.add_argument('-l', metavar='[int]', help="number of residues per line", type=int, default=80)
	parser.add_argument('-a', metavar='[int]', help="define the accuracy of the input data", type=int, default=1)
	parser.add_argument('-t', metavar='bool', help="get all information for the different comparisons in log file", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-tr', metavar='bool', help="get a tuple of target-residues", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-ck', metavar='bool', help="check reference method", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-rf', metavar='[str]', help="method to compare sequences ('defaul', 'pair', 'multiple')", type=str, default='default', choices=['default', 'pair', 'multiple'])
	parser.add_argument('-lco', metavar='bool', help="display a list of the available colors", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-ico', metavar='[int]', help="multiply the intensity of the color", type=int, default=1)
	parser.add_argument('-sco', metavar='[int]', help="sequence color (when 'ssc' descriptor is applied)", type=int, default=6)                 # Default: Red in colors.csv
	parser.add_argument('-mco', metavar='[int]', help="mutation color (when 'fmac' or 'pimc' descriptors are applied)", type=int, default=137)   # Default: Gray in colors.csv
	parser.add_argument('-tco', metavar='bool', help="display color at the level of 3D structure through pymol script", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-dop', metavar='bool', help="get a operation descriptors dictionary (only available for rf = 'multiple')", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-dco', metavar='bool', help="get a color descriptors dictionary (only available for rf = 'multiple')", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-dvi', metavar='bool', help="get a visualization descriptors dictionary (only available for rf = 'multiple')", default=False, type=str2bool, nargs='?', const=True)
	parser.add_argument('-dmx', metavar='bool', help="get a maximum restriction descriptors dictionary (only available for rf = 'multiple')", default=False, type=str2bool, nargs='?', const=True)

	args = parser.parse_args()

	# Parameters of visualization
	check_name_seq = args.ns             # Checking if it is introduced a tuple of sequences name
	ali_method = args.va                 # Method to visualize the comparisons separated (1) or not (2)
	join = args.j                        # To see different comparison together or not (useful just for 'pair' or 'multiple' comparison methods). Note: just for ali_method=2.
	conservation = args.hc
	turn_comparison = args.t             # Visualize comparisons in log file
	first_num = args.n                   # Number that should appear in the alignment at first
	res_per_col = args.l                 # Residues per line (default = 80). Note: Try to use multiples of 80

	# Parameters for input data
	check_tar_res = args.tr              # Checking if it is considered a specific set of residues
	accuracy = args.a                    # Number of decimal places used for eventual calculation

	# Parameters for comparisons
	check_reference = args.ck            # If it is desired to do comparisons
	method_reference = args.rf           # Method for comparison: 'default', 'pair' or 'multiple'

	# Parameters for operations
	feature_ope = 'operation'
	check_dict_ope = args.dop          # If it is desired to do operations

	# Parameters to define maximum restriction
	feature_mxr = 'maximum'
	check_dict_mxr = args.dmx

	# Parameters to define color (reference without color always)
	feature_col = 'color'
	factor_color = args.ico                # To improve the intensity of color (default = 1)
	check_dict_col = args.dco
	sequence_col = args.sco
	mutation_col = args.mco
	tridimensional_col = args.tco

	if args.lco:
		dir_run=os.path.dirname(os.path.abspath(__file__))
		color_dir=dir_run+"/colors.csv"
		name_color_file=color_dir
		colors=pd.read_csv(name_color_file)
		print(colors.to_string(index=1))
		quit()

	# Parameters to define visualization
	feature_vis = 'visualization'
	check_dict_vis = args.dvi

	# Output files
	ali_log=open('cprisma.log', 'w+')

	cover(ali_log)

	# Saving command-line
	ali_log.write("Command-line: ")
	ali_log.write(" ".join(sys.argv))
	ali_log.write("\n\n")
	print(f'Command-line: {" ".join(sys.argv)}\n')

	# Input parameters #

	target_residues, name_sequence, descriptor_ope, descriptor_col, descriptor_vis, descriptor_mxr, dict_ref, dict_ope, dict_col, dict_vis, dict_mxr=array_get()

	# Operation
	calculation=Operation(turn_comparison,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence,
		  feature_ope,descriptor_ope,dict_ope,check_dict_ope)

	# Maximum restriction
	maximum=Maximum(turn_comparison,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence,
	feature_mxr,descriptor_mxr,dict_mxr,check_dict_mxr)

	# Color
	color=Color(turn_comparison,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence,
	feature_col,descriptor_col,dict_col,check_dict_col)

	# Visualization
	visualization=Visualization(turn_comparison,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence,
				feature_vis,descriptor_vis,dict_vis,check_dict_vis)

	#### Loop to execute the classes ####

	x=0
	for feature in (calculation,maximum,color,visualization):

		feature.treatment(x)
		feature.reference(check_reference,dict_ref,method_reference,x)

		if x == 0:
			df=pd.DataFrame()
			# It is created a matrix with data input
			df=feature.matrix_operation(df)
			df_ope=copy.deepcopy(df)
		elif x == 1:
			# The new matrix created in Operation Class is passed to Maximum Class
			feature.matrix_operation(df_ope)
			df_mxr=feature.matrix_maximum()
		elif x == 2:
			# The new matrix created in Operation Class is passed to Color Class for pkac descriptor
			feature.matrix_operation(df_ope)
			df_col, df_col_array=feature.matrix_color(mutation_col,sequence_col)
		elif x == 3:
			dict_ref,df_vis=feature.matrix_visualization()

		x+=1

	###### Invoking aligment object ######

	aligment=Alignment(turn_comparison,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence,
	   method_reference, dict_ref, df_ope, df_mxr, df_col, df_col_array, df_vis)

	aligment.treatment(x)
	aligment.add_color(factor_color, first_num, join, conservation, res_per_col, tridimensional_col, check_dict_vis, ali_method)

	ali_log.close()
