import numpy as np
import sys
import os
#from math import exp
#from Bio.PDB import *


'''
For each FNM file -- add 0-vectors to so that 10 modes are avail

USAGE:
python add_zeros_to_fnm.py 
stdin/stdout

'''



def check_if_list_has_only_floats(list_to_check):
	status = True
	list_of_float_like_chars = ("-", ".", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
	for element in list_to_check:
		for charact in element:
			if charact not in list_of_float_like_chars:
				return False
	return status


num_modes = 0
avail_fnm_file_lines = list()
fnm_file_to_rd = sys.stdin
mode___2___num_vectors = {}
number_of_vectors_in_mode = 0  ## re-setting the counter
for line in fnm_file_to_rd:
	avail_fnm_file_lines.append(line)
	if line[0] is not "#":  ## ie, if we're not at file header
		ln_elmnts = list()
		ln_elmnts = line.split()
		if "BEGIN MODE" in line:
			num_modes += 1
			mode_indx = int(ln_elmnts[2])
			number_of_vectors_in_mode = 0  ## re-setting the counter
		if (len(ln_elmnts) == 3)  and  (check_if_list_has_only_floats(ln_elmnts)):
			number_of_vectors_in_mode += 1
		if line[0:3] == "END":
			mode___2___num_vectors[mode_indx] = number_of_vectors_in_mode
fnm_file_to_rd.close()
	#for mode in mode___2___num_vectors:
	#	print fnm_file + "    mode: " + str(mode) + "     num_modes: " + str(num_modes)

out_file_name_to_wrt = sys.stdout
for avail_line in avail_fnm_file_lines:
	out_file_name_to_wrt.write(avail_line)


num_remaining_modes_to_print = 10 - num_modes
starting_mode = 17 - num_remaining_modes_to_print
mode_to_print = starting_mode
while mode_to_print <= 16:
	out_file_name_to_wrt.write("BEGIN MODE " + str(mode_to_print) + "\n")
	vect = 0
	while vect < number_of_vectors_in_mode:
		out_file_name_to_wrt.write("0.0000000000 0.0000000000 0.0000000000\n")
		vect += 1
	out_file_name_to_wrt.write("END\n\n")
	mode_to_print += 1
out_file_name_to_wrt.close()


#print "\n\n --- run complete -- \n\n"

