import numpy as np
import sys
import os
import math
from decimal import *
#from math import exp
#from Bio.PDB import *


'''
This script calculated the correlated motions for pairs of residues which are in contact. The correlation 
between residues i and j is calculated as:

	Cov(i,j) / sqrt[  avg[(del_i)^2]   *   avg[(del_j)^2]  ]

	where Cov(i,j) = avg[  del_i * del_j  ]

	and where del_i corresponds to the VECTOR in residue i for a particular mode

	****** NOTE: all of the first 10 modes from the MMTK output are treated equally -- is this how it should be??


example usage:
###python correlate_motions.py /Users/admin/Desktop/rsch/allostery/networks/rosvall/infomap_undir/1n78.net 1N78_bio_anm_compat.fnm_t10 > 1N78_bio_corr_weighted.net
###python correlate_motions.py /Users/admin/Desktop/rsch/allostery/networks/make_contact_map/1FCK_out.net 1FCK.fnm_t10 > 1fck_o.txt
###python correlate_motions.py /Users/admin/Desktop/rsch/allostery/networks/rosvall/infomap_undir/1n78_tmp.net 1N78_bio_anm_compat.fnm_t10 > r.txt


python correlate_motions_ver2.py 2DUB_cont_map_BIN.txt 2DUB_CA.pdb 2DUB.fnm_t10 > 2DUB_WGH_neg_log_cont_map.txt

python correlate_motions_ver2.py 3KMW_cont_map_BIN.txt 3KMW_CA.pdb 3KMW.fnm_t10 > 3KMW_WGH_neg_log_cont_map.txt

python correlate_motions_ver2.py 1EDH_cont_map_BIN.txt 1EDH_CA.pdb 1EDH.fnm_t10 > 1EDH_WGH_neg_log_cont_map.txt



'''


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 3):
   sys.stderr.write("Usage: " + sys.argv[0] + " <heavy/ca> <distance_cutoff> \n")
   sys.stderr.write("where, for ex:\n")
   sys.stderr.write("2DUB_CA.pdb\n")
   sys.stderr.write("2DUB.fnm_t10\n")
   sys.exit()




##  Read in the CA pdb file to determine the number of residues
ca_pdb = sys.argv[1]
ca_pdb_to_rd = open(ca_pdb, "r")
num_residues_frm_CA_PDB_file = 0
for line in ca_pdb_to_rd:
	num_residues_frm_CA_PDB_file += 1
#print "num_residues_frm_CA_PDB_file:  " + str(num_residues_frm_CA_PDB_file)
ca_pdb_to_rd.close()




##  Read in data from ANM file:
anm_file = sys.argv[2]
anm_file_handle = open(anm_file, "r")
mode_index___2___num_NODES = {}
modes_data = []
for line in anm_file_handle:
	if "BEGIN MODE 7" in line: ##  Here we're at the very first mode, so make new mode struct
		mode = 0
		del_x = list()
		del_y = list()
		del_z = list()
		modes_data.append([])
	elif "BEGIN MODE" in line:  ##  Then we're at the start of a mode other than the first in the anm file
		##  First store all data for the prev mode:
		modes_data[mode].append(del_x)
		modes_data[mode].append(del_y)
		modes_data[mode].append(del_z)
		if (len(del_x) != len(del_y))  or  (len(del_x) != len(del_z))  or  (len(del_y) != len(del_z)):
			print "\nERROR: UNEQUAL DISTRIBUTION OF NUM NODES ! \n"
			sys.exit()
		else:
			mode_index___2___num_NODES[mode] = len(del_x)

		##  Now start building struct for new mode
		mode += 1
		del_x = list()
		del_y = list()
		del_z = list()
		modes_data.append([])
	if "." in line:   ##  if there is a "." in the line, then it contains a vector of values!
		# break up line and store vals for del_x, y, z
		vals = list()
		vals = line.split()
		if len(vals) != 3:
			print "\nERROR: INCOMPLETE VECTOR ! \n"
			sys.exit()
		del_x.append(float(vals[0]))
		del_y.append(float(vals[1]))
		del_z.append(float(vals[2]))

## At very end, save data for last mode:
modes_data[mode].append(del_x)
modes_data[mode].append(del_y)
modes_data[mode].append(del_z)
## Give error & exit if is something looks amiss
if (len(del_x) != len(del_y))  or  (len(del_x) != len(del_z))  or  (len(del_y) != len(del_z)):
	print "\nERROR: UNEQUAL DISTRIBUTION OF NUM NODES ! \n"
	sys.exit()
else:
	mode_index___2___num_NODES[mode] = len(del_x)

num_modes = mode + 1

##  For each mode, ensure that the number of NODES in the ANM file to matches num nodes in the CA_PDB file
num_nodes_in_fnm_network = len(del_x)
m = 0
while m < num_modes:
	if mode_index___2___num_NODES[m] != num_residues_frm_CA_PDB_file:
		print "\nERROR: num nodes in anm file does not correspond to num residues in the CA pdb file!\n"
		sys.exit()
	m += 1

'''
##  As a check, print data to ensure that it was stored properly:
print "\n\nStarting to print the fnm data structs here..."
m = 0
while m < len(modes_data):
	print "BEGIN MODE " + str(m+7)
	n = 0
	while n < mode_index___2___num_NODES[m]:
		#print "x,y,z : " + str(modes_data[m][0][n]) + " " + str(modes_data[m][1][n]) + " " + str(modes_data[m][2][n])
		#print str(modes_data[m][0][n]) + " " + str(modes_data[m][1][n]) + " " + str(modes_data[m][2][n])
		print str(format(modes_data[m][0][n], '.10f')) + " " + str(format(modes_data[m][1][n], '.10f')) + " " + str(format(modes_data[m][2][n], '.10f'))
		n += 1
	print "END\n"
	m += 1
print "\nFinished printing the fnm data structs here.\n\n"
## NOTE: After print this -- run diff on the output and the original fnm file to see if there are any inconsistencies!!
'''
anm_file_handle.close()



##  Read in the binary contact network
bin_net_to_rd = sys.stdin
#bin_net_to_rd = open(bin_net, "r")
map_ln_indx___2___ints = {}
i = 0
for line in bin_net_to_rd:
	#print line
	ln_elms = list()
	ln_elms = line.split()
	map_ln_indx___2___ints[i] = ln_elms
	i += 1
num_residues_in_contact_network = i

##  Note -- we shouldn't even be able to reach this state in the program's 
##  execution if this is this case, but let's just check in case:
if num_residues_frm_CA_PDB_file != num_residues_in_contact_network:
	print "\n ERROR: num_residues_frm_CA_PDB_file != num_residues_in_contact_network \n"
	print "\t\t num_residues_frm_CA_PDB_file:  " + str(num_residues_frm_CA_PDB_file)
	print "\t\t num_residues_in_contact_network:  " + str(num_residues_in_contact_network)
	sys.exit()
'''
##  Print the binary contact net to ensure all was stored properly:
i = 0
while i < num_residues_in_contact_network:
	for elm in map_ln_indx___2___ints[i]:
		print elm,
	print ""
	i += 1
## NOTE: After print this -- run diff on the output and the original binary contact net to see if there are any inconsistencies!!
'''
#bin_net_to_rd.close()






## For each residue i, determine avg[(del_i)^2] :
node = 0
avg_sqr_vect = list()
sqr_vect = list()
while node < num_nodes_in_fnm_network:
	mode = 0
	tot_sq_vct = 0.0
	while mode < num_modes:
		x_component = float(modes_data[mode][0][node])
		y_component = float(modes_data[mode][1][node])
		z_component = float(modes_data[mode][2][node])
		sq_vct = (x_component * x_component)  +  (y_component * y_component)  +  (z_component * z_component)
		#print "x,y,z for mode " + str(mode) + " and node " + str(node) + " :  " + str(modes_data[mode][0][node]) + " " + str(modes_data[mode][1][node]) + " " + str(modes_data[mode][2][node]) + "   sq_vct: " + str(sq_vct)
		tot_sq_vct = tot_sq_vct + sq_vct
		mode += 1
	avg_s_v = float(tot_sq_vct) / float(num_modes)
	avg_sqr_vect.append(avg_s_v)
	node += 1
'''
##  As a quality check, print some output:
print "\n\n\n"
node = 0
while node < num_nodes_in_fnm_network:
	print "avg for node " + str(node) + " is :  " + str(avg_sqr_vect[node])
	node += 1
'''




def determine_weight(i,j):
	node_i = i
	node_j = j
	#print "nodes in pairs: i,j:  " + str(node_i) + " " + str(node_j)

	mode = 0
	sum_of_dot_products_accross_modes = 0.0
	while mode < num_modes:
		dot_prd = float(   (modes_data[mode][0][node_i] * modes_data[mode][0][node_j])  +  (modes_data[mode][1][node_i] * modes_data[mode][1][node_j])  +  (modes_data[mode][2][node_i] * modes_data[mode][2][node_j])   )
		#print "dot prod in mode " + str(mode) + " for node i and j " + str(node_i) + " " + str(node_j) + " is:  " + str(dot_prd)
		sum_of_dot_products_accross_modes = float(sum_of_dot_products_accross_modes) + float(dot_prd)
		mode += 1
	avg_dot_prod_ij = float(sum_of_dot_products_accross_modes) / float(num_modes)
	#print "av in node i and j " + str(node_i) + " " + str(node_j) + " is:  " + str(avg_dot_prod_ij)
	##  determine C(i,j) :
	C_ij = avg_dot_prod_ij / math.sqrt(avg_sqr_vect[node_i] * avg_sqr_vect[node_j])
	#wght = -1.0 * math.log(abs(C_ij))
	abs_C_ij = abs(C_ij)
	wght = -1.0 * math.log(abs_C_ij)
	#print "C_ij_here:  " + str(C_ij) + "   node i,j: " + str(node_i) + " " + str(node_j)
	#print "wght:       " + str(wght) + "\n"
	return wght




i = 0
map_ln_indx___2___WEIGHTS = {}
while i < num_residues_in_contact_network:
	weights_for_this_line = list()
	j = 0
	while j < num_residues_in_contact_network:
		if map_ln_indx___2___ints[i][j] == "1":  ## then these residues are in contact! Weight by correlated motions
			#print "i,j:  " + str(i) + " " + str(j) 
			## get weight
			weight = determine_weight(i,j)
			#print "weight for nodes " + str(i) + " " + str(j) + " is:  " + str(weight)
			#weight = 1
		else:
			weight = 0
		weights_for_this_line.append(abs(weight))
		j += 1
	map_ln_indx___2___WEIGHTS[i] = weights_for_this_line
	i += 1


i = 0
while i < num_residues_in_contact_network:
	weights_for_this_line = map_ln_indx___2___WEIGHTS[i]
	for w in weights_for_this_line:
		print w,
		#print w
	print ""
	i += 1


