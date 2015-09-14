import numpy as np
import sys
import os
#from math import exp
#from Bio.PDB import *

'''
This script takes a PDB structure as input, and generates an output 
file which is a matrix to represent contacts between residues.
'''

##  This script is stored in the directory:
##  /Users/admin/Desktop/rsch/allostery/networks/make_contact_map/

##  Usage:
'''
python /Users/admin/Desktop/rsch/allostery/networks/make_contact_map/make_contact_map.py 1X6X_SIMPLIFIED.pdb heavy 4.5 1 1X6X_cont_map.txt 1X6X_CA.pdb
'''


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 6):
   sys.stderr.write("Usage: " + sys.argv[0] + " <input_pdb> <heavy/ca> <distance_cutoff> \n")
   sys.stderr.write("where:\n")
   sys.stderr.write("   <input_pdb> is a PDB file used to generate the contact map\n")
   sys.stderr.write("   <heavy> or <ca> indicates whether the contacts should be defined in terms of all heavy atoms or CA only atoms\n")
   sys.stderr.write("   <distance_cutoff> is the distance cutoff the define contact \n")
   sys.stderr.write("   <num_nearest_neighbors_to_exclude> is the number of nearest sequence neighbors (on each side) that should be ignored as trivial contacts in the calculations (think a val of 1 was used for Dynamical networks in tRNA:protein complexes (Sethi et al, 2009) ) \n")
   sys.stderr.write("   <CA_only_file> is the pdb f/above, but w/the CA-lines only -- this is is provided to enforce that the residues considred match those in the corresponding FNM files \n")
   sys.exit()



###  Assign input variable:
pdb_file = sys.argv[1]
pdb_file_to_read = open(pdb_file, "r")
at_types_to_consider = sys.argv[2]
distance_cutoff = float(sys.argv[3])
num_nearest_neighbors_to_exclude = int(sys.argv[4])
ca_file = sys.argv[5]
ca_file_to_read = open(ca_file, "r")



## Determine the residues which actually have alpha carbon atoms
list_of_residues_with_alpha_carbons = list()
for line in ca_file_to_read:
	ca_res_id = line[22:27]
	ca_chain = line[21]
	ca_res_name = line[17:20]
	residue_data =list()
	residue_data.append(ca_res_id)
	residue_data.append(ca_chain)
	residue_data.append(ca_res_name)
	list_of_residues_with_alpha_carbons.append(residue_data)
	#print line
	#print str(ca_res_id) + "    " + str(ca_chain) + "    " + str(ca_res_name)
	#print "\n\n\n"



###  Read and store coordinates
l = 0
residue_list = list()
res_id_list = list()
res_name_list = list()
chain_list = list()
x_coord = list()
y_coord = list()
z_coord = list()
list_of_possible_residue_names = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "PTR", "PSR", "SEP", "PTH", "PHY", "PHS", "PHT", "MSE", "MSO", "DAL", "NAL", "02K", "HR7", "2MR", "ASA", "PHD", "CSO", "SCH", "OCY", "SCH", "CMT", "CSS", "CAS", "CME", "CSX", "OCS", "CSW", "03Y", "PCA", "CGU", "EME", "B3E", "SAR", "56A", "MK8", "NLE", "MLE", "ALY", "LGY", "BTK", "PRK", "MLY", "LLP", "M3L", "MLZ", "PRS", "MEA", "0BN", "3PX", "HYP", "TPO", "TTQ", "HT7", "PTR", "TYS", "LE1", "UNK"]
for line in pdb_file_to_read:
	##  If this is a "HETATM" is that is simply a modified residue (such 
	##  as a phosphorylated TYR), then it's fine, so keep such a line
	if (line[0:4] == "ATOM")   or   (line[0:6] == "HETATM"):
		res_NAME = line[17:20]  ## ex: "ARG"
		atom_type = line[13:16] ## ex: "CA "
		is_a_residue = "no"
		if (res_NAME in list_of_possible_residue_names):
			## Be sure that there's a CA associated w/this residue!!
			rca = 0
			while rca < len(list_of_residues_with_alpha_carbons):
				if (   (line[22:27] == list_of_residues_with_alpha_carbons[rca][0])  and  (line[21] == list_of_residues_with_alpha_carbons[rca][1])  and  (line[17:20] == list_of_residues_with_alpha_carbons[rca][2])   ):
				#if (line[22:27] == list_of_residues_with_alpha_carbons[rca][0])  and  ():
					is_a_residue = "yes"
					rca = 1000000  # just a ghetto way to terminate the loop when the condition is met...
				rca += 1
		elif atom_type == "CA ":  ## If there's an alpha carbon atom, then this is in fact a residue -- add to list and make note of new residue type!
			is_a_residue = "yes"
			#print "\n\n New modified residue name:  " + res_NAME + "\n\n"
			#print line + "\n\n"
		if is_a_residue == "yes":
			if str(at_types_to_consider) == "heavy":
				##  First be sure resid looks normal -- Standard is something like "520". NON-standard is something such as "520A"
				res_id_pre = line[22:27]
				is_standard_res_number = "yes"
				for c in res_id_pre:
					if c.isalpha():
						is_standard_res_number = "no"
				
				##  Be sure res name looks normal -- Standard is something like "ARG". NON-standard is something 
				##  such as "AARG" or "BARG" -- ie: non-standard is anything w/more than 3 alphanumeric chars
				res_NAME_pre = line[16:20]
				is_normal_res_NAME = "yes"
				num_alpha_numeric_chars = 0
				for c in res_NAME_pre:
					if c.isalnum():
						num_alpha_numeric_chars += 1
				if num_alpha_numeric_chars > 3:  ## then we may have something such as "BARG"
					#print "NON-standard res_NAME_pre here:  |" + str(res_NAME_pre) + "|\n" + str(line)
					if res_NAME_pre[0:1] is not "A":  ## ie -- take only the lines w/"AARG", and not "BARG"
						is_normal_res_NAME = "no"
					else:  ## Then we have "AARG"
						res_NAME_rev = " " + str(res_NAME_pre[1:4])
						#print "res_NAME_rev:  |" + str(res_NAME_rev) + "|"
						res_NAME = res_NAME_rev
				else:
					res_NAME = res_NAME_pre
				if (is_standard_res_number == "yes")  and  (is_normal_res_NAME == "yes"):
					x = float(line[30:38])
					y = float(line[38:46])
					z = float(line[46:54])
					x_coord.append(x)
					y_coord.append(y)
					z_coord.append(z)
					res_id = int(line[22:27])
					res_id_list.append(res_id)
					res_name_list.append(res_NAME)
					chain = line[21]
					##print "chain is:  |" + str(chain) + "|"
					#if not chain.isalpha(): -- SOME CHAINS CAN BE NUMERIC
					#	print "\n\nThis chain is not a standard chain!  See line: "
					#	print line + "\n\n"
					#	quit()
					chain_list.append(chain)
					if len(res_id_list) > 1:
						if res_id_list[l] != res_id_list[l-1]:
							raw_res_index += 1
					else:
						raw_res_index = 0
					residue_list.append(raw_res_index)
					l += 1

##  Print data to see if it was stored as intended!
'''
print "\n\n\n\nStart new print segment here \n"
i = 0
while i < len(x_coord):
	print "residue: " + str(residue_list[i]) + "    res_id:" + str(res_id_list[i]) + "    chain: " + str(chain_list[i]) + "    x: " + str(x_coord[i]) + "    y: " + str(y_coord[i]) + "    z: " + str(z_coord[i]) + "\n"
	i += 1
'''





###  Build data structs for pdb
l = 0
res_data = []
while l < len(residue_list):  ## Really just iterating through all ATOM/HETATM lines w/residues here
	if l == 0:  ## Here, we're at the very first line, so make new residue struct
		rd = 0
		x_crd = list()
		y_crd = list()
		z_crd = list()
		res_data.append([])
	elif residue_list[l] != residue_list[l-1]: ## New residue reached, so make new residue struct
		## New residue reached, so first store all data from prev residue:
		res_data[rd].append(x_crd)
		res_data[rd].append(y_crd)
		res_data[rd].append(z_crd)
		res_data[rd].append(residue_list[l-1])
		res_data[rd].append(res_id_list[l-1])
		res_data[rd].append(res_name_list[l-1])
		res_data[rd].append(chain_list[l-1])

		# update rd index, empty coord lists, and append new [] to res_data:
		rd += 1
		x_crd = list()
		y_crd = list()
		z_crd = list()
		res_data.append([])

	x_crd.append(x_coord[l])
	y_crd.append(y_coord[l])
	z_crd.append(z_coord[l])

	l += 1
## At very end -- save data for last residue:
res_data[rd].append(x_crd)
res_data[rd].append(y_crd)
res_data[rd].append(z_crd)
res_data[rd].append(residue_list[l-1])
res_data[rd].append(res_id_list[l-1])
res_data[rd].append(res_name_list[l-1])
res_data[rd].append(chain_list[l-1])



## print statements to be sure that all was stored as intended
'''
print "\nprnt st here...\n"
r = 0
while r < len(res_data):
	print "residue: " + str(res_data[r][3])
	print "res_id: " + str(res_data[r][4])
	print "res_name: " + str(res_data[r][5])
	print "chain: " + str(res_data[r][6])
	print "x vals: " + str(res_data[r][0])
	print "y vals: " + str(res_data[r][1])
	print "z vals: " + str(res_data[r][2])
	e = 0
	while e < len(res_data[r][0]):
		print "x: " + str(res_data[r][0][e]) + "   y: " + str(res_data[r][1][e]) + "   z: " + str(res_data[r][2][e])
		e += 1
	print "\n\n"
	r += 1
'''


## Determine which pairs of residues are in contact:
contact_status = []
num_contacts = 0
i = 0
while i < len(res_data):
	contact_status.append([])
	j = 0
	while j < len(res_data):
		in_contact = 0  # by default, assume residues i & j are not in contact
		atom_i = 0
		while atom_i < len(res_data[i][0]):  ## For each atom in residue i
			atom_j = 0
			while atom_j < len(res_data[j][0]):  ## For each atom in residue j
				distance = (    (res_data[i][0][atom_i] - res_data[j][0][atom_j])**2  +  (res_data[i][1][atom_i] - res_data[j][1][atom_j])**2  +  (res_data[i][2][atom_i] - res_data[j][2][atom_j])**2    )**(0.5)
				#print "distance: " + str(distance) + "       cutoff: " + str(distance_cutoff) + "\n"
				diff_btwn_residues_in_raw_sequence = abs(res_data[i][4] - res_data[j][4])
				if (   (distance <= distance_cutoff)   and   (diff_btwn_residues_in_raw_sequence > num_nearest_neighbors_to_exclude)    ):
					in_contact = 1
					#print "HERE: distance = " + str(distance) + "   distance_cutoff: " + str(distance_cutoff) + "   diff_in_seq: " + str(diff_btwn_residues_in_raw_sequence) + "   res_i: " + str(res_data[i][3]) + "  atom_i:  " + str(atom_i) + "   res_j: " + str(res_data[j][3]) + "  atom_j:  " + str(atom_j) + "\n"
					#print "CONTACT HERE:  " + "res_i: " + str(res_data[i][3]) + "    res_j: " + str(res_data[j][3]) + "\n"
					#break  #########################################################################################################################################################################################################################################################################################################
				atom_j += 1
			atom_i += 1

		if in_contact == 1:
			num_contacts += 1
		contact_status[i].append(in_contact)
		j += 1
	i += 1

#print "\n\nlen cont statuts = num pairs = " + str(len(contact_status)) + "\n\n"



r = 0
while r < len(res_data):
	##  For each line in cont_map, the beginning of the line should contain info about the specific node.
	##  Thus, print the RESIDUE, RESID, res_NAME, and CHAIN (ex: "1235 336 PHE B")
	sys.stdout.write(str(res_data[r][3]) + " "  + str(res_data[r][4]) + " " + str(res_data[r][5]) + " " + str(res_data[r][6]) + "      ",)
	j = 0
	while j < len(res_data):
		sys.stdout.write(str(contact_status[r][j]) + " ",)
		j += 1
	sys.stdout.write("\n")
	r += 1



##  Close files
pdb_file_to_read.close()
ca_file_to_read.close()
