#!/usr/bin/env python

import os
import sys
import subprocess

'''
Given an input file of GN-critical residues w/the the following format 
(along w/the SAME PDB which was originally used by gncommunities to generate this file):

	###  Column 1: <com_pair_index> is an index designating a specific pair of interacting COMMUNITIES
	###  Column 2: <gn_com_pair_A> is COMMUNITY A (community numbering from GN script)
	###  Column 3: <gn_com_pair_B> is COMMUNITY B (community numbering from GN script)
	###  Column 4: <node_a> is the node id (RESIDUE) in community A
	###  Column 5: <node_b> is the node id (RESIDUE) in community B
	###  Column 6: <betweenness> is 'betweenness' value between the residues a and b
	0	1	3	2	230	1141.597168
	1	1	4	11	13	15647.082031
	2	1	5	8	97	7155.361816
	3	1	7	50	57	12536.102539
	4	1	9	68	70	4875.646484
	5	2	4	373	375	25041.375000
	6	2	10	335	380	9975.492188


generate an output file listing the GN critical residues with the following format:

	###  Column 1: RESIDUE
	###  Column 2: RESID
	###  Column 3: RESTYPE
	###  Column 4: CHAIN
	106 115 ASN A
	108 117 TYR A
	112 121 GLU A
	13 22 ALA A
	143 152 LEU A

USAGE:
python residues__2__resids_types_chains_wo_VMD.py ./crit/1edh_crit_GN.txt ./1edh_crit_mappings.txt ./1EDH_SIMPLIFIED.pdb

'''


###  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 3):
	sys.stderr.write("Usage: " + sys.argv[0] + " <outputDir>\n")
	sys.stderr.write("where:\n")
	sys.stderr.write("   <output_file (see 'new' format described above)>  \n")
	sys.stderr.write("   <SIMPLIFIED PDB file>  \n")
	sys.exit()



##  Read in the SIMPLIDIED PDB and populate the following dictionaries:
residue___2___resID = {}
residue___2___resType = {}
residue___2___chain = {}
simplified_pdb_file = sys.argv[2]
simplified_pdb_file_to_rd = open(simplified_pdb_file, "r")
residue = -1
ln_indx = 0
for line in simplified_pdb_file_to_rd:
	atom_type = line[13:15]
	#print "atom_type:  |" + atom_type + "|"
	if atom_type == "CA":
		residue += 1
		resID = int(line[22:26])
		#print "resID: |" + str(resID) + "|"
		resType = line[17:20]
		#print "resType: |" + resType + "|"
		chain = line[21]
		#print "chain: |" + chain + "|"
		residue___2___resID[residue] = resID
		residue___2___resType[residue] = resType
		residue___2___chain[residue] = chain
simplified_pdb_file_to_rd.close()



##  Read in the *_crit_GN.txt file to get a list of all the critical residues
unique_list_of_critical_RESIDUES = list()
crit_GN_file_to_rd = sys.stdin
for line in crit_GN_file_to_rd:
	if "###" not in line:
		line_elmnts = list()
		line_elmnts = line.split()
		crit_res_A = line_elmnts[3]
		crit_res_B = line_elmnts[4]
		if crit_res_A not in unique_list_of_critical_RESIDUES:
			unique_list_of_critical_RESIDUES.append(crit_res_A)
		if crit_res_B not in unique_list_of_critical_RESIDUES:
			unique_list_of_critical_RESIDUES.append(crit_res_B)
'''
for residue in unique_list_of_critical_RESIDUES:
	print str(residue)
'''
unique_list_of_critical_RESIDUES_sorted = list()
unique_list_of_critical_RESIDUES_sorted = sorted(unique_list_of_critical_RESIDUES)


## Print output 
out_file = sys.argv[1]
out_file_to_wrt = open(out_file, "w")
out_file_to_wrt.write("###  Column 1: RESIDUE \n")
out_file_to_wrt.write("###  Column 2: RESID \n")
out_file_to_wrt.write("###  Column 3: RESTYPE \n")
out_file_to_wrt.write("###  Column 4: CHAIN \n")
for residue in unique_list_of_critical_RESIDUES_sorted:
	string_to_print = str(int(residue)) + " " + str(residue___2___resID[int(residue)]) + " " + residue___2___resType[int(residue)] + " " + residue___2___chain[int(residue)] + "\n"
	out_file_to_wrt.write(string_to_print)
out_file_to_wrt.close()


