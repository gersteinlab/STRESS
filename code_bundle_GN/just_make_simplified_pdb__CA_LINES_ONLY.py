import numpy as np
import sys
import os
#from math import exp
#from Bio.PDB import *


'''
This script takes a PDB structure as input, and generates a simplified PDB
that does not contain residues such as "AGLY/BGLY", and it converts modified
residues to their appropriate names --- for ex:

	HETATM ... 	SEP ....

	 gets changed to:

	ATOM ...  	SER  ...

'''


##  This script is stored in the directory:
##  /Users/admin/Desktop/rsch/allostery/networks/make_contact_map/
##  
##  usage:
##  python just_make_simplified_pdb__CA_LINES_ONLY.py ./SIMPLIDIED_pdbs/3W00_SIMPLIFIED.pdb ./SIMPLIDIED_pdbs_CA_ONLY/3W00_CA.pdb


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 3):
   sys.stderr.write("Usage: " + sys.argv[0] + " <input_pdb> <output_pdb>  \n")
   sys.stderr.write("where:\n")
   sys.stderr.write("   <input_pdb> is a PDB file used to generate the contact map\n")
   sys.stderr.write("   <output_pdb> is the output simplified PDB name\n")
   sys.exit()


###  Assign input variable:
pdb_file = sys.argv[1]
pdb_file_to_read = open(pdb_file, "r")
simplified_pdb_out = sys.argv[2]
simplified_pdb_out_print = open(simplified_pdb_out, "w")



###  Read and store coordinates from the binding_leverage_file_to_read
l = 0
res_id_list = list()
raw_res_index_list = list()
list_of_possible_modified_residue_names = ["PTR", "PSR", "SEP", "PTH", "PHY", "PHS", "PHT", "MSE", "MSO", "DAL", "NAL", "02K", "HR7", "2MR", "ASA", "PHD", "CSO", "SCH", "OCY", "SCH", "CMT", "CSS", "CAS", "CME", "CSX", "OCS", "CSW", "03Y", "PCA", "CGU", "EME", "B3E", "SAR", "56A", "MK8", "NLE", "MLE", "ALY", "LGY", "BTK", "PRK", "MLY", "LLP", "M3L", "MLZ", "PRS", "MEA", "0BN", "3PX", "HYP", "TPO", "TTQ", "HT7", "PTR", "TYS", "LE1"]
for line in pdb_file_to_read:
	##  If this is a "HETATM" is that is simply a modified residue (such 
	##  as a phosphorylated TYR), then it's fine, so keep such a line
	atom_type = line[13:16]
	if line[0:6] == "HETATM":
		res_name = line[17:20]
		is_just_modified_residue = "no"
		if res_name in list_of_possible_modified_residue_names:
			is_just_modified_residue = "yes"
		if is_just_modified_residue == "yes":
			##  Be sure residue number looks normal -- Normal is something like "520". NON-standard is something such as "520A"
			res_id_pre = line[22:27]
			is_standard_res_number = "yes"
			for c in res_id_pre:
				if c.isalpha():
					is_standard_res_number = "no"
			
			##  Be sure res name looks normal -- Normal is something like "ARG". NON-standard is something 
			##  such as "AARG" or "BARG" -- ie: non-normal is anything w/more than 3 alphanumeric chars
			res_name_pre = line[16:20]
			is_normal_res_name = "yes"
			num_alpha_numeric_chars = 0
			for c in res_name_pre:
				if c.isalnum():
					num_alpha_numeric_chars += 1
			if num_alpha_numeric_chars > 3:  ## then we may have something such as "BARG"
				print "NON-standard res_name_pre here:  |" + str(res_name_pre) + "|\n" + str(line)
				#print "here:  " + res_name_pre[0:1] + "\n\n\n"
				if res_name_pre[0:1] is not "A":  ## ie -- take only the lines w/"AARG", and not "BARG"
					is_normal_res_name = "no"
				else: ## then we have "AARG"
					res_name_rev = " " + str(res_name_pre[1:4])
					print "res_name_rev:  |" + str(res_name_rev) + "|"
					res_name = res_name_rev
			else:
				res_name = res_name_pre

			if (is_standard_res_number == "yes")  and  (is_normal_res_name == "yes"):
				if atom_type == "CA ":
					simplified_pdb_out_print.write(str(line[0:16]) + str(res_name) + str(line[20:-1]) + "\n")
				l += 1

		else:
			if str(res_name) != "HOH":
				print "\n Modified residue name:  " + res_name
				print line

	elif line[0:4] == "ATOM":
		##  Be sure residue number looks normal -- Normal is something like "520". NON-standard is something such as "520A"
		res_id_pre = line[22:27]
		is_standard_res_number = "yes"
		for c in res_id_pre:
			if c.isalpha():
				is_standard_res_number = "no"
		
		##  Be sure res name looks normal -- Normal is something like "ARG". NON-standard is something 
		##  such as "AARG" or "BARG" -- ie: non-normal is anything w/more than 3 alphanumeric chars
		res_name_pre = line[16:20]
		is_normal_res_name = "yes"
		num_alpha_numeric_chars = 0
		for c in res_name_pre:
			if c.isalnum():
				num_alpha_numeric_chars += 1
		if num_alpha_numeric_chars > 3:  ## then we may have something such as "BARG"
			print "NON-standard res_name_pre here:  |" + str(res_name_pre) + "|\n" + str(line)
			#print "here:  " + res_name_pre[0:1] + "\n\n\n"
			if res_name_pre[0:1] is not "A":  ## ie -- take only the lines w/"AARG", and not "BARG"
				is_normal_res_name = "no"
			else: ## then we have "AARG"
				res_name_rev = " " + str(res_name_pre[1:4])
				print "res_name_rev:  |" + str(res_name_rev) + "|"
				res_name = res_name_rev
		else:
			res_name = res_name_pre

		if (is_standard_res_number == "yes")  and  (is_normal_res_name == "yes"):
			if atom_type == "CA ":
				simplified_pdb_out_print.write(str(line[0:16]) + str(res_name) + str(line[20:-1]) + "\n")
			l += 1
		else:
			print "NON-standard full line for pdb file:  " + str(pdb_file) + "   \n" + str(res_id_pre) + "\n" + str(line) + "\n"

##  Close files
pdb_file_to_read.close()
simplified_pdb_out_print.close()

