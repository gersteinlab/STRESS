#!/usr/bin/env python

import os
import sys
import subprocess

'''
Given a directory containing GN_output files -- ex: see directory 
	/Users/admin/Desktop/rsch/allostery/networks/networkTools/MASS_PROCESSING/crit/
this script parses through each file, and extracts the data from columns 4 and 5 (which
respresent critical RESIDUES) and then maps those residues to the associated RESIDs, 
RESTYPEs, and CHAINs. For each file, the format is something like:

				###  Column 1: RESIDUE
				###  Column 2: RESID
				###  Column 3: RESTYPE
				###  Column 4: CHAIN
				106 115 ASN A
				108 117 TYR A
				112 121 GLU A
				13 22 ALA A
				143 152 LEU A
'''

###   Example usage:
###   python residues__2__resids_types_chains.py /Users/admin/Desktop/rsch/allostery/networks/networkTools/MASS_PROCESSING/crit_TRIAL/ /Users/admin/Desktop/rsch/allostery/networks/conservation/MASS_PROCESSING/residues__2__resids_types_chains/ /Users/admin/Desktop/rsch/allostery/networks/conservation/MASS_PROCESSING/SIMPLIFIED_PDBs/ /Users/admin/Desktop/rsch/allostery/networks/conservation/MASS_PROCESSING/TEMP_file_with_mappings.txt > v.txt

###   python residues__2__resids_types_chains.py /Users/admin/Desktop/rsch/allostery/networks/networkTools/MASS_PROCESSING/crit/ /Users/admin/Desktop/rsch/allostery/networks/conservation/MASS_PROCESSING/residues__2__resids_types_chains/ /Users/admin/Desktop/rsch/allostery/networks/conservation/MASS_PROCESSING/SIMPLIFIED_PDBs/ /Users/admin/Desktop/rsch/allostery/networks/conservation/MASS_PROCESSING/TEMP_file_with_mappings.txt > v.txt


###  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 5):
   sys.stderr.write("Usage: " + sys.argv[0] + " <inputFile> <outputDir>\n")
   sys.stderr.write("where:\n")
   sys.stderr.write("   <input_directory_GN_files>  \n")
   sys.stderr.write("   <output_directory>  \n")
   sys.stderr.write("   <directory_w_pdbs>  \n")
   sys.stderr.write("   <temp_file_for_mappings>  ... \n")
   sys.exit()


###  Assign input variables
input_directory = sys.argv[1]
output_directory = sys.argv[2]
directory_w_pdbs = sys.argv[3]
temp_file_for_mappings = sys.argv[4]




##  Create a temp tcl script to laod the pdb file into VMD  -- ie, generates the needed TEMP_load_root_pdb.tcl file
def create_temp_pdb_load_tcl(pdb_file, temp_file_for_mappings, contiguous_string_list_of_residues):
   temp_file = "./" + "TEMP_load_root_pdb.tcl"
   temp_tcl_r = open(temp_file, "w")
   temp_tcl_r.write(" \n \
      proc load_root_pdb {} { \n \
         set root_file \"" + pdb_file + "\" \n \
         mol new $root_file \n \
         display resetview \n \
      } \n \
      load_root_pdb \n \
      set filename \"" + temp_file_for_mappings + "\" \n \
      set fileId [open $filename \"w\"] \n \
      set contiguous_string_list_of_residues \"" + contiguous_string_list_of_residues + "\" \n \
      puts $fileId \"[lsort -unique [[atomselect top \"residue $contiguous_string_list_of_residues\"] get {residue resid resname chain}]]\" \n \
      close $fileId \n \
      "
   )
   temp_tcl_r.close()


input_files = os.listdir(input_directory)
for f in input_files:
	full_file_name = str(input_directory) + str(f)
	print str(full_file_name) + "\n\n\n\n\n"
	crit_input_file = open(full_file_name, "r")
	list_of_residues = list()
	for line in crit_input_file:
		if "#" not in line:
			line_contents = list()
			line_contents = line.split()
			list_of_residues.append(line_contents[3])
			list_of_residues.append(line_contents[4])
	list_of_residues = list(set(list_of_residues))  ##  This is simply making list unique

	##  For each residue in list_of_residues -- Use VMD to get the associated info
	domain_id = f[0:4]  # first get the domain name
	full_pdb_file_name = str(directory_w_pdbs) + str(domain_id) + "_SIMPLIFIED.pdb"
	#print "full_pdb_file_name:  " + full_pdb_file_name  ########################################################

	contiguous_string_list_of_residues = ""
	for r in list_of_residues:
		contiguous_string_list_of_residues = str(contiguous_string_list_of_residues) +  " " + r

	create_temp_pdb_load_tcl(full_pdb_file_name, temp_file_for_mappings, contiguous_string_list_of_residues)
	os.system("/Applications/VMD\ 1.9.1.app/Contents/vmd/vmd_MACOSXX86 -dispdev text -e load_pdbs_and_make_fasta.tcl")

	##  We now have a temp file containing the mappings of crtitical RESIDUES 
	##  to RESIDS/CHAINS/RESTYPES.  Open this temp file and save the contents
	temp_map = "./" + "TEMP_file_with_mappings.txt"
	temp_map_r = open(temp_map, "r")
	raw_map_data = temp_map_r.readline()
	raw_map_data_trim = raw_map_data[1:-2]
	residue_maps = raw_map_data_trim.split("} {")
	map_file = str(output_directory) + str(domain_id) + "_crit_mappings.txt"
	map_file_w = open(map_file, "w")
	map_file_w.write("###  Column 1: RESIDUE \n")
	map_file_w.write("###  Column 2: RESID \n")
	map_file_w.write("###  Column 3: RESTYPE \n")
	map_file_w.write("###  Column 4: CHAIN \n")
	for rm in residue_maps:
		map_file_w.write(str(rm) + "\n")

	temp_map_r.close()
	map_file_w.close()
	os.system("rm TEMP_load_root_pdb.tcl")
	os.system("rm TEMP_file_with_mappings.txt")
	crit_input_file.close()



print "\n\n\n  ---  complete  --- \n\n\n "
