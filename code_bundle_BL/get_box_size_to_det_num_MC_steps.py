import sys
import numpy as np
import os
#from math import exp
#from Bio.PDB import *


'''
This script takes a PDB structure as input, and reports on the num of MC steps 
that should be set when running the BL calculations on that pdb according to:

Determine number of MC steps in BL simulations:
	From Berezovsky et al: 
	"We produced probe locations with the number of simulations set to 10 times the number of residues, 
	and the number of MC steps to 1000 times the size of the simulation box measured in Angstroms (the box 
	size is set according to the Methods section). The number of atoms in the probe (probe size) 
	was set to 4 universally...The boundary conditions are periodic and the size of the cubic simulation box is set to twice 
	the maximum size of the protein along any of the x, y or z-axes."

	--> Do 10x the number of MC steps previously used by Berezovsky et al.

'''

##  This script is stored in the directory:
##  /Users/admin/Desktop/rsch/allostery/networks/make_contact_map/
##  
##  usage:
##  python get_box_size_to_det_num_MC_steps.py /Users/admin/Desktop/rsch/allostery/stamp_aligns_seq/scripts/berezovsky_exs/whole_pdbs/clusts_2/1j3h_1atp/pdbs/raw_pdbs/1J3H_TRL.pdb
##  python get_box_size_to_det_num_MC_steps.py SIMPLIDIED_pdbs/106M_SIMPLIFIED.pdb


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 2):
   sys.stderr.write("Usage: " + sys.argv[0] + " <input_pdb> <output_pdb>  \n")
   sys.stderr.write("where:\n")
   sys.stderr.write("   <input_pdb> is a PDB file \n")
   sys.exit()



###  Assign input variable:
pdb_file = sys.argv[1]
pdb_file_to_read = open(pdb_file, "r")


x_coords = list()
y_coords = list()
z_coords = list()
for line in pdb_file_to_read:
	if line[0:4] == "ATOM":
		x_val = line[30:38]
		y_val = line[38:46]
		z_val = line[46:54]

		##  Perform a basic sanity check to be sure that these values look like real coordinates 
		##  (ie, that we don't have a weird PDB w/coordinates at the wrong indeces w/in a line!)
		##  	> Each coordinate must consist of ONE "."
		##  	> Each period must be followed by exactly 3 numeric values 
		##  	> There must be only one, two, or three numeric chars BEFORE the period
		if (   (x_val.count('.') != 1)  or  (y_val.count('.') != 1)  or  (z_val.count('.') != 1)   ):
			print "\n\n\n\nERROR IN RECORDING COORDINATE W/WEIRD PERIODS! -- SEE LINE: \n" + line
			sys.exit()
		else:
			##  Each period must be followed by exactly 3 numeric values:
			x_val_elems = list()
			x_val_elems = x_val.split('.')
			num_numer_chars = 0
			for ch in x_val_elems[1]:
				if (   (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					num_numer_chars += 1
			if num_numer_chars != 3:
				print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
				sys.exit()

			y_val_elems = list()
			y_val_elems = y_val.split('.')
			num_numer_chars = 0
			for ch in y_val_elems[1]:
				if (   (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					num_numer_chars += 1
			if num_numer_chars != 3:
				print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
				sys.exit()

			z_val_elems = list()
			z_val_elems = z_val.split('.')
			num_numer_chars = 0
			for ch in z_val_elems[1]:
				if (   (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					num_numer_chars += 1
			if num_numer_chars != 3:
				print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
				sys.exit()


			##  There must be only one, two, or three numeric chars BEFORE the period, and all chars before period must be either a space, a negative sign, or a numeric digit
			num_numer_chars = 0
			for ch in x_val_elems[0]:
				if (   (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					num_numer_chars += 1
				## must be either a numeric char or a negative sign or a space
				if not (   (ch == " ") or (ch == "-") or (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
					sys.exit()
			if not (  (num_numer_chars == 1) or (num_numer_chars == 2) or (num_numer_chars == 3)  ):
				print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
				sys.exit()

			num_numer_chars = 0
			for ch in y_val_elems[0]:
				if (   (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					num_numer_chars += 1
				## must be either a numeric char or a negative sign or a space
				if not (   (ch == " ") or (ch == "-") or (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
					sys.exit()
			if not (  (num_numer_chars == 1) or (num_numer_chars == 2) or (num_numer_chars == 3)  ):
				print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
				sys.exit()

			num_numer_chars = 0
			for ch in z_val_elems[0]:
				if (   (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					num_numer_chars += 1
				## must be either a numeric char or a negative sign or a space
				if not (   (ch == " ") or (ch == "-") or (ch == "1") or (ch == "2") or (ch == "3") or (ch == "4") or (ch == "5") or (ch == "6") or (ch == "7") or (ch == "8") or (ch == "9") or (ch == "0")   ):
					print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
					sys.exit()
			if not (  (num_numer_chars == 1) or (num_numer_chars == 2) or (num_numer_chars == 3)  ):
				print "\n\n\n\nERROR IN RECORDING COORDINATE! -- SEE LINE: \n" + line
				sys.exit()

		x_val_numeric = float(x_val)
		y_val_numeric = float(y_val)
		z_val_numeric = float(z_val)
		'''
		print line
		print "|" + str(x_val_numeric) + "|        |" + str(y_val_numeric) + "|        |" + str(z_val_numeric) + "|"
		print "\n\n\n\n"
		'''
		x_coords.append(x_val_numeric)
		y_coords.append(y_val_numeric)
		z_coords.append(z_val_numeric)


x_range = max(x_coords) - min(x_coords)
y_range = max(y_coords) - min(y_coords)
z_range = max(z_coords) - min(z_coords)

ranges = list()
ranges.append(x_range)
ranges.append(y_range)
ranges.append(z_range)

max_range = max(ranges)

##  Berezovsky: "...the size of the cubic simulation box is set to twice 
##  the maximum size of the protein along any of the x, y or z-axes"
max_range_by_2 = max_range * 2.0

final_num_mc_steps_to_use = int(max_range_by_2) * 1000 * 10 ##  Berezovsky: "...the number of MC steps to 1000 times the size of the simulation box measured in Angstroms"


'''
xmin = min(x_coords)
xmax = max(x_coords)

ymin = min(y_coords)
ymax = max(y_coords)

zmin = min(z_coords)
zmax = max(z_coords)

print "\nxmin:  " + str(xmin)
print "xmax:  " + str(xmax)
print "x_range:  " + str(x_range) + "\n"

print "\nymin:  " + str(ymin)
print "ymax:  " + str(ymax)
print "y_range:  " + str(y_range) + "\n"

print "\nzmin:  " + str(zmin)
print "zmax:  " + str(zmax)
print "z_range:  " + str(z_range) + "\n"

print "\nranges:  " + str(ranges) + "\n"
print "\nmax_range:  " + str(max_range) + "\n"
print "\nmax_range_by_2:  " + str(max_range_by_2) + "\n"

print "\nfinal_num_mc_steps_to_use:  " + str(final_num_mc_steps_to_use) + "\n"

print "\npdb_file:  " + str(pdb_file) + "\n"
'''

print str(pdb_file) + "\t" + str(final_num_mc_steps_to_use)


pdb_file_to_read.close()



