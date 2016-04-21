import sys
import numpy as np
import os
import random
# import matplotlib.pyplot as plt
# import os.path


##  Performs gap calculations: For each of the N candidate binding sites (for example, N = 958 for 2DW7),
##  calculate delta_BL(i)/DELTA_BL, where i is the LINE index in the *BL.dat file (not the BL site index)
##  and DELTA_BL is a constant which equals the (maximum BL score) - (minimum BL score) in the *BL.dat file
##  (note that, since the min score is generally 0.000000, DELTA_BL is generally just the BL score of the 
##  top-scoring BL site).
##  For example, according to 2DW7.txt, the raw trial index 150 hits the binding site that corresponds to 
##  the site in the *BL.dat file that has an index of 20. The BL site with index 20 occurs as the 278'th 
##  line in the *BL.dat file. The value delta_BL(i), with i=150 in this example, is:
##  [BL score associated with the 278th line in the *BL.dat file] - [BL score associated with the 279th lin in the *BL.dat file]
##  = 0.000002 - 0.000002 = 0. Thus, here, delta_BL(278) = 0. Also, DELTA_BL = 0.002004 - 0.000000 = 0.002004

'''
example usage:
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_2DW7__BL.dat 2DW7.txt 30 > 2DW7__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_4PIW__BL.dat 4PIW.txt 30 > 4PIW__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_1UF8__BL.dat 1UF8.txt 30 > 1UF8__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_3ZWC__BL.dat 3ZWC.txt 30 > 3ZWC__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_1DZK__BL.dat 1DZK.txt 30 > 1DZK__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_1AHH__BL.dat 1AHH.txt 30 > 1AHH__2.txt

python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_1URZ__BL.dat 1URZ.txt 30 > 1URZ__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_3RPK__BL.dat 3RPK.txt 30 > 3RPK__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_1GD8__BL.dat 1GD8.txt 30 > 1GD8__2.txt
python gap_theory_calculations.py ../../site__vs__freq_of_occurrence/input_files_annot/out_ba_1I0C__BL.dat 1I0C.txt 30 > 1I0C__2.txt

'''


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 3):
    print "\nYou provided " + str(len(sys.argv)-1) + " arguments even though 2 were expected!\n"
    sys.stderr.write("Usage: " + sys.argv[0] + " <binding_leverage_file> <num_bins>\n")
    sys.stderr.write("where:\n")
    sys.stderr.write("   <binding_leverage_file> is a *BL.dat file listing the processed binding sites, along w/their BL scores\n")
    sys.stderr.write("   <site_mappings_file> is the output file generated from map_raw_sites__to__BL_sites.py (from STDIN) \n")
    sys.stderr.write("   <num_bins> is a value specifying the number of bins in the histogram to be plotted in R or excel \n")
    sys.exit()

###  Assign input variables
binding_leverage_file = sys.argv[1]  ## the *BL.dat binding leverage output file listing all the candidate sites
binding_leverage_file_to_read = open(binding_leverage_file, "r")
#site_mappings_file = sys.argv[2]
#site_mappings_file_to_read = open(site_mappings_file, "r")
#num_bins = float(sys.argv[3])

site_mappings_file_to_read = sys.stdin
num_bins = float(sys.argv[2])

BL_dat_lines = list()
for line in binding_leverage_file_to_read:
    BL_dat_lines.append(line)

mappings_lines = list()
for line in site_mappings_file_to_read:
    mappings_lines.append(line)

binding_leverage_file_to_read.close()
site_mappings_file_to_read.close()




BL_index___2___BL_line_index = {}
BL_line_index___2___BL_score = {}
list_of_BL_scores = list()
i = 0
while i < len(BL_dat_lines):
    ln_elems = list()
    ln_elems = BL_dat_lines[i].split()
    list_of_BL_scores.append(float(ln_elems[1]))
    BL_index___2___BL_line_index[ln_elems[0]] = i
    BL_line_index___2___BL_score[i] = float(ln_elems[1])
    i += 1

DELTA_BL = float(max(list_of_BL_scores) - min(list_of_BL_scores))

'''
for key in BL_index___2___BL_line_index:
    print str(key) + "    " + str(BL_index___2___BL_line_index[key])

for key in BL_line_index___2___BL_score:
    print str(key) + "    " + str(BL_line_index___2___BL_score[key])
'''



i = 0
ratios = list()
while i < len(mappings_lines):  ## we're disregarding the last line in the 
    if "raw_trial_index" not in mappings_lines[i]:  ## ie, if this is not the header line
        ln_elems = list()
        ln_elems = mappings_lines[i].split()

        ## Get delta_BL(i)
        bl_index = ln_elems[1]
        bl_line_index = BL_index___2___BL_line_index[bl_index]
        if bl_line_index < (  len(BL_dat_lines) - 1  ): ## ie, if we're not at the very last line in *BL.dat
            delta_BL_i = BL_line_index___2___BL_score[bl_line_index] - BL_line_index___2___BL_score[bl_line_index+1]
            if DELTA_BL != 0:
                ratio_i = delta_BL_i / DELTA_BL
            else:
                ratio_i = 0
            ratios.append(ratio_i)
            #print "ratio  " + str(i) + "  " + str(ratio_i)
            #print "ratio  " + str(ln_elems[1]) + "  " + str(ratio_i)
            print "ratio  " + str(BL_index___2___BL_line_index[ln_elems[1]]) + "  " + str(ratio_i)

            '''
            print "mappings_line:   " + str(mappings_lines[i])
            print "bl_index:  " + str(bl_index)
            print "BL_line_index:  " + str(BL_index___2___BL_line_index[bl_index])
            print "BL_line_index___2___BL_score[bl_line_index]:  " + str(BL_line_index___2___BL_score[bl_line_index])
            print "BL_line_index___2___BL_score[bl_line_index+1]:  " + str(BL_line_index___2___BL_score[bl_line_index+1])
            print "delta_BL_i:  " + str(delta_BL_i)
            print "DELTA_BL:  " + str(DELTA_BL)
            print "ratio_i:  " + str(ratio_i)
            print "\n\n\n\n"
            '''
    i += 1







'''
##  Below is for generating the bins to make a histogram in Excel or R
print "\n\n\n\n\n"

ratios_range = max(ratios) - min(ratios)
incrmnt = float(ratios_range) / num_bins

brk_strng = "brk <- c("
i = 0.0
while i <= (ratios_range):
    if i == 0.0:
        brk_strng = brk_strng + str(i)
    else:
        brk_strng = brk_strng + ", " + str(i)
    i += incrmnt

print brk_strng
'''

