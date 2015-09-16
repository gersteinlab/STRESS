import sys
import numpy as np
import os
import random
# import matplotlib.pyplot as plt
# import os.path


##  For each candidate site, this script assigns its corresponding processed BL site 
##  (ie, for each candidate binding site w/at least 1 residue in the raw output file,
##  this file determines the *BL.dat site to which it corresponds).


'''
example usage:
python map_raw_sites__to__BL_sites.py out_ba_2DW7__BL.dat out_ba_2DW7.txt 0.00005 > 2DW7.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_4PIW__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_4PIW.txt 0.00005 > 4PIW.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1UF8__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1UF8.txt 0.00005 > 1UF8.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_3ZWC__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_3ZWC.txt 0.00005 > 3ZWC.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1DZK__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1DZK.txt 0.00005 > 1DZK.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1AHH__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1AHH.txt 0.00005 > 1AHH.txt

python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1URZ__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1URZ.txt 0.000005 > 1URZ.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_3RPK__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_3RPK.txt 0.000005 > 3RPK.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1GD8__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1GD8.txt 0.000005 > 1GD8.txt
python map_raw_sites__to__BL_sites.py /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1I0C__BL.dat /Users/admin/Desktop/rsch/allostery/binding_leverage/chains_analyses/main_plots/site__vs__freq_of_occurrence/input_files_annot/out_ba_1I0C.txt 0.000005 > 1I0C.txt

'''


#  If the user did not enter the correct number of args, notify the user and terminate the program
if (len(sys.argv) != 4):
    print "\nYou provided " + str(len(sys.argv)) + " arguments even though 3 were expected!\n"
    sys.stderr.write("Usage: " + sys.argv[0] + " <binding_leverage_file> <x>\n")
    sys.stderr.write("where:\n")
    sys.stderr.write("   <binding_leverage_file> is a *BL.dat file\n")
    sys.stderr.write("   <raw_out_file_from_BL> for ex: out_ba_2DW7.txt \n")
    sys.stderr.write("   <incr_val> is a float for setting the bin size when making the resultant histogram in R \n")
    sys.exit()

###  Assign input variables
binding_leverage_file = sys.argv[1]  ## the *BL.dat binding leverage output file listing all the candidate sites
binding_leverage_file_to_read = open(binding_leverage_file, "r")
raw_BL_out_file = sys.argv[2]
raw_BL_out_file_to_read = open(raw_BL_out_file, "r")
incrmnt = float(sys.argv[3])


BL_dat_lines = list()
for line in binding_leverage_file_to_read:
    BL_dat_lines.append(line)

raw_out_lines = list()
for line in raw_BL_out_file_to_read:
    raw_out_lines.append(line)

binding_leverage_file_to_read.close()
raw_BL_out_file_to_read.close()



BL_site___2___list_of_resIDs_in_dat = {}
BL_site___2___BL_score = {}
i = 0
while i < len(BL_dat_lines):
    ln_elems = list()
    ln_elems = BL_dat_lines[i].split()
    res_ids_w_chain = ln_elems[2:len(ln_elems)]
    list_of_resIds = list()
    for r in res_ids_w_chain:
        r_elems = list()
        r_elems = r.split("_")
        list_of_resIds.append(int(r_elems[0]))
    BL_site___2___list_of_resIDs_in_dat[int(ln_elems[0])] = list_of_resIds
    BL_site___2___BL_score[int(ln_elems[0])] = float(ln_elems[1])
    i += 1

'''
for key in BL_site___2___list_of_resIDs_in_dat:
    print str(key) + "     " + str(BL_site___2___list_of_resIDs_in_dat[key])
'''




##  Build dictionary:  trial_index___2___resIDs
i = 0
trial_index___2___resIDs = {}
list_of_trial_IDs = list()
while i < len(raw_out_lines):
    if ". Number of residues in binding site " in raw_out_lines[i]:
        #print raw_out_lines[i]
        line_elems = list()
        line_elems = raw_out_lines[i].split()
        num_res = int(line_elems[10])
        if num_res > 0:
            ## Determine the trial_index
            #curr_trial_run_line = raw_out_lines[i-12]
            curr_trial_run_line = raw_out_lines[i-11]
            #curr_trial_run_line = raw_out_lines[i-11]
            trl_run_ln_elems = list()
            trl_run_ln_elems = curr_trial_run_line.split()
            trial_index = int(trl_run_ln_elems[2])
            list_of_trial_IDs.append(trial_index)
            ## Determine the associated list of resIDs
            residues_line = raw_out_lines[i+1]
            rs_line_elems = list()
            rs_line_elems = residues_line.split()
            rs_line_elems_ints = list()
            for rs in rs_line_elems:
                rs_line_elems_ints.append(int(rs))
            trial_index___2___resIDs[trial_index] = rs_line_elems_ints
            #print str(trial_index) + "     " + str(trial_index___2___resIDs[trial_index])
    i += 1




##  Build dictionaries:
##      trial_index___2___list_of_merged_trials   ## if trial_index is subsidiary (ie, if trial_index is itself merged w/an earlier trial), provide flag indicating such
i = 0
list_of_merged_pairs = list()
conglomerate_list_of_all_merged_trials = list()
while i < len(raw_out_lines):
    line_elems = list()
    line_elems = raw_out_lines[i].split()
    if (   ("Merging sites " in raw_out_lines[i])  and  (" & " in raw_out_lines[i])  and  (len(line_elems) == 5)   ):
        #print raw_out_lines[i]
        root_trial_index = line_elems[2]
        subsidiary_trial_index = line_elems[4]
        pair_list = list()
        pair_list.append(int(root_trial_index))
        pair_list.append(int(subsidiary_trial_index))
        conglomerate_list_of_all_merged_trials.append(int(root_trial_index))
        conglomerate_list_of_all_merged_trials.append(int(subsidiary_trial_index))
        list_of_merged_pairs.append(pair_list)
        #print str(pair_list)
    i += 1



trial_index___2___list_of_merged_trials = {}
list_of_ALL_subsidiary_trials = list()
for key in trial_index___2___resIDs:
    list_of_subsidiary_trials = list()
    for pair in list_of_merged_pairs:
        if key == pair[0]:  ## then key is a root trial
            list_of_subsidiary_trials.append(pair[1])
            list_of_ALL_subsidiary_trials.append(pair[1])
    trial_index___2___list_of_merged_trials[key] = list_of_subsidiary_trials



BL_site___2___list_of_trial_indeces = {}
i = 0  ## note that i designates the index of the BL sites here
already_assigned = list()
for trl_indx in list_of_trial_IDs:
    if trl_indx in conglomerate_list_of_all_merged_trials:
        #print str(trl_indx) + "    " + str(trial_index___2___list_of_merged_trials[trl_indx])
        if trl_indx not in already_assigned:
            for pair in list_of_merged_pairs:
                if trl_indx == pair[0]:  ## then this is a root trial
                    tmp_list2 = list()
                    tmp_list2.append(trl_indx)
                    for t in trial_index___2___list_of_merged_trials[trl_indx]:
                        tmp_list2.append(t)
                    BL_site___2___list_of_trial_indeces[i] = tmp_list2
                    for t in tmp_list2:
                        already_assigned.append(t)
                    i += 1
                    break
    else:  ## then this trial is all by itself in determining a BL site -- assign a 1-membered list
        if trl_indx not in list_of_ALL_subsidiary_trials:
            tmp_list = list()
            tmp_list.append(trl_indx)
            BL_site___2___list_of_trial_indeces[i] = tmp_list
            i += 1

'''
for key in BL_site___2___list_of_trial_indeces:
    print "i: " + str(key) + "    " + str(BL_site___2___list_of_trial_indeces[key])
'''


##  Quality control checks: If all assignments of BL sites to trial runs are correct, 
##  then it should be the case that the BL site resIDs constitute a subset of the resIDs 
##  in the union of merged trials 
for key in BL_site___2___list_of_trial_indeces:
    trial_index___2___resIDs[trial_index] = rs_line_elems
    ## get "BL_list" = list of all resIDs in the BL site
    BL_list = BL_site___2___list_of_resIDs_in_dat[key]

    ## get "resIDs_of_merged_trials" = list of all resIDs in the respective set of merged trials
    resIDs_of_merged_trials = list()
    for trial in BL_site___2___list_of_trial_indeces[key]:
        for r in trial_index___2___resIDs[trial]:
            resIDs_of_merged_trials.append(int(r))
    #print str(key) + "    " + str(BL_site___2___list_of_trial_indeces[key]) + "     " + str(BL_site___2___list_of_resIDs_in_dat[key]) + "    " + str(resIDs_of_merged_trials)        

    ## check to confirm that BL_list is indeed a subset of resIDs_of_merged_trials -- if not, print error and exit
    for b in BL_list:
        if b not in resIDs_of_merged_trials:
            print "\n\nERROR:  The resID " + str(b) + " of BL site " + str(key) + " is not contained in any of the respective merge sites -- indicating a false mapping of this BL site to trial run sites!!\n\n"
            sys.exit()



print "raw_trial_index    associated_BL_site_index    BL_score"
for trial in list_of_trial_IDs:
    for key in BL_site___2___list_of_trial_indeces:
        if trial in BL_site___2___list_of_trial_indeces[key]:  ## then this trial correpsonds to BL site key
            print str(trial) + "    " + str(key) + "    " + str(BL_site___2___BL_score[key])



##  Loop through ALL trials (958 for 2DW7), and for each, get the 
##  corresponding BL score through BL_site___2___list_of_trial_indeces.
##  We want to plot BL score  vs  freq of occurrence of that BL score
##  Make dictionary:  BL_score___2___freq_of_score

##  First get a list of all the BL scores
list_of_BL_scores = list()
for key in BL_site___2___BL_score:
    list_of_BL_scores.append(float(BL_site___2___BL_score[key]))
list_of_BL_scores_u = list(set(list_of_BL_scores))
'''
for scr in list_of_BL_scores_u:
    print str(scr)
'''





BL_score___2___list_of_BL_indeces = {}
for scr in list_of_BL_scores_u:
    list_of_BL_indeces = list()
    for bl_indx in BL_site___2___BL_score:    
        if BL_site___2___BL_score[bl_indx] == scr:
            list_of_BL_indeces.append(bl_indx)
    BL_score___2___list_of_BL_indeces[scr] = list_of_BL_indeces
'''
for scr in BL_score___2___list_of_BL_indeces:
    print str(scr) + "    " + str(len(BL_score___2___list_of_BL_indeces[scr]))
'''



BL_score___2___freq_of_score = {}
for scr in list_of_BL_scores_u:
    num_matches = 0
    for bl_indx in BL_score___2___list_of_BL_indeces[scr]:
        for trial in list_of_trial_IDs:
            if trial in BL_site___2___list_of_trial_indeces[bl_indx]:
                num_matches += 1
    BL_score___2___freq_of_score[scr] = num_matches
'''
for scr in BL_score___2___freq_of_score:
    print str(scr) + "    " + str(BL_score___2___freq_of_score[scr])
'''






##  Below is for generating the bins to make a histogram in Excel or R
'''
print "\n\n\n\n\n"
brk_strng = "brk <- c("
i = 0.0
while i <= (max(list_of_BL_scores_u) + incrmnt):
    if i == 0.0:
        brk_strng = brk_strng + str(i)
    else:
        brk_strng = brk_strng + ", " + str(i)
    i += incrmnt

print brk_strng
'''


