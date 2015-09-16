#!/bin/bash
## This is the pipeline for the BL and GN modules
## Notice it will generate tmp files 
## Make sure the file names are unique
## Usage: pipeline_both.sh [INPUTPDB] [OUTPUTDIR] [JOBID] [BLNEEDED] [GNNEEDED]
## BL final output in $ID__top_sites_BL.dat
## GN final output in $ID_crit_mappings.txt

## EXIT CODE (binary)
## 00000 => Success
## 00001 => FNM number < 10, protein too small
## 00010 => PDB file checking failed
## 00100 => BL pipeline failed: PDB file too large, suggest user to run locally
## 01000 => BL pipeline failed on step 4/5 (main C prog && post-analysis)
## 10000 => GN pipeline failed on step 4/5 (main C prog && post-analysis)
EXIT_CODE=0

############################################################
### Script paths

# Top level dirs
CODE_BUNDLE_BL="./surface_residue"
CODE_BUNDLE_GN="./interior_residue"
PYTHONPATH="python"
PERLPATH="perl"

# Shared code used in both BL and GN pipelines
SIMPPRG=$CODE_BUNDLE_BL"/just_make_simplified_pdb_and_rm_HETATM.py"
RMCAPRG=$CODE_BUNDLE_BL"/just_make_simplified_pdb__CA_LINES_ONLY.py"
CALPHAPRG=$CODE_BUNDLE_BL"/calpha_modes.py"
TOPNMPRG=$CODE_BUNDLE_BL"/top_nm.pl"

# BL specific code
BOXSIZEPRG=$CODE_BUNDLE_BL"/get_box_size_to_det_num_MC_steps.py"
MAINCPRG=$CODE_BUNDLE_BL"/code_pack/bindingSites" 
MAPPRG=$CODE_BUNDLE_BL"/map_raw_sites__to__BL_sites.py"
GAPPRG=$CODE_BUNDLE_BL"/gap_theory_calculations.py"
HICPRG=$CODE_BUNDLE_BL"/identify_high_confidence_BL_sites.py"

# GN specific code
CMAPPRG=$CODE_BUNDLE_GN"/make_contact_map.py"
CMAPPRGGN=$CODE_BUNDLE_GN"/mod_cont_map_for_gn.py"
CORRMPRG=$CODE_BUNDLE_GN"/correlate_motions_ver2.py"
GNCOMMPRG=$CODE_BUNDLE_GN"/networkTools/gncommunities"
HIBTWNSSPRG=$CODE_BUNDLE_GN"/get_highest_btwnss_edges_frm_gn_output.py"
R2RPRG=$CODE_BUNDLE_GN"/residues__2__resids_types_chains_wo_VMD.py"


############################################################
### Output file paths

# Shared files
while (( "$#" )); do
key="$1"

case $key in 
	-h | --help)
	echo "To run STRESS, simple do:
		./stress -I INPUTFILE -O OUTPUTDIR -L LOGFILE [-s|--surface] [-i|--interior] [-h|--help]
		-I | --input:		followed by the inputfile (pdb file)
		-O | --output:		followed by the output directory
		-L | --logfile:		followed by the logfile
		-s | --surface:		run the surface-critical module(optional)
		-i | --interior:	run the interior-critical module(optonal)
		-h | --help:		for help(optional)
		To run the program, you have to at least choose one of surface/interior modules, or both."
	exit
	;;
	-I | --input)
	INPUT="$2"
	shift
	;;
	-O | --output)
	OUTPUTDIR="$2"
	shift
	;;
	-L | --logfile)
	LOGFILE=$2
	shift
	;;
	-s | --surface)
	BLNEEDED="true"
	;;
	-i | --interior)
	GNNEEDED="true"
	;;
	*)
	echo "Unknown option. Please use -h or --help to see help information."
	exit
	;;
esac

shift
done

{
	{ [ -n "$INPUT" ] && [ -n "$OUTPUTDIR" ] && [ -n "$LOGFILE" ] 
	} && { 
		[ -n "$BLNEEDED" ] || [ -n "$GNNEEDED" ] 
	}
} || { 
	echo "Missing argument. Please use -h or --help to see help information." 
	exit 
}

if [[ -z "`ls $INPUT`" ]]; then
	echo "Input file doesn't exit. Please try again."
	exit
fi


INPUTID=$(echo $INPUT | sed 's/.pdb[1]?//' | sed 's/.pdb//' | sed 's/.*\///')
SIMPPDB=$OUTPUTDIR$INPUTID"_SIMP.pdb"
CAONLYPDB=$OUTPUTDIR$INPUTID"_CA.pdb"
TOPFNM=$OUTPUTDIR$INPUTID".fnm_t10"
LOGTAG="" # Place holder for job IDs

# BL specific files
BLFILENAME=$INPUTID"_BL.dat" # Based on name convention in C program
BLFILE=$OUTPUTDIR$BLFILENAME
OUTFILE=$OUTPUTDIR$INPUTID".out"

# GN specific files
WGHCONTMAP=$OUTPUTDIR$INPUTID"_WGH_cont_map.txt"
COMMOUT=$OUTPUTDIR$INPUTID"_o_wght.out"
COMMLOG=$OUTPUTDIR$INPUTID"_o_wght.log"
OUTPUT=$OUTPUTDIR$INPUTID"_crit_mapping.txt"



############################################################
### Execute scripts

## Shared pipeline
echo "`date -u`: ${LOGTAG} Starting pipeline" >> $LOGFILE
echo "`date -u`: ${LOGTAG} INPUTID: ${INPUTID}, BL: ${BLNEEDED}, GN: ${GNNEEDED}" >> $LOGFILE

# step I: generate CA only PDB
echo "`date -u`: ${LOGTAG} Starting step I: generate CA only PDB" >> $LOGFILE
{
`$PYTHONPATH $SIMPPRG $INPUT $SIMPPDB` &&
$PYTHONPATH $RMCAPRG $SIMPPDB $CAONLYPDB
} || {
	echo "`date -u`: ${LOGTAG} PDB file checking failed!" >> $LOGFILE
	exit 2
}
echo "`date -u`: ${LOGTAG} Finished step I: generate CA only PDB" >> $LOGFILE

# step II: generate fnm_t10 file
echo "`date -u`: ${LOGTAG} Starting step II: generate fnm_t10 file" >> $LOGFILE
`$PYTHONPATH $CALPHAPRG $CAONLYPDB f | $PERLPATH $TOPNMPRG 10 > $TOPFNM`
FNMSIZE=$(grep "BEGIN" $TOPFNM | wc -l)
if [ $FNMSIZE -ne 10 ]; then
	echo "`date -u`: ${LOGTAG} Not enough FNM(10), protein too small!" >> $LOGFILE
	exit 1
fi
echo "`date -u`: ${LOGTAG} Finished step II: generate fnm_t10 file" >> $LOGFILE


## BL specific pipeline
[ "$BLNEEDED" == "true" ] &&
{
	echo "`date -u`: ${LOGTAG} Running BL module" >> $LOGFILE
	MCSTEP=$($PYTHONPATH $BOXSIZEPRG $SIMPPDB | cut -f2)
	SSTEP="$(((`wc -l $CAONLYPDB | awk '{print $1}'`) * 10))"

#	[ $SSTEP -gt 20000 ] && 
#	{
#		echo "`date -u`: ${LOGTAG} BL: PDB file too large!" >> $LOGFILE
#		let "EXIT_CODE += 4"
#	} ||
	{
		# For dev., we only use 50 simulation steps here
		# SSTEP=50
		echo "`date -u`: ${LOGTAG} BL: Finished step III: determine BL run parameters" >> $LOGFILE

		# Params should be passed in variables...but let's leave them like these now
		{
		`$MAINCPRG $CAONLYPDB 4 $INPUTID $SSTEP $MCSTEP $SIMPPDB $TOPFNM > $OUTFILE` &&
		mv $BLFILENAME $BLFILE
		$PYTHONPATH $MAPPRG $BLFILE $OUTFILE 0.00005 |
			$PYTHONPATH $GAPPRG $BLFILE 30 | 
			$PYTHONPATH $HICPRG $BLFILE 0.055 0.70 0.33 $INPUTID $OUTPUTDIR
		} && {
			echo "`date -u`: ${LOGTAG} BL: Finished steps IV and V: generate final output file" >> $LOGFILE
			echo "`date -u`: ${LOGTAG} BL module SUCCESS!" >> $LOGFILE
		} || {
			echo "`date -u`: ${LOGTAG} BL: Step 4/5 Failed!" >> $LOGFILE
			let "EXIT_CODE += 8"
		}
	}
} || echo "`date -u`: ${LOGTAG} Skipping BL module" >> $LOGFILE


## GN specific pipeline
[ "$GNNEEDED" == "true" ] &&
{
	echo "`date -u`: ${LOGTAG} Running GN module" >> $LOGFILE
	{
		`$PYTHONPATH $CMAPPRG $SIMPPDB heavy 4.5 1 $CAONLYPDB | 
		$PYTHONPATH $CMAPPRGGN |
	    $PYTHONPATH $CORRMPRG $CAONLYPDB $TOPFNM > $WGHCONTMAP` && 
		`$GNCOMMPRG $WGHCONTMAP $COMMOUT > $COMMLOG` &&
		$PYTHONPATH $HIBTWNSSPRG $COMMOUT |
		$PYTHONPATH $R2RPRG $OUTPUT $SIMPPDB
	} && {
		echo "`date -u`: ${LOGTAG} GN: Finished steps III and IV: generate final output file" >> $LOGFILE
		echo "`date -u`: ${LOGTAG} GN module SUCCESS!" >> $LOGFILE
	} || {
		echo "`date -u`: ${LOGTAG} GN: Step 3/4 Failed!" >> $LOGFILE
		let "EXIT_CODE += 16"
	}
} || echo "`date -u`: ${LOGTAG} Skipping GN module" >> $LOGFILE

echo "`date -u`: ${LOGTAG} Pipeline finished. Exiting with exit code: ${EXIT_CODE}" >> $LOGFILE
exit $EXIT_CODE
