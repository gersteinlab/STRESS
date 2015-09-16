#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_rng.h>

#include "pdbFile.h"
#include "surfaceProbe.h"
#include "network.h"
#include "transformations.h"
#include "modes.h"
#include "bindingSiteMeasures.h"


// sample command to run script:
//		./bindingSites ./4KGR_CA_ba1.pdb 4 out_ba_4KGR_ 200 40000 ./4KGR_SIMPLIFIED_BA1.pdb ./4KGR.fnm_t10 > out_ba_4KGR.txt



int main(int argc, char *argv[])
{
	pdbFile *macromolecule; //contains all information about the macromolecule in pdbInput
	pdbFile *macromolecule_w_heavy; // this line was not in orig
	FILE *pdbInput; //pdb file (input)
	FILE *pdbInput_heavy; //this line was not in orig
	box *pbc; //Contains all the information about the PBC
	ligand *lig; //Contains all the information about the ligand
	int Ntrials, NstepsPerTrial = 0; //Specifications of MC simulations
	// FILE *output1, *output2, *output3; //output files
	FILE *output2;
	Graph *network; //network
	int *residues, residueCount;
	int i, j;
	bindingSite *allSites, *mergedSites;
	int numMergedSites;
	modePtr normalModes;
	FILE *fnm_t10;
	char outputFileName[1000];
	
	/*
	// below *was* not in orig
	gsl_rng* rng;
		//rng = gsl_rng_alloc(gsl_rng_ran3);
	    //init random number generator
    long seed = time(NULL);
    //gsl_rng* rng = gsl_rng_alloc(gsl_rng_ran3);
    printf("Initial random number seed = %lu\n", seed);
    gsl_rng_set(rng, seed);
    */
	
	//Checking the arguments.
	if ((argc < 4) || (argc > 8)) // was:  if ((argc < 6) || (argc > 8))
	{
		printf("Error in number of arguments. Number of arguments should be between 3-5. It is %d. Usage:\n", (argc-1));
		printf("./bindingSites <pdbFile> <numberAtomsLigand> <outputFile> [<numberOfTrials> [<numberOfStepsPerTrial>]]\n");
		printf("where:\n");
		printf("   pdbFile refers to the input file -- the pdb file must be the same as that used to genereate the anm output\n");
		printf("   numberAtomsLigand refers to the number of atoms in the ligand\n");
		printf("   outputPrefix refers to the prefix for all output files\n");
		printf("   numberOfTrials refers to number of Monte Carlo surface probe simulations. Default value = 1000\n");
		printf("   numberOfStepsPerTrial refers to number of steps per surface probe simulation. Default value = 10000\n");
		printf("   pdb_w_heavy_atoms refers to the input file WITH HEAVY ATOMS\n");  // this line was not in orig
		printf("   fnm_t10 file refers to file with the MMTK output file \n");
		exit(1);
	}
	
	//Creating memory to save and set the macromolecule coordinates
	macromolecule = (pdbFile *) malloc(sizeof(pdbFile));
	if (macromolecule == NULL)
	{
		printf("ERROR: No memory available to store macromolecule information\n");
		exit(1);
	}

	// below was not in orig
	macromolecule_w_heavy = (pdbFile *) malloc(sizeof(pdbFile));
	if (macromolecule_w_heavy == NULL)
	{
		printf("ERROR: No memory available to store macromolecule_w_heavy information\n");
		exit(1);
	}
	
	if ((pdbInput = fopen(argv[1], "r")) == NULL)
	{
		printf("pdb file %s does not exist or open\n", argv[1]);
		exit(1);
	}
	
	// below was not in orig
	if ((pdbInput_heavy = fopen(argv[6], "r")) == NULL)
	{
		printf("pdb (heavy) file %s does not exist or open\n", argv[6]);
		exit(1);
	}
	
	macromolecule->nAtoms = findNumberOfAtomsPdbFile(pdbInput);
	printf("Number of atoms in file is %d\n", macromolecule->nAtoms);
	rewind(pdbInput);
	
	if ((macromolecule->atom = (atomStr *) calloc(macromolecule->nAtoms, sizeof(atomStr))) == NULL)
	{
		printf("No space to save the coordinates of all atoms\n");
		exit(1);
	}
	readPdbFile(macromolecule, pdbInput);
	fclose(pdbInput);	


	// below was not in orig
	macromolecule_w_heavy->nAtoms = findNumberOfAtomsPdbFile(pdbInput_heavy);
	printf("Number of atoms in (heavy) file is %d\n", macromolecule_w_heavy->nAtoms);
	rewind(pdbInput_heavy);	
	if ((macromolecule_w_heavy->atom = (atomStr *) calloc(macromolecule_w_heavy->nAtoms, sizeof(atomStr))) == NULL)
	{
		printf("No space to save the coordinates of all atoms in heavy file\n");
		exit(1);
	}
	readPdbFile(macromolecule_w_heavy, pdbInput_heavy);
	fclose(pdbInput_heavy);	

	//Setting network
	if ((network = (Graph *) malloc(sizeof(Graph))) == NULL)
	{
		printf("ERROR: No space to make network\n");
		exit(1);
	}
	
	constructNetwork(network, macromolecule);

	if ((residues = (int *)malloc(network->nres*sizeof(int))) == NULL)
	{
		printf("ERROR: no space for saving network information\n");
		exit(1);
	}
	residueCount = macromolecule->nAtoms;
	for (i=0; i<residueCount; i++)
	{
		residues[i] = i;
	}
	
	floydWarshall(network, residues, residueCount);
	calculateLC(network, residues, residueCount);
	//printf("Got back\n");
	
	//Setting up ligand
	if ((lig = (ligand *) malloc(sizeof(ligand))) == NULL)
	{
		printf("ERROR:No space to save the coordinates of ligand\n");
		exit(1);
	}
	
	lig->nAtoms = atoi(argv[2]);
	printf("Number of atoms in ligand is %d\n", lig->nAtoms);
	if ((lig->nAtoms < 2) || (lig->nAtoms > 8)) 
	{
		printf("ERROR: ligand has to be made up of 2-8 atoms\n");
		exit(1);
	}
	if ((lig->atom = calloc(lig->nAtoms, sizeof(coord))) == NULL)
	{
		printf("ERROR: No space to save the coordinates of ligand\n");
		exit(1);
	}
	
	//Setting PBC
	if ((pbc = (box *) malloc(sizeof(box))) == NULL) 
	{
		printf("ERROR: No memory available to store box information\n");
		exit(1);
	}
	
	calculatePBCbox(macromolecule, pbc);
	printf("The center of the box is located at %f %f %f\n", pbc->center[0], pbc->center[1], pbc->center[2]);
	printf("The size of the box is %f %f %f\n", pbc->size[0], pbc->size[1], pbc->size[2]);

	if ((fnm_t10 = fopen(argv[7], "r")) == NULL)
	{
		printf("ERROR: Cannot open fnm_t10 %s\n", argv[7]);
		exit(1);
	}

/*
	strcpy(outputFileName, argv[3]);
    strcat(outputFileName, "_bindingSites.dat");
	if ((output1 = fopen(outputFileName, "w")) == NULL)
	{
		printf("ERROR:Cannot open output file %s\n", outputFileName);
		exit(1);
	}
*/	
	strcpy(outputFileName, argv[3]);
    strcat(outputFileName, "_BL.dat");
	if ((output2 = fopen(outputFileName, "w")) == NULL)
	{
		printf("ERROR:Cannot open output file %s\n", outputFileName);
		exit(1);
	}
/*
	strcpy(outputFileName, argv[3]);
    strcat(outputFileName, "_LC.dat");
	if ((output3 = fopen(outputFileName, "w")) == NULL)
	{
		printf("ERROR:Cannot open output file %s\n", outputFileName);
		exit(1);
	}
*/	
	if (argc > 5)
	{
		Ntrials = atoi(argv[4]);
		if (Ntrials < 1) {
			printf("ERROR: number of trials = %d and no trials possible\n", Ntrials);
			exit(1);
		}
	} else 
	{
		Ntrials = 1000;
	}
	printf("Number of trials are %d\n", Ntrials);
	
	if (argc >= 6) // *was*  if (argc == 8)
	{
		printf("argc >= 8, and Num trials is: %d    NstepsPerTrial is: %d \n", Ntrials, NstepsPerTrial);
		NstepsPerTrial = atoi(argv[5]);
		if (NstepsPerTrial < 1) {
			printf("ERROR: number of steps per trial = %d and no trials possible\n", NstepsPerTrial);
			exit(1);
		}
	} else 
	{
		NstepsPerTrial = 100000;
	}
	printf("Number of steps per trial are %d\n", NstepsPerTrial);

	if ((allSites = (bindingSite *) malloc(Ntrials*sizeof(bindingSite))) == NULL)
	{
		printf("Error: No space for all Sites\n");
		exit(1);
	}
	//Doing the surface probe simulation
	runSurfaceProbeSimulation(macromolecule, lig, pbc, Ntrials, NstepsPerTrial, network->LC, allSites, macromolecule_w_heavy);  //   was: runSurfaceProbeSimulation(macromolecule, lig, pbc, Ntrials, NstepsPerTrial, network->LC, allSites);
	//runSurfaceProbeSimulation(macromolecule, lig, pbc, Ntrials, NstepsPerTrial, network->LC, allSites);   // NOTE: this was the orig
	
	//printf("I got back here %d\n", numMergedSites);
	//for (i=0; i<numMergedSites; i++)
	//{
		//printf("nRes %d\n", mergedSites[i].nRes);
	//}
	numMergedSites = getNumberFinalSites(allSites, Ntrials);
	if ((mergedSites = (bindingSite*)malloc(numMergedSites*sizeof(bindingSite))) == NULL)
	{
		printf("Error: No space to save merged sites\n");
		exit(1);
	}
	copyFinalSites(allSites, mergedSites, Ntrials);
	
	// outputFinalSites(mergedSites, numMergedSites, output1);

	/*
	for (i=0; i<numMergedSites; i++)
	{
		printf("Number of residues in site %d are %d\n", i, mergedSites[i].nRes);
		printf("They are:\n");
		for (j=0; j<mergedSites[i].nRes; j++)
		{
			printf("%d ", mergedSites[i].resid[j]);  // was:  printf("%d ", mergedSites[i].residue[j]);  
		}
		printf("\n");
	}
	*/

	// below for block was not in orig
	for (i=0; i<numMergedSites; i++)
	{
		printf("The resid's are:\n");
		for (j=0; j<mergedSites[i].nRes; j++)
		{
			//printf("%d  ", mergedSites[i].resid[j]);
			printf("%d_%s  ", mergedSites[i].resid[j], macromolecule->atom[mergedSites[i].residue[j]].chain);
		}
		printf("\n");
	}
	
	if ((normalModes = (modePtr)malloc(sizeof(modeStr))) == NULL)
	{
		printf("ERROR: No space to save normal modes\n");
		exit(1);
	}




	////////////    Allocate space for modes 1 - 10  ////////////
	// Mode 1:
	if ((normalModes->fnm_t10_data_mode1 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode1\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode1[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode1\n");
			exit(1);
		}
	}

	// Mode 2:
	if ((normalModes->fnm_t10_data_mode2 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode2\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode2[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode2\n");
			exit(1);
		}
	}

	// Mode 3:
	if ((normalModes->fnm_t10_data_mode3 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode3\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode3[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode3\n");
			exit(1);
		}
	}

	// Mode 4:
	if ((normalModes->fnm_t10_data_mode4 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode4\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode4[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode4\n");
			exit(1);
		}
	}

	// Mode 5:
	if ((normalModes->fnm_t10_data_mode5 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode5\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode5[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode5\n");
			exit(1);
		}
	}

	// Mode 6:
	if ((normalModes->fnm_t10_data_mode6 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode6\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode6[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode6\n");
			exit(1);
		}
	}

	// Mode 7:
	if ((normalModes->fnm_t10_data_mode7 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode7\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode7[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode7\n");
			exit(1);
		}
	}

	// Mode 8:
	if ((normalModes->fnm_t10_data_mode8 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode8\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode8[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode8\n");
			exit(1);
		}
	}

	// Mode 9:
	if ((normalModes->fnm_t10_data_mode9 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode9\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode9[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode9\n");
			exit(1);
		}
	}

	// Mode 10:
	if ((normalModes->fnm_t10_data_mode10 = (float **)malloc(3*macromolecule->nAtoms*sizeof(float *))) == NULL)
	{
		printf("ERROR: No space to save fnm_t10_data_mode10\n");
		exit(1);
	}
	for (i=0; i<3*macromolecule->nAtoms; i++)
	{
		if ((normalModes->fnm_t10_data_mode10[i] = (float *)malloc(3*macromolecule->nAtoms*sizeof(float))) == NULL)
		{
			printf("ERROR: No space to save fnm_t10_data_mode10\n");
			exit(1);
		}
	}

	printf("macromolecule->nAtoms  %d\n", macromolecule->nAtoms);
	read_fnm_t10(normalModes->fnm_t10_data_mode1, normalModes->fnm_t10_data_mode2, normalModes->fnm_t10_data_mode3, normalModes->fnm_t10_data_mode4, normalModes->fnm_t10_data_mode5, normalModes->fnm_t10_data_mode6, normalModes->fnm_t10_data_mode7, normalModes->fnm_t10_data_mode8, normalModes->fnm_t10_data_mode9, normalModes->fnm_t10_data_mode10, fnm_t10, macromolecule->nAtoms);

	/*
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode1 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode1[i][0], normalModes->fnm_t10_data_mode1[i][1], normalModes->fnm_t10_data_mode1[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode2 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode2[i][0], normalModes->fnm_t10_data_mode2[i][1], normalModes->fnm_t10_data_mode2[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode3 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode3[i][0], normalModes->fnm_t10_data_mode3[i][1], normalModes->fnm_t10_data_mode3[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode4 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode4[i][0], normalModes->fnm_t10_data_mode4[i][1], normalModes->fnm_t10_data_mode4[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode5 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode5[i][0], normalModes->fnm_t10_data_mode5[i][1], normalModes->fnm_t10_data_mode5[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode6 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode6[i][0], normalModes->fnm_t10_data_mode6[i][1], normalModes->fnm_t10_data_mode6[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode7 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode7[i][0], normalModes->fnm_t10_data_mode7[i][1], normalModes->fnm_t10_data_mode7[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode8 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode8[i][0], normalModes->fnm_t10_data_mode8[i][1], normalModes->fnm_t10_data_mode8[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode9 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode9[i][0], normalModes->fnm_t10_data_mode9[i][1], normalModes->fnm_t10_data_mode9[i][2]);
	}
	printf("\n\n");
	for (i = 0; i<(macromolecule->nAtoms); i++)
	{
		printf("reading data from fnm_t10_data_mode10 as:    %f   %f   %f \n", normalModes->fnm_t10_data_mode10[i][0], normalModes->fnm_t10_data_mode10[i][1], normalModes->fnm_t10_data_mode10[i][2]);
	}
	*/


	calculateBindingLeverage(macromolecule, normalModes, normalModes->fnm_t10_data_mode1, normalModes->fnm_t10_data_mode2, normalModes->fnm_t10_data_mode3, normalModes->fnm_t10_data_mode4, normalModes->fnm_t10_data_mode5, normalModes->fnm_t10_data_mode6, normalModes->fnm_t10_data_mode7, normalModes->fnm_t10_data_mode8, normalModes->fnm_t10_data_mode9, normalModes->fnm_t10_data_mode10, 10, mergedSites, numMergedSites, output2, macromolecule_w_heavy);
	//calculateBindingLeverage(macromolecule, normalModes, 10, mergedSites, numMergedSites, output2, macromolecule_w_heavy);
	//calculateBindingLeverage(macromolecule, normalModes, 10, mergedSites, numMergedSites, output2, macromolecule_w_heavy);


	// calculateLCmergedSite(macromolecule, mergedSites, numMergedSites, network, output3);
	
	//Freeing memory
	// D14 Ad
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode1[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode2[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode3[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode4[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode5[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode6[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode7[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode8[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode9[i]);
	for (i=0; i<3*macromolecule->nAtoms; i++)
		free(normalModes->fnm_t10_data_mode10[i]);

	free(normalModes->fnm_t10_data_mode1);
	free(normalModes->fnm_t10_data_mode2);
	free(normalModes->fnm_t10_data_mode3);
	free(normalModes->fnm_t10_data_mode4);
	free(normalModes->fnm_t10_data_mode5);
	free(normalModes->fnm_t10_data_mode6);
	free(normalModes->fnm_t10_data_mode7);
	free(normalModes->fnm_t10_data_mode8);
	free(normalModes->fnm_t10_data_mode9);
	free(normalModes->fnm_t10_data_mode10);

	free(normalModes);
	free(macromolecule->atom);
	free(macromolecule);
	free(macromolecule_w_heavy->atom);	// was not in orig
	free(macromolecule_w_heavy);		//  was not in orig
	free(pbc);
	free(lig->atom);
	free(lig);


	for (i=0; i<Ntrials; i++)
	{
		if (allSites[i].nRes > 0)
		{
			//printf("Got here i %d\n", i);
			free(allSites[i].resid);
			free(allSites[i].residue);
			//printf("%d\n", allSites[i].numSprings);
			for (j=0; j<allSites[i].numSprings; j++)
			{
				//printf("j %d\n", j);
				free(allSites[i].spring[j]);
			}
			if (allSites[i].numSprings > 0) free(allSites[i].spring);
		}
	}
	free(allSites);


	for (i=0; i<numMergedSites; i++)
	{
		//printf("i %d\n", i);
		free(mergedSites[i].residue);
		free(mergedSites[i].resid);
		for (j=0; j<mergedSites[i].numSprings; j++)
		{
			//printf("j %d\n", j);
			free(mergedSites[i].spring[j]);
		}
		if (mergedSites[i].numSprings > 0) free(mergedSites[i].spring);
	}
	free(mergedSites);

	for (i=0; i<network->nedges; i++)
	{
		free(network->edge[i]);
	}
	for (i=0; i<network->nres; i++)
	{
		free(network->dis[i]);
		free(network->shortDis[i]);
	}
	free(network->dis);
	free(network->shortDis);
	free(network->edge);
	free(network->LC);
	free(network);
	free(residues);

	//fclose(output1);
	fclose(output2);
	//fclose(output3);
	fclose(fnm_t10);

	printf("\n\nfin\n\n\n");

	return 0;
}
