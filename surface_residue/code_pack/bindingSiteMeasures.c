#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "modes.h"
#include "pdbFile.h"
#include "network.h"
#include "surfaceProbe.h"
#include "bindingSiteMeasures.h"



void calculateBindingLeverage(pdbFile *macromolecule, modePtr normalModes, float **fnm_t10_data_mode1, float **fnm_t10_data_mode2, float **fnm_t10_data_mode3, float **fnm_t10_data_mode4, float **fnm_t10_data_mode5, float **fnm_t10_data_mode6, float **fnm_t10_data_mode7, float **fnm_t10_data_mode8, float **fnm_t10_data_mode9, float **fnm_t10_data_mode10, int numModes, bindingSite *mergedSites, int numMergedSites, FILE *output, pdbFile *macromolecule_w_heavy)
{
	int i, j;
	float BL[numMergedSites];
	float bl = 0.;
	int order[numMergedSites];
	int nRes = 0;
	
	for (i=0; i<numMergedSites; i++)
	{
		order[i] = i;
	}
	//animateNormalModes(macromolecule, R, normalModes, numModes);

	for (i = 0; i<numMergedSites; i++)
	{
		BL[i] = 0.0;
	}

	for (i = 0; i<numMergedSites; i++)
	{
		bl = 0.;
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode1);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode2);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode3);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode4);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode5);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode6);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode7);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode8);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode9);
		bl = bl + calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, fnm_t10_data_mode10);

		//printf("\n\nbl val is: %f  \n", bl);
		mergedSites[i].BL = bl;
		BL[i] = mergedSites[i].BL;
	}


	/*
	for (i = 0; i<numMergedSites; i++)
	{
		   //mergedSites[i].BL = calculateBindingLeverageSite(mergedSites, i, macromolecule, R, macromolecule_w_heavy);     // was mergedSites[i].BL = calculateBindingLeverageSite(mergedSites, i, macromolecule, R);
		//mergedSites[i].BL = calculateBindingLeverageSite(mergedSites, i, macromolecule, R);
		mergedSites[i].BL = calculateBindingLeverageSite__AddingModes(mergedSites, i, macromolecule, R, normalModes->fnm_t10_data_mode1);

		BL[i] = mergedSites[i].BL;

		//printf("BL is: %f \n", BL[i]);
	}
	*/





	sortArray(BL, numMergedSites, order);
	
	printf("->printing Binding Leverage to output file\n");
	
	for (i=0; i<numMergedSites; i++)
	{
		fprintf(output, "%6d  %4f  ", order[i], mergedSites[order[i]].BL);
		if (mergedSites[order[i]].nRes > 10)
		{
			nRes = 10;
		}
		else
		{
			nRes = mergedSites[order[i]].nRes;
		}
		for (j=0; j<nRes; j++)
		{
			fprintf(output, "%4d_%s ", mergedSites[order[i]].resid[j], macromolecule->atom[mergedSites[order[i]].residue[j]].chain);  // was:  fprintf(output, "%4d ", mergedSites[order[i]].residue[j]);   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			
			//fprintf(output, "%4d ", mergedSites[order[i]].resid[j]);  // was:  fprintf(output, "%4d ", mergedSites[order[i]].residue[j]);   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
		fprintf(output, "\n");
	}
	
	printf("<-printing Binding Leverage to output file\n");
	
	return;
}




//// ***** Below is an extremely inefficient implementation for calculating binding leverage:
float calculateBindingLeverageSite__AddingModes(bindingSite *site, int siteNum, pdbFile *macromolecule, float **fnm_t10_data_modeX)
{
	float R[macromolecule->nAtoms*3]; // NOTE (DC): R[] will represent the new atom locations (after animating w/normal modes)
	int i, j, res1, res2;
	float r0 = 0., r = 0., Dr = 0.;
	float BL = 0.;
	//printf("Got here %d\n", site[siteNum].numSprings);

	for (i=0; i<macromolecule->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			R[i*3 + j] = 0;
		}
	}

	// printf("\n\n unperturbed based on fnm_t10_data_modeX: \n");
	for (i=0; i<macromolecule->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			//R[i*3 + j] += normalModes->fnm_t10_data_modeX[i][j];
			//R[i*3 + j] = &(fnm_t10_data_modeX[i][j]);
			R[i*3 + j] = fnm_t10_data_modeX[i][j];
		}
		R[i*3] += macromolecule->atom[i].x;
		R[i*3 + 1] += macromolecule->atom[i].y;
		R[i*3 + 2] += macromolecule->atom[i].z;
		//printf("i %d  fnm_x   %f    fnm_y   %f    fnm_z   %f  \n", i, fnm_t10_data_modeX[i][0], fnm_t10_data_modeX[i][1], fnm_t10_data_modeX[i][2]);
		//printf("i %d R[i*3]   %f  R[i*3+1]  %f  R[i*3+2]  %f  \n", i, R[i*3], R[i*3+1], R[i*3+2]);
		//printf("i %d unp x    %f    unp y   %f    unp z   %f  \n", i, macromolecule->atom[i].x, macromolecule->atom[i].y, macromolecule->atom[i].z);
		//printf("\n");
	}
	
	for (i=0; i<site[siteNum].numSprings; i++)
	{
		//printf("i %d numSprings %d\n", i, site[siteNum].numSprings);
		res1 = site[siteNum].spring[i]->res1;
		res2 = site[siteNum].spring[i]->res2;

		//printf("macromolecule->atom[res1].x: %f", macromolecule->atom[res1].x);
		r0 = sqrt((macromolecule->atom[res1].x - macromolecule->atom[res2].x)*(macromolecule->atom[res1].x - macromolecule->atom[res2].x) + (macromolecule->atom[res1].y - macromolecule->atom[res2].y)*(macromolecule->atom[res1].y - macromolecule->atom[res2].y) + (macromolecule->atom[res1].z - macromolecule->atom[res2].z)*(macromolecule->atom[res1].z - macromolecule->atom[res2].z));

		//printf("R[3*res1], R[3*res1 + 1], R[3*res1 + 2]:  %f  %f  %f  \n", R[3*res1], R[3*res1 + 1], R[3*res1 + 2]);
		r = sqrt((R[3*res1] - R[3*res2])*(R[3*res1] - R[3*res2]) + (R[3*res1 + 1] - R[3*res2 + 1])*(R[3*res1 + 1] - R[3*res2 + 1]) + (R[3*res1 + 2] - R[3*res2 + 2])*(R[3*res1 + 2] - R[3*res2 + 2]));
		//printf("r0 %f r %f\n", r0, r);
		Dr = r - r0;
		// take absolute value of Dr (the function "abs" was outputting all Dr vals as 0.000000 for some reason)
		if (Dr < 0.0)
		{
			Dr = -1.0 * Dr;
		}
		//printf("Dr is:  %f   r0 is: %f   r is: %f   \n", Dr, r0, r);
		BL += Dr*Dr;
	}
	//printf("siteNum = %d BL = %f\n", siteNum, BL);
	// printf("bl val here is:  %f   \n", BL);	
	return BL;
}





float calculateBindingLeverageSite(bindingSite *site, int siteNum, pdbFile *macromolecule, float *R)
{
	int i, res1, res2;
	float r0, r, Dr;
	float BL = 0.;
	//printf("Got here %d\n", site[siteNum].numSprings);
	
	for (i=0; i<site[siteNum].numSprings; i++)
	{
		//printf("i %d numSprings %d\n", i, site[siteNum].numSprings);
		res1 = site[siteNum].spring[i]->res1;
		res2 = site[siteNum].spring[i]->res2;

		//printf("macromolecule->atom[res1].x: %f", macromolecule->atom[res1].x);
		r0 = sqrt((macromolecule->atom[res1].x - macromolecule->atom[res2].x)*(macromolecule->atom[res1].x - macromolecule->atom[res2].x) + (macromolecule->atom[res1].y - macromolecule->atom[res2].y)*(macromolecule->atom[res1].y - macromolecule->atom[res2].y) + (macromolecule->atom[res1].z - macromolecule->atom[res2].z)*(macromolecule->atom[res1].z - macromolecule->atom[res2].z));

		r = sqrt((R[3*res1] - R[3*res2])*(R[3*res1] - R[3*res2]) + (R[3*res1 + 1] - R[3*res2 + 1])*(R[3*res1 + 1] - R[3*res2 + 1]) + (R[3*res1 + 2] - R[3*res2 + 2])*(R[3*res1 + 2] - R[3*res2 + 2]));
		//printf("r0 %f r %f\n", r0, r);
		Dr = fabsf(r - r0);
		BL += Dr*Dr;
	}
	//printf("siteNum = %d BL = %f\n", siteNum, BL);
	return BL;
}




void calculateLCmergedSite(pdbFile *macromolecule, bindingSite *mergedSites, int numMergedSites, Graph *network, FILE *output)
{
	int i, j;
	int order[numMergedSites];
	float LC[numMergedSites];
	int nRes;

	for (i=0; i<numMergedSites; i++)
        {
                order[i] = i;
        }

	
	for (i=0; i<numMergedSites; i++)
	{
		mergedSites[i].LC = calculateLCsite(mergedSites, i, network);
		LC[i] = mergedSites[i].LC;
	}
	sortArray(LC, numMergedSites, order);
	
	for (i=0; i<numMergedSites; i++)
	{
		fprintf(output, "%6d  %4f  ", order[i], mergedSites[order[i]].LC);
		if (mergedSites[order[i]].nRes > 10)
		{
			nRes = 10;
		}
		else
		{
			nRes = mergedSites[order[i]].nRes;
		}
		for (j=0; j<nRes; j++)
		{
			fprintf(output, "%4d ", mergedSites[order[i]].resid[j]);  // was:  fprintf(output, "%4d ", mergedSites[order[i]].residue[j]);
		}
		fprintf(output, "\n");
	}
	
	return;
}

float calculateLCsite(bindingSite *mergedSites, int siteNum, Graph *network)
{
	int i, nRes;
	float LC = 0.;
	
	if (mergedSites[siteNum].nRes >= 10)
	{
		nRes = 10;
	} else
	{
		nRes = mergedSites[siteNum].nRes;
	}
	
	for (i=0; i<nRes; i++)
	{
		LC += network->LC[mergedSites[siteNum].residue[i]];
	}
	LC /= nRes;
	
	//printf("siteNum = %d LC = %f\n", siteNum, LC);
	return LC;
}

// TODO: BAD sorting 
// valgrind complains about uninitialised A[order[j]]
void sortArray(float *A, int numMergedSites, int *order)
{
	int i, j;
	int swap;
	
	for (i=0; i<numMergedSites; i++)
	{
		for (j=(i+1); j<numMergedSites; j++)
		{
			if (A[order[j]] > A[order[i]])
			{
				
				swap = order[i];
				order[i] = order[j];
				order[j] = swap;
			}
		}
		//printf("%d %d\n", i, order[i]);
	}
	
	return;
}
