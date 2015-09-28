/*
 * University of Illinois Open Source License
 * Copyright 2009, 2011, 2012 Luthey-Schulten Group, 
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *                           University of Illinois at Urbana-Champaign
 *                           http://www.scs.uiuc.edu/~schulten
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimers in the documentation 
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Anurag Sethi, John Eargle
 * Modifications by: Shantao Li, Declan Clarke (lab of Mark Gerstein, 
 * Yale University, 2015)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "pdbFile.h"
#include "surfaceProbe.h"
#include "transformations.h"

const float change_in_bondangle_range = 3.14159/3.;
const float update_pivot_range = 3.14159/6.;
const float square_thresh_for_dist_to_springs = 9.0;	 // ****NOTE 2: In paper, this dist is reported as 3.5, but in Berezovsky's code, a value of 3.0 is used  --  see Berezovsky's analysis.c

// STL: some globals
inode **Hash;
int TraverseOrder[3] = { 0, -1, 1 };
int GridSize[3];
float minR[3];

void runSurfaceProbeSimulation(pdbFile *macromolecule, ligand *lig, box *pbc, int Ntrials, int NstepsPerTrial, float *LC, bindingSite *allSites, pdbFile *macromolecule_w_heavy)  //  was:  void runSurfaceProbeSimulation(pdbFile *macromolecule, ligand *lig, box *pbc, int Ntrials, int NstepsPerTrial, float *LC, bindingSite *allSites)
{
	int i, j, lig_atom, lig_res, lig_atom_to_print, tem;
	inode *temptr = NULL;
	float v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v1, v2, v1_v2_dot_prod, bond_angle; // these are for bond angle calculations
	float maxR[3];
	float U[lig->nAtoms + 1];	// Energy array; the last slot stores the sum
	// bool V[lig->nAtoms*2 - 1];	// Validation array;the last slot stores the bond angle resutls (y/n)
	bool V[lig->nAtoms];
	// float U[NstepsPerTrial];
	ligand *ligNew, *ligTemp, *ligTempOld;
	int stableSites = 0;

	// srand(12345);
	srand(time(NULL));
	printf("->runSurfaceProbeSimulation\n");
	if ((ligNew = (ligand *) malloc(sizeof(ligand))) == NULL)
	{
		printf("Not enough space for performing surface probe simulation\n");
		exit(1);
	}
	
	ligNew->nAtoms = lig->nAtoms;
	if ((ligNew->atom = (coord *) calloc(ligNew->nAtoms, sizeof(coord))) == NULL)
	{
		printf("Not enough space for performing surface probe simulation\n");
		exit(1);
	}
	
	if ((ligTemp = (ligand *) malloc(sizeof(ligand))) == NULL || (ligTempOld = (ligand *) malloc(sizeof(ligand))) == NULL)
        {
                printf("Not enough space for performing surface probe simulation\n");
                exit(1);
        }

        ligTemp->nAtoms = lig->nAtoms;
	ligTempOld->nAtoms = lig->nAtoms;
        if ((ligTemp->atom = (coord *) calloc(ligTemp->nAtoms, sizeof(coord))) == NULL \
		|| (ligTempOld->atom = (coord *) calloc(ligTempOld->nAtoms, sizeof(coord))) == NULL)
        {
                printf("Not enough space for performing surface probe simulation\n");
                exit(1);
        }


	for (i=0; i<3; i++)
	{
		minR[i] = pbc->center[i] - pbc->size[i]/2.;
		maxR[i] = pbc->center[i] + pbc->size[i]/2.;
		GridSize[i] = (int)(pbc->size[i]/5 + 0.5);
		//printf("%d\n", GridSize[i]);
		// Check PBC in energy calculation
		if (pbc->size[i]/2 < 10.1) {
			printf("Space is too small, you have to consider PBC in energy calculation!\n");
			exit(1);
		}
	}
	
	// GridSize = ceil(pbc->size[0]/5); // STL: Just want to be safe
	Hash = (inode**)malloc(sizeof(inode*)*GridSize[0]*GridSize[1]*GridSize[2]);
	for (i=0; i<GridSize[0]*GridSize[1]*GridSize[2]; i++)
		Hash[i] = NULL;

	// STL: Hash the macromolecule atoms
	for (int k=0; k<macromolecule_w_heavy->nAtoms; k++) {
		tem = (int)((macromolecule_w_heavy->atom[k].x - minR[0])/5) + GridSize[0] * ((int)((macromolecule_w_heavy->atom[k].y - minR[1])/5) \
				+ (int)(((macromolecule_w_heavy->atom[k].z - minR[2])/5))*GridSize[1]);
		//printf("%d, %d, %d, %d, %d\n", ((int)((macromolecule_w_heavy->atom[k].x - minR[0])/5)), (int)((macromolecule_w_heavy->atom[k].y - minR[1])/5), (int)((macromolecule_w_heavy->atom[k].z - minR[2])/5), tem, k);
		temptr = Hash[tem];
		if (!temptr) {
			temptr = malloc(sizeof(inode));
			Hash[tem] = temptr;
			temptr->poz = k;
			temptr->next = NULL;
		}
		else {	
			while (temptr->next) 
				temptr = temptr->next; 
			temptr->next = malloc(sizeof(inode));
			temptr->next->poz = k;
			temptr->next->next = NULL;
		}
	}
	
	for (i=0; i<Ntrials; i++)
	{
		printf("Trial run %d\n", i);
		//Resetting random number generator
		// srand(time(NULL));
		// Initialize U and V
		U[lig->nAtoms] = 10000;

		while (U[lig->nAtoms] > 1000) {
			//Guessing probe initial coordinates (atom by atom)
			// printf("Guessing initial coordinates for ligand\n");
			guessFirstAtomCoordinates(lig, minR, maxR);
			guessSecondAtomCoordinates(lig);
			for (j=2; j<lig->nAtoms; j++) guessAnotherAtomCoordinates(lig, j);
			//Making sure ligand is in box
			transformCoordinatesToSameBox(lig, minR, maxR, pbc);
			
			//Initial potential energy
			for (j = 0; j<lig->nAtoms; j++) V[j] = false;
			calculateEnergyCAmodel(macromolecule_w_heavy, lig, U, V, pbc, ligTempOld);
		}
		
		/*while (U[lig->nAtoms] > 1000) {
			performOneMonteCarloStep(macromolecule_w_heavy, lig, ligNew, pbc, minR, maxR, U, V);
		}*/
		//U[0] = calculateEnergyCAmodel(macromolecule, lig);   <---- this used for CA model
		// printf(" U is %f\n", U[0]);
		
		//Perform Simulations
		for (j=1; j<NstepsPerTrial; j++)
		{
			performOneMonteCarloStep(macromolecule_w_heavy, lig, ligNew, pbc, minR, maxR, U, V, &ligTemp, &ligTempOld);
			//performOneMonteCarloStep(macromolecule, lig, ligNew, pbc, minR, maxR, U, j);
			//printf("step %d U %f\n", j, U[j]);
		}

		printf("Ufinal is %f\n", U[lig->nAtoms]);

		float del_x_sqrd;
		float del_y_sqrd;
		float del_z_sqrd;
		float CA_CA_bond_dist_in_lig;
		if (U[lig->nAtoms] < 0.)  // NOTE (DC) : ie, if the last U val from the NstepsPerTrial MC steps is favorable (ie, negative) -- the ligand's position and orientation at that particular site are good (binding site is good candidate)
		{
			stableSites++;
			// Since this is a stable ligand site, let's print it's coordinates and view the structure in Pymol (does it hit the protein interior?)

			makeMoleculeWhole(lig, ligNew, pbc);  // *nec* to keep this??

			lig_res = 330;
			lig_atom_to_print = 0;
			for (lig_atom=0; lig_atom < ligNew->nAtoms; lig_atom++)
			{
				lig_atom_to_print = lig_atom_to_print + 1;
				printf("ATOM    %d  CA  ALA B   %d      %.3f  %.3f  %.3f                       C\n", lig_res, lig_atom_to_print, ligNew->atom[lig_atom].r[0], ligNew->atom[lig_atom].r[1], ligNew->atom[lig_atom].r[2]);  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Print a record of all CA-CA bond distances in the final ligand configuration:
				if (lig_atom>0)
				{
					del_x_sqrd = (ligNew->atom[lig_atom].r[0]-ligNew->atom[lig_atom-1].r[0]) * (ligNew->atom[lig_atom].r[0]-ligNew->atom[lig_atom-1].r[0]);
					del_y_sqrd = (ligNew->atom[lig_atom].r[1]-ligNew->atom[lig_atom-1].r[1]) * (ligNew->atom[lig_atom].r[1]-ligNew->atom[lig_atom-1].r[1]);
					del_z_sqrd = (ligNew->atom[lig_atom].r[2]-ligNew->atom[lig_atom-1].r[2]) * (ligNew->atom[lig_atom].r[2]-ligNew->atom[lig_atom-1].r[2]);
					CA_CA_bond_dist_in_lig = sqrt(del_x_sqrd + del_y_sqrd + del_z_sqrd);
					printf("CA_CA_bond_dist_in_lig: %f \n", CA_CA_bond_dist_in_lig);  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}


				// Print a record of all bond angles:
				if (lig_atom>0 && lig_atom<ligNew->nAtoms-1)
				{
					// To calculate the bond angle, use the following formalism:
					//     a · b = ax × bx + ay × by  + az × bz
					//     a · b = |a| × |b| × cos(θ)
					//     θ = acos[(a · b) / (|a| × |b|)]
			    	v1_x = ligNew->atom[lig_atom-1].r[0] - ligNew->atom[lig_atom].r[0];
			    	v1_y = ligNew->atom[lig_atom-1].r[1] - ligNew->atom[lig_atom].r[1];
			    	v1_z = ligNew->atom[lig_atom-1].r[2] - ligNew->atom[lig_atom].r[2];

			    	v2_x = ligNew->atom[lig_atom + 1].r[0] - ligNew->atom[lig_atom].r[0];
			    	v2_y = ligNew->atom[lig_atom + 1].r[1] - ligNew->atom[lig_atom].r[1];
			    	v2_z = ligNew->atom[lig_atom + 1].r[2] - ligNew->atom[lig_atom].r[2];

					v1 = sqrt((v1_x * v1_x)  +  (v1_y * v1_y)  +  (v1_z * v1_z));
					v2 = sqrt((v2_x * v2_x)  +  (v2_y * v2_y)  +  (v2_z * v2_z));

					v1_v2_dot_prod = (v1_x * v2_x)  +  (v1_y * v2_y)  +  (v1_z * v2_z);
					bond_angle = acos(v1_v2_dot_prod / (v1 * v2));
					printf("%f BOND_ANGLE  =  %f \n", v1_v2_dot_prod, bond_angle);
				}
				lig_res += 1;
			}
			findBindingSite(macromolecule, lig, allSites, i);
		} else  // NOTE (DC): then the num of residues in the favorable binding site are 0
		{
			allSites[i].nRes = 0;
		}
		
		printf("Binding site simulation %d. Number of residues in binding site %d\n", i, allSites[i].nRes);
		for (j=0; j<allSites[i].nRes; j++)
		{
			printf("%d ", allSites[i].resid[j]);
		}
		printf("\n");
		
		if (allSites[i].nRes > 10)
		{
			sort(allSites, LC, i);
		}	
	}
	printf("number of stable sites %d\n", stableSites);
	mergeBindingSites(allSites, Ntrials, LC);	
	free(ligNew->atom);
	free(ligTemp->atom);
	free(ligTempOld->atom);
	free(ligNew);
	free(ligTemp);
	free(ligTempOld);
	for (i=0; i<GridSize[0]*GridSize[1]*GridSize[2]; i++) {
		while (Hash[i]) {
			temptr = Hash[i];
			Hash[i] = Hash[i]->next;
			free(temptr);
		}
	}
	free(Hash);
	printf("<-runSurfaceProbeSimulation\n");	
	return;
}


void guessFirstAtomCoordinates(ligand *lig, float *minR, float *maxR) 
{
	float randNum;
	int i;

	//Choosing random x,y,z coordinate within box for 1st atom
	for (i=0; i<3; i++)
	{
		randNum = (float)rand()/RAND_MAX;
		lig->atom[0].r[i] = (maxR[i] - minR[i])*randNum + minR[i];
	}
	//printf("Ufinal ligand first atom guess: %f %f %f\n", lig->atom[0].r[0], lig->atom[0].r[1], lig->atom[0].r[2]);
	return;	
}

void guessSecondAtomCoordinates(ligand *lig) 
{
	float phi, theta;
	
	//Choosing a random point on sphere that is 3.8 Angstroms from 1st atom
	phi = 2*3.14159*((float)rand()/RAND_MAX);
	theta = acos(2.*((float)rand()/RAND_MAX) - 1.);
	lig->atom[1].r[0] = lig->atom[0].r[0] + 3.8 * sin(theta) * cos(phi);
	lig->atom[1].r[1] = lig->atom[0].r[1] + 3.8 * sin(theta) * sin(phi);
	lig->atom[1].r[2] = lig->atom[0].r[2] + 3.8 * cos(theta);
	//printf("ligand second atom guess: %f %f %f\n", lig->atom[1].r[0], lig->atom[1].r[1], lig->atom[1].r[2]);
	//sleep(1);
	
	return;
}

void guessAnotherAtomCoordinates(ligand *lig, int i)
{
	float phi, theta;
	float r[3],r1[3], r2[3];
	int j;
	float T[9];
	//float angle;
	
	//Finding vector that is transformed to z-axis
	for (j=0; j<3; j++)
	{
		r[j] = lig->atom[i-2].r[j] - lig->atom[i-1].r[j];
		//printf("%f ", r[j]);
	}
	//printf("\n");
	
	calculateTransformationMatrixZaxis(r, T);

	//Choosing a random point on quarter-sphere that is 3.8 Angstrom from 2nd atom and z-axis
	phi = 2*3.14159*((float)rand()/RAND_MAX);
	theta = acos(((float)rand()/RAND_MAX) - 1.);
	r1[0] = 3.8 * sin(theta) * cos(phi);
	r1[1] = 3.8 * sin(theta) * sin(phi);
	r1[2] = 3.8 * cos(theta);
	//printf("theta is %f\n", theta*180/3.1416);
	
	//Reverse transformation and get it back to same coordinate system
	r2[0] = r1[0]*T[0] + r1[1]*T[3] + r1[2]*T[6];
	r2[1] = r1[0]*T[1] + r1[1]*T[4] + r1[2]*T[7];
	r2[2] = r1[0]*T[2] + r1[1]*T[5] + r1[2]*T[8];
	//printf("%f\n", sqrt(test2[0]*test2[0] + test2[1]*test2[1] + test2[2]*test2[2]));
	//angle = acos((r2[0]*r[0] + r2[1]*r[1] + r2[2]*r[2])/(3.8*3.8))*180/3.1416;
	//printf("angle is %f\n", angle);
	
	lig->atom[i].r[0] = lig->atom[i-1].r[0] + r2[0];
	lig->atom[i].r[1] = lig->atom[i-1].r[1] + r2[1];
	lig->atom[i].r[2] = lig->atom[i-1].r[2] + r2[2];
	//printf("ligand %d atom guess: %f %f %f\n", i+1, lig->atom[i].r[0], lig->atom[i].r[1], lig->atom[i].r[2]);
	
	//sleep(1);
	return;
}



void calculateEnergyCAmodel(pdbFile *macromolecule, ligand *lig, float *U, bool *V, box *pbc, ligand *ligNew)
{
	float r;
	int my_x, my_y, my_z;
	int j, k, checkkey = -1;
	float v1_x, v1_y, v1_z, v2_x, v2_y, v2_z; //v1, v2, v1_v2_dot_prod, bond_angle; // these are for bond angle calculations
	inode *temptr = NULL;

	// First: calculate the energy associated w/the internal deg. of freedom in the ligand (infinite energy if unfavorable bond or torsion angles!!)
	// Berezovsky: range for favorable bond anlges = 90 - 180 deg
	makeMoleculeWhole(lig, ligNew, pbc);
	for (j=1; j<ligNew->nAtoms-1; j++)
	{
		if (V[j-1] && V[j] && V[j+1]) continue;
		// To calculate the bond angle, use the following formalism:
		//     a · b = ax × bx + ay × by  + az × bz
		//     a · b = |a| × |b| × cos(θ)
		//     θ = acos[(a · b) / (|a| × |b|)]
		v1_x = ligNew->atom[j-1].r[0] - ligNew->atom[j].r[0];
		v1_y = ligNew->atom[j-1].r[1] - ligNew->atom[j].r[1];
		v1_z = ligNew->atom[j-1].r[2] - ligNew->atom[j].r[2];

		v2_x = ligNew->atom[j + 1].r[0] - ligNew->atom[j].r[0];
		v2_y = ligNew->atom[j + 1].r[1] - ligNew->atom[j].r[1];
		v2_z = ligNew->atom[j + 1].r[2] - ligNew->atom[j].r[2];

		// v1 = (v1_x * v1_x)  +  (v1_y * v1_y)  +  (v1_z * v1_z);
		// v2 = (v2_x * v2_x)  +  (v2_y * v2_y)  +  (v2_z * v2_z);

		// v1_v2_dot_prod = (v1_x * v2_x)  +  (v1_y * v2_y)  +  (v1_z * v2_z);
		// bond_angle = acos(v1_v2_dot_prod / (v1 * v2));
		// bond_angle = v1_v2_dot_prod / sqrt(v1 * v2); // delayed sqrt
		if (((v1_x * v2_x)  +  (v1_y * v2_y)  +  (v1_z * v2_z)) >= 0)
		{
				U[lig->nAtoms] += 10000.;
				// Then invalid the angle bits
				// V[lig->atoms - 1 + j] = false;
				return;
		}
		//printf("AA:%f\n", (v1_x * v2_x)  +  (v1_y * v2_y)  +  (v1_z * v2_z));
	}
	// Now -- calculate the energy associated w/the ligand-protein interactions
	// Decide the order first
	int end = lig->nAtoms;
	if (V[end-1]) {
		for (j = lig->nAtoms - 1; j > 0 && V[j]; j--) 
			if (U[j] > 1000.) return;
		end = j + 1;
	}
	// printf("%d !!!<<<<\n", end);
	for (j = 0; j<end; j++)
	{
		if (V[j]) {
			if (U[j] > 1000.) return;
			continue; 
		}
		U[j] = 0.;
		V[j] = true;  // Will be valid;
		my_x = (int)((lig->atom[j].r[0] - minR[0])/5);
		my_y = (int)((lig->atom[j].r[1] - minR[1])/5);
		my_z = (int)((lig->atom[j].r[2] - minR[2])/5);
		// printf("%d   %d    %d\n", (int)((lig->atom[j].r[0] - minR[0])/5), (int)((lig->atom[j].r[1] - minR[1])/5), (int)((lig->atom[j].r[2] - minR[2])/5));
	//	mykey  = (int)((lig->atom[j].r[0] - minR[0])/5) + GridSize[0]*((int)((lig->atom[j].r[1]- minR[1])/5) + GridSize[1]*((int)((lig->atom[j].r[2] - minR[2])/5)));
		for (int k_z=0; k_z<3; k_z++) {
			if (TraverseOrder[k_z] + my_z < 0 || TraverseOrder[k_z] + my_z >= GridSize[2]) continue;
			for (int k_y=0; k_y<3; k_y++) {
				if (TraverseOrder[k_y] + my_y < 0 || TraverseOrder[k_y] + my_y >= GridSize[1]) continue;
				for (int k_x=0; k_x<3; k_x++) {
					if (TraverseOrder[k_x] + my_x < 0 || TraverseOrder[k_x] + my_x >= GridSize[0]) continue;
					checkkey = my_x + TraverseOrder[k_x] + GridSize[0] \
							* (my_y + TraverseOrder[k_y] + (my_z + TraverseOrder[k_z]) * GridSize[1]);
					// if ((checkkey = mykey + TraverseOrder[k_x] + GridSize[0] * (TraverseOrder[k_y] + TraverseOrder[k_z] * GridSize[1])) < 0 || checkkey >= GridSize[0]*GridSize[1]*GridSize[2]) continue;
		/*			printf("%f   %f    %f\n", lig->atom[j].r[0], lig->atom[j].r[1],lig->atom[j].r[2]);*/
					//printf("%d   %d    %d\n", (int)((lig->atom[j].r[0] - minR[0])/5), (int)((lig->atom[j].r[1] - minR[1])/5), (int)((lig->atom[j].r[2] - minR[2])/5));
		//			printf("checkkey: %d\n",checkkey);
					temptr = Hash[checkkey];
					while (temptr) { 			
						k = temptr->poz;
						temptr = temptr->next;
			//		printf("%f   %f    %f    %d\n", macromolecule->atom[k].x, macromolecule->atom[k].y, macromolecule->atom[k].z, checkkey); 
						// NOTE (DC): r is the distance btwn a ligand atom and a residue atom
						r = ((lig->atom[j].r[0] - macromolecule->atom[k].x) \
							* (lig->atom[j].r[0] - macromolecule->atom[k].x)) \
							+ ((lig->atom[j].r[1] - macromolecule->atom[k].y) \
							* (lig->atom[j].r[1] - macromolecule->atom[k].y)) \
							+ ((lig->atom[j].r[2] - macromolecule->atom[k].z) \
							* (lig->atom[j].r[2] - macromolecule->atom[k].z));
						// printf("j %d k %d r %f\n", j, k, r);
	
						if (r < 25.)  // was 8. in orig CA model
						{
							//printf("j %d k %d r %f\n", j, k, r);
							//printf("lig %f %f %f\n", lig->atom[j].r[0], lig->atom[j].r[1], lig->atom[j].r[2]);
							//printf("atom %f %f %f\n", macromolecule->atom[k].x, macromolecule->atom[k].y, macromolecule->atom[k].z);
							if (r >= 12.25)  // was 5.5 in orig CA model
							{
								U[j] -= 0.35;
		//						printf("NEW: %d\n", k);
							} else if (r >= 9.0)  // was 4.5 in orig CA model
							{
								U[j] += 10.;
							} else
							{
								U[j] += 10000.;
								U[lig->nAtoms] = 10000;
								// printf("j %d k %d r %f\n", j, k, r);
								return;
							}
						}
					}
				}
			}
		}
	}

	// Update the total U
	U[lig->nAtoms] = 0.;
	for (int i = 0; i<lig->nAtoms; i++) U[lig->nAtoms] += U[i];
	U[lig->nAtoms] *= 0.75;  /// this *was* in Berezovsky code
}

/*
float calculateEnergyCAmodelold(pdbFile *macromolecule, ligand *lig)
{
        float r;
        int j, k;
        float v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v1, v2, v1_v2_dot_prod, bond_angle; // these are for bond angle calculations
        float U;
        U = 0.;

        // First: calculate the energy associated w/the internal deg. of freedom in the ligand (infinite energy if unfavorable bond or torsion angles!!)
        // Berezovsky: range for favorable bond anlges = 90 - 180 deg
        for (j=0; j<lig->nAtoms; j++)
        {
                if (j>0 && j<lig->nAtoms-1)
                {
                        // To calculate the bond angle, use the following formalism:
                        //     a · b = ax × bx + ay × by  + az × bz
                        //     a · b = |a| × |b| × cos(θ)
                        //     θ = acos[(a · b) / (|a| × |b|)]
                        v1_x = lig->atom[j-1].r[0] - lig->atom[j].r[0];
                        v1_y = lig->atom[j-1].r[1] - lig->atom[j].r[1];
                        v1_z = lig->atom[j-1].r[2] - lig->atom[j].r[2];

                        v2_x = lig->atom[j + 1].r[0] - lig->atom[j].r[0];
                        v2_y = lig->atom[j + 1].r[1] - lig->atom[j].r[1];
                        v2_z = lig->atom[j + 1].r[2] - lig->atom[j].r[2];

                        v1 = sqrt((v1_x * v1_x)  +  (v1_y * v1_y)  +  (v1_z * v1_z));
                        v2 = sqrt((v2_x * v2_x)  +  (v2_y * v2_y)  +  (v2_z * v2_z));

                        v1_v2_dot_prod = (v1_x * v2_x)  +  (v1_y * v2_y)  +  (v1_z * v2_z);
                        bond_angle = acos(v1_v2_dot_prod / (v1 * v2));
                        if (  (bond_angle <= 1.57079) || (bond_angle >= 3.14159)  )
                        {
                                   	U += 10000.;
					break;
                        }
                }
        }

        // Now -- calculate the energy associated w/the ligand-protein interactions
        for (j=0; j<lig->nAtoms; j++)
        {
                for (k=0; k<macromolecule->nAtoms; k++)
                {
                        // NOTE (DC): r is the distance btwn a ligand atom and a residue atom
                        r = ((lig->atom[j].r[0] - macromolecule->atom[k].x) * (lig->atom[j].r[0] - macromolecule->atom[k].x)) + ((lig->atom[j].r[1] - macromolecule->atom[k].y) * (lig->atom[j].r[1] - macromolecule->atom[k].y)) + ((lig->atom[j].r[2] - macromolecule->atom[k].z) * (lig->atom[j].r[2] - macromolecule->atom[k].z));
                        //printf("j %d k %d r %f\n", j, k, r);
                        if (r < 25.)  // was 8. in orig CA model
                        {
                                //printf("j %d k %d r %f\n", j, k, r);
                                //printf("lig %f %f %f\n", lig->atom[j].r[0], lig->atom[j].r[1], lig->atom[j].r[2]);
                                //printf("atom %f %f %f\n", macromolecule->atom[k].x, macromolecule->atom[k].y, macromolecule->atom[k].z);
                                if (r >= 12.25)  // was 5.5 in orig CA model
                                {
                                        U -= 0.35;*/
			/*		printf("OLD: %f  %f  %f  %d\n", macromolecule->atom[k].x, macromolecule->atom[k].y, macromolecule->atom[k].z, k);
					printf("OLD: %f  %f  %f\n", lig->atom[j].r[0], lig->atom[j].r[1], lig->atom[j].r[2]);
					printf("LOC: %d   %d    %d\n", (int)((lig->atom[j].r[0] - minR[0])/5), (int)((lig->atom[j].r[1] - minR[1])/5), (int)((lig->atom[j].r[2] - minR[2])/5)); */
  /*                              } else if (r >= 9.0)  // was 4.5 in orig CA model
                                {
                                        U += 10.;
                                } else
                                {
                                        U += 10000.;
					return U;
                                }
                        }
                }
        }
        U *= 0.75;  /// this *was* in Berezovsky code
        return U;
}*/

void performOneMonteCarloStep(pdbFile *macromolecule, ligand *lig, ligand *ligNew, box *pbc, \
	 	float *minR, float *maxR, float *U, bool *V, ligand **ligTemp, ligand **ligTempOld)
{
	float randNum;
	float Unew[lig->nAtoms + 1];
	bool Vnew[lig->nAtoms + 1];
	affected result;

	result.start = 0;
	result.end = lig->nAtoms;
	
	// ligNew = ligTempOld;

	for (int i = 0; i < lig->nAtoms; i++) {
		Unew[i] = U[i];
		Vnew[i] = V[i];
		for (int j = 0; j < 3; j++) 
			ligNew->atom[i].r[j] = (*ligTempOld)->atom[i].r[j];
	}
	Unew[lig->nAtoms] = U[lig->nAtoms];
	// *nec* --> ***** will need to modify the below to act based on the appropriate number of atoms for a ligand moleucle!! *****
	//Choose type of step - translation, rotation, or internal degrees of freedom
	randNum = (float)rand()/RAND_MAX;

	if (lig->nAtoms >= 4)  // then there are 4 potential deg of freedom that can be changed: transl, rot, bond angle, or torsion angle
	{
		if (randNum <= 0.25)   // then perform a rigid translation of the ligand
		{
			translateCoordinates(lig, ligNew, pbc, 2.0, minR, maxR);
		} 
		else if (randNum <= 0.50)  // then perform a rigid rotation of the ligand
		{
			rotateCoordinates(lig, ligNew, pbc, 3.14159/7.0, minR, maxR);
		}else if (randNum <= 0.75)  // then change one bond angle
		{
			//guessAngles(lig, ligNew, pbc, minR, maxR);   // *was* in original
			result = changeOneBondAngle(lig, ligNew, pbc, minR, maxR);  // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / 
		}
		else // then change one torsion angle
		{
			result = changeOneTorsionAngle(lig, ligNew, pbc, minR, maxR);
		}
	}
	else if (lig->nAtoms == 3) // then there are 3 deg of freedom that may be explored: transl, rot, and bond angle:
	{
		if (randNum <= 0.33333)   // then perform a rigid translation of the ligand
		{
			translateCoordinates(lig, ligNew, pbc, 2.0, minR, maxR);
		} 
		else if (randNum <= 0.66666)  // then perform a rigid rotation of the ligand
		{
			rotateCoordinates(lig, ligNew, pbc, 3.14159/7.0, minR, maxR);
		}else // then change one bond angle
		{
			//guessAngles(lig, ligNew, pbc, minR, maxR);   // *was* in original
			changeOneBondAngle(lig, ligNew, pbc, minR, maxR);  // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / 
		}		
	}
	else if (lig->nAtoms == 2) // then there are 2 deg of freedom that may be explored: transl & rot
	{
		if (randNum <= 0.5)   // then perform a rigid translation of the ligand
		{
			translateCoordinates(lig, ligNew, pbc, 2.0, minR, maxR);
		} 
		else  // then perform a rigid rotation of the ligand
		{
			rotateCoordinates(lig, ligNew, pbc, 3.14159/7.0, minR, maxR);
		}
	}
	else  // then the lig has one atom only (never the case in Berezovsky's paper, but let's just be general and assume that such is possible) -- only rigid translations allowed
	{
		translateCoordinates(lig, ligNew, pbc, 2.0, minR, maxR);		
	}

	//printf("SE: %d   %d\n", result.start, result.end);
	for (int i=result.start; i < result.end; i++) Vnew[i] = false;

	calculateEnergyCAmodel(macromolecule, ligNew, Unew, Vnew, pbc, *ligTemp);
	//float Uold = calculateEnergyCAmodelold(macromolecule, ligNew);
	
	//printf("U: %f	%f\n", Unew[lig->nAtoms], Uold);
	if (Unew[lig->nAtoms] < 1000.) doMetroPolis(U, Unew, lig, ligNew, V, Vnew, ligTemp, ligTempOld);
	// printf("OLDTEMP: %f\n", (*ligTempOld)->atom[0].r[0]);
	return;
}


affected changeOneTorsionAngle(ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR)
{
	makeMoleculeWhole(lig, ligNew, pbc);

	// Initialize some variables
	float v1_x, v1_y, v1_z; // this will represent the vector about which the torsion angle should be changed
	float origin_x, origin_y, origin_z; // the atom that places the vector in 3d space
	affected result;
	result.start = 0;
	result.end = lig->nAtoms;
	
	//Start by selecting a random torsion angle that should be changed (torsion_index_to_be_changed)  
	float randNum = (float)rand()/RAND_MAX;  // *nec* to *change*??
	int num_available_torsion_angles = ligNew->nAtoms - 3;
	int torsion_angle_index_to_be_changed = 0;
	if (num_available_torsion_angles > 1) {
		float bin = 1.0 / num_available_torsion_angles;
		float j = 1.0;
		for (int torsion_angle_index=0; torsion_angle_index < num_available_torsion_angles; torsion_angle_index++)
		{
			if (randNum <= (j * bin))
			{
				torsion_angle_index_to_be_changed = torsion_angle_index;
				break;
			}
			j += 1.0;
		}
		randNum = (float)rand()/RAND_MAX;
	}
	//printf("from changeOneTorsionAngle:  chosen torsion_angle_index_to_be_changed  =  %d \n", torsion_angle_index_to_be_changed);

	// Select a random torsion angle change, which (in accordance w/Berezovsky), is  
	// chosen from a uniform distribution between -PI/12 to +PI/12
    // float delta_angle = update_pivot_range * ((float)rand()/RAND_MAX - 0.5);  // *nec* to change rand number?? *******
    float delta_angle = update_pivot_range * (randNum - 0.5);
    //printf("from changeOneTorsionAngle:  delta_angle= %f \n", delta_angle);

	// Now that the specific torsion angle and delta_angle are determined, implement the angle change
	// by rotating the ligand on one side of the the torsion vector by the amount delta_angle
    int middle_torsion_angle_index = num_available_torsion_angles / 2;
    //printf("from changeOneTorsionAngle:  middle_torsion_angle_index = %d \n", middle_torsion_angle_index);

    // Note that the axis of rotation about the torsion_angle_index_to_be_changed is
    // just the atom-atom bond vector defining the torsion angle;
    v1_x = ligNew->atom[torsion_angle_index_to_be_changed + 2].r[0] - ligNew->atom[torsion_angle_index_to_be_changed + 1].r[0];
    v1_y = ligNew->atom[torsion_angle_index_to_be_changed + 2].r[1] - ligNew->atom[torsion_angle_index_to_be_changed + 1].r[1];
    v1_z = ligNew->atom[torsion_angle_index_to_be_changed + 2].r[2] - ligNew->atom[torsion_angle_index_to_be_changed + 1].r[2];

    // The origin is just the first atom defining this vector (ie, it places the 
    // vector in 3D space):
	origin_x = ligNew->atom[torsion_angle_index_to_be_changed + 1].r[0];
	origin_y = ligNew->atom[torsion_angle_index_to_be_changed + 1].r[1];
	origin_z = ligNew->atom[torsion_angle_index_to_be_changed + 1].r[2];

    if (torsion_angle_index_to_be_changed < middle_torsion_angle_index || (torsion_angle_index_to_be_changed == middle_torsion_angle_index \
		&& (float)rand()/RAND_MAX <= 0.5)) // Then rotate the ligand atoms to the "LEFT" of the torsion_angle_index_to_be_changed by delta_angle radians
    {
		int atm;
		for (atm=0; atm < (torsion_angle_index_to_be_changed + 1); atm++)
		{
			// rotate atom atm by delta_angle about the central torsion vector
			float new_x = x_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, v1_x, v1_y, v1_z, delta_angle);
			float new_y = y_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, v1_x, v1_y, v1_z, delta_angle);
			float new_z = z_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, v1_x, v1_y, v1_z, delta_angle);

			//  Below are print statements to see that new_* were calculated correctly 
			//    Note that results have been checked against an online calculator:
			//    http://twist-and-shout.appspot.com/ 
			/*
			printf("\n\nTOR_ROTATIONS \n");
			printf("TOR_ROTATIONS: ORIGIN:       x = %f    y = %f    z = %f   \n", origin_x, origin_y, origin_z);
			printf("TOR_ROTATIONS: CROSS:        x = %f    y = %f    z = %f   \n", v1_x, v1_y, v1_z);
			printf("TOR_ROTATIONS: ORIGINAL PT:  x = %f    y = %f    z = %f   \n", ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2]);
			printf("TOR_ROTATIONS: DELTA_ANGLE:  %f  \n", delta_angle);
			printf("TOR_ROTATIONS: NEW x = %f  \n", new_x);
			printf("TOR_ROTATIONS: NEW y = %f  \n", new_y);
			printf("TOR_ROTATIONS: NEW z = %f  \n", new_z);
			printf("TOR_ROTATIONS: DIPLACEMENT: =  %f  \n", displacement_frm_tor_bond_angle_change);
			printf("TOR_ROTATIONS \n\n");
			*/
			ligNew->atom[atm].r[0] = new_x;
			ligNew->atom[atm].r[1] = new_y;
			ligNew->atom[atm].r[2] = new_z;
		}
		result.end = atm;
    }
    else // rotate the ligand atoms to the "RIGHT" of the torsion_angle_index_to_be_changed by delta_angle radians
    {
		int atm;
		for (atm=(torsion_angle_index_to_be_changed+3); atm < ligNew->nAtoms; atm++)
		{
			// rotate atom atm by delta_angle about the central torsion vector
			float new_x = x_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, v1_x, v1_y, v1_z, delta_angle);
			float new_y = y_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, v1_x, v1_y, v1_z, delta_angle);
			float new_z = z_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, v1_x, v1_y, v1_z, delta_angle);

			//  Below are print statements to see that new_* were calculated correctly 
			//  Note that results have been checked against an online calculator:
			//  http://twist-and-shout.appspot.com/
			/*
			printf("\n\nTOR_ROTATIONS \n");
			printf("TOR_ROTATIONS: ORIGIN:       x = %f    y = %f    z = %f   \n", origin_x, origin_y, origin_z);
			printf("TOR_ROTATIONS: CROSS:        x = %f    y = %f    z = %f   \n", v1_x, v1_y, v1_z);
			printf("TOR_ROTATIONS: ORIGINAL PT:  x = %f    y = %f    z = %f   \n", ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2]);
			printf("TOR_ROTATIONS: DELTA_ANGLE:  %f  \n", delta_angle);
			printf("TOR_ROTATIONS: NEW x = %f  \n", new_x);
			printf("TOR_ROTATIONS: NEW y = %f  \n", new_y);
			printf("TOR_ROTATIONS: NEW z = %f  \n", new_z);
			printf("TOR_ROTATIONS: DIPLACEMENT: =  %f  \n", displacement_frm_tor_bond_angle_change);
			printf("TOR_ROTATIONS \n\n");
			*/
			ligNew->atom[atm].r[0] = new_x;
			ligNew->atom[atm].r[1] = new_y;
			ligNew->atom[atm].r[2] = new_z;
		}
		result.start = atm - 1;
    }
	transformCoordinatesToSameBox(ligNew, minR, maxR, pbc);
	return result;
}


affected changeOneBondAngle(ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR)
{
	
	/*float a = ligNew->atom[0].r[0];
	float b = ligNew->atom[3].r[1];
	
	makeMoleculeWhole(lig, ligNew, pbc);

	if (a != ligNew->atom[0].r[0] || b != ligNew->atom[3].r[1]) 
		printf("WARNING! %f:%f	%f:%f\n", a, ligNew->atom[0].r[0], b, ligNew->atom[3].r[1]);*/

	// initialize some variables 
	float v1_x, v1_y, v1_z;
	//float v1[3];
	float v2_x, v2_y, v2_z;
	//float v2[3];
	float cross_x, cross_y, cross_z;
	float origin_x, origin_y, origin_z;
	affected result;
	result.start = 0;
	result.end = lig->nAtoms;


	//Start by selecting a random bond angle that should be changed (angle_index_to_be_changed)  
	float randNum = (float)rand()/RAND_MAX;  // *nec* to *change*??
	int num_available_angles = ligNew->nAtoms - 2;
	int angle_index_to_be_changed = 0;
	if (num_available_angles > 1 ) {
		float bin = 1.0 / num_available_angles;
		float j = 1.0;
		for (int angle_index=0; angle_index < num_available_angles; angle_index++)
		{
			if (randNum <= (j * bin))
			{
				angle_index_to_be_changed = angle_index;
				break;
			}
			j += 1.0;
		}
		randNum = (float)rand()/RAND_MAX;
	}
	//printf("from changeOneBondAngle:  chosen angle index  =  %d \n", angle_index_to_be_changed);


	// Select a random angle change, which (in accordance w/Berezovsky), is  
	// chosen from a uniform distribution between -PI/6 to +PI/6
    float delta_angle = change_in_bondangle_range * (randNum - 0.5);  // *nec* to change rand number?? *******
    //printf("from changeOneBondAngle:  change_in_bondangle_range = %f \n", delta_angle);


	// Now that the specific angle and delta_angle are determined, implement the angle change
	// by rotating the ligand on one side of the the angle by the amount delta_angle
    int middle_angle_index = num_available_angles / 2;
    //printf("from changeOneBondAngle:  middle_angle_index = %d \n", middle_angle_index);
    if (angle_index_to_be_changed < middle_angle_index 
	|| (angle_index_to_be_changed == middle_angle_index && rand()/RAND_MAX < 0.5 ))
	// Then rotate the ligand atoms to the "LEFT" of the angle_index_to_be_changed by delta_angle radians
    {
    	// In order to determine the axis of rotation about the angle_index_to_be_changed, take 
    	// the cross product of the vectors that define the angle_index_to_be_changed:
    	v1_x = ligNew->atom[angle_index_to_be_changed].r[0] - ligNew->atom[angle_index_to_be_changed + 1].r[0];
    	v1_y = ligNew->atom[angle_index_to_be_changed].r[1] - ligNew->atom[angle_index_to_be_changed + 1].r[1];
    	v1_z = ligNew->atom[angle_index_to_be_changed].r[2] - ligNew->atom[angle_index_to_be_changed + 1].r[2];

    	v2_x = ligNew->atom[angle_index_to_be_changed + 2].r[0] - ligNew->atom[angle_index_to_be_changed + 1].r[0];
    	v2_y = ligNew->atom[angle_index_to_be_changed + 2].r[1] - ligNew->atom[angle_index_to_be_changed + 1].r[1];
    	v2_z = ligNew->atom[angle_index_to_be_changed + 2].r[2] - ligNew->atom[angle_index_to_be_changed + 1].r[2];

		cross_x = v1_y*v2_z - v1_z*v2_y;
		cross_y = v1_z*v2_x - v1_x*v2_z;
		cross_z = v1_x*v2_y - v1_y*v2_x;
		//printf("crosses: %f   %f   %f  \n", cross_x, cross_y, cross_z);
		origin_x = ligNew->atom[angle_index_to_be_changed + 1].r[0];
		origin_y = ligNew->atom[angle_index_to_be_changed + 1].r[1];
		origin_z = ligNew->atom[angle_index_to_be_changed + 1].r[2];

		//  Now, rotate the ligand atoms to the "LEFT" of the angle_index_to_be_changed 
		//  by delta_angle radians (perform rotation around line defined by cross_*)
		int atm;
		for (atm=0; atm < (angle_index_to_be_changed + 1); atm++)
		{
			// rotate atom atm by delta_angle about the line cross_*
			float new_x = x_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, cross_x, cross_y, cross_z, delta_angle);
			float new_y = y_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, cross_x, cross_y, cross_z, delta_angle);
			float new_z = z_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, cross_x, cross_y, cross_z, delta_angle);

			//  Below are print statements to see that new_* were calculated correctly 
			//    Note that results have been checked against an online calculator:
			//    http://twist-and-shout.appspot.com/ 
			/*
			printf("\n\nROTATIONS \n");
			printf("ROTATIONS: ORIGIN:       x = %f    y = %f    z = %f   \n", origin_x, origin_y, origin_z);
			printf("ROTATIONS: CROSS:        x = %f    y = %f    z = %f   \n", cross_x, cross_y, cross_z);
			printf("ROTATIONS: ORIGINAL PT:  x = %f    y = %f    z = %f   \n", ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2]);
			printf("ROTATIONS: DELTA_ANGLE:  %f  \n", delta_angle);
			printf("ROTATIONS: NEW x = %f  \n", new_x);
			printf("ROTATIONS: NEW y = %f  \n", new_y);
			printf("ROTATIONS: NEW z = %f  \n", new_z);
			printf("ROTATIONS: DIPLACEMENT: =  %f  \n", displacement_frm_bond_angle_change);
			printf("ROTATIONS \n\n");
			*/
			ligNew->atom[atm].r[0] = new_x;
			ligNew->atom[atm].r[1] = new_y;
			ligNew->atom[atm].r[2] = new_z;
		}
		result.end = atm;
    }
    else // rotate the ligand atoms to the "RIGHT" of the angle_index_to_be_changed by delta_angle radians
    {
    	// In order to determine the axis of rotation about the angle_index_to_be_changed, take 
    	// the cross product of the vectors that define the angle_index_to_be_changed:
    	v1_x = ligNew->atom[angle_index_to_be_changed].r[0] - ligNew->atom[angle_index_to_be_changed + 1].r[0];
    	v1_y = ligNew->atom[angle_index_to_be_changed].r[1] - ligNew->atom[angle_index_to_be_changed + 1].r[1];
    	v1_z = ligNew->atom[angle_index_to_be_changed].r[2] - ligNew->atom[angle_index_to_be_changed + 1].r[2];

    	v2_x = ligNew->atom[angle_index_to_be_changed + 2].r[0] - ligNew->atom[angle_index_to_be_changed + 1].r[0];
    	v2_y = ligNew->atom[angle_index_to_be_changed + 2].r[1] - ligNew->atom[angle_index_to_be_changed + 1].r[1];
    	v2_z = ligNew->atom[angle_index_to_be_changed + 2].r[2] - ligNew->atom[angle_index_to_be_changed + 1].r[2];

		cross_x = v1_y*v2_z - v1_z*v2_y;
		cross_y = v1_z*v2_x - v1_x*v2_z;
		cross_z = v1_x*v2_y - v1_y*v2_x;
		//printf("crosses: %f   %f   %f  \n", cross_x, cross_y, cross_z);
		origin_x = ligNew->atom[angle_index_to_be_changed + 1].r[0];
		origin_y = ligNew->atom[angle_index_to_be_changed + 1].r[1];
		origin_z = ligNew->atom[angle_index_to_be_changed + 1].r[2];

		//  Now, rotate the ligand atoms to the "RIGHT" of the angle_index_to_be_changed 
		//  by delta_angle radians (perform rotation around line defined by cross_*)
		int atm;
		for (atm=(angle_index_to_be_changed+1); atm < ligNew->nAtoms; atm++)
		{
			// rotate atom atm by delta_angle about the line cross_*
			float new_x = x_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, cross_x, cross_y, cross_z, delta_angle);
			float new_y = y_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, cross_x, cross_y, cross_z, delta_angle);
			float new_z = z_AfterRotation(ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2], origin_x, origin_y, origin_z, cross_x, cross_y, cross_z, delta_angle);
			//  Below are print statements to see that new_* were calculated correctly 
			//  Note that results have been checked against an online calculator:
			//  http://twist-and-shout.appspot.com/
			
			/*
			printf("\n\nROTATIONS \n");
			printf("ROTATIONS: ORIGIN:       x = %f    y = %f    z = %f   \n", origin_x, origin_y, origin_z);
			printf("ROTATIONS: CROSS:        x = %f    y = %f    z = %f   \n", cross_x, cross_y, cross_z);
			printf("ROTATIONS: ORIGINAL PT:  x = %f    y = %f    z = %f   \n", ligNew->atom[atm].r[0], ligNew->atom[atm].r[1], ligNew->atom[atm].r[2]);
			printf("ROTATIONS: DELTA_ANGLE:  %f  \n", delta_angle);
			printf("ROTATIONS: NEW x = %f  \n", new_x);
			printf("ROTATIONS: NEW y = %f  \n", new_y);
			printf("ROTATIONS: NEW z = %f  \n", new_z);
			printf("ROTATIONS \n\n");
			*/
			ligNew->atom[atm].r[0] = new_x;
			ligNew->atom[atm].r[1] = new_y;
			ligNew->atom[atm].r[2] = new_z;
		}
		result.start = atm - 1;
    }
	transformCoordinatesToSameBox(ligNew, minR, maxR, pbc);
	return result;
}


void guessAngles(ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR)
{
	int j;
	
	makeMoleculeWhole(lig, ligNew, pbc);
	for (j=2; j<lig->nAtoms; j++)
	{
		guessAnotherAtomCoordinates(ligNew, j);
	}
	
	return;
}


//  NOTE: formulas for all this are taken from:  http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
//  NOTE: this was checked using: http://twist-and-shout.appspot.com/
float x_AfterRotation(float x_orig, float y_orig, float z_orig, float origin_x, float origin_y, float origin_z, float cross_x, float cross_y, float cross_z, float delta_angle)
{
	float a = origin_x;
	float b = origin_y;
	float c = origin_z;

	float u = cross_x;
	float v = cross_y;
	float w = cross_z;

	float x = x_orig;
	float y = y_orig;
	float z = z_orig;

	float L = (u*u) + (v*v) + (w*w);

	float x_new = (  (a * ((v*v) + (w*w))  -  u * (b*v + c*w - u*x - v*y - w*z)) * (1 - cos(delta_angle)) + L*x*cos(delta_angle) + sqrt(L)*( -(c*v) + b*w - w*y + v*z )*sin(delta_angle)  ) / L; 
	return x_new;
}

//  NOTE: formulas for all this are taken from:  http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
//  NOTE: this was checked using: http://twist-and-shout.appspot.com/
float y_AfterRotation(float x_orig, float y_orig, float z_orig, float origin_x, float origin_y, float origin_z, float cross_x, float cross_y, float cross_z, float delta_angle)
{
	float a = origin_x;
	float b = origin_y;
	float c = origin_z;

	float u = cross_x;
	float v = cross_y;
	float w = cross_z;

	float x = x_orig;
	float y = y_orig;
	float z = z_orig;

	float L = (u*u) + (v*v) + (w*w);

	float y_new = (  (b * ((u*u) + (w*w))  -  v * (a*u + c*w - u*x - v*y - w*z)) * (1 - cos(delta_angle)) + L*y*cos(delta_angle) + sqrt(L)*( c*u - a*w + w*x - u*z )*sin(delta_angle)  ) / L; 
	return y_new;
}

//  NOTE: formulas for all this are taken from:  http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
//  NOTE: this was checked using: http://twist-and-shout.appspot.com/
float z_AfterRotation(float x_orig, float y_orig, float z_orig, float origin_x, float origin_y, float origin_z, float cross_x, float cross_y, float cross_z, float delta_angle)
{
	float a = origin_x;
	float b = origin_y;
	float c = origin_z;

	float u = cross_x;
	float v = cross_y;
	float w = cross_z;

	float x = x_orig;
	float y = y_orig;
	float z = z_orig;

	float L = (u*u) + (v*v) + (w*w);

	float z_new = (  (c * ((u*u) + (v*v))  -  w * (a*u + b*v - u*x - v*y - w*z)) * (1 - cos(delta_angle)) + L*z*cos(delta_angle) + sqrt(L)*( -(b*u) + a*v - v*x + u*y )*sin(delta_angle)  ) / L;
	return z_new;
}


void doMetroPolis(float *U, float *Unew, ligand *lig, ligand *ligNew, bool *V, bool *Vnew, ligand **ligTemp, ligand **ligTempOld)
{
	int j, k;
	float randNum;
	
	// printf("AA: %f %f\n", lig->atom[0].r[0], (*ligTempOld)->atom[0].r[0]);
	randNum = (float)rand()/RAND_MAX;	
	if ((Unew[lig->nAtoms] <= U[lig->nAtoms]) || (randNum < exp(-(Unew[lig->nAtoms] - U[lig->nAtoms]))))
	{
		//Accept step
		//printf("Unew %f %f %f < %f\n", Unew, U[i-1], randNum, exp(-(Unew - U[i-1])));
		for (j=0; j<lig->nAtoms; j++)
		{
			V[j] = Vnew[j];
			U[j] = Unew[j];
			for (k=0; k<3; k++)
			{
				lig->atom[j].r[k] = ligNew->atom[j].r[k];
			}
		}
		U[lig->nAtoms] = Unew[lig->nAtoms];
		ligand *swap;
		swap = *ligTempOld;
		*ligTempOld = *ligTemp;
		*ligTemp = swap;
	} 
	return;
}


void calculateCenterOfMass(ligand *lig, float *Rcom)
{
	int i, j;
	
	for (j=0; j<3; j++)
	{
		Rcom[j] = 0.;
	}
	
	for (i=0; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			Rcom[j] += lig->atom[i].r[j];
			//printf("j %d Rcom[j] %f lig->atom[i].r[j] %f\n", j, Rcom[j], lig->atom[i].r[j]);
		}
	}
	
	for (j=0; j<3; j++)
	{
		Rcom[j] = Rcom[j]/lig->nAtoms;
		//printf("%d %f ", lig->nAtoms, Rcom[j]);
	}
	//printf("\n");
	
	return;
}



void findBindingSite(pdbFile *macromolecule, ligand *lig, bindingSite *allSites, int n)
{
	int i, j, space;
	// float r[macromolecule->nAtoms];

	allSites[n].nRes = 0;
	space = 8; // Just a prediction
	if ((allSites[n].resid = (int *) malloc(space*sizeof(int))) == NULL)
        {
                printf("No space to save binding sites\n");
                exit(1);
        }
        if ((allSites[n].residue = (int *) malloc(space*sizeof(int))) == NULL)
        {
                printf("No space to save binding sites\n");
                exit(1);
        }

	
	for (i = 0; i<macromolecule->nAtoms; i++)
	{
		for (j=0; j<lig->nAtoms; j++)
		{
			if (((macromolecule->atom[i].x - lig->atom[j].r[0])*(macromolecule->atom[i].x - lig->atom[j].r[0]) \
				+ (macromolecule->atom[i].y - lig->atom[j].r[1])*(macromolecule->atom[i].y - lig->atom[j].r[1]) \
				+ (macromolecule->atom[i].z - lig->atom[j].r[2])*(macromolecule->atom[i].z - lig->atom[j].r[2])) \
				< 64.) { // was:  if (r[i][j] < 8.)
				if (++allSites[n].nRes > space) {
					space = space*2; 
					if ((allSites[n].resid = (int *) realloc(allSites[n].resid, space*sizeof(int))) == NULL) {
						printf("No space to save binding sites\n");
						exit(1);
					}
					if ((allSites[n].residue = (int *) realloc(allSites[n].residue, space*sizeof(int))) == NULL) {
						printf("No space to save binding sites\n");
						exit(1);
					}
				}
				allSites[n].resid[allSites[n].nRes - 1] = macromolecule->atom[i].resID;
				allSites[n].residue[allSites[n].nRes - 1] = i;
				// allSites[n].nRes++;
				break;
				//j = lig->nAtoms; // NOTE (DC): ok, we found that @ least one residue is in contact w/residue -- so terminate j loop  --- count this residue as being in contact w/ligand (add to residue count)
			}
		}
	}

	if (allSites[n].nRes == 0) {
		free(allSites[n].resid);
		free(allSites[n].residue);
		return;
	}
	else {		
		allSites[n].resid = (int *) realloc(allSites[n].resid, allSites[n].nRes*sizeof(int));
		allSites[n].residue = (int *) realloc(allSites[n].residue, allSites[n].nRes*sizeof(int));
	}

	/*if ((allSites[n].resid = (int *) malloc(allSites[n].nRes*sizeof(int))) == NULL)
	{
		printf("No space to save binding sites\n");
		exit(1);
	}
	if ((allSites[n].residue = (int *) malloc(allSites[n].nRes*sizeof(int))) == NULL)
	{
		printf("No space to save binding sites\n");
		exit(1);
	}
	
	k = 0;
	for (i = 0; i<macromolecule->nAtoms; i++)
	{
	//	for (j=0; j<lig->nAtoms; j++)
	//	{
			if (r[i] < 64.)  // was:  if (r[i][j] < 8.)
			{
				allSites[n].resid[k] = macromolecule->atom[i].resID;  	// NOTE (DC): "resid" refers to the actual residue index as it is assigned in the PDB file (doesn't start @ 0)
				allSites[n].residue[k] = i; 							// NOTE (DC): "residue" just refers to an artificial index to keep track of the residues within this program (starts at 0)
	//			j = lig->nAtoms;
				k++;
			}
	//	}
	}*/
	
	// allSites[n].numSprings = calculateNumberOfSprings(macromolecule, lig, allSites, n);
	// printf("Number of springs in site %d are %d\n", n, allSites[n].numSprings);

	allSites[n].numSprings = 0;
	if (allSites[n].nRes < 2) return;
	int numSprings = allSites[n].nRes * (allSites[n].nRes - 1) / 2; // Could use Prediction! Can be optimized!
	if ((allSites[n].spring = (edgePtr *)malloc(numSprings*sizeof(edgePtr))) == NULL)
	{
		printf("Error: No space to save the spring in binding sites\n");
		exit(1);
	}
	
	setSpringsInBindingSite(macromolecule, lig, allSites, n);
	
	// Resize!
	if (allSites[n].numSprings == 0) 
		free(allSites[n].spring);
	else 
		allSites[n].spring = (edgePtr*)realloc(allSites[n].spring, allSites[n].numSprings*sizeof(edgePtr));
	return;
}

void sort(bindingSite *allSites, float *LC, int nSim)
{
	int i, j, swap;
	
	for (i=0; i<allSites[nSim].nRes; i++) 
	{
		for (j=(i+1); j<allSites[nSim].nRes; j++) 
		{
			if (LC[allSites[nSim].residue[i]] < LC[allSites[nSim].residue[j]])
			{
				swap = allSites[nSim].residue[i];
				allSites[nSim].residue[i] = allSites[nSim].residue[j];
				allSites[nSim].residue[j] = swap;
				swap = allSites[nSim].resid[i];
				allSites[nSim].resid[i] = allSites[nSim].resid[j];
				allSites[nSim].resid[j] = swap;
			}
		}
	}
	
	printf("Sorted LC\n");
	for (i=0; i<allSites[nSim].nRes; i++)
		printf("%d ", allSites[nSim].residue[i]);
	printf("\n");
	
	return;		
}

void mergeBindingSites(bindingSite *allSites, int Ntrials, float *LC)
{
	int i, j, k;
	// float Sim[Ntrials][Ntrials];
	float **Sim;
	float max = 0.;
	int mergingSites;
	
	Sim = calloc(sizeof(float*), Ntrials);
	for (i=0; i<Ntrials; i++) 
		Sim[i] = calloc(sizeof(float), Ntrials);

	printf("->merging Sites\n");
	
	

	k = 0;
	for (i=0; i<Ntrials; i++)
	{
		if (allSites[i].nRes > 0)
		{
			mergingSites = 1;
			while (mergingSites == 1)
			{
				for (j=(i+1); j<Ntrials; j++)
				{
					//printf("i %d j %d\n", i, j);
					//Calculate Similarities
					if (allSites[j].nRes > 0)
					{
						//printf("j %d nres %d\n", j, allSites[j].nRes);
						Sim[i][j] = calculateSimilarity(allSites, i, j);
						Sim[j][i] = Sim[i][j];
					} else
					{
						//printf("j %d nres %d\n", j, allSites[j].nRes);
						Sim[i][j] = 0.;
						Sim[j][i] = 0.;
					}
				}
				max = 0.;
				for (j=(i+1); j<Ntrials; j++)
				{
					if (max < Sim[i][j]) 
					{
						//printf("j %d Sim %f\n", j, Sim[i][j]);
						max = Sim[i][j];
						k = j;
					}
				}
				if (max > 0.7)
				{
					//printf("k %d Sim %f\n", k, Sim[i][k]);
					mergeSites(allSites, LC, i, k);
					Sim[i][k] = 0.;
					Sim[k][i] = 0.;
				}
				else
				{
					mergingSites = 0;
				}
			}	 
		} else 
		{
			for (j=(i+1); j<Ntrials; j++)
			{
				Sim[i][j] = 0.;
				Sim[j][i] = 0.;
			}
		}
		Sim[i][i] = 0.;
	}
	
	printf("<-merging Sites\n");

	for (i=0; i<Ntrials; i++)
		free(Sim[i]);
	free(Sim);

	return;
}				

float calculateSimilarity(bindingSite *allSites, int site1, int site2)
{
	int i, j;
	int n1, n2;
	int intElements = 0;
	float sim;
	
	if (allSites[site1].nRes > 10)
	{
		n1 = 10;
	} else
	{
		n1 = allSites[site1].nRes;
	}
	
	if (allSites[site2].nRes > 10)
	{
		n2 = 10;
	} else
	{
		n2 = allSites[site2].nRes;
	}
	
	for (i=0; i<n1; i++)
	{
		for (j=0; j<n2; j++)
		{
			if (allSites[site1].residue[i] == allSites[site2].residue[j])
			{
				intElements++;
				j = n2;
			}
		}
	}
	
	//printf("intElements %d\n", intElements);
	sim = (float)intElements/(n1 + n2 - intElements);
	
	return sim;
}

void mergeSites(bindingSite *allSites, float *LC, int i, int j)
{
	int k, l, n;
	int *residues;
	int *resIDs;
	int intElements = 0;
	
	printf("Merging sites %d & %d \n",i ,j);
	for (k=0; k<allSites[i].nRes; k++)
	{
		for (l=0; l<allSites[j].nRes; l++)
		{
			if (allSites[i].residue[k] == allSites[j].residue[l])
			{
				intElements++;
				l = allSites[j].nRes;
			}
		}
	}
	
	int totalElements = allSites[i].nRes + allSites[j].nRes - intElements;
	
	if ((resIDs = (int *) malloc(totalElements*sizeof(int))) == NULL)
	{
		printf("Error: Couldnt merge elements\n");
		exit(1);
	}
	if ((residues = (int *) malloc(totalElements*sizeof(int))) == NULL)
	{
		printf("Error: Couldnt merge elements\n");
		exit(1);
	}
	
	for (k=0; k<allSites[i].nRes; k++)
	{
		resIDs[k] = allSites[i].resid[k];
		residues[k] = allSites[i].residue[k];
	}
	n = allSites[i].nRes;
	for (k=0; k<allSites[j].nRes; k++)
	{
		intElements = 0;
		for (l=0; l<allSites[i].nRes; l++)
		{
			if (allSites[j].residue[k] == residues[l])
			{
					intElements = 1;
			}
		}
		if (intElements == 0)
		{
			residues[n] = allSites[j].residue[k];
			resIDs[n] = allSites[j].resid[k];
			n++;
		}
	}
	//printf("totalElements = %d n = %d\n", totalElements, n);
	
	allSites[i].nRes = totalElements;
	free(allSites[i].resid);
	free(allSites[i].residue);
	//printf("Got here\n");
	allSites[j].nRes = 0;
	free(allSites[j].resid);
	free(allSites[j].residue);
	if ((allSites[i].resid = (int *) malloc(totalElements*sizeof(int))) == NULL)
	{
		printf("Error: no space to merge sites\n");
		exit(1);
	}
	if ((allSites[i].residue = (int *) malloc(totalElements*sizeof(int))) == NULL)
	{
		printf("Error: no space to merge sites\n");
		exit(1);
	}
	for (k=0; k<totalElements; k++)
	{
		allSites[i].resid[k] = resIDs[k];
		allSites[i].residue[k] = residues[k];
	}
	sort(allSites, LC, i);
	
	mergeSprings(allSites, i, j);
	free(resIDs);
	free(residues);
	
	return;
}

int getNumberFinalSites(bindingSite *allSites, int Ntrials)
{
	int numSites, i;
	
	numSites = 0;
	for (i=0; i<Ntrials; i++)
	{
		if (allSites[i].nRes != 0)
		{
			numSites++;
		}
	}
	printf("Total number of merged sites are %d\n", numSites);
	return numSites;
}

void copyFinalSites(bindingSite *allSites, bindingSite *mergedSites, int Ntrials)
{
	int i, j, k;
	
	j = 0;
	for (i=0; i<Ntrials; i++)
	{
		if (allSites[i].nRes != 0)
		{
			mergedSites[j].nRes = allSites[i].nRes;
			//printf("nRes %d %d\n", mergedSites[j].nRes, allSites[i].nRes);
			if ((mergedSites[j].resid = (int*) malloc(allSites[i].nRes*sizeof(int))) == NULL)
			{
				printf("Error: No space to save merged sites\n");
				exit(1);
			}
			if ((mergedSites[j].residue = (int*) malloc(allSites[i].nRes*sizeof(int))) == NULL)
			{
				printf("Error: No space to save merged sites\n");
				exit(1);
			}
			for (k=0; k<allSites[i].nRes; k++)
			{
				mergedSites[j].residue[k] = allSites[i].residue[k];
				mergedSites[j].resid[k] = allSites[i].resid[k];
				// printf("Merged sites resid and resname are: %d %s   ",  allSites[i].resid[k]);     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}			
			//printf("mergedSites nRes %d %d\n", mergedSites[j].nRes, allSites[i].nRes);
			
			mergedSites[j].numSprings = allSites[i].numSprings;
			//printf("mergedSites nRes %d\n", mergedSites[j].nRes);
			if (mergedSites[j].numSprings > 0 && (mergedSites[j].spring = (edgePtr*) malloc(allSites[i].numSprings*sizeof(edgePtr))) == NULL)
			{
				printf("Error: No space to save merged sites\n");
				exit(1);
			}
			
			for (k=0; k<allSites[i].numSprings; k++)
			{
				mergedSites[j].spring[k] = (edgePtr)malloc(sizeof(EdgeStr));
				mergedSites[j].spring[k]->res1 = allSites[i].spring[k]->res1;
				mergedSites[j].spring[k]->res2 = allSites[i].spring[k]->res2;
			}
			j++;
		}
	}
	
	//printf("j %d \n", j);

	/*for (i=0; i<j; i++)
	{
		printf("Number of residues in site %d are %d\n", i, mergedSites[i].nRes);
		printf("They are:\n");
		for (k=0; k<mergedSites[i].nRes; k++)
		{
			printf("%d ", mergedSites[i].residue[k]);
		}
		printf("\n");
	}*/
	
	return;
}


/*
int calculateNumberOfSprings(pdbFile *macromolecule, ligand *lig, bindingSite *allSites, int n)
{
	int i, j, k;
	int numSprings = 0;
	float line[3], r, line2[3], dot, r2, dis;
	
	for (i=0; i<allSites[n].nRes; i++)
	{
		//printf("residue %d r %f %f %f\n", allSites[n].residue[i], macromolecule->atom[allSites[n].residue[i]].x, macromolecule->atom[allSites[n].residue[i]].y, macromolecule->atom[allSites[n].residue[i]].z);
		for (j=(i+1); j<allSites[n].nRes; j++)
		{
			//printf("residue %d r %f %f %f\n", allSites[n].residue[j], macromolecule->atom[allSites[n].residue[j]].x, macromolecule->atom[allSites[n].residue[j]].y, macromolecule->atom[allSites[n].residue[j]].z);
			line[0] = macromolecule->atom[allSites[n].residue[j]].x - macromolecule->atom[allSites[n].residue[i]].x;
			line[1] = macromolecule->atom[allSites[n].residue[j]].y - macromolecule->atom[allSites[n].residue[i]].y;
			line[2] = macromolecule->atom[allSites[n].residue[j]].z - macromolecule->atom[allSites[n].residue[i]].z;
			r = sqrt(line[0]*line[0] + line[1]*line[1] + line[2]*line[2]);
			for (k=0; k<3; k++) {
				line[k] = line[k]/r;
			}
			//printf("line %f %f %f\n", line[0], line[1], line[2]);
			for (k=0; k<lig->nAtoms; k++)
			{
				//printf("k %d r %f %f %f\n",k, lig->atom[k].r[0], lig->atom[k].r[1], lig->atom[k].r[2]); 
				line2[0] = macromolecule->atom[allSites[n].residue[j]].x - lig->atom[k].r[0];
				line2[1] = macromolecule->atom[allSites[n].residue[j]].y - lig->atom[k].r[1];
				line2[2] = macromolecule->atom[allSites[n].residue[j]].z - lig->atom[k].r[2];
				r2 = sqrt(line2[0]*line2[0] + line2[1]*line2[1] + line2[2]*line2[2]);
				//printf("line2 %f  %f %f %f\n", line2[0], line2[1], line[2], r2);
				dot = line[0]*line2[0] + line[1]*line2[1] + line[2]*line2[2];
				//printf("j %d k %d r2 %f dot %f\n", j, k, r2, dot);
				dis = r2*r2 - dot*dot;

				//printf("dis = %f \n", dis);
				if (dis <= square_thresh_for_dist_to_springs)  // ****NOTE1 (DC): may need to change this value for all-heavy atom model?
								 // From paper: "Between each pair of Ca atoms i and j, whose connecting line passes within
								 // 3.5 A ̊ of any probe atom, we place a spring of length dij (dashed lines in the figure)"
								 // ****NOTE 2: In paper, this dist is reported as 3.5, but in Berezovsky's code, a value of 
								 // 3.0 is used  --  see Berezovsky's analysis.c
				{
					//dist_btwn_resids = sqrt(  (macromolecule->atom[allSites[n].residue[i]].x - macromolecule->atom[allSites[n].residue[j]].x)*(macromolecule->atom[allSites[n].residue[i]].x - macromolecule->atom[allSites[n].residue[j]].x)   +   (macromolecule->atom[allSites[n].residue[i]].y - macromolecule->atom[allSites[n].residue[j]].y)*(macromolecule->atom[allSites[n].residue[i]].y - macromolecule->atom[allSites[n].residue[j]].y)   +   (macromolecule->atom[allSites[n].residue[i]].z - macromolecule->atom[allSites[n].residue[j]].z)*(macromolecule->atom[allSites[n].residue[i]].z - macromolecule->atom[allSites[n].residue[j]].z)  );  // was not in orig
					//printf("x_c_1: %f y_c_1: %f z_c_1: %f     x_c_2: %f y_c_2: %f z_c_2: %f   _dst_: %f \n", macromolecule->atom[allSites[n].residue[i]].x, macromolecule->atom[allSites[n].residue[i]].y, macromolecule->atom[allSites[n].residue[i]].z, macromolecule->atom[allSites[n].residue[j]].x, macromolecule->atom[allSites[n].residue[j]].y, macromolecule->atom[allSites[n].residue[j]].z, dist_btwn_resids); // was not in orig
					numSprings++;
					//printf("Got in here, numSprings = %d \n", numSprings);
					k = lig->nAtoms;
				}
			}
		}
	}
	
	//printf("Number of springs = %d\n", numSprings);
	return numSprings;
}
*/


void setSpringsInBindingSite(pdbFile *macromolecule, ligand *lig, bindingSite *allSites, int n)
{
	int i, j, k;
	int numSprings = 0;
	float dist_btwn_resids;
	float r, dot, r2, dis;
	float line[3] =  { 0., 0., 0. };
	float line2[3] =  { 0., 0., 0. };
	
	for (i=0; i<allSites[n].nRes; i++)
	{
		for (j=(i+1); j<allSites[n].nRes; j++)
		{
			line[0] = macromolecule->atom[allSites[n].residue[j]].x - macromolecule->atom[allSites[n].residue[i]].x;
			line[1] = macromolecule->atom[allSites[n].residue[j]].y - macromolecule->atom[allSites[n].residue[i]].y;
			line[2] = macromolecule->atom[allSites[n].residue[j]].z - macromolecule->atom[allSites[n].residue[i]].z;
			r = sqrt(line[0]*line[0] + line[1]*line[1] + line[2]*line[2]);
			for (k=0; k<3; k++) {
				line[k] = line[k]/r;
			}
			for (k=0; k<lig->nAtoms; k++)
			{
				line2[0] = macromolecule->atom[allSites[n].residue[j]].x - lig->atom[k].r[0];
				line2[1] = macromolecule->atom[allSites[n].residue[j]].y - lig->atom[k].r[1];
				line2[2] = macromolecule->atom[allSites[n].residue[j]].z - lig->atom[k].r[2];
				r2 = line2[0]*line2[0] + line2[1]*line2[1] + line2[2]*line2[2]; // sqrt removed
				dot = line[0]*line2[0] + line[1]*line2[1] + line[2]*line2[2];
				dis = r2 - dot*dot;
				if (dis <= square_thresh_for_dist_to_springs)  // NOTE (DC): may need to change this value for all-heavy atom model?
								 // From paper: "Between each pair of Ca atoms i and j, whose connecting line passes within
								 // 3.5 A ̊ of any probe atom, we place a spring of length dij (dashed lines in the figure)"
				{
					allSites[n].spring[numSprings] = (edgePtr)malloc(sizeof(EdgeStr));
					allSites[n].spring[numSprings]->res1 = allSites[n].residue[i];  // *was* change to the following?:  allSites[n].spring[numSprings]->res1 = allSites[n].resid[i];
					//printf("HERE: allSites[n].residue[i] = %d    allSites[n].resid[i] = %d    \n", allSites[n].residue[j], allSites[n].resid[j]);  // NOTE (DC): "residue" just refers to an artificial index to keep track of the residues within this program (starts at 0)
					allSites[n].spring[numSprings]->res2 = allSites[n].residue[j];  // *was* change to the following?:  allSites[n].spring[numSprings]->res2 = allSites[n].resid[j];
					//printf("spring_ %d res1 %d res2 %d\n", n, allSites[n].spring[numSprings]->res1, allSites[n].spring[numSprings]->res2);
					dist_btwn_resids = sqrt(  (macromolecule->atom[allSites[n].residue[i]].x - macromolecule->atom[allSites[n].residue[j]].x)*(macromolecule->atom[allSites[n].residue[i]].x - macromolecule->atom[allSites[n].residue[j]].x)   +   (macromolecule->atom[allSites[n].residue[i]].y - macromolecule->atom[allSites[n].residue[j]].y)*(macromolecule->atom[allSites[n].residue[i]].y - macromolecule->atom[allSites[n].residue[j]].y)   +   (macromolecule->atom[allSites[n].residue[i]].z - macromolecule->atom[allSites[n].residue[j]].z)*(macromolecule->atom[allSites[n].residue[i]].z - macromolecule->atom[allSites[n].residue[j]].z)  );  // was not in orig
					//printf("res1: %d  x_c_1: %f y_c_1: %f z_c_1: %f      res2: %d x_c_2: %f y_c_2: %f z_c_2: %f   _dst_: %f \n", allSites[n].spring[numSprings]->res1, macromolecule->atom[allSites[n].residue[i]].x, macromolecule->atom[allSites[n].residue[i]].y, macromolecule->atom[allSites[n].residue[i]].z,          allSites[n].spring[numSprings]->res2, macromolecule->atom[allSites[n].residue[j]].x, macromolecule->atom[allSites[n].residue[j]].y, macromolecule->atom[allSites[n].residue[j]].z, dist_btwn_resids); // was not in orig
					//printf("numSprings_ %d\n", numSprings);
					numSprings++;
					break;
				}
			}
		}
	}
	allSites[n].numSprings = numSprings;
	//printf("Number of springs are %d\n", numSprings);
	return;
}

void mergeSprings(bindingSite *allSites, int i, int j)
{
	int k, l, n;
	int intSprings = 0;
	edgePtr *springs;
	
	//Finding number of common springs
	for (k=0; k<allSites[i].numSprings; k++)
	{
		for (l=0; l<allSites[j].numSprings; l++)
		{
			if ((allSites[i].spring[k]->res1 == allSites[j].spring[l]->res1) && (allSites[i].spring[k]->res2 == allSites[j].spring[l]->res2))
			{
				intSprings++;
				l = allSites[j].numSprings;
			}
		}
	}
	//printf("intSprings = %d\n", intSprings);
	int totalSprings = allSites[i].numSprings + allSites[j].numSprings - intSprings;
	printf("totalSprings = %d\n", totalSprings);
	
	if ((springs = (edgePtr *) malloc(totalSprings*sizeof(edgePtr))) == NULL)
	{
		printf("Error: Couldnt merge elements\n");
		exit(1);
	}
	
	//Merging springs on to springs structure
	for (k=0; k<allSites[i].numSprings; k++)
	{
		springs[k] = (edgePtr)malloc(sizeof(EdgeStr));
		springs[k]->res1 = allSites[i].spring[k]->res1;
		springs[k]->res2 = allSites[i].spring[k]->res2;
		//printf("k %d spring %d %d\n", k, springs[k]->res1, springs[k]->res2);
	}
	n = allSites[i].numSprings;
	for (k=0; k<allSites[j].numSprings; k++)
	{
		intSprings = 0;
		for (l=0; l<allSites[i].numSprings; l++)
		{
			if ((allSites[j].spring[k]->res1 == allSites[i].spring[l]->res1) && (allSites[j].spring[k]->res2 == allSites[i].spring[l]->res2))
			{
					intSprings = 1;
					l = allSites[i].numSprings;
			}
		}
		if (intSprings == 0)
		{
			springs[n] = (edgePtr)malloc(sizeof(EdgeStr));
			springs[n]->res1 = allSites[j].spring[k]->res1;
			springs[n]->res2 = allSites[j].spring[k]->res2;
			//printf("n %d spring %d %d\n", n, springs[n]->res1, springs[n]->res2);
			n++;
		}
	}
	
	
	//Freeing old springs structure
	for (k=0; k<allSites[i].numSprings; k++)
	{
		free(allSites[i].spring[k]);
	}
	if (allSites[i].numSprings > 0) free(allSites[i].spring);
	allSites[i].numSprings = totalSprings;
	//printf("Got here\n");
	for (k=0; k<allSites[j].numSprings; k++)
	{
		free(allSites[j].spring[k]);
	}
	if (allSites[j].numSprings > 0) free(allSites[j].spring);
	allSites[j].numSprings = 0;
	
	
	//Copying springs into merged structure
	if (totalSprings > 0 && (allSites[i].spring = (edgePtr *) malloc(totalSprings*sizeof(edgePtr))) == NULL)
	{
		printf("Error: no space to merge sites\n");
		exit(1);
	}
	
	for (k=0; k<totalSprings; k++)
	{
		allSites[i].spring[k] = (edgePtr)malloc(sizeof(EdgeStr));
		allSites[i].spring[k]->res1 = springs[k]->res1;
		allSites[i].spring[k]->res2 = springs[k]->res2;
		free(springs[k]);
	}
	free(springs);
	
	return;
}

void outputFinalSites(bindingSite *mergedSites, int numMergedSites, FILE *output)
{
	int i, j;
	
	printf("->printing Merged Sites to output file\n");
	
	for (i=0; i<numMergedSites; i++)
	{
		fprintf(output, "%6d  ", i);
		for (j=0; j<mergedSites[i].nRes; j++)
		{
			fprintf(output, "%4d ", mergedSites[i].residue[j]);
		}
		fprintf(output, "\n");
	}
	
	printf("<-printing Merged Sites to output file\n");

	return;
}
