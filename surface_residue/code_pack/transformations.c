#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <gsl/gsl_rng.h>

#include "pdbFile.h"
#include "surfaceProbe.h"
#include "transformations.h"

void calculatePBCbox(pdbFile *macromolecule,box *pbc)
{
	int i;
	float min[3], max[3];
	
	printf("->calculating size of box\n");
	max[0] = -10000.; max[1] = -10000.; max[2] = -10000.;
	min[0] = 10000.; min[1] = 10000.; min[2] = 10000.;
	//Calculating extremities of macromolecule in each direction
	for (i=0; i<macromolecule->nAtoms; i++)
	{
		if (max[0] < macromolecule->atom[i].x) {
			//printf(" max[0] %f macromolecule->atom[i].x %f i %d\n", max[0], macromolecule->atom[i].x, i);
			max[0] = macromolecule->atom[i].x;
		}
		if (max[1] < macromolecule->atom[i].y) {
			max[1] = macromolecule->atom[i].y;
		}
		if (max[2] < macromolecule->atom[i].z) {
			max[2] = macromolecule->atom[i].z;
		}
		if (min[0] > macromolecule->atom[i].x) {
			min[0] = macromolecule->atom[i].x;
		}
		if (min[1] > macromolecule->atom[i].y) {
			min[1] = macromolecule->atom[i].y;
		}
		if (min[2] > macromolecule->atom[i].z) {
			min[2] = macromolecule->atom[i].z;
		}
	}
	//printf("max %f %f %f min %f %f %f\n", max[0], max[1], max[2], min[0], min[1], min[2]);
	for (i=0; i<3; i++) {
		pbc->center[i] = (max[i] + min[i])/2.; //Center is the mean of the extremities of the macromolecule
		pbc->size[i] = (max[i] - min[i])*2.00; //Size is twice the size of the macromolecule.
	}
	printf("<-calculating size of box\n");
	
	return;
}

void transformCoordinatesToSameBox(ligand *lig, float *minR, float *maxR, box *pbc)
{
	int i, j;
	
	
	//Transforming coordinates back into primary box
	for (i=0; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			while (lig->atom[i].r[j] < minR[j])
			{
				lig->atom[i].r[j] += pbc->size[j];
			} 
	/*	}
		
		for (j=0; j<3; j++)
		{*/
			while (lig->atom[i].r[j] > maxR[j])
			{
				lig->atom[i].r[j] -= pbc->size[j];
			}
		}
		//printf("%f %f %f\n", lig->atom[i].r[0], lig->atom[i].r[1], lig->atom[i].r[2]);
	}
	
	return;
}


void calculateTransformationMatrixZaxis(float *r, float *T)
{
	float Txz[9], Tz[9];
	
	//Transforming r vector to z-axis (http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html)
	Txz[0] = r[0]/sqrt(r[0]*r[0] + r[1]*r[1]);
	Txz[1] = r[1]/sqrt(r[0]*r[0] + r[1]*r[1]);
	Txz[2] = 0.;
	Txz[3] = -r[1]/sqrt(r[0]*r[0] + r[1]*r[1]);
	Txz[4] = r[0]/sqrt(r[0]*r[0] + r[1]*r[1]);
	Txz[5] = 0.;
	Txz[6] = 0.;
	Txz[7] = 0.;
	Txz[8] = 1.;
	Tz[0] = r[2]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	Tz[1] = 0.;
	Tz[2] = -sqrt(r[0]*r[0] + r[1]*r[1])/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	Tz[3] = 0.;
	Tz[4] = 1.;
	Tz[5] = 0.;
	Tz[6] = sqrt(r[0]*r[0] + r[1]*r[1])/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	Tz[7] = 0.;
	Tz[8] = r[2]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

	T[0] = Tz[0]*Txz[0] + Tz[1]*Txz[3] + Tz[2]*Txz[6];
	T[1] = Tz[0]*Txz[1] + Tz[1]*Txz[4] + Tz[2]*Txz[7];
	T[2] = Tz[0]*Txz[2] + Tz[1]*Txz[5] + Tz[2]*Txz[8];
	T[3] = Tz[3]*Txz[0] + Tz[4]*Txz[3] + Tz[5]*Txz[6];
	T[4] = Tz[3]*Txz[1] + Tz[4]*Txz[4] + Tz[5]*Txz[7];
	T[5] = Tz[3]*Txz[2] + Tz[4]*Txz[5] + Tz[5]*Txz[8];
	T[6] = Tz[6]*Txz[0] + Tz[7]*Txz[3] + Tz[8]*Txz[6];
	T[7] = Tz[6]*Txz[1] + Tz[7]*Txz[4] + Tz[8]*Txz[7];
	T[8] = Tz[6]*Txz[2] + Tz[7]*Txz[5] + Tz[8]*Txz[8];

	return;
}




void makeMoleculeWhole(ligand *lig, ligand *ligNew, box *pbc)
{
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA at start of makeMoleculeWhole i %d lig->atom[i].r %f %f %f\n", t, lig->atom[t].r[0], lig->atom[t].r[1], lig->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/

	int i, j;
	float dr;

	//for (j=0; j<lig->nAtoms; j++)
	for (i=0; i<3; i++)
		ligNew->atom[0].r[i] = lig->atom[0].r[i];
	
	//printf("i %d ligNew->atom[i] %f %f %f\n", 0, lig->atom[0].r[0], lig->atom[0].r[1], lig->atom[0].r[2]);	
	for (i=1; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			dr = lig->atom[i].r[j] - ligNew->atom[i-1].r[j];
			if (dr > 4)  // *was*:  if (dr[j] > 3.8)
			{
				ligNew->atom[i].r[j] = lig->atom[i].r[j] - pbc->size[j];
			} else if (dr < -4)  // *was*:  } else if (dr[j] < -3.8) 
			{
				ligNew->atom[i].r[j] = lig->atom[i].r[j] + pbc->size[j];
				// printf("%f, %f, %f\n", ligNew->atom[i].r[j], lig->atom[i].r[j], dr);
			} else
			{
				ligNew->atom[i].r[j] = lig->atom[i].r[j];
			}
		}
		//printf("i %d ligNew->atom[i] %f %f %f\n", i, lig->atom[i].r[0], lig->atom[i].r[1], lig->atom[i].r[2]);
	}
	
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA at end of makeMoleculeWhole i %d ligNew->atom[i].r %f %f %f\n", t, ligNew->atom[t].r[0], ligNew->atom[t].r[1], ligNew->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/

	return;
}



void translateCoordinates(ligand *lig, ligand *ligNew, box *pbc, float stepSize, float *minR, float *maxR)
{
	float t[3];
	int i, j;


	// below code is original (local):
	/*
	r = stepSize*(float)rand()/RAND_MAX;
	phi = 2*3.14159*((float)rand()/RAND_MAX);
	theta = acos(((float)rand()/RAND_MAX) - 1.);
	
	t[0] = r * sin(theta) * cos(phi);
	t[1] = r * sin(theta) * sin(phi);
	t[2] = r * cos(theta);
	//printf("%f %f %f %f\n", t[0], t[1], t[2], r);
	
	for (i=0; i<lig->nAtoms; i++)
	{
		//printf("atom %d ", i);
		for (j=0; j<3; j++)
		{
			ligNew->atom[i].r[j] = lig->atom[i].r[j] + t[j];
			//printf("%f ", lig->atom[i].r[j]);
		}
		//printf("\n");
	}
	*/


	// Below algorithm is the same approach used by Berezovsky:

	const double update_trans_range = 3.;
	// the vector "t" represents the translation vector
    t[0] = update_trans_range*(((float)rand()/RAND_MAX)-0.5);
    t[1] = update_trans_range*(((float)rand()/RAND_MAX)-0.5);
    t[2] = update_trans_range*(((float)rand()/RAND_MAX)-0.5);
    // makeMoleculeWhole(lig, ligNew, pbc);
    for (i=0; i < lig->nAtoms; i++)
    {
    	for (j=0; j < 3; j++)
    	{
    		lig->atom[i].r[j] = ligNew->atom[i].r[j];
    	}
    }

	for (i=0; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			ligNew->atom[i].r[j] = lig->atom[i].r[j] + t[j];
			//printf("%f ", lig->atom[i].r[j]);
		}
		//printf("\n");
	}
	
	transformCoordinatesToSameBox(ligNew, minR, maxR, pbc);
	//makeMoleculeWhole(lig, ligNew, pbc);

	return;	
}





void rotateCoordinates(ligand *lig, ligand *ligNew, box *pbc, float stepSize, float *minR, float *maxR)
{
	/*
	printf("CA_CA -- in rotateCoordinates:   minR: %f %f %f     maxR: %f %f %f  \n", minR[0], minR[1], minR[2], maxR[0], maxR[1], maxR[2]);
	int t;
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA at start of rotateCoordinates i %d lig->atom[i].r %f %f %f\n", t, lig->atom[t].r[0], lig->atom[t].r[1], lig->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/


	float phi, psi, theta;
	float Tx[9], Ty[9], Tz[9], Txy[9], T[9];
	float Rcom[3];
	int i, j;
	
	theta = stepSize*((float)rand()/RAND_MAX);
	phi = stepSize*((float)rand()/RAND_MAX);
	psi = stepSize*((float)rand()/RAND_MAX);

	calculateTransformationMatrixRotateXaxis(Tx, theta);
	calculateTransformationMatrixRotateYaxis(Ty, phi);
	calculateTransformationMatrixRotateZaxis(Tz, psi);
	multiplyMatrices(Ty, Tx, Txy, 3, 3, 3);
	multiplyMatrices(Tz, Txy, T, 3, 3, 3);
	
	//for (i=0; i<3; i++)
			//printf("%f %f %f\n", T[3*i], T[3*i+1], T[3*i+2]);
	
	
	// makeMoleculeWhole(lig, ligNew, pbc);
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA after makeMoleculeWhole %d ligNew->atom[i].r %f %f %f\n", t, ligNew->atom[t].r[0], ligNew->atom[t].r[1], ligNew->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/



	/*
		// If molecule is NOT actually whole (ie, if unrealistic bond  
		// lengths), then there's a problem, so exit out of program
		int lig_res = 330;
		int lig_atom;
		int lig_atom_to_print = 0;
		double CA_CA_bond_dist_in_lig;
		double del_x_sqrd;
		double del_y_sqrd;
		double del_z_sqrd;

			for (lig_atom=0; lig_atom < ligNew->nAtoms; lig_atom++)
			{
				// Print a record of all CA-CA bond distances in the final ligand configuration:
				if (lig_atom>0)
				{
					del_x_sqrd = (ligNew->atom[lig_atom].r[0]-ligNew->atom[lig_atom-1].r[0]) * (ligNew->atom[lig_atom].r[0]-ligNew->atom[lig_atom-1].r[0]);
					del_y_sqrd = (ligNew->atom[lig_atom].r[1]-ligNew->atom[lig_atom-1].r[1]) * (ligNew->atom[lig_atom].r[1]-ligNew->atom[lig_atom-1].r[1]);
					del_z_sqrd = (ligNew->atom[lig_atom].r[2]-ligNew->atom[lig_atom-1].r[2]) * (ligNew->atom[lig_atom].r[2]-ligNew->atom[lig_atom-1].r[2]);
					CA_CA_bond_dist_in_lig = sqrt(del_x_sqrd + del_y_sqrd + del_z_sqrd);
					printf("CA_CA_bond_dist_in_lig: %f \n", CA_CA_bond_dist_in_lig);  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//if (CA_CA_bond_dist_in_lig < 3.6  ||  CA_CA_bond_dist_in_lig > 4.0)
					if (CA_CA_bond_dist_in_lig > 4.0)
					{
						printf("CA_CA extreme bond dist\n");
						exit("\nEXIT\n");
					}
				}
				lig_res += 1;
			}
		printf("CA_CA \n\n\n", t, lig->atom[t].r[0], lig->atom[t].r[1], lig->atom[t].r[2]);
	*/
	


	calculateCenterOfMass(ligNew, Rcom);
	//printf("ligNew\n");
	//for (i=0; i<lig->nAtoms; i++) {
		//for (j=0; j<3; j++) {
			//printf("%f ", ligNew->atom[i].r[j]);
		//}
		//printf("\n");
	//}
	//printf("Rcom - %f %f %f\n", Rcom[0], Rcom[1], Rcom[2]);
	for (i=0; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			ligNew->atom[i].r[j] = ligNew->atom[i].r[j] - Rcom[j];
		}
		//printf("i %d ligNew->atom[i].r %f %f %f\n", i, ligNew->atom[i].r[0], ligNew->atom[i].r[1], ligNew->atom[i].r[2]);
	}
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA after subtracting Rcoms %d ligNew->atom[i].r %f %f %f\n", t, ligNew->atom[t].r[0], ligNew->atom[t].r[1], ligNew->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/

	rotateCoordinatesTransformationMatrix(ligNew, T);
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA after rotateCoordinatesTransformationMatrix %d ligNew->atom[i].r %f %f %f\n", t, ligNew->atom[t].r[0], ligNew->atom[t].r[1], ligNew->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/

	
	for (i=0; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			ligNew->atom[i].r[j] = ligNew->atom[i].r[j] + Rcom[j];
		}
	}
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA after adding Rcoms %d ligNew->atom[i].r %f %f %f\n", t, ligNew->atom[t].r[0], ligNew->atom[t].r[1], ligNew->atom[t].r[2]);
	}
	printf("CA_CA \n");
	*/


	transformCoordinatesToSameBox(ligNew, minR, maxR, pbc);
	/*
	for (t=0; t<lig->nAtoms; t++)
	{
		printf("CA_CA after transformCoordinatesToSameBox %d ligNew->atom[i].r %f %f %f\n", t, ligNew->atom[t].r[0], ligNew->atom[t].r[1], ligNew->atom[t].r[2]);
	}
	*/

	return;
}

void calculateTransformationMatrixRotateYaxis(float *Ty, float theta)
{
	Ty[0] = cos(theta);
	Ty[1] = 0.;
	Ty[2] = sin(theta);
	Ty[3] = 0.;
	Ty[4] = 1.;
	Ty[5] = 0.;
	Ty[6] = -sin(theta);
	Ty[7] = 0.;
	Ty[8] = cos(theta);
	
	return;
}

void calculateTransformationMatrixRotateXaxis(float *Tx, float phi)
{
	Tx[0] = 1.;
	Tx[1] = 0.;
	Tx[2] = 0.;
	Tx[3] = 0.;
	Tx[4] = cos(phi);
	Tx[5] = -sin(phi);
	Tx[6] = 0.;
	Tx[7] = sin(phi);
	Tx[8] = cos(phi);
	
	return;
}

void calculateTransformationMatrixRotateZaxis(float *Tz, float psi)
{
	Tz[0] = cos(psi);
	Tz[1] = -sin(psi);
	Tz[2] = 0.;
	Tz[3] = sin(psi);
	Tz[4] = cos(psi);
	Tz[5] = 0.;
	Tz[6] = 0.;
	Tz[7] = 0.;
	Tz[8] = 1.;

	return;
}

void multiplyMatrices(float *A, float *B, float *C, int l, int m, int n)
{
	int i, j, k;
	float sum;
	
	//Calculating C = A*B; B - m*n matrix and A - l*m matrix
	for (i=0; i<l; i++)
	{
		for (j=0; j<n; j++)
		{
			sum = 0.;
			for (k=0; k<m; k++)
			{
				sum += A[i*m + k] * B[k*n + j];
			}
			C[i*n + j] = sum;			
		}
	}
}

void rotateCoordinatesTransformationMatrix(ligand *lig, float *T)
{
	int i, j, k;
	float sum[3];
	
	for (i=0; i<lig->nAtoms; i++)
	{
		for (j=0; j<3; j++)
		{
			sum[j] = 0.;
			for (k=0; k<3; k++)
			{
				sum[j] += T[j*3 + k]*lig->atom[i].r[k];
			}
		}
		for (j=0; j<3; j++) 
			lig->atom[i].r[j] = sum[j];
		//printf("i %d lig->atom[i] %f %f %f\n", i, lig->atom[i].r[0], lig->atom[i].r[1], lig->atom[i].r[2]);
	}
	
	return;
}
