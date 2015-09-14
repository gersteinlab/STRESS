#ifndef _surfaceProbe_h
#define _surfaceProbe_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "pdbFile.h"
#include "network.h"

//Contains information about the box
typedef struct 
{
	float center[3]; //The coordinates of the center of the box
	float size[3]; //The size of the box which is twice the dimensions of the macromolecule in all three axis of the coordinate frame
} box;

//Contains coordinates of an atom
typedef struct
{
	float r[3]; //x, y, and z coordinates of the atom
} coord; 


//Contains information about the ligand
typedef struct
{
	int nAtoms; //Number of atoms in the ligand
	coord *atom; //Array of coordinates of each atom in the ligand
} ligand;

//Contains binding site residues
typedef struct
{
	int nRes; //Number of residues in the binding site
	int *resid; //ResID of residues in the binding site
	int *residue; //Residue numbers of residues in the binding site
	edgePtr *spring;
	int numSprings;
	float LC, BL; //Local centrality and Binding Leverage of binding site
} bindingSite;


/*run MC simulation*/
// NOTE: the below line is the orig, and the line following it is the revised:
//void runSurfaceProbeSimulation(pdbFile *macromolecule, ligand *lig, box *pbc, int Ntrials, int NstepsPerTrial, float *LC, bindingSite *allSites);

void runSurfaceProbeSimulation(pdbFile *macromolecule, ligand *lig, box *pbc, int Ntrials, int NstepsPerTrial, float *LC, bindingSite *allSites, pdbFile *macromolecule_w_heavy);


//*param: macromolecule - contains macromolecule coordinate information
//*param: ligand - contains ligand information
//*param: pbc - periodic boundary conditions for box
//*param: Ntrials - number of MC simulations
//*param: NstepsPerTrial - number of steps per simulation
//*param: output - file for parsing output*/

//guess coordinates for first atom within box
void guessFirstAtomCoordinates(ligand *lig, float *minR, float *maxR);
//*param: lig - contains ligand information
//*param: minR - coordinates for one end of box
//*param: maxR - coordinates for other end of box*/

//guess coordinates for second atom of ligand
void guessSecondAtomCoordinates(ligand *lig);
//*param: lig - contains ligand information after 1st atom coordinates are known*/

//guess coordinates for another atom of ligand
void guessAnotherAtomCoordinates(ligand *lig, int i);
//*param: lig - contains ligand information after at least two atom coordinates are known
//*param: i - atom index for which the coordinates are being guessed*/

//Calculate the potential energy of the system in the C-alpha model (original paper)
void calculateEnergyCAmodel(pdbFile *macromolecule, ligand *lig, float *U, bool *V, box *pbc, ligand *ligTemp);
//*param: macromolecule - contains all the macromolecule information
//*param: lig - contains all the ligand information
//*returns: the potential energy of the system given ligand and macromolecule information*/

//perform one MC step
void performOneMonteCarloStep(pdbFile *macromolecule, ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR, \
					float *U, bool *V, ligand **ligTemp, ligand **ligTempOld);
//*param: macromolecule - contains all the macromolecule information
//*param: lig - contains all the ligand information
//*param: ligNew - contains trial MC step information
//*param: pbc - contains periodic boundary conditions
//*param: minR - the minimum coordinates of the box for each axis 
//*param: maxR - the maximum coordinates of the box for each axis 
//*param: U - potential enegy of system
//*param: i - step of MC simulation*/
// STL: removed i; add U as the energy array and V as the validation array

void calculateCenterOfMass(ligand *lig, float *Rcom);

void doMetroPolis(float *U, float *Unew, ligand *lig, ligand *ligNew, bool *V, bool *Vnew, ligand **ligTemp, ligand **ligTempOld);

void guessAngles(ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR);

void findBindingSite(pdbFile *macromolecule, ligand *lig, bindingSite *allSites, int n);

void sort(bindingSite *allSites, float *LC, int i);
/**
 * Sorts the number in residues according to each residue's LC score
 * @param residues is array to be sorted
 * @param LC is score for all residues
 * @param nres is number of elements in residues
 */
 
void mergeBindingSites(bindingSite *allSites, int Ntrials, float *LC);

float calculateSimilarity(bindingSite *allSites, int i, int j);

void mergeSites(bindingSite *allSites, float *LC, int i, int j);

void copyFinalSites(bindingSite *allSites, bindingSite *mergedSites, int Ntrials);

int getNumberFinalSites(bindingSite *allSites, int Ntrials);

int calculateNumberOfSprings(pdbFile *macromolecule, ligand *lig, bindingSite *allSites, int n);

void setSpringsInBindingSite(pdbFile *macromolecule, ligand *lig, bindingSite *allSites, int n);

void mergeSprings(bindingSite *allSites, int i, int j);

void outputFinalSites(bindingSite *mergedSites, int numMergedSites, FILE *output);



// STL
typedef struct Inode{
	int poz;
	struct Inode *next;
} inode;

typedef struct Affected{
	int start;
	int end;
} affected;

// *NEW*
affected changeOneBondAngle(ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR);
affected changeOneTorsionAngle(ligand *lig, ligand *ligNew, box *pbc, float *minR, float *maxR);
float x_AfterRotation(float x_orig, float y_orig, float z_orig, float origin_x, float origin_y, float origin_z, float cross_x, float cross_y, float cross_z, float delta_angle);
float y_AfterRotation(float x_orig, float y_orig, float z_orig, float origin_x, float origin_y, float origin_z, float cross_x, float cross_y, float cross_z, float delta_angle);
float z_AfterRotation(float x_orig, float y_orig, float z_orig, float origin_x, float origin_y, float origin_z, float cross_x, float cross_y, float cross_z, float delta_angle);


// float calculateEnergyCAmodelold(pdbFile *macromolecule, ligand *lig);

#endif
