#ifndef _transformations_h
#define _transformations_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pdbFile.h"
#include "network.h"

//calculating the size and center of the box for surface probe calculations
void calculatePBCbox(pdbFile *macromolecule,box *pbc);
//* param: macromolecule - which contains the coordinates of all the atoms in the macromolecule
//* pbc: contains information about the box

//calculate matrix for transformation of vector to z-axis*/
void calculateTransformationMatrixZaxis(float *r, float *T);
//*param: r - vector for which transformation matrix is calculated
//*param: T - transformation matrix that is returned to the program*/

//Transform all coordinates to same box according to the PBC
void transformCoordinatesToSameBox(ligand *lig, float *minR, float *maxR, box *pbc);
//*param: lig - contains ligand information and coordinates
//*param: minR - the minimum coordinates of the box for each axis 
//*param: maxR - the maximum coordinates of the box for each axis 
//*param: pbc - the periodic boundary condition information*/

//remove periodic boundary conditions
void makeMoleculeWhole(ligand *lig, ligand *ligNew, box *pbc);
//*param: lig - contains ligand information and coordinates
//*param: ligNew - returned ligand in which the PBC are removed
//*param: pbc - the periodic boundary condition information*/

//translate coordinates
void translateCoordinates(ligand *lig,  ligand *ligNew, box *pbc, float stepSize, float *minR, float *maxR);
//*param: lig -  contains all the ligand information
//*param: ligNew - contains trial MC step information
//*param: pbc - contains periodic boundary conditions
//*param: stepSize - maximum difference in distance of center of mass between lig and ligNew in trial run
//*param: minR - the minimum coordinates of the box for each axis 
//*param: maxR - the maximum coordinates of the box for each axis */

//rotate coordinates
void rotateCoordinates(ligand *lig, ligand *ligNew, box *pbc, float stepSize, float *minR, float *maxR);
//*param: lig -  contains all the ligand information
//*param: ligNew - contains trial MC step information
//*param: pbc - contains periodic boundary conditions
//*param: stepSize - maximum difference in angle along x,y, and z axes between lig and ligNew in trial run
//*param: minR - the minimum coordinates of the box for each axis 
//*param: maxR - the maximum coordinates of the box for each axis */

//calculate the transformation matrix for rotation around y-axes
void calculateTransformationMatrixRotateYaxis(float *Ty, float theta);
//*param: Ty - transformation matrix
//*param: theta - angle to rotate by*/

//calculate the transformation matrix for rotation around x-axes
void calculateTransformationMatrixRotateXaxis(float *Tx, float theta);
//*param: Tx - transformation matrix
//*param: theta - angle to rotate by*/

//calculate the transformation matrix for rotation around z-axes
void calculateTransformationMatrixRotateZaxis(float *Tz, float theta);
//*param: Tz - transformation matrix
//*param: theta - angle to rotate by*/

void multiplyMatrices(float *A, float *B, float *C, int l, int m, int n);

void rotateCoordinatesTransformationMatrix(ligand *lig, float *T);

#endif