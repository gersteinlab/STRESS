#ifndef _pdb_h
#define _pdb_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct 
{
	char name[4]; //name of atom
	char resName[4]; // residue name of atom
	int resID; //residue number of atom
	char chain[2]; //residue number of atom
	float x, y, z; //x,y, and z coordinates of atom
	float beta; //beta of atom
	float occupancy; //occupancy of atom
	char element[3]; //element of atom
} atomStr;


typedef struct 
{
   int nAtoms; //Number of atoms 
   atomStr *atom; //Information about each atom
} pdbFile;

//This function calculates the number of macromolecular atoms (not HETATM) in the pdb file
int findNumberOfAtomsPdbFile(FILE *pdbInput); 
//* param: pdbInput is the input pdb file for which the number of atoms are calculated.
//* returns value of number of atoms in the pdb file

//This function reads the pdb file and initializes the whole pdbInput data structure
void readPdbFile(pdbFile *macromolecule, FILE *pdbInput);
//* param: macromolecule is the data structure which contains the information in the pdb file
//* param: pdbInput is the input pdb file from which the information is read

#endif
