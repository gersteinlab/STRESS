#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pdbFile.h"

int findNumberOfAtomsPdbFile(FILE *pdbInput)
{
	char line[10000];
	int posn;
	int nAtoms = 0;
	
	printf("-> finding number of atoms in pdb file\n");
	while (fgets(line, 10000, pdbInput) != NULL)
	{
		posn = strncmp(line, "ATOM", 4);
		if (posn == 0) {
			nAtoms++;
			//printf("line %s", line);
			//printf("posn %d nAtoms %d\n", posn, nAtoms);
		} 
	}
	
	printf("<- finding number of atoms in pdb file\n");
	
	return nAtoms;
}

void readPdbFile(pdbFile *macromolecule, FILE *pdbInput)
{
	char line[10000];
	int i = 0;
	int posn;
	char resID[5], x[8], y[8], z[8], O[5], beta[5];
	
	printf("-> reading macromolecule pdb\n");
	
	//Note: Reading information from lines that start with ATOM alone (macromolecule)
	while (fgets(line, 10000, pdbInput) != NULL)
	{
		posn = strncmp(line, "ATOM", 4);
		if (posn == 0) {
			strncpy(macromolecule->atom[i].name, line+13, 3);
			macromolecule->atom[i].name[3] = '\0';
			strncpy(macromolecule->atom[i].resName, line+17, 3);
			macromolecule->atom[i].resName[3] = '\0';  // was not commented out before
			strncpy(resID, line+22, 4);
			resID[4] = '\0';
			macromolecule->atom[i].resID = atoi(resID);
			resID[4] = '\0';
			macromolecule->atom[i].resID = atoi(resID);

			///  Newly-introduced (dc)
			strncpy(macromolecule->atom[i].chain, line+21, 1);
			macromolecule->atom[i].chain[1] = '\0';

			strncpy(x, line+30, 8);  // line was originally:   strncpy(x, line+31, 7);  ########################################################
			x[7] = '\0';
			macromolecule->atom[i].x = atof(x);
			strncpy(y, line+38, 8);  // line was originally:   strncpy(y, line+39, 7);  ########################################################
			y[7] = '\0';
			macromolecule->atom[i].y = atof(y);
			strncpy(z, line+46, 8);  //  line was originally:   strncpy(z, line+47, 7);  #######################################################
			z[7] = '\0';
			macromolecule->atom[i].z = atof(z);
			strncpy(O, line+56, 4);
			O[4] = '\0';
			macromolecule->atom[i].occupancy = atof(O);
			strncpy(beta, line+62, 4);
			beta[4] = '\0';
			macromolecule->atom[i].beta = atof(beta);
			strncpy(macromolecule->atom[i].element, line+77, 2);
			macromolecule->atom[i].element[2] = '\0';
			//printf("line %s%s %s %d %f %f %f %f %f %s   \n  chain: |%s| \n\n", line, macromolecule->atom[i].name, macromolecule->atom[i].resName, macromolecule->atom[i].resID, macromolecule->atom[i].x, macromolecule->atom[i].y, macromolecule->atom[i].z,  macromolecule->atom[i].occupancy,  macromolecule->atom[i].beta, macromolecule->atom[i].element, macromolecule->atom[i].chain);
			i++;
		}
		//printf("i %d\n", i);
	}
	
	printf("<- reading macromolecule pdb\n");
	
	return;
}

