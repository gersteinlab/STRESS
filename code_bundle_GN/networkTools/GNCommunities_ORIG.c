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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Graph.h"

int main (int argc, char *argv[])
{
  int i;
  int residueCount;
  int *residues;
  Graph protein;
  Graph *proteinPtr;
  FILE *input;
  FILE *output;
  FILE *logfile;
  FILE *edgeConnectivityMatrixFile;

  if (argc != 3) {
    /*printf("Error in number of arguments. \n");*/
    printf("The command line should be of the format: \n");
    printf(">./gncommunities <contactMap> <output>\n");
    printf(" where:\n");
    printf("   contactMap is file with contact map.\n");
    printf("   output is prefix for output file.\n");
    exit(1);
  }

  /* Check that input file exists. */
  printf("Input file: %s\n",argv[1]);
  if ( (input=fopen(argv[1], "r")) == NULL ) {
    printf("Input file (contactMap) does not exist.\n");
    exit(1);
  }

  output = fopen(argv[2], "w");
  logfile = fopen("output.log", "w");
  edgeConnectivityMatrixFile = fopen("betweenness.dat", "w");
  proteinPtr = &protein;
  get(proteinPtr, input);
  getEdges(proteinPtr);

  /* Initialize residues to contain all of the residues in tempGraph. */
  residueCount = proteinPtr->nres;
  residues = (int *) calloc(proteinPtr->nres, sizeof(int));
  for (i=0; i<residueCount; i++) {
    residues[i] = i;
  }
  floydWarshall(proteinPtr, residues, residueCount);
  
  edgeConnectivity(proteinPtr, residues, residueCount);
  printEdgeConnectivities(proteinPtr, edgeConnectivityMatrixFile);
  fclose(edgeConnectivityMatrixFile);

  protein.Community=NULL;

  /* REACTIVATE */
  gnewman(&protein, output, logfile);
  /* */
  
  fclose(output);

  /* For reading in previously generated community information.
     Need to turn gnewman off if this is on. */
  /*getComm(protein, input);*/
  
  /* REACTIVATE */  
  getFlowComm(protein, logfile);
  getIntercommunityFlow(protein, logfile);
  /* */

  /*nodeConnectivity(proteinPtr);*/
  /*characteristicPathLengthNodes(protein, logfile);*/ 

  /* REACTIVATE */
  //  graphProperties(protein, logfile);
  /* */

  fclose(logfile);

  return 0;
}
