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
  /* XXX - need to calculate mean edge distance to set default value for edgeLengthOffset */
  int edgeLengthOffset = 0;
  int source;
  int target;
  int residueCount;
  int *residues;
  //  int getShortestPath = 0;
  char outputFileName[1024];
  char betweennessFileName[1024];
  char shortestDistanceFileName[1024];
  char predecessorFileName[1024];
  Graph protein;
  Graph *proteinPtr;
  //  FILE *res1file = NULL;
  //  FILE *res2file = NULL;
  FILE *input;
  FILE *output;
  FILE *logfile;
  FILE *betweennessFile;
  FILE *shortestDistanceFile;
  FILE *predecessorFile;

  /* Turn off buffering for stdout so that it always prints to the screen. */
  //setbuf(stdout, NULL);

  if (argc != 6 && argc != 7) {
    printf("The command line should be of the format: \n");
    printf(">subopt <contactMap> <output> <lengthOffset> <source> <target> [<prevOutput>]\n");
    printf(" where:\n");
    printf("   contactMap - file with contact map.\n");
    printf("   output - prefix for output files (.out, .betweenness, .shortDist, .pred).\n");
    printf("   lengthOffset - the amount added to the optimal path length to get the maximum suboptimal path length.\n");
    printf("   source - the node index for the source residue.\n");
    printf("   target - the node index for the target residue.\n");
    printf("   prevOutput - the file with the shortest distance matrix from a previous subopt run.\n");
    exit(1);
  }

  /* Check that input file exists. */
  if ( (input=fopen(argv[1], "r")) == NULL ) {
    printf("Input file (contactMap) does not exist.\n");
    exit(1);
  }

  /* Get lengthOffset, source index, and target index. */
  edgeLengthOffset = atoi(argv[3]);
  source = atoi(argv[4]);
  target = atoi(argv[5]);

  output=fopen(strcat(strcpy(outputFileName,argv[2]),".out"), "w");
  logfile=fopen("output.log", "w");
  proteinPtr = &protein;
  get(proteinPtr, input);
  getEdges(proteinPtr);

  /* Initialize residues to contain all of the residues in tempGraph. */
  residueCount = proteinPtr->nres;
  residues = (int *) calloc(residueCount, sizeof(int));
  for (i=0; i<residueCount; i++) {
    residues[i] = i;
  }
  if (argc == 6) {
    floydWarshall(proteinPtr, residues, residueCount);

    /* Print shortest distance matrix. */
    strcat(strcpy(shortestDistanceFileName,argv[2]),".shortDist");
    printf("Print shortest distance matrix: %s.\n",shortestDistanceFileName);
    shortestDistanceFile = fopen(shortestDistanceFileName, "w");
    printShortestDistances(proteinPtr, shortestDistanceFile);
    fclose(shortestDistanceFile);
    printf("Done\n");

    /* Print predecessor list matrix. */
    strcat(strcpy(predecessorFileName,argv[2]),".pred");
    printf("Print predecessor list matrix: %s.\n",predecessorFileName);
    predecessorFile = fopen(predecessorFileName, "w");
    printPredecessors(proteinPtr, predecessorFile);
    fclose(predecessorFile);
    printf("Done\n");

    /* Print betweenness matrix. */
    strcat(strcpy(betweennessFileName,argv[2]),".betweenness");
    printf("Print betweenness matrix: %s.\n",betweennessFileName);
    betweennessFile = fopen(betweennessFileName, "w");
    edgeConnectivity(proteinPtr, residues, residueCount);
    printEdgeConnectivities(proteinPtr, betweennessFile);
    fclose(betweennessFile);
    printf("Done\n");
  } else {
    /* Read in shortest distance matrix from file. */
    strcat(strcpy(shortestDistanceFileName,argv[6]),".shortDist");
    printf("Read shortest distance matrix: %s.\n", shortestDistanceFileName);
    readShortestDistances(proteinPtr, shortestDistanceFileName);
    printf("Done\n");
    
    /* Read in shortest path predecessors from file. */
    strcat(strcpy(predecessorFileName,argv[6]),".pred");
    printf("Read predecessor matrix: %s.\n", predecessorFileName);
    readPredecessors(proteinPtr, predecessorFileName);    
    printf("Done\n");
  }

  getRes(&protein.Flow, source, target);

  printf("Calculating shortest distance between selections.\n");
  //shortestPath(protein, &protein.Flow);
  shortestPath(&protein);
  subOpt(protein, edgeLengthOffset, &protein.Flow, output, logfile);
  
  fclose(logfile);

  free(proteinPtr->Pred); // JULY_1

  return 0;
}
