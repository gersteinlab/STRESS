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
    printf(">./graphNetwork <contactMap> <output>\n");
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

  output=fopen(argv[2], "w");
  logfile=fopen("output.log", "w");
  edgeConnectivityMatrixFile=fopen("betweenness.dat", "w");
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
  printEdgeConnectivities(protein, edgeConnectivityMatrixFile);
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
  graphProperties(protein, logfile);
  /* */

  fclose(logfile);

  free(proteinPtr->Pred);  // JULY_1

  return 0;
}
