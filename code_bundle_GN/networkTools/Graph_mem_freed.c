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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Graph.h"
#define BIG 1000000

edgePtr copyEdge(edgePtr edge)
{
  /*printf("->copyEdge\n");*/
  edgePtr newEdge = (edgePtr) malloc(sizeof(EdgeStr));

  newEdge->res1 = edge->res1;
  newEdge->res2 = edge->res2;
  newEdge->edgeConn = edge->edgeConn;

  /*printf("<-copyEdge\n");*/
  return newEdge;
}



Graph* copyGraph(Graph *graph)
{
  printf("->copyGraph\n");
  int i, j=0;
  Graph *newGraph = (Graph*) malloc(sizeof(Graph));

  newGraph->nres = graph->nres;
  newGraph->nedges = graph->nedges;
  
  if ((newGraph->dis = (int **) calloc(newGraph->nres, sizeof(int*))) == NULL)
    printf("No memory space allocatable for finding communities - newGraph->dis.\n");
  for (i=0; i<newGraph->nres; i++) {
    if ((newGraph->dis[i] = (int *) calloc(newGraph->nres, sizeof(int))) == NULL)
      printf("No memory space allocatable for reading Contact Map - newGraph->dis[%d].\n",i);
  }
  for (i=0; i<newGraph->nres; i++) {
    for (j = 0; j<newGraph->nres; j++) {
      newGraph->dis[i][j] = graph->dis[i][j];
    }
  }
	
  if ((newGraph->shortDis = (int **) calloc(newGraph->nres, sizeof(int*))) == NULL)
    printf("No memory space allocatable for making Shortest Distance Matrix.\n");
  for (i=0; i<newGraph->nres; i++) {
    if ((newGraph->shortDis[i] = (int *) calloc(newGraph->nres, sizeof(int))) == NULL)
      printf("No memory space allocatable for making Shortest Distance Matrix. \n");
  }
  for (i=0; i<newGraph->nres; i++) {
    for (j = 0; j<newGraph->nres; j++) {
      newGraph->shortDis[i][j] = graph->shortDis[i][j];
    }
  }

  if ((newGraph->Pred = (nodePtr **) calloc(newGraph->nres, sizeof(nodePtr*))) == NULL)
    printf("No memory space allocatable for making Shortest Path Matrix.\n");
  for (i=0; i<newGraph->nres; i++)
    if ((newGraph->Pred[i] = (nodePtr *) calloc(newGraph->nres, sizeof(nodePtr))) == NULL)
      printf("No memory space allocatable for making Shortest Path Matrix. \n");
  for (i=0; i<newGraph->nres; i++) {
    for (j=0; j<newGraph->nres; j++) {
      newGraph->Pred[i][j] = copyList(graph->Pred[i][j]);
    }
  }

  newGraph->nodeConn = NULL;

  if ((newGraph->edge = (edgePtr *) calloc(newGraph->nedges, sizeof(edgePtr))) == NULL)
    printf("No memory space allocatable for calculating node Connectivity.\n");
  for (i=0; i<newGraph->nedges; i++) {
    newGraph->edge[i] = copyEdge(graph->edge[i]);
  }

  newGraph->Community = NULL;
  newGraph->Flow = NULL;

  printf("<-copyGraph\n");
  return newGraph;
}


void get(Graph *prot, FILE *input)
{
  int i, j;
  char line[BIG];
  //float **dis;
  float dis = 0;
 
  prot->nres = 0;
  printf("Reading Distance matrix file. \n");
  while (fgets(line, BIG, input) != NULL) {
    prot->nres++;
  }
  printf("Number of residues is %d.\n", prot->nres);

  /*if ((dis = (float **) calloc(prot->nres, sizeof(float*))) == NULL)
    printf("No memory space allocatable for reading Contact Map.\n");
  //printf("hey1\n");
  for (i=0; i<prot->nres; i++) {
    if ((dis[i] = (float *) calloc(prot->nres, sizeof(float))) == NULL) {
      printf("No memory space allocatable for reading Contact Map. \n");
    }
    }*/

  /*#printf("hey2\n");*/
  if ((prot->dis = (int **) calloc(prot->nres, sizeof(int*))) == NULL) {
    printf("No memory space allocatable for reading Contact Map.\n");
  }

  for (i=0; i<prot->nres; i++) {
    if ((prot->dis[i] = (int *) calloc(prot->nres, sizeof(int))) == NULL) {
      printf("No memory space allocatable for reading Contact Map. \n");
    }
  }
  rewind(input);

  // Check for '0' in first place of first row
  // Used with new adjacencyMatrix.tcl which stores the atomselection string
  // in the first line of the matrix file.
  fgets(line, BIG, input);
  if (line[0] == '0') {
    rewind(input);
  } else {
    strcpy(prot->line,line);
    prot->nres--;
  }

  for (i=0; i<prot->nres; i++) {
    for (j=0; j<prot->nres; j++) {
      //fscanf(input, "%f ", &dis[i][j]);
      //prot->dis[i][j] = floor(dis[i][j]*100 + 0.5);
      fscanf(input, "%f ", &dis);
      prot->dis[i][j] = floor(dis*100 + 0.5);
    }
    fscanf(input,"\n");
  }
  fclose(input);

  for (i=0; i<prot->nres; i++) {
    for (j=(i+1); j<prot->nres; j++) {
      if (prot->dis[i][j] != prot->dis[j][i])
	printf("Warning: There is an error in the distance matrix:  i %d j %d dis[i][j] %d dis[j][i] %d\n", i, j, prot->dis[i][j], prot->dis[j][i]);
    }
  }
      
  if ((prot->shortDis = (int **) calloc(prot->nres, sizeof(int*))) == NULL)
    printf("No memory space allocatable for making Shortest Distance Matrix.\n");

  for (i=0; i<prot->nres; i++)
    if ((prot->shortDis[i] = (int *) calloc(prot->nres, sizeof(int))) == NULL)
      printf("No memory space allocatable for making Shortest Distance Matrix. \n");
  if ((prot->Pred = (nodePtr **) calloc(prot->nres, sizeof(nodePtr*))) == NULL)
    printf("No memory space allocatable for making Shortest Path Matrix.\n");
  for (i=0; i<prot->nres; i++)
    if ((prot->Pred[i] = (nodePtr *) calloc(prot->nres, sizeof(nodePtr))) == NULL)
      printf("No memory space allocatable for making Shortest Path Matrix. \n");
	
  prot->nedges = 0;
  for (i=0; i<prot->nres;i++)
    for (j=(i+1); j<prot->nres; j++)
      if (prot->dis[i][j] != 0)
	prot->nedges++;
	
  printf("The number of edges are %d\n", prot->nedges);
  if ((prot->edge = (edgePtr *) calloc(prot->nedges, sizeof(edgePtr))) == NULL)
    printf("No memory space allocatable for calculating node Connectivity.\n");

  printf("<-get\n");

  return;
}


void removeEdge(Graph *graph, int edgeIndex)
{
  int i=0;
  edgePtr *newEdges;

  /* Set distance to 0 for pair of nodes belonging to the deleted edge. */
  graph->dis[graph->edge[edgeIndex]->res1][graph->edge[edgeIndex]->res2] = 0;
  graph->dis[graph->edge[edgeIndex]->res2][graph->edge[edgeIndex]->res1] = 0;

  if ((newEdges = (edgePtr *) calloc(graph->nedges-1, sizeof(edgePtr))) == NULL)
    printf("No memory space allocatable for edges.\n");
  for (i=0; i<graph->nedges; i++) {
    if (i<edgeIndex) {
      newEdges[i] = graph->edge[i];
    } else if (i>edgeIndex) {
      newEdges[i-1] = graph->edge[i];
    } else {
      free(graph->edge[i]);
    }
  }

  free(graph->edge);
  graph->edge = newEdges;
  graph->nedges--;

  return;
}


void floydWarshall(Graph *prot, int *residues, int residueCount)
{
  printf("->floydWarshall\n");
  printf("  residueCount = %d\n",residueCount);
  int i, j, k, newDistance;

  /* Clean Pred linked lists. */
  for (i=0; i<residueCount; i++) {
    for (j=i; j<residueCount; j++) {
      if (prot->Pred[residues[i]][residues[j]] != NULL) {
	deleteList(&(prot->Pred[residues[i]][residues[j]]));
	//prot->Pred[residues[i]][residues[j]] = NULL;
	if (prot->Pred[residues[i]][residues[j]] != NULL) {
	  printf("*** NOT NULL 1 ***\n");
	}
      }
      /* Clean [j][i] as well as [i][j], but only if i!=j. */
      if (i != j) {
	if (prot->Pred[residues[j]][residues[i]] != NULL) {
	  deleteList(&(prot->Pred[residues[j]][residues[i]]));
	  //prot->Pred[residues[j]][residues[i]] = NULL;
	  if (prot->Pred[residues[i]][residues[j]] != NULL) {
	    printf("*** NOT NULL 2 ***\n");
	  }
	}
      }
    }
  }

  /* Initialization of shortest distance and shortest path matrices. */
  for (i=0; i<residueCount; i++) {
    prot->shortDis[residues[i]][residues[i]] = 0;
    for (j=i+1; j<residueCount; j++) {
      if (prot->dis[residues[i]][residues[j]] != 0) {
	prot->shortDis[residues[i]][residues[j]] = prot->dis[residues[i]][residues[j]];
	prot->shortDis[residues[j]][residues[i]] = prot->dis[residues[i]][residues[j]];
      } else {
	prot->shortDis[residues[i]][residues[j]] = BIG;  /* Arbitrary large number */
	prot->shortDis[residues[j]][residues[i]] = BIG;  /* Arbitrary large number */
      }
    }
  }
  printf("  Initialization done.\n");

  /* Dynamic programming algorithm. */
  printf("  Dynamic programming to calculate pairwise shortest distances.\n");
  for (k=0; k < residueCount; k++) {
    printf("%d ", k);
    for (i=0; i < residueCount; i++) {
      //printf("%d ", i);
      for (j=i; j < residueCount; j++) {
	/*printf("%d %d %d \n", k, j, i);*/
	newDistance = prot->shortDis[residues[i]][residues[k]] + prot->shortDis[residues[k]][residues[j]];
	if (prot->shortDis[residues[i]][residues[j]] > newDistance) {
	  prot->shortDis[residues[i]][residues[j]] = newDistance;
	  prot->shortDis[residues[j]][residues[i]] = newDistance;
	}
      }
    }
  }
  printf("\n");
  printf("  Distances calculated.\n");

  /* Build lists containing all residues for each path from i to j. */
  for (i=0; i<residueCount; i++) {
    printf("%d ", i);
    for (k=0; k<residueCount; k++) {
      for (j=0; j<residueCount; j++) {
	if ( (prot->dis[residues[k]][residues[j]] != 0) &&
	     (prot->shortDis[residues[i]][residues[j]] ==
	      (prot->shortDis[residues[i]][residues[k]] + prot->shortDis[residues[k]][residues[j]])) ) {
	  pushNode(&prot->Pred[residues[i]][residues[j]], residues[k]);
	  /*printf("i %d j %d k %d \n", i, j, k);*/
	}
      }
    }
  }
  printf("\n");

  printf("<-floydWarshall\n");
  return;
}


//void shortestPath(Graph protein, flowPtr *Flow)
void shortestPath(Graph *protein)
{
  int i, j;
  nodePtr currPtr1, currPtr2, PathPtr;
  flowPtr Flow = protein->Flow;
  Flow->shortestDis = BIG; /* Arbitrary large number */
  for (i=0; i <Flow->nres1; i++) {
    for (j=0; j<Flow->nres2; j++) {
      if (protein->shortDis[Flow->source[i]][Flow->target[j]] < Flow->shortestDis) {
	Flow->shortestDis = protein->shortDis[Flow->source[i]][Flow->target[j]];
      }
    }
  }
  
  if (Flow->shortestDis == BIG) {
	  printf("There are no paths between Source(s) and Target(s).\n");
	  exit(1);
  }

  printf("The shortest distance is %d\n", Flow->shortestDis);
  Flow->sources = NULL;
  Flow->targets = NULL;
  for (i=0; i<Flow->nres1; i++) {
    for (j=0; j<Flow->nres2; j++) {
      if (protein->shortDis[Flow->source[i]][Flow->target[j]] == Flow->shortestDis) {
	pushNode(&Flow->sources, (Flow->source[i]));
	pushNode(&Flow->targets, (Flow->target[j]));
      }
    }
  }

  currPtr1 = Flow->sources;
  currPtr2 = Flow->targets;
  while (currPtr1 != NULL) {
    PathPtr = NULL;
    allPaths(protein, currPtr1->residue, currPtr2->residue, &PathPtr);
    currPtr1 = currPtr1->next;
    currPtr2 = currPtr2->next;
  }

  return;
}


void subOpt(Graph protein, int edgeLengthOffset, flowPtr *subOptPaths, FILE *output3, FILE *logfile)
{
  printf("->subOpt\n");
  int i, j;
  int jPred;         // Predecessor node to j in (target->source) optimal path
  int jSucc;         // Successor node to j in (target->source) optimal path
  int currGraphNode;
  int pathCount = 0;
  int source = 0;
  int target;
  int max;
  int pathDist;
  int suboptPathDist;
  nodePtr currPtr;
  nodePtr prevPtr;
  nodePtr triedPathsPtr = NULL;  // Holds all tested paths that have been less than cutoff
  nodePtr path = NULL;           // Holds source side of (source->target) shortest path, but backwards (so, l->source)
  nodePtr tempPath = NULL;       // Holds temporary path information before path is verified as good
  nodePtr allowedPaths = NULL;
  nodePtr reversePath = NULL;
  
  (*subOptPaths)->shortestDis = BIG;
  for (i=0; i <(*subOptPaths)->nres1; i++) {
    for (j=0; j<(*subOptPaths)->nres2; j++) {
      if (protein.shortDis[(*subOptPaths)->source[i]][(*subOptPaths)->target[j]] < (*subOptPaths)->shortestDis) {
	(*subOptPaths)->shortestDis = protein.shortDis[(*subOptPaths)->source[i]][(*subOptPaths)->target[j]];
      }
    }
  }
  
  fprintf(output3, "The shortest distance is here %d\n", (*subOptPaths)->shortestDis);
  if ((*subOptPaths)->shortestDis == BIG) {
    printf("There are no paths between the Source(s) and Target(s).\n");
    //exit(1);
    return;
  }
  printf("The shortest distance is %d\n", (*subOptPaths)->shortestDis);
  printf("The edge length offset is %d\n", edgeLengthOffset);
  fprintf(output3,"The edge length offset is %d\n", edgeLengthOffset);
  printf("Suboptimal paths will be no longer than %d\n", (*subOptPaths)->shortestDis + edgeLengthOffset);
  fprintf(output3,"Suboptimal paths will be no longer than %d\n", (*subOptPaths)->shortestDis + edgeLengthOffset);
  
  /* Looking for all possible sources and targets with shortest distance smaller than (*subOptPaths)->shortestDis + edgeLengthOffset */
  (*subOptPaths)->sources = NULL;
  (*subOptPaths)->targets = NULL;
  for (i=0; i<(*subOptPaths)->nres1; i++) {
    for (j=0; j<(*subOptPaths)->nres2; j++) {
      if (protein.shortDis[(*subOptPaths)->source[i]][(*subOptPaths)->target[j]] <= (edgeLengthOffset + (*subOptPaths)->shortestDis)) {
	pushNode(&(*subOptPaths)->sources, (*subOptPaths)->source[i]);
	pushNode(&(*subOptPaths)->targets, (*subOptPaths)->target[j]);
      }
    }
  }

  /* Print out lists of the sources and targets */
  printList((*subOptPaths)->sources);
  printf("The sources are: \n");
  printf("%d ", (*subOptPaths)->sources->residue);
  fprintf(output3, "The sources are: ");
  fprintf(output3, "%d ", (*subOptPaths)->sources->residue);
  
  printf("\nThe targets are: ");
  printf("%d ", (*subOptPaths)->targets->residue);
  printf("\n");
  fprintf(output3, "\nThe targets are: ");
  fprintf(output3,"%d ", (*subOptPaths)->targets->residue);
  fprintf(output3, "\n");
  (*subOptPaths)->subOptPtr = NULL;

  fprintf(output3, "The final paths are:\n");

  /* Regardless of which is identified as target or source, use the node with smaller ID as source. */
  if ((*subOptPaths)->targets->residue > (*subOptPaths)->sources->residue) {
    target = (*subOptPaths)->targets->residue;
    source = (*subOptPaths)->sources->residue;
  } else {
    target = (*subOptPaths)->sources->residue;
    source = (*subOptPaths)->targets->residue;
  }

  //path=NULL;
  /*if (triedPathsPtr != NULL) {
    deleteList(&triedPathsPtr);
    }*/
  //allowedPaths = NULL;
  //triedPathsPtr=NULL;

  i = target;
  pushNode(&path, i);
  pushNode(&triedPathsPtr, i);
  pushNode(&allowedPaths, i);
  //fprintf(output3, "%d, ", i);
  
  /* Load shortest path nodes backwards (target->source) into path, triedPathsPtr, and allowedPaths. */
  while (protein.Pred[source][i]->residue != source) {
    //fprintf(output3, "%d, ", protein.Pred[source][i]->residue);
    pushNode(&path, protein.Pred[source][i]->residue);
    pushNode(&triedPathsPtr, protein.Pred[source][i]->residue);
    pushNode(&allowedPaths, protein.Pred[source][i]->residue);
    i = protein.Pred[source][i]->residue;
  }
  //fprintf(output3, "%d, \n", source);
  pushNode(&path, source);
  pushNode(&triedPathsPtr, source);
  pushNode(&allowedPaths, source);
  //fprintf(output3, "  distance: %d\n", protein.shortDis[(*subOptPaths)->sources->residue][(*subOptPaths)->targets->residue]);
  reversePath = reverse(path);
  fprintList(output3,reversePath);
  fprintf(output3, "(%d)\n",(*subOptPaths)->shortestDis);
  //fprintf(output3, "\n");
  deleteList(&reversePath);
  pathCount++;
  
  prevPtr = NULL;
  
  /* Looking for suboptimal paths by tracing the path backwards. */
  while (path != NULL) {

    /* Set currPtr to source. */
    currPtr = path;
    while (currPtr->next != NULL) {	
      prevPtr = currPtr;
      currPtr = currPtr->next;
    }
    
    /* Moving from source to target removing the last member of path. */
    /* Save i as node in path closest to source, j as next closest. */
    if (currPtr->residue != target) {
      jSucc = currPtr->residue;
      j = prevPtr->residue;
      delete(&path);
      
      /* Add distances to know the shortest distance of the residue until
	 pathDist = prot->shortDis[target][j] (as it could be a suboptimal path). */
      //printf("jSucc: %d, j: %d\n", jSucc,j);
      pathDist = 0;
      currPtr=path;
      prevPtr=path;
      while (currPtr!=NULL) {
	pathDist += protein.dis[prevPtr->residue][currPtr->residue];
	prevPtr = currPtr;
	currPtr = currPtr->next;
      }
      /* Might have to add something for currGraphNode<jSucc too (looking for suboptimal paths when the first optimal path is found). */
      
      /* Go through all graph nodes. */
      for (currGraphNode=0; currGraphNode<protein.nres; currGraphNode++) {
	if (j == target) {
	  jPred = target;
	} else {
	  jPred = protein.Pred[target][j]->residue;
	}

	suboptPathDist = pathDist + protein.dis[j][currGraphNode] + protein.shortDis[currGraphNode][source];
	/* j and currGraphNode have an edge between them */
	/* currGraphNode isn't the predecessor of j */
	/* currGraphNode isn't the successor of j */
	/* if the suboptimal path defined by path, j, and currGraphNode is shorter than the cutoff */
	if ( (protein.dis[j][currGraphNode] != 0) &&
	     (jPred != currGraphNode) &&
	     (currGraphNode != jSucc) &&
	     (suboptPathDist <= (edgeLengthOffset + (*subOptPaths)->shortestDis)) ) {
	  
	  /*fprintf(output3,"\n");
	    fprintf(output3,"  -- pathDist: %d; protein.dis[%d][%d]: %d; protein.shortDis[%d][%d]: %d; currGraphNode: %d\n",
	    pathDist, j, currGraphNode, protein.dis[j][currGraphNode], currGraphNode, source, protein.shortDis[currGraphNode][source], currGraphNode);
	    fprintf(output3,"  -- distance: (%d + %d + %d) = %d; currGraphNode: %d\n",
	    pathDist, protein.dis[j][currGraphNode], protein.shortDis[currGraphNode][source],
	    pathDist + protein.dis[j][currGraphNode] + protein.shortDis[currGraphNode][source], currGraphNode);*/
	  
	  /* Load shortest path (currGraphNode->source) into tempPath where currGraphNode is a node off the optimal path between target and source. */
	  deleteList(&tempPath);
	  tempPath = copyList(path);
	  i = currGraphNode;
	  pushNode(&tempPath, i);
	  if (i != source) {
	    while (protein.Pred[source][i]->residue != source) {
	      pushNode(&tempPath, protein.Pred[source][i]->residue);
	      i = protein.Pred[source][i]->residue;
	    }
	    pushNode(&tempPath, source);
	  }
	  printList(tempPath);
	  printf("\n");
	  
	  /* if path has no matches to previously tried paths and has no cycles,
	     append it onto allowedPaths, and print it to the output file */
	  if (match(triedPathsPtr,tempPath) == 0 &&
	      hasCycle(path) == 0) {
	    /* Add nodes to triedPathsPtr so that this path is not tried again. */
	    deleteList(&path);
	    path = copyList(tempPath);
	    append(&triedPathsPtr,path);
	    
	    append(&allowedPaths,path);
	    reversePath = reverse(path);
	    fprintList(output3,reverse(path));
	    deleteList(&reversePath);
	    fprintf(output3, "(%d)\n",suboptPathDist);
	    //fprintf(output3,"\n");
	    pathCount++;
	    
	    currGraphNode = protein.nres;
	  }
	}
      }
      prevPtr = NULL;
    } else {
      /* Reached if (currPtr->residue == target) */
      /* Deleting path removes the last node in the list (closest to source),
	 and if there is only one node left it knocks the program out of the
	 enclosing while loop. */
      delete(&path);
    }
  }
  
  currPtr = (*subOptPaths)->subOptPtr;
  if (currPtr == NULL) {
    (*subOptPaths)->subOptPtr = allowedPaths;
  } else {	
    while (currPtr->next != NULL) {
      currPtr = currPtr->next;
    }
    currPtr->next = allowedPaths;
  }
  triedPathsPtr = NULL;
	  
  fprintf(output3, "Number of paths is %d\n", pathCount);
  printf("Number of paths is %d\n", pathCount);
  if (pathCount != 0) {
    fprintf(logfile, "%d %d \n", source, pathCount);
    printf("%d %d \n", source, pathCount);
  }
  
  /* Path matching to check which nodes and edges occur the most often in suboptimal paths */
  if (((*subOptPaths)->nodeMatch = (int *) calloc(protein.nres, sizeof(int))) == NULL)
    printf("No memory space allocatable for finding communities.\n");
  for (i=0; i<protein.nres; i++)
    (*subOptPaths)->nodeMatch[i] = 0;
  currPtr = (*subOptPaths)->subOptPtr;
  while (currPtr != NULL) {
    (*subOptPaths)->nodeMatch[currPtr->residue] = (*subOptPaths)->nodeMatch[currPtr->residue] + 1;
    currPtr = currPtr->next;
  }
  max = 0;
  for (i=0; i<protein.nres; i++) 
    if (max < (*subOptPaths)->nodeMatch[i])
      max = (*subOptPaths)->nodeMatch[i];
  fprintf(output3, "The residues that occur the most in the suboptimal paths are: \n");
  printf("The residues that occur the most in the suboptimal paths are: \n");
  for (i=0; i<protein.nres; i++) {
    if ((*subOptPaths)->nodeMatch[i] >= (0.5 * max)) {
      fprintf(output3, "%d %d\n", i, ((*subOptPaths)->nodeMatch[i]));
      printf("%d %d\n", i, ((*subOptPaths)->nodeMatch[i]));
    }
  }

  if (((*subOptPaths)->edgeMatch = (int *) calloc(protein.nedges, sizeof(int))) == NULL)
    printf("No memory space allocatable for finding suboptimal paths.\n");

  prevPtr = (*subOptPaths)->subOptPtr;
  currPtr = (*subOptPaths)->subOptPtr->next;
  while (currPtr != NULL) {
    for (i=0; i<protein.nedges; i++)
      if (((protein.edge[i]->res1 ==prevPtr->residue) && (protein.edge[i]->res2 == currPtr->residue)) ||
	  ((protein.edge[i]->res2 ==prevPtr->residue) && (protein.edge[i]->res1 == currPtr->residue)))
	(*subOptPaths)->edgeMatch[i] = (*subOptPaths)->edgeMatch[i] + 1;
    prevPtr = currPtr;
    currPtr = currPtr->next;
  }

  max = 0;
  for (i=0; i<protein.nedges; i++) {
    if (max < (*subOptPaths)->edgeMatch[i]) {
      max = (*subOptPaths)->edgeMatch[i];
    }
  }
  fprintf(output3, "The edges that occur the most in the suboptimal paths are:\n");
  for (i=0; i<protein.nedges; i++) {
    if ((*subOptPaths)->edgeMatch[i] >= (0.75 * max)) {
      fprintf(output3, "%d %d %d\n", (protein.edge[i]->res1), (protein.edge[i]->res2), ((*subOptPaths)->edgeMatch[i]));
    }
  }

  printf("<-subOpt\n");

  return;
}


void allPaths(Graph *prot, int source, int target, nodePtr *PathPtr)
{
  int i, j, k;
  nodePtr currPtr, prevPtr, path, altPtr;
	
  i = target;
  path = NULL;
  pushNode(&path, target);
  pushNode(&(*PathPtr), target);
  while (prot->Pred[source][i]->residue != source) {
    pushNode(&path, prot->Pred[source][i]->residue);
    pushNode(&*PathPtr, prot->Pred[source][i]->residue);
    i = prot->Pred[source][i]->residue;
  }
  pushNode(&path, source);
  pushNode(&*PathPtr, source);
  currPtr = path;
  prevPtr = NULL;
	
  k = 1;
  /*Looking for alternative paths by tracing the path backwards and looking for alternative paths*/
  while (path != NULL) {
    k = k + 1;
    while (currPtr->next != NULL) {
      prevPtr = currPtr;
      currPtr = currPtr->next;
    }
		
    if (currPtr->residue !=target) {
      i = currPtr->residue;
      j = prevPtr->residue;
      delete(&path);
			
      currPtr = prot->Pred[source][j];
      while (currPtr->residue != i) {
	currPtr = currPtr->next;
      }
      if (currPtr->next!=NULL) {
	currPtr = currPtr->next;
	altPtr = path;
	while (altPtr != NULL) {
	  pushNode(&*PathPtr, altPtr->residue);
	  altPtr = altPtr->next;
	}
	i = currPtr->residue;
	pushNode(&*PathPtr, i);
	pushNode(&path, i);
	while (prot->Pred[source][i]->residue != source) {
	  pushNode(&path, prot->Pred[source][i]->residue);
	  pushNode(&*PathPtr, prot->Pred[source][i]->residue);
	  i = prot->Pred[source][i]->residue;
	}
	pushNode(&path, source);
	pushNode(&*PathPtr, source);
      }
      currPtr = path;
      prevPtr = NULL;
    } else {
      delete(&path);
      currPtr = path;
    }
    /*free(&prevPtr);*/
  }
	
  /*Printing all paths*/
  currPtr = *PathPtr;
  printf("The path is:\n");
  while (currPtr != NULL) {
    printf("%d ", (currPtr->residue));
    currPtr = currPtr->next;
    if ((currPtr != NULL) && (currPtr->residue == target))
      printf("\n");
    /*free(&currPtr);*/
  }
  printf("\n");
  return;
}


void nodeConnectivity(Graph *prot)
{
  int i,j, k, *done, numPred;
  float max, *conn;
  nodePtr predPtr;
	
  k=0;
  if ((prot->nodeConn = (float *) calloc(prot->nres, sizeof(float))) == NULL)
    printf("No memory space allocatable for calculating node Connectivity.\n");
  if ((conn = (float *) calloc(prot->nres, sizeof(float))) == NULL)
    printf("No memory space allocatable for calculating node Connectivity.\n");
  if ((done = (int *) calloc(prot->nres, sizeof(int))) == NULL)
    printf("No memory space allocatable for calculating node Connectivity.\n");
	
  for (i=0; i<prot->nres; i++)
    prot->nodeConn[i] = 0;
	
  for (i=0; i<prot->nres; i++) {
    max = 0;
    for (j=0; j<prot->nres; j++) {
      /*Parsing out which paths need to be calculated*/
      if ((prot->shortDis[i][j] > 0) && (prot->shortDis[i][j] < BIG)) {
	done[j] = 0;
	conn[j] = 1;
      } else {
	done[j] = 1;
	conn[j] = 0;
      }
			
      /*Finding longest distance and corresponding node from node i*/
      if ((prot->shortDis[i][j] > max) && (prot->shortDis[i][j] < BIG)) {
	max = prot->shortDis[i][j];
	k = j;
      }
    }
    //allDone = 0;
		
    while (max !=0) {
      //succScore = conn[k];
      done[k] = 1;
      predPtr = prot->Pred[i][k];
      numPred = 0;
      while (predPtr !=NULL) {
	predPtr = predPtr->next;
	numPred = numPred + 1;
      }
      predPtr = prot->Pred[i][k];
      while (predPtr != NULL) {
	conn[predPtr->residue] = conn[predPtr->residue] + conn[k]/numPred;
	predPtr = predPtr->next;
      }
      max = 0;
      for (j=0; j<prot->nres; j++) {
	if ((prot->shortDis[i][j] > max) && (done[j]!=1)) {
	  max = prot->shortDis[i][j];
	  k = j;
	}
      }
    }
		
    for (j=0; j<prot->nres; j++) {
      prot->nodeConn[j] = prot->nodeConn[j] + conn[j];
    }
  }
  for (i=0;i<prot->nres; i++)
    printf("Node %d Connecitivity = %f\n", i, prot->nodeConn[i]);

  return;
}


void getEdges(Graph *prot)
{
  printf("->getEdges\n");
  int i, j, k;

  /* Clear existing edges from prot. */
  for (i=0; i<prot->nedges; i++) {
    if(prot->edge[i] != NULL) {
      free(prot->edge[i]);
    }
  }
  prot->nedges = 0;

  k = 0;
  for (i=0; i<prot->nres;i++) {
    for (j=(i+1); j<prot->nres; j++) {
      //printf("prot->dis[%d][%d] = %d\n",i,j,prot->dis[i][j]);
      if (prot->dis[i][j] != 0) {
	prot->edge[k] = (edgePtr)malloc(sizeof(EdgeStr));
	prot->edge[k]->res1 = i;
	prot->edge[k]->res2 = j;
	printf("res1 %d res2 %d k %d\n", prot->edge[k]->res1, prot->edge[k]->res2, k);
	k++;
      }
    }
  }

  prot->nedges = k;
  printf("<-getEdges\n");
  return;
}


void edgeConnectivity(Graph *proteinGraph, int *residues, int residueCount)
{
  printf("->edgeConnectivity\n");

  int i, j, k, l, numPred, reverseIndex1, reverseIndex2;
  int *done, *reverseIndices;
  float max, tempConnectivity;
  float *nodeConnectivity;
  float **conn;
  nodePtr predPtr;
  k=0;

  if ((done = (int *) calloc(residueCount, sizeof(int))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");
  if ((reverseIndices = (int *) calloc(proteinGraph->nres, sizeof(int))) == NULL)
    printf("No memory space allocatable for calculating reverseIndices.\n");
  for (i=0; i<proteinGraph->nres; i++) {
    reverseIndices[i] = -1;
  }
  for (i=0; i<residueCount; i++) {
    reverseIndices[residues[i]] = i;
  }
  if ((nodeConnectivity = (float *) calloc(residueCount, sizeof(float))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");

  
  if ((conn = (float **) calloc(residueCount, sizeof(float*))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");
  for (i=0; i<residueCount; i++)
    if ((conn[i] = (float *) calloc(residueCount, sizeof(float))) == NULL)
      printf("No memory space allocatable for calculating edge Connectivity.\n");
  
  printf("  proteinGraph->nedges = %d\n",proteinGraph->nedges);

  /* Clear existing edge connectivity data from edges connecting specified residues. */
  for (i=0; i<proteinGraph->nedges; i++) {
    while (i < proteinGraph->nedges &&
	   (reverseIndices[proteinGraph->edge[i]->res1] == -1 ||
	    reverseIndices[proteinGraph->edge[i]->res2] == -1)) {
      i++;
    }
    if (i < proteinGraph->nedges) {
      proteinGraph->edge[i]->edgeConn = 0;
    }
  }

  /* */
  for (i=0; i<residueCount; i++) {
    printf("%d ",i);
    for (j=0; j<residueCount; j++) {
      for (l=j; l<residueCount; l++) {
	conn[j][l] = 0;
	conn[l][j] = 0;
      }
    }

    max = 0;
    /*for (j=0; j<residueCount; j++) {*/
    for (j=0; j<residueCount; j++) {
      /* Determine which paths need to be calculated. */
      if ((proteinGraph->shortDis[residues[i]][residues[j]] > 0) && (proteinGraph->shortDis[residues[i]][residues[j]] < BIG)) {
	done[j] = 0;
	nodeConnectivity[j] = 1;

	/* Find longest distance and corresponding node from node i. */
	if (proteinGraph->shortDis[residues[i]][residues[j]] > max) {
	  max = proteinGraph->shortDis[residues[i]][residues[j]];
	  k = j;
	}
      } else {
	done[j] = 1;
	nodeConnectivity[j] = 0;
      }
    }
    
    /*printf("  hey4: %d\n", i);*/
    /* */
    while (max != 0) {
      done[k] = 1;
      predPtr = proteinGraph->Pred[residues[i]][residues[k]];
      numPred = 0;
      while (predPtr !=NULL) {
	predPtr = predPtr->next;
	numPred++;
      }
      predPtr = proteinGraph->Pred[residues[i]][residues[k]];
      while (predPtr != NULL) {
	/* Get index of predPtr->residue in residues. */
	//m = 0;
	//flag = 0;
	reverseIndex1 = reverseIndices[predPtr->residue];
	/*reverseIndex1 = -1;
	while (reverseIndex1 == -1 && flag == 0) {
	  if (residues[m] == predPtr->residue) {
	    reverseIndex1 = m;
	  }
	  if (m >= residueCount) {
	    printf("  Error: didn't find reverseIndex1; m=%d, predPtr->residue=%d\n",m,predPtr->residue);
	    flag = 1;
	  }
	  m++;
	  }*/
	tempConnectivity = nodeConnectivity[k]/numPred;
	nodeConnectivity[reverseIndex1] += tempConnectivity;
	conn[reverseIndex1][k] = tempConnectivity;
	conn[k][reverseIndex1] = tempConnectivity;
	predPtr = predPtr->next;
      }
      max = 0;
      for (j=0; j<residueCount; j++) {
	if ((proteinGraph->shortDis[residues[i]][residues[j]] > max) && (done[j] != 1)) {
	  max = proteinGraph->shortDis[residues[i]][residues[j]];
	  k = j;
	}
      }
    }
    
    for (j=0; j<proteinGraph->nedges; j++) {
      reverseIndex1 = reverseIndices[proteinGraph->edge[j]->res1];
      reverseIndex2 = reverseIndices[proteinGraph->edge[j]->res2];
      if (reverseIndex1 != -1 && reverseIndex2 != -1) {
	proteinGraph->edge[j]->edgeConn += conn[reverseIndex1][reverseIndex2];
      }
    }
  }
  
  free(done);
  free(reverseIndices);
  free(nodeConnectivity);
  for (i=0; i<residueCount; i++)
    free(conn[i]);
  free(conn);
  
  printf("<-edgeConnectivity\n");
  return;
}


void addNodeToCommunityTree(Graph *proteinGraph, Graph *tempGraph, int numComm, int numCommPrev) {

  printf("->addNodeToCommunityTree\n");
  int i, j, addComm;
  float *a, *e;
  nodePtr newPtr;
  commPtr comm = proteinGraph->Community;
  addComm = 0;
  
  if ((e = (float *) calloc(tempGraph->nres, sizeof(float))) == NULL)
    printf("No memory space allocatable for calculating community modularity - e.\n");
  if ((a = (float *) calloc(tempGraph->nres, sizeof(float))) == NULL)
    printf("No memory space allocatable for calculating community modularity - a.\n");

  while (numCommPrev < numComm) {
    printf("numComm = %d, numCommPrev = %d\n",numComm,numCommPrev);
    for (i=0; i<tempGraph->nres; i++) {
      for (j = (i + 1); j<tempGraph->nres; j++) {
	if ((comm->leaf[i] == comm->leaf[j]) &&
	    (comm->community[i] != comm->community[j])) {
	  addComm = comm->community[j];
	  /* jth community will have a different predecessor now */
	  i = tempGraph->nres;
	  j = tempGraph->nres;
	}
      }
    }

    numCommPrev++;
    newPtr = malloc(sizeof(struct Node));
    for (i=0; i <tempGraph->nres; i++) {
      if (comm->community[i] == addComm) {
	newPtr->next = comm->leaf[i];
	newPtr->residue = tempGraph->nres + numCommPrev - 1;
	comm->leaf[i] = newPtr;
      }
    }

    /* Calculation of modularity. Uses full original protein graph. */
    for (i=0; i<numCommPrev; i++) {
      e[i] = 0;
      a[i] = 0;
    }

    /* All this business with subtracting (proteinGraph->nres) to get the right
       array index is horribleness based on the fact that the community indices
       crazily start at (proteinGraph->nres + 1).  This seriously needs to be
       fixed.  Also, funnily enough, the "residue" field in the leaf structure
       actually refers to the community index, not a residue index. */
    for (i=0; i<proteinGraph->nres; i++) {
      for (j=(i+1); j<proteinGraph->nres; j++) {
	/*printf("  i=%d, j=%d\n",i,j);*/
	if (proteinGraph->dis[i][j] != 0) {
	  /*a[(comm->leaf[i]->residue)-proteinGraph->nres] = a[(comm->leaf[i]->residue)-proteinGraph->nres] + 1;*/
	  /*a[(comm->leaf[j]->residue)-proteinGraph->nres] = a[(comm->leaf[j]->residue)-proteinGraph->nres] + 1;*/
	  /*printf(" a[%d], a[%d]\n",(comm->leaf[i]->residue)-proteinGraph->nres,(comm->leaf[j]->residue)-proteinGraph->nres);*/
	  a[(comm->leaf[i]->residue)-proteinGraph->nres]++;
	  a[(comm->leaf[j]->residue)-proteinGraph->nres]++;
	  /*printf(" i %d j %d community %d\n", i, j, ((comm->leaf[i]->residue)-proteinGraph->nres));*/
	  if (comm->leaf[j]->residue == comm->leaf[i]->residue) {
	    /*e[(comm->leaf[i]->residue)-proteinGraph->nres] = e[(comm->leaf[i]->residue)-proteinGraph->nres] + 1;*/
	    e[(comm->leaf[i]->residue)-proteinGraph->nres]++;
	  }
	}
      }
    }

    comm->Q[numCommPrev-1] = 0;
    for (i=0; i<numCommPrev; i++) {
      e[i] = e[i]/(proteinGraph->nedges);
      a[i] = a[i]/(2*proteinGraph->nedges);
      comm->Q[numCommPrev-1] += e[i] - (a[i] * a[i]);
      /*printf("i=%d, e[i]=%f a[i]=%f Q=%f\n", i, e[i], a[i], comm->Q[numCommPrev-1]);*/
    }
  }

  free(e);
  free(a);

  printf("<-addNodeToCommunityTree\n");
  return;
}


void gnewman(Graph *protein, FILE *output3, FILE *output4, char *community_tcl_fl)
{
  printf("->gnewman\n");
  Graph *tempGraph = copyGraph(protein); /* Copy of protein; protein stays in original form. */
  int i, j, i1;
  int edgeMax, residueCount, numComm, numCommPrev, addComm;
  int *residues;
  float max;
  float **edgeComm;
  nodePtr newPtr;

  /* Initialize residues to contain all of the residues in tempGraph. */
  residueCount = tempGraph->nres;
  //printf("residueCount = %d\n",residueCount);
  residues = (int *) calloc(tempGraph->nres, sizeof(int));
  for (i=0; i<residueCount; i++) {
    residues[i] = i;
  }

  /* Initialize community datastructures. */
  /* The community structure will be held in a tree structure with each node
     pointing towards their predecessor node.
     Initially assume one community with all leaves pointing towards root. */
  protein->Community = malloc(sizeof(commStr));
  protein->Community->leaf = (nodePtr *) calloc(tempGraph->nres, sizeof(nodePtr));
  if ((protein->Community->community = (int *) calloc(tempGraph->nres, sizeof(int))) == NULL)
    printf("No memory space allocatable for finding communities.\n");
  if ((protein->Community->Q = (float *) calloc(tempGraph->nres, sizeof(float))) == NULL)
    printf("No memory space allocatable for calculating community modularity.\n");
  for (i=0; i<tempGraph->nres; i++)
    protein->Community->Q[i] = 0;

  /* Set all leaf pointers to the initial community including all residues. */
  newPtr = malloc(sizeof(struct Node));
  for (i=0; i<tempGraph->nres; i++)
    protein->Community->leaf[i] = newPtr;
  newPtr->residue = tempGraph->nres;
  newPtr->next = NULL;
	  
  for (i=0; i<tempGraph->nres; i++)
    protein->Community->community[i] = 0;

  numComm = 1;
  numCommPrev = 1;
  addComm = 0;

  /*
    Local variables used within WHILE (removable '*'):
    tempGraph
    addComm
    protein
    numComm
    numCommPrev
    *max
    *edgeMax
    *residueCount
    *residues
    *newPtr
  */

  /*********
   * WHILE *
   *********/
  //printf("tempGraph->nedges = %d\n",tempGraph->nedges);
  while (tempGraph->nedges != 0) {
    /* Check for number of communities. */
    for (i=0; i<tempGraph->nres; i++) {
      addComm = 0;
      for (j=i+1; j<tempGraph->nres; j++) {
	if ( tempGraph->shortDis[i][j] == BIG &&
	     protein->Community->community[i] == protein->Community->community[j] ) {
	  protein->Community->community[j] = numComm;
	  addComm = 1;
	}
      }
      if (addComm == 1)
	numComm++;
    }
    printf("Number of communities %d, communities in prev iteration %d and edges %d\n", numComm, numCommPrev, tempGraph->nedges);

    /* if number of communities has increased after edge was removed, then we
     * need to add nodes to the tree.  The community number is stored in the
     * residue field of the node (node->residue = nres + communitynumber - 1).
     */
    if (numComm != numCommPrev) {
      addNodeToCommunityTree(protein, tempGraph, numComm, numCommPrev);
      numCommPrev = numComm;
    }

    edgeConnectivity(tempGraph, residues, residueCount);

    /* Identify the edge with the maximum edge connectivity. */
    max = 0;
    edgeMax = 0;
    for (i=0; i<tempGraph->nedges; i++) {
      if (max < tempGraph->edge[i]->edgeConn) {
	max = tempGraph->edge[i]->edgeConn;
	edgeMax = i;
      }
    }
    printf("Maximum edge Connectivity %f between residues %d and %d\n", max, tempGraph->edge[edgeMax]->res1+1, tempGraph->edge[edgeMax]->res2+1);
    
    /* Find distance only for those nodes in same community as deleted edge. */
    residueCount = 0;
    /*printf("Community %d", community[tempGraph->edge[edgeMax]->res1]);*/
    for (i1=0; i1<tempGraph->nres; i1++) {
      if (protein->Community->community[i1] == protein->Community->community[tempGraph->edge[edgeMax]->res1]) {
	residues[residueCount] = i1;
	residueCount++;
      }
    }

    /* Remove edge with maximum edge connectivity. */
    removeEdge(tempGraph, edgeMax);

    floydWarshall(tempGraph, residues, residueCount);
    printf("Shortest distance and shortest path matrices calculated.\n");
    
    /* Special case where last edge is removed from the network. */
    if (tempGraph->nedges == 1) {
      numCommPrev = numCommPrev + 1;
      newPtr = malloc(sizeof(struct Node));
      newPtr->next = protein->Community->leaf[tempGraph->edge[edgeMax]->res2];
      newPtr->residue = tempGraph->nres + numCommPrev - 1;
      protein->Community->leaf[tempGraph->edge[edgeMax]->res2] = newPtr;
    }
  }
  /*************
   * END WHILE *
   *************/


  /* Find community with optimum Q value. */
  protein->Community->Q[0] = 0;
  protein->Community->Q[tempGraph->nres-1] = -1;
  protein->Community->optQ = 0;
  numComm = 1;
  for (i=0; i<tempGraph->nres; i++) {
    /*fprintf(output3, "Community Q[%d]: %f\n", i, protein->Community->Q[i]);*/
    if (protein->Community->optQ < protein->Community->Q[i]) { 
      protein->Community->optQ = protein->Community->Q[i];
      numComm = i;
    }
  }

  /* numComm needs to have the number of communities, not the associated index
     of the 0-indexed array. */
  numComm++;
  fprintf(output3, "The optimum number of communities is %d and modularity value is %f\n", numComm, protein->Community->optQ);
  fprintf(output3, "The optimum number of communities is %d and modularity value is %f\n", numComm, protein->Community->optQ);
  protein->Community->communityNumber = numComm;

  /* Print out communities with the residues they contain. */
  for (i=0; i<numComm; i++) {
    fprintf(output3, "The residues in community %d are: ", (i+1));
    printf("The residues in community %d are: ", (i+1));
    for (j=0; j<tempGraph->nres; j++) {
      newPtr = protein->Community->leaf[j];
      /*printf("j=%d ", j);
	printf("newPtr->residue = %d\n",newPtr->residue);*/
      while ((newPtr->residue) >= (tempGraph->nres + numComm)) {
	newPtr = newPtr->next;
	printf("newPtr->residue = %d\n",newPtr->residue);
      }
      if (newPtr->residue == (tempGraph->nres + i)) {
	/*printf("^^^\n");*/
	fprintf(output3, "%d ", j);
	protein->Community->community[j] = i;
      }
    }
    fprintf(output3, "\n");
  }
	
  /* Find edges with the most edge Connectivity for communication within each
     community and for communication between communities */
  if ((edgeComm = (float **) calloc(numComm, sizeof(float*))) == NULL)
    printf("No memory space allocatable for finding communities.\n");
  for (i=0; i<numComm; i++)
    if ((edgeComm[i] = (float *) calloc(numComm, sizeof(float))) == NULL)
      printf("No memory space allocatable for finding communities. \n");
  for (i=0; i<numComm; i++)
    for (j=0;j<numComm; j++)
      edgeComm[i][j] = 0;

  for (i1=0; i1<protein->nedges; i1++) {
    if (edgeComm[protein->Community->community[protein->edge[i1]->res1]][protein->Community->community[protein->edge[i1]->res2]] < protein->edge[i1]->edgeConn) {
      edgeComm[protein->Community->community[protein->edge[i1]->res1]][protein->Community->community[protein->edge[i1]->res2]] = protein->edge[i1]->edgeConn;
      edgeComm[protein->Community->community[protein->edge[i1]->res2]][protein->Community->community[protein->edge[i1]->res1]] = protein->edge[i1]->edgeConn;
    }
  }
	
  /* Setting max to be the normalization for radius of cylinder connecting main information centers for edges */
  max = 0;
  for (i=0; i<numComm; i++)
    for (j=(i+1); j<numComm; j++) 
      if (edgeComm[i][j] > max) 
	max = edgeComm[i][j];
	
  /* Print Community.tcl for displaying communities in VMD. */
  FILE *output2 = fopen(community_tcl_fl, "w");
  for (i=0; i<numComm; i++) {
    fprintf(output2, "mol representation NewCartoon 0.300000 6.000000 4.100000 0\n");
    fprintf(output2, "mol color ColorID %d\n", i);
    fprintf(output2, "mol selection {residue ");
    for (j=0; j<tempGraph->nres; j++) {
      if (protein->Community->community[j] == i) {
	fprintf(output2, "%d ", j);
      }
    }
    fprintf(output2, "}\n");
    fprintf(output2, "mol material Opaque\n");
    fprintf(output2, "mol addrep top\n");
  }
	
  fprintf(output2, "graphics top color black\n");
  for (i=0; i<numComm; i++) {
    for (j=(i+1); j<numComm; j++) {
      if (edgeComm[i][j] != 0) {
	fprintf(output3, "The highest score in edge connectivities between communities %d and %d with highest edge connectivity is %f\n", (i+1), (j+1), edgeComm[i][j]);
	for (i1=0; i1<protein->nedges; i1++) {
	  if ((((protein->Community->community[protein->edge[i1]->res1] == i) &&
		(protein->Community->community[protein->edge[i1]->res2] == j)) ||
	       ((protein->Community->community[protein->edge[i1]->res1] == j) &&
		(protein->Community->community[protein->edge[i1]->res2] == i))) &&
	      (protein->edge[i1]->edgeConn >= 0.75 * edgeComm[protein->Community->community[protein->edge[i1]->res1]][protein->Community->community[protein->edge[i1]->res2]])) {
	    fprintf(output3, "%d %d %f \n", (protein->edge[i1]->res1), (protein->edge[i1]->res2), protein->edge[i1]->edgeConn);
	    fprintf(output2, "mol representation VDW 1.000000 8.000000\n");
	    fprintf(output2, "mol color ColorID %d\n", protein->Community->community[protein->edge[i1]->res1]);
	    fprintf(output2, "mol selection {residue %d and name CA P}\n", protein->edge[i1]->res1);
	    fprintf(output2, "mol material Opaque\n");
	    fprintf(output2, "mol addrep top\n");
	    fprintf(output2, "mol representation VDW 1.000000 8.000000\n");
	    fprintf(output2, "mol color ColorID %d\n", protein->Community->community[protein->edge[i1]->res2]);
	    fprintf(output2, "mol selection {residue %d and name CA P}\n", protein->edge[i1]->res2);
	    fprintf(output2, "mol material Opaque\n");
	    fprintf(output2, "mol addrep top\n");
	    fprintf(output2, "set sel1 [atomselect top \"residue %d and name CA P\"]\n", protein->edge[i1]->res1);
	    fprintf(output2, "set sel2 [atomselect top \"residue %d and name CA P\"]\n", protein->edge[i1]->res2);
	    fprintf(output2, "set coord1 [lindex [$sel1 get {x y z}] 0]\n");
	    fprintf(output2, "set coord2 [lindex [$sel2 get {x y z}] 0]\n");
	    fprintf(output2, "graphics top cylinder $coord1 $coord2 radius %f\n", edgeComm[i][j]/max);
	    fprintf(output2, "$sel1 delete\n");
	    fprintf(output2, "$sel2 delete\n");
	  }
	}
      }
    }
  }
  fclose(output2);

  /* Print out community tree. */
  fprintf(output3, "Community Tree\n");
  printf("Community Tree\n");
  treePtr* treeArray = malloc((numComm+1)*sizeof(treePtr));
  for (i=0; i<=numComm; i++) {
    treeArray[i] = NULL;
  }
  treeArray[0] = newTree(NULL,1);
  treeArray[1] = treeArray[0];
  for (i=0; i<numComm; i++) {
    int currentCommunity = 1;
    int parentCommunity = 1;
    fprintf(output3, "Community list for community %d is: ", (i+1));
    printf("Community list for community %d is: ", (i+1));
    for (j=0; j<tempGraph->nres; j++) {
      newPtr = protein->Community->leaf[j];
      //printf("j=%d ", j);
      //printf("newPtr->residue = %d\n",newPtr->residue);
      while ((newPtr->residue) >= (tempGraph->nres + numComm)) {
	newPtr = newPtr->next;
	//printf("newPtr->residue = %d\n",newPtr->residue);
      }
      if (newPtr->residue == (tempGraph->nres + i)) {
	fprintf(output3, "%d: ", j);
	printf("%d: ", j);
	while (newPtr->residue > tempGraph->nres + 1) {
	  //printf("^^^\n");
	  if (currentCommunity == 1) {
	    currentCommunity = newPtr->residue - tempGraph->nres + 1;
	  } else if (parentCommunity == 1) {
	    parentCommunity = newPtr->residue - tempGraph->nres + 1;
	  }
	  fprintf(output3, "%d ", newPtr->residue - tempGraph->nres + 1);
	  printf("%d ", newPtr->residue - tempGraph->nres + 1);
	  newPtr = newPtr->next;
	}
	if (currentCommunity == 1) {
	  currentCommunity = newPtr->residue - tempGraph->nres + 1;
	} else if (parentCommunity == 1) {
	  parentCommunity = newPtr->residue - tempGraph->nres + 1;
	}
	fprintf(output3, "%d ", newPtr->residue - tempGraph->nres + 1);
	printf("%d ", newPtr->residue - tempGraph->nres + 1);
	/*fprintf(output3, "currentCommunity: %d, parentCommunity: %d\n",currentCommunity,parentCommunity);
	printf("currentCommunity: %d, parentCommunity: %d\n",currentCommunity,parentCommunity);
	fprintTree(output3, treeArray[0]);
	printTree(treeArray[0]);
	fprintf(output3, "\n");
	printf("\n");*/
	if (currentCommunity > 1) {
	  treePtr child1 = addChild(treeArray[parentCommunity],1,parentCommunity);
	  treeArray[currentCommunity] = addChild(treeArray[parentCommunity],2,currentCommunity);
	  treeArray[parentCommunity] = child1;
	}
	break;
      }
    }
    fprintf(output3, "\n");
    printf("\n");
  }
  //fprintf(output3, "Printing Community Tree\n");
  //printf("Printing Community Tree\n");
  fprintTree(output3, treeArray[0]);
  printTree(treeArray[0]);
  deleteTree(treeArray[0]);
  free(treeArray);
  
  fprintf(output4, "Community Membership\n");
  for (i=0; i<protein->nres; i++) {
    fprintf(output4, "%d %d\n", i, protein->Community->community[i]);
  }
  
  fprintf(output4, "Modularity Score\n");
  for (i=0; i<protein->nres; i++) {
    fprintf(output4, "%d %.4f\n", i, protein->Community->Q[i]);
  }

  free(residues);
  for (i=0; i<numComm; i++)
    free(edgeComm[i]);
  free(edgeComm);

  printf("<-gnewman\n");
  return;
}


void getComm(Graph protein, FILE *input)
{
  printf("->getComm\n");
  FILE *ip;
  int i;

  protein.Community =malloc(sizeof(commStr));
  if ((protein.Community->community = (int *) calloc(protein.nres, sizeof(int))) == NULL)
    printf("No memory space allocatable for finding communities.\n");
  ip = fopen("Community45.log", "r");
  for (i=0; i<protein.nres; i++) {
    fscanf(input, "%d\n", &(protein.Community->community[i]));
  }
  protein.Community->communityNumber = 0;
  for (i=0; i<protein.nres; i++) {
    if (protein.Community->community[i] > protein.Community->communityNumber)
      protein.Community->communityNumber = protein.Community->community[i];
  }
  protein.Community->communityNumber = protein.Community->communityNumber + 1;
  printf("Number of communities is %d\n", protein.Community->communityNumber);
  fclose(ip);

  printf("->getComm\n");
  return;
}


void getFlowComm(Graph protein, FILE *output4)
{
  printf("->getFlowComm\n");
  int i, j, k, l, numPairs, max, m, n, numPred;
  int *done;
  int **CommLink;
  nodePtr predPtr;
  float max1, TotalFlow, DirectFlow, percentFlow;
  float *edgeConn, *nodeConn;
  float **conn, **CommChange, **TotalCommChange;
  
  max = 0;
  m = 1;
  
  printf("number of communities = %d\n", protein.Community->communityNumber);
  if ((edgeConn = (float *) calloc(protein.nedges, sizeof(float))) == NULL)
    printf("No memory space allocatable for finding flow between communities.\n");
  if ((nodeConn = (float *) calloc(protein.nres, sizeof(float))) == NULL)
    printf("No memory space allocatable for finding flow between communities.\n");
  if ((done = (int *) calloc(protein.nres, sizeof(int))) == NULL)
    printf("No memory space allocatable for finding flow between communities.\n");

  if ((conn = (float **) calloc(protein.nres, sizeof(float*))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");
  for (i=0; i<protein.nres; i++)
    if ((conn[i] = (float *) calloc(protein.nres, sizeof(float))) == NULL)
      printf("No memory space allocatable for calculating flow between communities.\n");

  if ((CommChange = (float **) calloc(protein.Community->communityNumber, sizeof(float*))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");
  for (i=0; i<protein.Community->communityNumber; i++)
    if ((CommChange[i] = (float *) calloc(protein.Community->communityNumber, sizeof(float))) == NULL)
      printf("No memory space allocatable for calculating flow between communities.\n");

  if ((TotalCommChange = (float **) calloc(protein.Community->communityNumber, sizeof(float*))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");
  for (i=0; i<protein.Community->communityNumber; i++)
    if ((TotalCommChange[i] = (float *) calloc(protein.Community->communityNumber, sizeof(float))) == NULL)
      printf("No memory space allocatable for calculating flow between communities.\n");

  if ((CommLink = (int **) calloc(protein.Community->communityNumber, sizeof(int*))) == NULL)
    printf("No memory space allocatable for calculating edge Connectivity.\n");
  for (i=0; i<protein.Community->communityNumber; i++)
    if ((CommLink[i] = (int *) calloc(protein.Community->communityNumber, sizeof(int))) == NULL)
      printf("No memory space allocatable for calculating flow between communities.\n");

  for (i=0; i<protein.Community->communityNumber; i++)
    for (j=0; j<protein.Community->communityNumber; j++) {
      CommLink[i][j] = 0;
      TotalCommChange[i][j] = 0;
      CommChange[i][j] = 0;
    }

  /* Loop over all pairs of communities. */
  for (i=0; i<protein.Community->communityNumber; i++) {
    for (j=(i+1); j<protein.Community->communityNumber; j++) {
    //for (j=protein.Community->communityNumber-1; j>i; j--) {
      //fprintf(output4,"i=%d, j=%d\n",i,j);
      numPairs = 0;
      for (k=0; k<protein.nedges; k++)
	edgeConn[k] = 0;

      /* For each pair of communities, find out how many of the paths go through
	 other communities. */
      for (k=0; k<protein.Community->communityNumber; k++)
	for (l=0; l<protein.Community->communityNumber; l++)
	  CommChange[k][l] = 0;
      
      for (k=0; k<protein.nres; k++) {

	for (l=0; l<protein.nres; l++)
	  for (n=0; n<protein.nres; n++)
	    conn[l][n] = 0;

	/* if we're dealing with a residue from community i,
	   loop through all residues and identify edges between this residue
	   and residues in community j. */
	if (protein.Community->community[k] == i) {
	  //fprintf(output4," k=%d, protein.Community->community[%d] = %d\n",k,k,i);
	  max = 0;
	  for (l=0; l <protein.nres; l++) {
	    if (protein.Community->community[l] == j) {
	      numPairs++;
	      
	      /* Got to trace back all paths from residue k and l in the same
		 way as edge Connectivity was calculated.
		 Mark direct connections from residue k in community i to
		 residue l in community j. */
	      if ((protein.shortDis[k][l] > 0) && (protein.shortDis[k][l] < BIG)) {
		done[l] = 0;
		nodeConn[l] = 1;
	      } else {
		done[l] = 1;
		nodeConn[l] = 0;
	      }
	      
	      /*if ((protein.shortDis[k][l] > max) && (protein.shortDis[k][l] < BIG)) {
		max = protein.shortDis[k][l];
		m = l;
		} */
	    } else {
	      /* Mark all node-node (k to l) connections between i and !j. */
	      done[l] = 0;
	      nodeConn[l] = 0;
	    }

	    /* Find longest node-node (k to l) distance from community i to any
	       other community. */
	    if ((protein.shortDis[k][l] > max) && (protein.shortDis[k][l] < BIG)) {
	      max = protein.shortDis[k][l];
	      m = l;
	    }
	  }
	}
	
	//fprintf(output4,"-----WHILE-----\n");
	/*  */
	while (max != 0) {
	  done[m] = 1;   /* Do not use this node as an endpoint in the future. */
	  predPtr = protein.Pred[k][m];
	  numPred = 0;   /* Number of shortest path branches from m to k. */
	  //fprintf(output4,"%d",m);
	  while (predPtr != NULL) {
	    //fprintf(output4,"->%d",predPtr->residue);
	    predPtr = predPtr->next;
	    numPred++;
	  }
	  predPtr = protein.Pred[k][m];
	  //fprintf(output4, "\n  numPred = %d\n",numPred);
	  
	  /* Add up fraction of path that is not within the current community. */
	  while (predPtr != NULL) {
	    nodeConn[predPtr->residue] += nodeConn[m]/numPred;
	    conn[predPtr->residue][m] = nodeConn[m]/numPred;
	    conn[m][predPtr->residue] = nodeConn[m]/numPred;
	    /*fprintf(output4, "    m: nodeConn[%d] = %f\n",m,nodeConn[m]);
	      fprintf(output4, "    p: nodeConn[%d] = %f\n",predPtr->residue,nodeConn[predPtr->residue]);
	      for(z=0; z<protein.nres; z++) {
	      fprintf(output4, "%1.2f ",nodeConn[z]);
	      }
	      fprintf(output4, "\n");*/
	    if (protein.Community->community[predPtr->residue] != protein.Community->community[m]) {
	      CommChange[protein.Community->community[m]][protein.Community->community[predPtr->residue]] += conn[m][predPtr->residue];
	      //fprintf(output4,"      conn[%d][%d] = %f\n",m,predPtr->residue,conn[m][predPtr->residue]);
	    }
	    predPtr = predPtr->next;
	  }

	  /* Pick new max path between residue k and farthest residue l in the molecule */
	  max = 0;
	  for (l=0; l<protein.nres; l++) {
	    if ((protein.shortDis[k][l] > max) && (done[l] != 1)) {
	      max = protein.shortDis[k][l];
	      m = l;
	    }
	  }
	}
	
	/* Set edge connectivities for all edges. */
	for (l=0; l<protein.nedges; l++) {
	  edgeConn[l] += conn[protein.edge[l]->res1][protein.edge[l]->res2];
	}
      }

      /* Determine whether communities are connected. */
      max1 = 0;
      for (l=0; l<protein.nedges; l++) 
	if (max1 < edgeConn[l])
	  max1 = edgeConn[l];

      /* Only calculate flow between communities that are directly connected. */
      if (max1 != 0) {
	TotalFlow=0;
	/*fprintf(output4, "Flow between communities %d and %d\n", i, j);
	  for (l=0; l<protein.nedges; l++) {
	  fprintf(output4, "%d %d %f\n",
	  protein.edge[l]->res1, protein.edge[l]->res2, edgeConn[l]);
	  }
	  fprintf(output4, "\n");
	  for (k=0; k<protein.Community->communityNumber; k++) {
	  for (l=0; l<protein.Community->communityNumber; l++) {
	  fprintf(output4, "%f ", CommChange[k][l]);
	  }
	  fprintf(output4, "\n");
	  }*/
	for (k=0; k<protein.Community->communityNumber; k++) {
	  TotalFlow += fabs(CommChange[i][k] - CommChange[k][i]);
	  /*fprintf(output4, "k %d Flow %f\n", k, TotalFlow);*/
	}
	DirectFlow = fabs(CommChange[i][j] - CommChange[j][i]);
	/*fprintf(output4, "i=%d, Flow %f\n", i, TotalFlow);
	  fprintf(output4, "i=%d, j=%d, Direct Flow %f\n", i, j, DirectFlow);
	  fprintf(output4, "          CommChange[%d][%d] = %f\n",i,j,CommChange[i][j]);
	  fprintf(output4, "          CommChange[%d][%d] = %f\n",j,i,CommChange[j][i]);*/
	percentFlow = DirectFlow/TotalFlow * 100;

	/* CommLink is a matrix of flags for community flows responsible for >15%. */
	if (percentFlow >= 15) {
	  CommLink[i][j] = 1;
	  CommLink[j][i] = 1;
	} else {
	  CommLink[i][j] = 0;
	  CommLink[j][i] = 0;
	}
	/*fprintf(output4, "Total Flow = %f, Direct Flow = %f, percentFlow = %f CommLink = %d \n",
	  TotalFlow, DirectFlow, percentFlow, CommLink[i][j]);*/

	/*fprintf(output4,"***COMMCHANGE***\n");
	  for (k=0; k<protein.Community->communityNumber; k++) {
	  for (l=0; l<protein.Community->communityNumber; l++) {
	  fprintf(output4, "%f ", CommChange[k][l]);
	  }
	  fprintf(output4, "\n");
	  }
	  fprintf(output4, "\n");*/

	for (k=0; k<protein.Community->communityNumber; k++) {
	  for (l=0; l<protein.Community->communityNumber; l++) {
	    TotalCommChange[k][l] += fabs(CommChange[k][l] - CommChange[l][k]);
	    fprintf(output4, "%f ", TotalCommChange[k][l]);
	  }
	  fprintf(output4, "\n");
	}
	fprintf(output4, "\n");
      }
    }
  }
  
  fprintf(output4, "Total Community Flow is: \n");
  for (i=0; i<protein.Community->communityNumber; i++) {
    for (j=0; j<protein.Community->communityNumber; j++) {
      if (CommLink[i][j] == 0)
	TotalCommChange[i][j] = 0;
      fprintf(output4, "%f ", TotalCommChange[i][j]);
    }
    fprintf(output4, "\n");
  }
  
  free(edgeConn);
  free(nodeConn);
  free(done);

  for (i=0; i<protein.nres; i++)
    free(conn[i]);
  free(conn);

  for (i=0; i<protein.Community->communityNumber; i++)
    free(CommChange[i]);
  free(CommChange);


  for (i=0; i<protein.Community->communityNumber; i++)
    free(CommLink[i]);
  free(CommLink);

  for (i=0; i<protein.Community->communityNumber; i++)
    free(TotalCommChange[i]);
  free(TotalCommChange);

  printf("<-getFlowComm\n");
  return;
}


void getIntercommunityFlow(Graph protein, FILE *output)
{
  printf("->getIntercommunityFlow\n");
  int i, j, k, l;
  float edgeFlow = 0.0;
  float **intercommunityFlow;
  
  printf("number of communities = %d\n", protein.Community->communityNumber);

  if ((intercommunityFlow = (float **) calloc(protein.Community->communityNumber, sizeof(float*))) == NULL)
    printf("No memory space allocatable for calculating flow between communities.\n");
  for (i=0; i<protein.Community->communityNumber; i++)
    if ((intercommunityFlow[i] = (float *) calloc(protein.Community->communityNumber, sizeof(float))) == NULL)
      printf("No memory space allocatable for calculating flow between communities.\n");

  for (i=0; i<protein.Community->communityNumber; i++)
    for (j=0; j<protein.Community->communityNumber; j++) {
      intercommunityFlow[i][j] = 0;
    }


  /* Loop over all pairs of communities. */
  for (i=0; i<protein.Community->communityNumber; i++) {
    for (j=(i+1); j<protein.Community->communityNumber; j++) {
      /* For each pair of communities, sum edge weight for edges
	 between the communities. */      
      for (k=0; k<protein.nres; k++) {

	/* if we're dealing with a residue from the current community,
	   loop through all residues and identify edges between this residue
	   and residues in community j. */
	if (protein.Community->community[k] == i) {
	  for (l=0; l <protein.nres; l++) {
	    if (protein.Community->community[l] == j) {
	      
	      /* Connections from residue k in community i to residue l in
		 community j contribute to the flow between i & j. */
	      // XXX - issues if there actually is a node pair with perfect correlation
	      if ((protein.dis[k][l] > 0) && (protein.dis[k][l] < BIG)) {
		//printf("protein.dis[%d][%d] = %d\n",k,l,protein.dis[k][l]);
		// Distance values are multiplied by 100 when read in.
		// Distances are calculated by adjacencyMatrix.tcl as
		//   -log(corrVal) where corrVal is from the carma correlation
		//   matrix
		edgeFlow = exp( (-(float)protein.dis[k][l]) / 100.0 );
		intercommunityFlow[i][j] += edgeFlow;
		intercommunityFlow[j][i] += edgeFlow;
	      }
	    }
	  }
	}
      }
    }
  }

  fprintf(output, "Intercommunity Flow is: \n");
  for (i=0; i<protein.Community->communityNumber; i++) {
    for (j=0; j<protein.Community->communityNumber; j++) {
      fprintf(output, "%f ", intercommunityFlow[i][j]);
    }
    fprintf(output, "\n");
  }

  for(i=0; i<protein.Community->communityNumber; i++)
    free(intercommunityFlow[i]);
  free(intercommunityFlow);

  printf("<-getIntercommunityFlow\n");
  return;
}


void characteristicPathLengthNodes(Graph protein, FILE *output)
{
  /*Graph dummy;*/
  	Graph *dummyPtr = copyGraph(&protein);
  	int i, j, k, l, m, currentResIndex;
	int startNode, midNode, endNode;
	float charPathLength;
	float *charPathLengthNode, *diffCharPathLengthNode;
	int *residues;
	int residueCount;
	
	startNode = 0;
	midNode = 74;
	endNode = protein.nres;
	
	/*Calculating characteristic path length as the average path length between nodes that are connected by a path.*/
	printf("->characteristicPathLength\n");
	
	k = 0;
	charPathLength = 0;
	for (i=0; i<protein.nres; i++) {
		for (j=0; j<protein.nres; j++) {
			if (protein.shortDis[i][j] < BIG) {
				charPathLength = charPathLength + protein.shortDis[i][j];
				k = k + 1;
			}
		}
	}
	k = k - protein.nres; /* Removing diagonal terms in shortDis matrix*/
	charPathLength = charPathLength/k;
	fprintf(output, "Characteristic path length of network is %f\n", charPathLength);
	
	residueCount = dummyPtr->nres - 1;
	residues = (int *) calloc(dummyPtr->nres, sizeof(int));
	for (i=0; i<dummyPtr->nres; i++) {
		residues[i] = i;
	}
	charPathLengthNode = (float *) calloc((endNode-startNode+1), sizeof(float));
	diffCharPathLengthNode = (float *) calloc((endNode-startNode+1), sizeof(float));
	m = 0;
	fprintf(output, "Difference in characteristic lengths measured by removing nodes \n"); 
	fprintf(output, "Node   Difference in Characteristic path length\n");	
	for (l=startNode; l<endNode; l++) {
	  currentResIndex = 0;
		for (i=0; i<protein.nres; i++) {
			for (j=0; j<protein.nres; j++) {
				if ((i != l) && (j!=l)) {				
					dummyPtr->dis[i][j] = protein.dis[i][j];
				} else {
					dummyPtr->dis[i][j] = 0;
				}
			}
			if (i != l) {
			  residues[i] = currentResIndex;
			  currentResIndex++;
			}
		} 
		floydWarshall(dummyPtr, residues, residueCount);
		/*printf("I got here\n");*/
		k = 0;
		charPathLengthNode[m] = 0;
		for (i=0; i<dummyPtr->nres; i++) {
			for (j=0; j<dummyPtr->nres; j++) {
				if (dummyPtr->shortDis[i][j] < BIG) {
					charPathLengthNode[m] = charPathLengthNode[m] + dummyPtr->shortDis[i][j];
					k = k + 1;
				}
			}
		}
		k = k - protein.nres; /* Removing diagonal terms in shortDis matrix*/
		charPathLengthNode[m] = charPathLengthNode[m]/k;
		printf("Node %d\n", l);
		diffCharPathLengthNode[m] = charPathLengthNode[m] - charPathLength;
		fprintf(output,"  %d       %f\n", l, diffCharPathLengthNode[m]);
		m = m + 1;
	}
	
	m = 0;
	fprintf(output, "Difference in characteristic lengths measured by removing edges \n"); 
	fprintf(output, "Node   Difference in Characteristic path length\n");
	for (l=startNode; l<endNode; l++) {
		for (i=0; i<protein.nres; i++) {
			for (j=0; j<protein.nres; j++) {
				if ((i != l) && (j!=l)) {				
					dummyPtr->dis[i][j] = protein.dis[i][j];
				} 
				else if ((i==l) && (i <= midNode) && (j <= midNode) && (j >= startNode)) {
					dummyPtr->dis[i][j] = protein.dis[i][j];
				} else if ((j==l) && (j <= midNode) && (i <= midNode) && (i >= startNode)) {
					dummyPtr->dis[i][j] = protein.dis[i][j];
				} else if ((i==l) && (i > midNode) && (j > midNode) && (j < endNode)) {
					dummyPtr->dis[i][j] = protein.dis[i][j];
				} else if ((j==l) && (j > midNode) && (i > midNode) && (i >= startNode)) {
					dummyPtr->dis[i][j] = protein.dis[i][j];
				} else {
					dummyPtr->dis[i][j] = 0;
				}
			}
		}
				
		floydWarshall(dummyPtr, residues, residueCount);
		/*printf("I got here\n");*/
		
		k = 0;
		charPathLengthNode[m] = 0;
		for (i=0; i<dummyPtr->nres; i++) {
			for (j=0; j<dummyPtr->nres; j++) {
				if (dummyPtr->shortDis[i][j] < BIG) {
					charPathLengthNode[m] = charPathLengthNode[m] + dummyPtr->shortDis[i][j];
					k = k + 1;
				}
			}
		}
		k = k - protein.nres; 

		/* Removing diagonal terms in shortDis matrix*/
		charPathLengthNode[m] = charPathLengthNode[m]/k;
		printf("Node %d\n", l);
		diffCharPathLengthNode[m] = charPathLengthNode[m] - charPathLength;
		fprintf(output,"  %d       %f\n", l, diffCharPathLengthNode[m]);
		m = m + 1;
	}
	printf("<-characteristicPathLength\n");
}


void graphProperties(Graph protein, FILE *output)
{
	int i,j;
	int startRnaNode, endRnaNode;
	float distanceRna, degreeRna, distanceProt, degreeProt, distanceInt, degreeInt;
	
	startRnaNode = 0;
	endRnaNode = 74;
	distanceRna = 0;
	degreeRna = 0;
	for (i=startRnaNode; i<=endRnaNode; i++) {
		for (j=(i + 1); j<=endRnaNode; j++) {
			if (protein.dis[i][j] != 0) {
				distanceRna = distanceRna + protein.dis[i][j];
				degreeRna = degreeRna + 1;
			}
		}
	}
	distanceRna = distanceRna/degreeRna;
	degreeRna = degreeRna/(endRnaNode - startRnaNode + 1);
	printf("Average distance in RNA = %f \n", distanceRna);
	printf("Average degree in RNA = %f \n", degreeRna);

	distanceProt = 0;
	degreeProt = 0;
	for (i=(endRnaNode + 1); i < protein.nres; i++) {
		for (j=(i + 1); j < protein.nres; j++) {
			if (protein.dis[i][j] != 0) {
				distanceProt = distanceProt + protein.dis[i][j];
				degreeProt = degreeProt + 1;
			}
		}
	}
	distanceProt = distanceProt/degreeProt;
	degreeProt = degreeProt/(protein.nres - (endRnaNode - startRnaNode + 1));
	printf("Average distance in Protein = %f \n", distanceProt);
	printf("Average degree in Protein = %f \n", degreeProt);
	
	distanceInt = 0;
	degreeInt = 0;
	for (i=startRnaNode; i<=endRnaNode; i++) {
		for (j=(endRnaNode + 1); j<protein.nres; j++) {
			if (protein.dis[i][j] != 0) {
				distanceInt = distanceInt + protein.dis[i][j];
				degreeInt = degreeInt + 1;
			}
		}
	}
	distanceInt = distanceInt/degreeInt;
	degreeInt = degreeInt/(endRnaNode - startRnaNode + 1);
	printf("Average distance in Protein-RNA interface = %f \n", distanceInt);
	printf("Average degree in Protein-RNA interface per tRNA nucleotide = %f \n", degreeInt);

	return;
}


void printEdgeConnectivities(Graph *protein, FILE *output) {

  float** edgeConMatrix;
  int i=0;
  int j=0;

  if ((edgeConMatrix = (float **) calloc(protein->nres, sizeof(float*))) == NULL)
    printf("Error::printEdgeConnectivities: could not allocate memory for edgeConMatrix.\n");
  for (i=0; i<protein->nres; i++)
    if ((edgeConMatrix[i] = (float *) calloc(protein->nres, sizeof(float))) == NULL)
      printf("Error::printEdgeConnectivities: could not allocate memory for edgeConMatrix.\n");

  for(i=0; i<protein->nres; i++) {
    for(j=0; j<protein->nres; j++) {
      edgeConMatrix[i][j] = 0.0;
    }
  }

  /* Loop over all edges and put betweenness information in edgeConMatrix. */
  for (i=0; i<protein->nedges; i++) {
    if (protein->edge[i]->edgeConn != 0.0) {
      edgeConMatrix[protein->edge[i]->res1][protein->edge[i]->res2] = protein->edge[i]->edgeConn;
      edgeConMatrix[protein->edge[i]->res2][protein->edge[i]->res1] = protein->edge[i]->edgeConn;
    }
  }

  /* Print the matrix to output. */
  if (protein->line != NULL) {
    fprintf(output,"%s",protein->line);
  }
  for(i=0; i<protein->nres; i++) {
    for(j=0; j<protein->nres; j++) {
      fprintf(output,"%f ", edgeConMatrix[i][j]);
    }
    fprintf(output,"\n");
  }
  
  return;
}


void printShortestDistances(Graph *protein, FILE *output) {

  int i, j;

  /* Print shortest distances matrix. */
  for (i=0; i<protein->nres; i++) {
    for (j=0; j<protein->nres; j++) {
      fprintf(output,"%d ", protein->shortDis[i][j]);
    }
    fprintf(output,"\n");
  }
  
  return;
}


void printPredecessors(Graph *protein, FILE *output) {

  int i, j;
  nodePtr tempNode;

  /* Print predecessor list matrix. */
  /* Build lists containing all residues for each path from i to j. */
  for (i=0; i<protein->nres; i++) {
    for (j=0; j<protein->nres; j++) {
      tempNode = protein->Pred[i][j];
      if (tempNode != NULL) {
	fprintf(output,"%d",tempNode->residue);
	//printf("%d",tempNode->residue);
	tempNode = next(tempNode);
      } else {
	fprintf(output,"-1");
	//printf("-1");
      }
      while (tempNode != NULL) {
	fprintf(output,",%d",tempNode->residue);
	//printf(",%d",tempNode->residue);
	tempNode = next(tempNode);
      }
      fprintf(output," ");
      //printf(" ");
    }
    fprintf(output,"\n");
    //printf("\n");
  }

  return;
}


void readShortestDistances(Graph *protein, char *inputFileName) {

  int i,j,distance;
  FILE *input;

  /* protein->shortDis is already allocated in get() */

  /* Check that input file exists. */
  if ( (input=fopen(inputFileName, "r")) == NULL ) {
    printf("Input file (%s) does not exist.\n",inputFileName);
    exit(1);
  }

  for (i=0; i<protein->nres; i++) {
    for (j=0; j<protein->nres; j++) {
      fscanf(input, "%d ", &distance);
      //printf("%d ", distance);
      /*printf("%f\n", &dis[i][j]);*/
      protein->shortDis[i][j] = distance;
    }
    /*printf("%d,%d\n", i,j);*/
    fscanf(input,"\n");
    //printf("\n");
  }
  /*printf("hey4-3\n");*/
  fclose(input);
      
  return;
}


void readPredecessors(Graph *protein, char *inputFileName) {

  int i,j,predecessor;
  int strLen = 64 * protein->nres;
  FILE *input;
  char *entries[protein->nres];
  char line[strLen];
  char *predString;

  /* protein->Pred is already allocated in get() */

  /* Check that input file exists. */
  if ( (input=fopen(inputFileName, "r")) == NULL ) {
    printf("Input file (%s) does not exist.\n",inputFileName);
    exit(1);
  }

  for (i=0; i<protein->nres; i++) {
    fgets(line,strLen,input);
    entries[0] = strtok(line," ");
    for (j=1; j<protein->nres; j++) {
      //fscanf(input, "%d ", &distance);
      //printf("%f\n", &dis[i][j]);
      //protein->dis[i][j] = distance;
      entries[j] = strtok(NULL," ");
      //printf("%s ",entries[j]);
    }
    for (j=0; j<protein->nres; j++) {
      /* Get indices possibly separated by commas. */
      predString = strtok(entries[j],",");
      while (predString != NULL) {
	sscanf(predString,"%d",&predecessor);
	pushNode(&protein->Pred[i][j], predecessor);
	predString = strtok(NULL,",");
	//printf("hey2\n");
	//printf("%d-",predecessor);
      }
      //printf(" ");
    }
    //printf("\n");
    /*printf("%d,%d\n", i,j);*/
    //fscanf(input,"\n");
  }
  /*printf("hey4-3\n");*/
  fclose(input);
      
  return;
}

