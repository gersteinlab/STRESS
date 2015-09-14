#ifndef _network_h
#define _network_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pdbFile.h"

/**
 * This is a structure that indicates that there is an edge between res1 and 
 * res2. The edge connectivity of this edge is stored in the variable edgeConn.
 */
typedef struct {
  int res1, res2;
  //float edgeConn;
} EdgeStr;

typedef EdgeStr *edgePtr;

/**
 * Contains all graph information needed for calculation of pairwise path,
 * betweenness, community structure, and flow.
 */
typedef struct {
  int nres;          /**< Number of nodes or residues */
  int nedges;        /**< Number of edges */
  int **dis;         /**< Matrix of edge weights */
  int **shortDis;    /**< Matrix of shortest path distances */
  edgePtr *edge;     /**< Structures having information edges including its connectivity */
  float *LC;		 /**< Local clustering quantity */
} Graph;

void constructNetwork(Graph *network, pdbFile *macromolecule);
/**
 * Constructs the network given the macromolecule
 * returns all edges within network
 */

void floydWarshall(Graph *network, int *residues, int residueCount);

/**
 * Finds the shortest path between the source and target residues.
 * The network is input to this function.
 * @param network graph for shortest path calculation
 * @param residue is array of residues for which the shortest path calculation is done
 * @param residueCount is number of elementsi n residue
 */

void calculateLC(Graph *network, int *residue, int residueCount);

/**
 * Calculates LC or Local Clustering defined in Mitternacht and Berezovsky, Prot Engg, Des & Sel, pp 1-5, 2010
 * @param network graph for shortest path calculation
 * @param residue is array of residues for which the LC calculation is done
 * @param residueCount is number of elementsi n residue
 */
 
#endif