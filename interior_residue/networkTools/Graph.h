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

#ifndef GRAPH_H
#define GRAPH_H

#include "Community.h"
#include "List.h"
#include "Flow.h"
#include "Tree.h"

/**
 * This is a structure that indicates that there is an edge between res1 and 
 * res2. The edge connectivity of this edge is stored in the variable edgeConn.
 */
typedef struct {
  int res1, res2;
  float edgeConn;
} EdgeStr;

typedef EdgeStr *edgePtr;


/**
 * Copy constructor.
 * Copy all of the member data into another EdgeStr.
 * @param edge EdgeStr to copy.
 * @return newly created copy.
 */
edgePtr copyEdge(edgePtr edge);


/**
 * Contains all graph information needed for calculation of pairwise path,
 * betweenness, community structure, and flow.
 */
typedef struct {
  char line[10000];  /**< Atomselection string needed at top of betweenness file */
  int nres;          /**< Number of nodes or residues */
  int nedges;        /**< Number of edges */
  int **dis;         /**< Matrix of edge weights */
  int **shortDis;    /**< Matrix of shortest path distances */
  nodePtr (**Pred);  /**< Predecessor for each node in shortest path. For path from i to j, Pred[i][j] gives node connected to j. */
  float *nodeConn;   /**< Connectivity of each node */
  edgePtr *edge;     /**< Structures having information edges including its connectivity */
  commPtr Community; /**< Information about community structure in the protein */
  flowPtr Flow;
} Graph;


/**
 * Copy constructor.
 * Copy all of the member data into another Graph.
 * @param graph Graph to copy.
 * @return newly created copy.
 */
Graph* copyGraph(Graph *graph);

/**
 * Read in the input network file.
 * This is generated using a tcl file 
 * normally.  The input should have connectivity information of the nodes and 
 * the weights of the connected nodes in a matrix form.  These inputs are 
 * stored in the network structure and a pointer to this network is the 
 * input to this function.
 * @param prot graph to be generated.
 * @param input
 */
void get(Graph *prot, FILE *input);

/**
 * Remove an edge from a graph.
 * @param graph Graph from which to remove the edge.
 * @param edgeIndex index of edge to remove.
 */
void removeEdge(Graph *graph, int edgeIndex);

/**
 * Computes all pairwise shortest distance and shortest path matrices
 * using the Floyd-Warshall algorithm. These are stored in the network and a 
 * pointer to the network structure is input to this function.
 * @param prot graph for which the shortest pairwise distances are found
 * @param residues array of residue IDs for pairwise distance calculation.
 * @param residueCount number of residue IDs in residues array; may be smaller than the number of nodes in prot.
 */
void floydWarshall(Graph *prot, int *residues, int residueCount);

/**
 * Finds the shortest path between the source and target residues.
 * The pointer to flow is input because the shortest path is stored 
 * in the flow structure.  The network is also input to this function.
 * @param protein graph for shortest path calculation
 * @param Flow stores shortest path information
 */
void shortestPath(Graph *protein);
//void shortestPath(Graph protein, flowPtr *Flow);

/**
 * Find all suboptimal paths connecting any of the source 
 * residues to any one of the target residues.  The suboptimal paths are
 * placed in the structure Flow (hence pointer to it is input to function).  The
 * graph is input to the function in order to find the shortest distance and 
 * predecessor of each node along the path.
 * @param protein graph for suboptimal path calculation
 * @param edgeLengthOffset additional length threshold for suboptimal path definition
 * @param Flow stores shortest path information
 * @param output3
 * @param logfile
 */
void subOpt(Graph protein, int edgeLengthOffset, flowPtr *Flow, FILE *output3, FILE *output);

/**
 * Enumerate shortest paths between the source and target residues.
 * That is all the paths with the shortest distance connecting any 
 * one of the source residues to any one of the target residues. Pointers
 * to the path and the network itself are input arguments to this function.  
 * In addition, the source and target for which all shortest paths have to be 
 * found is also input to this function.
 * @param prot graph for shortest path calculation
 * @param source index of source
 * @param target index of target
 * @param PathPtr path
 */
void allPaths(Graph *prot, int source, int target, nodePtr *PathPtr);

/**
 * Find the node connectivity of all the nodes in the network.
 * These values are stored in the network structure. Hence, a pointer to the 
 * network is input to this function.
 * @param prot graph for connectivity calculation
 */
void nodeConnectivity(Graph *prot);

/**
 * Find all edges in the network and place them in an array of 
 * the variable type edgeStr. The edges information is stored in the array of 
 * edge structures.  A pointer to the network structure is hence input to this 
 * function.
 * @param prot graph for edge extraction
 */
void getEdges(Graph *prot);

/**
 * Find the edge connectivities of all edges in a given subnetwork.
 * The edge connectivities are stored in the array of edge structures in the 
 * graph structure.  Hence a pointer to the graph is input to this function.
 * @param prot graph for edge connectivity calculation
 * @param residues array of residue IDs for pairwise distance calculation.
 * @param residueCount number of residue IDs in residues array.
 */
void edgeConnectivity(Graph *proteinGraph, int *residues, int residueCount);

/**
 * Add nodes to the community tree.
 * When an edge is deleted, one or more new communities may be produced.  For
 * each new community, add a node to the community tree pointing to the
 * community node from whence it came.
 * @param proteinGraph
 * @param tempGraph
 * @param numComm
 * @param numCommPrev
 */
void addNodeToCommunityTree(Graph *proteinGraph, Graph *tempGraph, int numComm, int numCommPrev);

/**
 * Apply the Girvan Newman algorithm to find the optimal community structure
 * of the network.  The community information is stored in the 
 * community structure (hence pointer to it is input).  The network is also 
 * input to this function.
 * @param protein graph for optimal community structure calculation
 * @param output3
 * @param output4
 */
void gnewman(Graph *protein, FILE *output3, FILE *output4, char *community_tcl_fl);

/**
 * Get the community information if intermediate files are present for doing
 * flow analysis across communities easily.
 * In this function, the number of communities and the community
 * for each node is read as input and a pointer to the community
 * structure in the graph is one of the arguments.
 * @param protein graph for retrieval of community information
 * @param input
 */
void getComm(Graph protein, FILE *input);

/**
 * Calculate the amount of communication between different communities in the
 * network.
 * This information helps in calculating the width of
 * edges in the community diagram for the paper. For this subroutine, only the 
 * network is input in the form of Graph.
 * @param protein graph for community communication calculation
 * @param output4
 */
void getFlowComm(Graph protein, FILE *output4);

/**
 * Calculate the amount of communication between different communities in the
 * network.
 * This information helps in calculating the width of
 * edges in the community diagram for the paper. For this subroutine, only the 
 * network is input in the form of Graph.
 * @param protein graph for community communication calculation
 * @param output4
 */
void getIntercommunityFlow(Graph protein, FILE *output);

/**
 * This function is used to find the nodes from the tRNA that affect the 
 * characteristic path length the most. 
 * This is similar to the Ruth Nussinov measure for most important nodes in the 
 * network.
 * @param protein graph for characteristic path length calculations.
 * @param output4
 */
void characteristicPathLengthNodes(Graph protein, FILE *output);

/**
 * This function is used to find some properties of the graph such as degree and 
 * average distance between nearest neighbors in two separate macromolecules in 
 * the graph.
 * @param protein graph for calculating graph properties.
 * @param output4
 */
void graphProperties(Graph protein, FILE *output);

/**
 * Print a file with the matrix of pairwise edge connectivities.
 * @param protein
 * @param output
 */
void printEdgeConnectivities(Graph *protein, FILE *output);

/**
 * Print a file with the matrix of shortest path distances.
 * @param protein
 * @param output
 */
void printShortestDistances(Graph *protein, FILE *output);

/**
 * Print a file with the matrix of predecessor lists.
 * @param protein
 * @param output
 */
void printPredecessors(Graph *protein, FILE *output);

/**
 * Read a file with the matrix of shortest path distances.
 * @param protein
 * @param output
 */
void readShortestDistances(Graph *protein, char *inputFileName);

/**
 * Read a file with the matrix of predecessor lists.
 * @param protein
 * @param output
 */
void readPredecessors(Graph *protein, char *inputFileName);



#endif
