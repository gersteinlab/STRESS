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
 * Author(s): Anurag Sethi
 */

#ifndef NETWORKCOMPARISON_H
#define NETWORKCOMPARISON_H

#include "Community.h"
#include "List.h"
#include "Flow.h"
#include "Graph.h"


/**
 * Read a sequence alignment
 * @param Aln alignment for which the networks have to be read.
 * @param input input file from which to read the file names for all the networks
 */
void readNetworkFiles(Alignment *Aln, FILE *input);

/**
 * Find the number of edges in the alignment
 * @param Aln alignment to find number of edges for.
 */
int numberEdgesAlignment(Alignment *Aln);

/**
 * Create the edges in the alignment
 * @param Aln alignment to find number of edges for.
 */
void createAlnEdges(Alignment *Aln);

/**
 * Conservation of a node in an alignment
 * @param Aln alignment for which the conservation has to be calculated
 */
void nodeConservation(Alignment *Aln);

/**
 * Conservation of an edge in the alignment
 * @param Aln alignment for which the conservation has to be calculated
 */
void edgeConservation(Alignment *Aln);

/**
 * Number of edge matches between pair of sequence in the alignment
 * @param Aln alignment for which the identity has to be calculated
 */
void numberEdgeMatches(Alignment *Aln);

/**
 * Identity of edge matches between pair of sequences in the alignment
 * @param Aln alignment for which the identity has to be calculated
 * @param output output file to which edge identity distance matrix is written
 */
void calculateEdgeIdentity(Alignment *Aln, FILE *output);

/**
 * calculate shortest distance of all networks in the alignment
 * @param Aln alignment for which the shorest distance of all networks are calculated
 */
void calculateShortestDistance(Alignment *Aln);

/**
 * calculate the average of the difference in shortest path length between all pairs of aligned columns in two sequences
 * @param Aln alignment for which the calculations are done
 * @param output output file to which the average distance matrix is written
 */
void diffShortestDistanceAligned(Alignment *Aln, FILE *output);
 

#endif
