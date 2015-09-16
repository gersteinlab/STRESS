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

#ifndef FLOW_H
#define FLOW_H

#include "List.h"

/**
 * The main structure of flow in the protein between defined sets of residues.
 * Source and Target is an array of integers containing the 
 * information from the user of the source and target for the allosteric signal. 
 * nres1 and nres2 are the number of residues in source and target respectively.
 * shortestDis is the shortest distance between any of the source residues and 
 * any of the the target residues.  This may not be a unique path.  If there are 
 * multiple residues sources and targets that lead to the same distance as 
 * shortestDis, all the sources and targets with the same distance are saved in
 * a linked list to which sources and targets point. The subOptPtr will be a
 * linked list containing a particular suboptimal path from one of the targets
 * to one of the source residues. The allPathPtr is a pointer pointing to a
 * linked list with all the suboptimal paths from one of the targets to one of
 * the source residues.
 */
typedef struct {
  int *source;
  int *target;
  int nres1;
  int nres2;
  int shortestDis;
  nodePtr allPathPtr;
  nodePtr subOptPtr;
  nodePtr sources;  /** source nodes leading to shortest distances */
  nodePtr targets;  /** target nodes for shortest distances */
  int *nodeMatch;
  int *edgeMatch;
} flowStr;

typedef flowStr *flowPtr;



/**
 * Read in the 2nd and 3rd input files.
 * These files should contain 
 * the residue numbers for the source and target of the allosteric signal. These 
 * variables are stored in the structure Flow and a pointer to Flow is input to 
 * this function.
 * @param Flow flow structure to be created
 * @param source node index of source
 * @param target node index of target
 */
//void getRes(flowPtr *Flow, FILE *res1file, FILE *res2file);
void getRes(flowPtr *Flow, int source, int target);

#endif
