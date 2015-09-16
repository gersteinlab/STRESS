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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "Community.h"
#include "List.h"
#include "Flow.h"
#include "Graph.h"

typedef struct {
	char *alignedSequence; /**< The aligned sequence */
	char name[1000];       /**< Name of the biomolecule */
	char pdbName[1000];    /**< Name of the pdb file */
	char networkName[1000];/**< Name of the network file*/ 
	int *mapping;          /**< Mapping of the unaligned sequence to the aligned sequence*/
	int *revMapping;       /**< Reverse mapping of the columns in the aligned sequence to the unaligned sequence;*/
	int Nres;	       /**< Number of residues in unaligned sequence */
	Graph network;         /**< Network for sequence */
	char fileName[1000];          /**< All output of conservation related to this sequence is output in this file*/
} Sequence;

typedef Sequence *seqPtr;

typedef struct {
	int Naln;         		/**< Number of columns in alignment */
	int Nedges;       		/**< Number of edges (all possible) in alignment */
	int Nseq;         		/**< Number of sequences in alignment */
	int *edgeCol1, *edgeCol2; 	/**< Edges are formed between edgeCol1 and edgeCol2 */
	seqPtr *Sequence; 		/**< Structures have information about each sequence in the alignment */
	float *nodeCons;  		/**< Conservation of each column in the alignment */
	float *edgeCons;                /**< Conservation of each edge in the alignment */
	int **numMatches;               /**< Number of matches of edges between aligned residues in each pair of sequences */
	float **edgeIdentity;           /**< Percentage of edge matches for all pairs of sequences*/
	float **avgDiffShortestDistance;     /**< Average difference of shortest distance between all aligned residues for each pair of sequences */
} Alignment;

typedef Alignment *alignmentPtr;

void readSequenceAlignment(Alignment *Aln, FILE *input);

/**
 * Read a sequence alignment
 * @param Aln alignment to be read in.
 * @param input input file from which to read the alignment
 */

void readPdbFileNames(Alignment *Aln, FILE *input);

/**
 * Read the list of pdb file names for each protein in the alignment
 * @param Aln alignment to be read in.
 * @param input input file from which to read the pdb file names
 */
 
void mappingSequenceAlignment(Alignment *Aln);

/**
 * Mapping the unaligned sequence to the aligned sequence.
 * @param Aln alignment for which sequences are mapped.
 */
 
int numberResiduesSequence(seqPtr Seq);

/**
 * Finding the number of current Sequence
 * @param Seq current sequence for which the number of residues (without gaps) are found
 */
 
void deleteAlignment(Alignment *Aln);

/**
 * Delete the alignment
 * @param Aln alignment to be deleted
 */
 
void deleteSequence(seqPtr Seq);

/**
 * Delete a sequence
 * @param Seq sequence to be deleted
 */
 
 
#endif
