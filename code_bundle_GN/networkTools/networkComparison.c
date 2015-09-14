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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Graph.h"
#include "alignment.h"
#include "networkComparison.h"

#define BIG1 1000000

void readNetworkFiles(Alignment *Aln, FILE *input)
{
	int maxLineSize = 65534;
	char line[maxLineSize+1];
	int last, i;
	Graph *networkPtr;
	FILE *networkFile;
	
	printf("Reading networks\n");
	i = 0;
	while (!feof(input)) {
		
		// Get the next line
		if (fgets(line, maxLineSize, input) == NULL) 
			break;
		
		// Remove any trailing characters.
		last = strlen(line) - 1;
		while (line[last] == '\r' || line[last] == '\n') {
			line[last--] = '\0';
		}
		
		//printf("file is %s \n", line);
		//Reading network file for sequence i
		networkPtr = &(Aln->Sequence[i]->network);
		if ( (networkFile = fopen(line, "r")) == NULL ) {
			printf("Input file (contactMap) %s does not exist.\n", line);
			exit(1);
		}
		Aln->Sequence[i]->networkName[0] = '\0';
		strcat(Aln->Sequence[i]->networkName, line);
		get(networkPtr, networkFile);
		getEdges(networkPtr);
		
		//Checking whether the number of nodes in the network is equivalent to the number of residues on the protein/RNA
		if (networkPtr->nres != Aln->Sequence[i]->Nres) {
			printf("There is a mismatch in the number of nodes and the number of residues for sequence %d\n", i);
			printf("Number of residues in alignment %d.  Number of nodes in network %d\n", Aln->Sequence[i]->Nres, networkPtr->nres);
		}
		i = i + 1;
		fclose(networkFile);
	}
	free(networkPtr->Pred); // JULY_1
	return;
}

int numberEdgesAlignment(Alignment *Aln)
{
	int nedges;
	int alnIndex1, alnIndex2, seqIndex;
	int res1, res2;
	
	nedges = 0;
	
	// All unique edges between alignment columns are counted.
	for (alnIndex1=0; alnIndex1<Aln->Naln; alnIndex1++) {	
		for (alnIndex2=(alnIndex1+1); alnIndex2<Aln->Naln; alnIndex2++) {
			for (seqIndex=0; seqIndex<Aln->Nseq; seqIndex++) {
				res1 = Aln->Sequence[seqIndex]->revMapping[alnIndex1];
				res2 = Aln->Sequence[seqIndex]->revMapping[alnIndex2];
				if ((res1 != -1) && (res2 != -1)) 
					if (Aln->Sequence[seqIndex]->network.dis[res1][res2] != 0) {
						nedges = nedges + 1;
						seqIndex = Aln->Nseq;
					}
			}
		}
	}
	printf("Number of edges in alignment are %d\n", nedges);
	
	return nedges;
}
	
void createAlnEdges(Alignment *Aln)
{
	int alnIndex1, alnIndex2, seqIndex;
	int res1, res2;
	int edgeIndex;
	
	if ((Aln->edgeCol1 = (int *) calloc(Aln->Nedges, sizeof(int))) == NULL)
			printf("No memory space allocatable for saving conservation of node.\n");
	if ((Aln->edgeCol2 = (int *) calloc(Aln->Nedges, sizeof(int))) == NULL)
			printf("No memory space allocatable for saving conservation of node.\n");

	edgeIndex=0;
	for (alnIndex1=0; alnIndex1<Aln->Naln; alnIndex1++) {	
		for (alnIndex2=(alnIndex1+1); alnIndex2<Aln->Naln; alnIndex2++) {
			for (seqIndex=0; seqIndex<Aln->Nseq; seqIndex++) {
				res1 = Aln->Sequence[seqIndex]->revMapping[alnIndex1];
				res2 = Aln->Sequence[seqIndex]->revMapping[alnIndex2];
				if ((res1 != -1) && (res2 != -1)) 
					if (Aln->Sequence[seqIndex]->network.dis[res1][res2] != 0) {
						Aln->edgeCol1[edgeIndex] = alnIndex1;
						Aln->edgeCol2[edgeIndex] = alnIndex2;
						edgeIndex++;
						seqIndex = Aln->Nseq;
					}
			}
		}
	}
	/*
	for (edgeIndex=0; edgeIndex<Aln->Nedges; edgeIndex++) 
		printf("Edge %d between column %d and %d\n", edgeIndex, Aln->edgeCol1[edgeIndex], Aln->edgeCol2[edgeIndex]);*/
	return;
}
	
	
void nodeConservation(Alignment *Aln)
{
	int alnIndex, seqIndex, resIndex;
	FILE *output;
	
	printf("->nodeConservation\n");
	if ((Aln->nodeCons = (float *) calloc(Aln->Naln, sizeof(float))) == NULL)
			printf("No memory space allocatable for saving conservation of node.\n");
	
	printf("Measuring Node Conservation\n");
	for (alnIndex=0; alnIndex < Aln->Naln; alnIndex++) {
		Aln->nodeCons[alnIndex] = 0;
		for (seqIndex=0; seqIndex < Aln->Nseq; seqIndex++) {
			if ((Aln->Sequence[seqIndex]->alignedSequence[alnIndex] != '-') && (Aln->Sequence[seqIndex]->alignedSequence[alnIndex] != ' '))
				Aln->nodeCons[alnIndex]++;
		}
		Aln->nodeCons[alnIndex] = Aln->nodeCons[alnIndex]/Aln->Nseq;
		//printf("%d %f\n", alnIndex, Aln->nodeCons[alnIndex]);
	}
	
	for (seqIndex=0; seqIndex < Aln->Nseq; seqIndex++) {
		output = fopen(Aln->Sequence[seqIndex]->fileName, "a");
		fprintf(output,"Node Conservation\n");
		for (resIndex=0; resIndex < Aln->Sequence[seqIndex]->Nres; resIndex++) {
			fprintf(output,"%d %f\n", resIndex, Aln->nodeCons[Aln->Sequence[seqIndex]->mapping[resIndex]]);
		}	
		fclose(output);
	}
	
	printf("<-nodeConservation\n");
	return;
}

void edgeConservation(Alignment *Aln)
{
	int seqIndex, alnEdgeIndex, seqEdgeIndex;
	int alnIndex1, alnIndex2;
	int res1, res2;
	FILE *output;
	
	printf("->edgeConservation\n");
	printf("Measuring Edge Conservation\n");
	if ((Aln->edgeCons = (float *) calloc(Aln->Nedges, sizeof(float))) == NULL)
			printf("No memory space allocatable for saving conservation of edge.\n");
	
	//printf("Edge Conservation is:\n");
	for (alnEdgeIndex=0; alnEdgeIndex < Aln->Nedges; alnEdgeIndex++) {
		Aln->edgeCons[alnEdgeIndex] = 0;
		alnIndex1 = Aln->edgeCol1[alnEdgeIndex];
		alnIndex2 = Aln->edgeCol2[alnEdgeIndex];
		for (seqIndex=0; seqIndex < Aln->Nseq; seqIndex++) {
			res1 = Aln->Sequence[seqIndex]->revMapping[alnIndex1];
			res2 = Aln->Sequence[seqIndex]->revMapping[alnIndex2];
			if ((res1 != -1) && (res2 != -1))
				if (Aln->Sequence[seqIndex]->network.dis[res1][res2] != 0)
					Aln->edgeCons[alnEdgeIndex]++;
		}
		Aln->edgeCons[alnEdgeIndex] = Aln->edgeCons[alnEdgeIndex]/Aln->Nseq;
		//printf("%d %f\n", alnEdgeIndex, Aln->edgeCons[alnEdgeIndex]);
	}
	
	for (seqIndex=0; seqIndex < Aln->Nseq; seqIndex++) {
		output = fopen(Aln->Sequence[seqIndex]->fileName, "a");
		fprintf(output,"Edge Conservation\n");
		for (seqEdgeIndex=0; seqEdgeIndex < Aln->Sequence[seqIndex]->network.nedges; seqEdgeIndex++) {
			//printf("%d %d\n", Aln->Sequence[seqIndex]->network.edge[edgeIndex]->res1, Aln->Sequence[seqIndex]->network.edge[edgeIndex]->res2);
			for (alnEdgeIndex=0; alnEdgeIndex < Aln->Nedges; alnEdgeIndex++) {
				alnIndex1 = Aln->edgeCol1[alnEdgeIndex];
				alnIndex2 = Aln->edgeCol2[alnEdgeIndex];
				res1 = Aln->Sequence[seqIndex]->network.edge[seqEdgeIndex]->res1;
				res2 = Aln->Sequence[seqIndex]->network.edge[seqEdgeIndex]->res2;
				//alnIndex1 = Aln->Sequence[seqIndex]->mapping[res1];
				
				if (((Aln->Sequence[seqIndex]->mapping[res1] == alnIndex1) && (Aln->Sequence[seqIndex]->mapping[res2] == alnIndex2)) || ((Aln->Sequence[seqIndex]->mapping[res1] == alnIndex2) && (Aln->Sequence[seqIndex]->mapping[res2] == alnIndex1))) {
					fprintf(output, "%d %d %f\n", res1, res2, Aln->edgeCons[alnEdgeIndex]);
					alnEdgeIndex = Aln->Nedges;
				}
			}
		}
		fclose(output);
	}
	
	printf("<-edgeConservation\n");
	return;
}

void numberEdgeMatches(Alignment *Aln)
{
	int seqIndex1, seqIndex2;
	int alnIndex1, alnIndex2;
	int res1, res2, res3, res4;
	
	if ((Aln->numMatches = (int **) calloc(Aln->Nseq, sizeof(int*))) == NULL)
		printf("No memory space allocatable for calculating number of matches.\n");
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++)
		if ((Aln->numMatches[seqIndex1] = (int *) calloc(Aln->Nseq, sizeof(int))) == NULL)
			printf("No memory space allocatable for calculating number of matches. \n");
	
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		for (seqIndex2=(seqIndex1+1); seqIndex2<Aln->Nseq; seqIndex2++) {
			for (alnIndex1=0; alnIndex1<Aln->Naln; alnIndex1++) {
				res1 = Aln->Sequence[seqIndex1]->revMapping[alnIndex1];
				res2 = Aln->Sequence[seqIndex2]->revMapping[alnIndex1];
				if ((res1!=-1) && (res2!=-1)) {
					for (alnIndex2=(alnIndex1+1); alnIndex2<Aln->Naln; alnIndex2++) {
						res3 = Aln->Sequence[seqIndex1]->revMapping[alnIndex2];
						res4 = Aln->Sequence[seqIndex2]->revMapping[alnIndex2];
						if ((res3!=-1) && (res4!=-1)) {
							if ((Aln->Sequence[seqIndex1]->network.dis[res1][res3] != 0) && (Aln->Sequence[seqIndex2]->network.dis[res2][res4] != 0)) {
								Aln->numMatches[seqIndex1][seqIndex2]++;
							}
						}
					}
				}
			}
			Aln->numMatches[seqIndex2][seqIndex1] = Aln->numMatches[seqIndex1][seqIndex2];
		}
	}

	/*
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		for (seqIndex2=0; seqIndex2<Aln->Nseq; seqIndex2++)
			printf("%d  ", Aln->numMatches[seqIndex1][seqIndex2]);
		printf("\n");
	}*/
	return;
}

void calculateEdgeIdentity(Alignment *Aln, FILE *output)
{
	int seqIndex1, seqIndex2;
	int norm;
	
	printf("->calculateEdgeIdentity\n");
	if ((Aln->edgeIdentity = (float **) calloc(Aln->Nseq, sizeof(float*))) == NULL)
		printf("No memory space allocatable for calculating identity of edge matches.\n");
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++)
		if ((Aln->edgeIdentity[seqIndex1] = (float *) calloc(Aln->Nseq, sizeof(float))) == NULL)
			printf("No memory space allocatable for calculating identity of edge matches. \n");

	//printf("The edge Identity score is:\n");
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		Aln->edgeIdentity[seqIndex1][seqIndex1] = 1.0;
		for (seqIndex2=(seqIndex1+1); seqIndex2<Aln->Nseq; seqIndex2++) {
			if (Aln->Sequence[seqIndex1]->network.nedges<=Aln->Sequence[seqIndex2]->network.nedges)
				norm = Aln->Sequence[seqIndex1]->network.nedges;
			else 
				norm = Aln->Sequence[seqIndex2]->network.nedges;
			//printf("seqIndex1 %d seqIndex2 %d numMatches %d norm %d\n", seqIndex1, seqIndex2, Aln->numMatches[seqIndex1][seqIndex2], norm);
			Aln->edgeIdentity[seqIndex1][seqIndex2] = (float) Aln->numMatches[seqIndex1][seqIndex2]/norm;
			Aln->edgeIdentity[seqIndex2][seqIndex1] = Aln->edgeIdentity[seqIndex1][seqIndex2];
		}
	}
	
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		fprintf(output, "%s ", Aln->Sequence[seqIndex1]->name);
		for (seqIndex2=0; seqIndex2<Aln->Nseq; seqIndex2++)
			fprintf(output, "%f  ", (1 - Aln->edgeIdentity[seqIndex1][seqIndex2]));
		fprintf(output, "\n");
	}
	printf("<-calculateEdgeIdentity\n");
	return;
}

void calculateShortestDistance(Alignment *Aln)
{
	int residueCount, seqIndex1, i;
	int *residues;
	Graph *networkPtr;
	
	printf("->calculateShortestDistance\n");
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		residueCount = Aln->Sequence[seqIndex1]->Nres;
		residues = (int *) calloc(Aln->Sequence[seqIndex1]->Nres, sizeof(int));
		for (i=0; i<residueCount; i++) {
			residues[i] = i;
		}
		networkPtr = &(Aln->Sequence[seqIndex1]->network);
		floydWarshall(networkPtr, residues, residueCount);
		free(residues);
	}
	printf("<-calculateShortestDistance\n");
	return;
}

void diffShortestDistanceAligned(Alignment *Aln, FILE *output)
{
	int seqIndex1, seqIndex2, resIndex1, resIndex2, resIndex3, resIndex4, ctr;
	Graph *networkPtr1, *networkPtr2;
	float temp;
	
	printf("->diffShortestDistanceAligned\n");
	if ((Aln->avgDiffShortestDistance = (float **) calloc(Aln->Nseq, sizeof(float*))) == NULL)
		printf("No memory space allocatable for calculating identity of edge matches.\n");
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++)
		if ((Aln->avgDiffShortestDistance[seqIndex1] = (float *) calloc(Aln->Nseq, sizeof(float))) == NULL)
			printf("No memory space allocatable for calculating identity of edge matches. \n");
	
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		networkPtr1 = &(Aln->Sequence[seqIndex1]->network);
		for (seqIndex2=(seqIndex1+1); seqIndex2<Aln->Nseq; seqIndex2++) {
			Aln->avgDiffShortestDistance[seqIndex1][seqIndex2] = 0;
			ctr = 0;
			networkPtr2 = &(Aln->Sequence[seqIndex2]->network);
			for (resIndex1=0; resIndex1<Aln->Sequence[seqIndex1]->Nres; resIndex1++) {
				// resIndex3 (resIndex4) is the residue in sequence 2 aligned to resIndex1 (resIndex3) of sequence 1. 
				// If there is no residue aligned in sequence 2 in that column, then resIndex3 (or resIndex4) = - 1
				resIndex3 = Aln->Sequence[seqIndex2]->revMapping[Aln->Sequence[seqIndex1]->mapping[resIndex1]];
				if (resIndex3 != -1) {
					for (resIndex2=(resIndex1+1); resIndex2<Aln->Sequence[seqIndex1]->Nres; resIndex2++) {
						resIndex4 = Aln->Sequence[seqIndex2]->revMapping[Aln->Sequence[seqIndex1]->mapping[resIndex2]];
						if (resIndex4 != -1) {
							//printf("%d %d %d %d %d %d\n", resIndex1, resIndex2, resIndex3, resIndex4, seqIndex1, seqIndex2);
							if ((networkPtr1->shortDis[resIndex1][resIndex2] != BIG1) && (networkPtr2->shortDis[resIndex3][resIndex4] != BIG1)) {
								temp = abs(networkPtr1->shortDis[resIndex1][resIndex2] - networkPtr2->shortDis[resIndex3][resIndex4]);
								//printf("seqIndex1 %d seqIndex2 %d temp is %f\n", seqIndex1, seqIndex2, temp);
								Aln->avgDiffShortestDistance[seqIndex1][seqIndex2] = Aln->avgDiffShortestDistance[seqIndex1][seqIndex2] + temp;
								ctr = ctr + 1;
							}
						}
					}
				}
			}
			//printf("%f %d\n", Aln->avgDiffShortestDistance[seqIndex1][seqIndex2], ctr);
			Aln->avgDiffShortestDistance[seqIndex1][seqIndex2] = Aln->avgDiffShortestDistance[seqIndex1][seqIndex2]/(100*ctr);
			Aln->avgDiffShortestDistance[seqIndex2][seqIndex1] = Aln->avgDiffShortestDistance[seqIndex1][seqIndex2];
		}
	}
	
	for (seqIndex1=0; seqIndex1<Aln->Nseq; seqIndex1++) {
		fprintf(output, "%s ", Aln->Sequence[seqIndex1]->name);
		for (seqIndex2=0; seqIndex2<Aln->Nseq; seqIndex2++)
			fprintf(output, "%f  ",  Aln->avgDiffShortestDistance[seqIndex1][seqIndex2]);
		fprintf(output, "\n");
	}
	printf("<-diffShortestDistanceAligned\n");
	return;
}
