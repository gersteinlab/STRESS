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


void readSequenceAlignment(Alignment *Aln, FILE *input)
{
	int maxLineSize = 65534;
	char line[maxLineSize+1], sequence[100000];
	char* firstSpace;
	char* firstBar;
	int last, i, j;
	
	printf("Reading sequence alignment\n");
	
	/*Finding the number of columns in the alignment and the number of sequences in the alignment*/
	Aln->Nseq = 0;
	sequence[0] = '\0';
	if (fgets(line,maxLineSize, input) == NULL) 
		return;
	
	if (line[0] == '>')
		Aln->Nseq = Aln->Nseq + 1;
	fgets(line,maxLineSize, input);
	
	while (!feof(input) && (line[0]!='>')) {
		/* Removing any trailing characters*/
		last = strlen(line) - 1;
		while ((line[last] == '\r') || (line[last] == '\n')) {
			line[last] = '\0';
			last = last - 1;
		}
		strcat(sequence,line);
		fgets(line,maxLineSize, input);
	}
	
	Aln->Naln = strlen(sequence);
	
	if (line[0] == '>')
		Aln->Nseq = Aln->Nseq + 1;
	while (!feof(input)) {
		if (fgets(line,maxLineSize, input) == NULL) 
			break;
		if (line[0] == '>')
			Aln->Nseq = Aln->Nseq + 1;
		
	}
	printf("Number of sequences are %d\n", Aln->Nseq);
	printf("Number of columns in alignment are %d\n", Aln->Naln);
	rewind(input);
	
	if ((Aln->Sequence = (seqPtr *) calloc(Aln->Nseq, sizeof(seqPtr))) == NULL)
		printf("No memory space allocatable for saving alignment.\n");
	i = -1;
	
	while(!feof(input)) {
		
		// Get the next line
		if (fgets(line, maxLineSize, input) == NULL) 
			break;
		
		// Remove any trailing characters.
		last = strlen(line) - 1;
		while (line[last] == '\r' || line[last] == '\n') {
			line[last--] = '\0';
		}
		
		//If this is a new sequence
		if (line[0] == '>') {
			i = i + 1;
			
			//Creating memory for new sequence
			Aln->Sequence[i] = (seqPtr)malloc(sizeof(Sequence));
			if ((Aln->Sequence[i]->alignedSequence = (char *) calloc(Aln->Naln, sizeof(char))) == NULL)
				printf("No memory space allocatable for saving alignment.\n");
			Aln->Sequence[i]->name[0] = '\0';
			
			//Removing the characters after first Space in name
			firstSpace = strstr(line, " ");
			if (firstSpace != NULL) *firstSpace = '\0';
			
			//Using the characters between first Bar and second Bar for the name
			firstBar = strstr(line, "|");
			j = 0;
			if (firstBar != NULL) {
				while ((*(firstBar + 1) != '|') && (*(firstBar + 1) != '\0')) {
					line[j] = *(firstBar + 1);
					firstBar = firstBar + 1;
					j = j + 1;
				}
				line[j] = '\0';
			}
			strcat(Aln->Sequence[i]->name, line);
			
			//printf("Name of sequence %d is %s\n", i, Aln->Sequence[i]->name);
		} else {
			//Adding line to current sequence
			strcat(Aln->Sequence[i]->alignedSequence, line);
			/*printf("%d %s\n", i, Aln->Sequence[i]->alignedSequence);*/
		}
	}
	/*
	for (i=0; i<Aln->Nseq; i++) {
		printf("Sequence %d is \n", i);
		printf("%s\n", Aln->Sequence[i]->alignedSequence);
	}
	*/
	return;	
}

void readPdbFileNames(Alignment *Aln, FILE *input)
{
	int maxLineSize = 65534;
	char line[maxLineSize+1];
	int i;
	
	i = 0;
	printf("Reading pdb file names\n");
	while (!feof(input)) {
		fgets(line,maxLineSize, input);
		// Setting pdb file name for each sequence as each line in the pdbFileNameList file
		if (i < Aln->Nseq) {
			Aln->Sequence[i]->pdbName[0] = '\0';
			strcat(Aln->Sequence[i]->pdbName,line);
		}
		//printf("%d %s", i, line);
		if (line[0] != '\0')
			i = i + 1;
		line[0] = '\0';
	}
	
	// Checking for consistency in number of sequences
	if (i != Aln->Nseq) {
		printf("Mismatch in number of sequences and number of pdb files i %d Nseq %d\n", i, Aln->Nseq);
		exit(1);
	}
	return;
}


void mappingSequenceAlignment(Alignment *Aln) 
{
	int i, seqIndex, alnIndex;
	seqPtr Seq;
	
	for (i=0; i<Aln->Nseq; i++) {
		
		// Finding number of residues in sequence without the gaps
		Seq = Aln->Sequence[i];
		Seq->Nres = numberResiduesSequence(Seq);
		printf("Number of residues in sequence %s is %d\n", Seq->name, Seq->Nres);
		
		// Mapping residues onto columns of alignment and reverse mapping (aln column->seqResidueNumber with gaps mapped to -1)
		if ((Seq->mapping = (int *) calloc(Seq->Nres, sizeof(int))) == NULL)
			printf("No memory space allocatable for saving mapping in alignment.\n");
		if ((Seq->revMapping = (int *) calloc(Aln->Naln, sizeof(int))) == NULL)
			printf("No memory space allocatable for saving mapping in alignment.\n");
		for (seqIndex=0; seqIndex<Aln->Naln; seqIndex++)
			Seq->revMapping[seqIndex] = -1;
		seqIndex=0;
		for (alnIndex=0; alnIndex<Aln->Naln; alnIndex++) {
			if ((Seq->alignedSequence[alnIndex]!='-') && (Seq->alignedSequence[alnIndex] != ' ')) {
				Seq->mapping[seqIndex] = alnIndex;
				Seq->revMapping[alnIndex] = seqIndex;
				seqIndex++;
			}
		}
		/*
		printf("Mapping for sequence %d is :", i);
		for (seqIndex=0; seqIndex<Seq->Nres; seqIndex++) {
			printf("%d ", Seq->mapping[seqIndex]);
		}
		printf("\n");
		*/	
	}
	return;
}


int numberResiduesSequence(seqPtr Seq)
{
	int i, Naln, Nres;
	
	Naln = strlen(Seq->alignedSequence);
	// printf("%s\n", Seq->alignedSequence);
	Nres = 0;
	// printf("%s\n", Seq->name);
	for (i=0; i<Naln; i++) {
		/*printf("%c ", Seq->alignedSequence[i]);*/
		if ((Seq->alignedSequence[i] != '-') && (Seq->alignedSequence[i] != ' ')) {
			Nres = Nres + 1;
			//printf("%d %c %d\n", i, Seq->alignedSequence[i], Nres);
		}
	}
	// printf("\n");
	return Nres;
}

void deleteAlignment(Alignment *Aln)
{
	int seqIndex;
	
	if (Aln->nodeCons != NULL)
		free(Aln->nodeCons);
	if (Aln->edgeCons != NULL)
		free(Aln->edgeCons);
	if (Aln->edgeCol1 != NULL)
		free(Aln->edgeCol1);
	if (Aln->edgeCol2 != NULL)
		free(Aln->edgeCol2);

	for (seqIndex=0; seqIndex<Aln->Nseq; seqIndex++) {
		deleteSequence(Aln->Sequence[seqIndex]);
	}
	
	return;
}

void deleteSequence(seqPtr Seq)
{
	if (Seq->alignedSequence != NULL)
		free(Seq->alignedSequence);
	if (Seq->mapping != NULL)
		free(Seq->mapping);
	if (Seq->revMapping != NULL)
		free(Seq->revMapping);
	if (Seq->network.dis != NULL)
		free(Seq->network.dis);
		
}
