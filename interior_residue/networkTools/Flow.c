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
#include "Flow.h"

//void getRes(flowPtr *Flow, FILE *res1file, FILE *res2file)
void getRes(flowPtr *Flow, int source, int target)
{
  //int i;
  //char line[100];

  *Flow = malloc(sizeof(flowStr));
  //(*Flow)->nres1 = 0;
  (*Flow)->nres1 = 1;

  /*while (fgets(line, 100, res1file) != NULL)
    (*Flow)->nres1 = (*Flow)->nres1 + 1;
    printf("Number of residues at one end of signal transmitted is %d.\n", (*Flow)->nres1);*/
  //(*Flow)->nres2 = 0;
  (*Flow)->nres2 = 1;

  /*while (fgets(line, 100, res2file) != NULL)
    (*Flow)->nres2 = (*Flow)->nres2 + 1;
    printf("Number of residues at other end of signal transmitted is %d.\n", (*Flow)->nres2);*/
  
  if (((*Flow)->source = (int *) calloc((*Flow)->nres1, sizeof(int))) == NULL)
    printf("No memory space allocatable for reading Residue 1 selection.\n");
  if (((*Flow)->target = (int *) calloc((*Flow)->nres2, sizeof(int))) == NULL)
    printf("No memory space allocatable for reading Residue 2 selection.\n");

  /*rewind(res1file);
    rewind(res2file);
    for (i=0; i<(*Flow)->nres1; i++) {                                                   
    fscanf(res1file, "%d", &((*Flow)->source[i]));
    }
    for (i=0; i<(*Flow)->nres2; i++) {
    fscanf(res2file, "%d", &((*Flow)->target[i]));
    }
    fclose(res1file);
    fclose(res2file);*/

  (*Flow)->source[0] = source;
  (*Flow)->target[0] = target;

  return;
}
