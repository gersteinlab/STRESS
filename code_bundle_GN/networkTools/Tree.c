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
#include "Tree.h"


treePtr newTree(treePtr parent, int val)
{
  treePtr newTree = (treePtr) malloc(sizeof(struct Tree));
  newTree->value = val;
  newTree->parent = parent;
  newTree->child1 = NULL;  
  newTree->child2 = NULL;  

  return newTree;
}


treePtr copyTree(treePtr tree)
{
  //printf("->copyTree\n");
  if (tree == NULL) {
    return NULL;
  }

  treePtr newTree = (treePtr) malloc(sizeof(struct Tree));
  newTree->value = tree->value;
  newTree->parent = NULL;
  newTree->child1 = copyTree(tree->child1);
  newTree->child1->parent = newTree;
  newTree->child2 = copyTree(tree->child2);
  newTree->child2->parent = newTree;

  //printf("<-copyTree\n");
  return newTree;
}


treePtr addChild(treePtr tree, int child, int value)
{
  treePtr success = NULL;

  if (child == 1 &&
      tree->child1 == NULL) {
    tree->child1 = newTree(tree,value);
    success = tree->child1;
  } else if (child == 2 &&
	     tree->child2 == NULL) {
    tree->child2 = newTree(tree,value);
    success = tree->child2;
  }

  return success;
}


void deleteTree(treePtr tree)
{

  if (tree != NULL) {
    deleteTree(tree->child1);
    deleteTree(tree->child2);
    free(tree);
  }

  return;
}


int matchTree(treePtr tree, treePtr matchTree)
{
  int match = 0;
  //treePtr treeNode = tree;
  //treePtr matchTreeNode = matchTree;

  return match;
}


void printTree(treePtr tree)
{

  //printf(">printTree\n");
  if (tree->child1 != NULL &&
      tree->child2 != NULL) {
    printf("(");
    printTree(tree->child1);
    printf(",");
    printTree(tree->child2);
    printf(")");
  } else if (tree->child1 != NULL) {
    printf("(");
    printTree(tree->child1);
    printf(")");
  } else if (tree->child2 != NULL) {
    printf("(");
    printTree(tree->child2);
    printf(")");
  }
  printf("%d",tree->value);

  //printf("<printTree\n");
  return;
}


void fprintTree(FILE *file, treePtr tree)
{
  //printf(">fprintTree\n");
  if (tree->child1 != NULL &&
      tree->child2 != NULL) {
    fprintf(file,"(");
    fprintTree(file,tree->child1);
    fprintf(file,",");
    fprintTree(file,tree->child2);
    fprintf(file,")");
  } else if (tree->child1 != NULL) {
    fprintf(file,"(");
    fprintTree(file,tree->child1);
    fprintf(file,")");
  } else if (tree->child2 != NULL) {
    fprintf(file,"(");
    fprintTree(file,tree->child2);
    fprintf(file,")");
  }
  fprintf(file,"%d",tree->value);

  //printf("<fprintTree\n");
  return;
}
