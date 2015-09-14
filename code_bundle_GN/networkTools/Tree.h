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

#ifndef TREE_H
#define TREE_H

/**
 * Linked list with pointers to both head and tail nodes to facilitate
 * both pushing and popping.
 */
struct Tree {
  int value;
  struct Tree *parent;
  struct Tree *child1;
  struct Tree *child2;
};

typedef struct Tree *treePtr;


/**
 * Constructor.
 * Create a new tree.
 * @param parent Parent for this tree, NULL if none.
 * @param val Value for this tree.
 * @return Newly created tree.
 */
treePtr newTree(treePtr parent, int val);

/**
 * Copy constructor.
 * Copy all of the member data into another Tree.
 * @param tree Tree to copy.
 * @return Newly created copy.
 */
treePtr copyTree(treePtr tree);

/**
 * Add child if there is none present already.
 * @param tree Tree to add the child to.
 * @param child Which child to add (1 or 2).
 * @param value Value for the child.
 * @return Pointer to child, NULL if unsuccessful.
 */
treePtr addChild(treePtr tree, int child, int value);

/**
 * Delete an entire linked list.
 * As the pointer may be deleted in this function, a pointer to the first
 * pointer is input to this function.
 * @param tree Pointer to head of linked list
 */
void deleteTree(treePtr tree);

/**
 * Determine if a tree contains a second tree.
 * The first tree is searched for an exact copy of the second tree
 * @param tree Pointer to root of tree.
 * @param matchTree Pointer to head of linked list to search for
 * @return 0 if there are no matches, 1 if at least one match exists
 */
int matchTree(treePtr tree, treePtr matchTree);

/**
 * Print out the linked list residues.
 * Prints out the residues in the whole linked list in the 
 * order they are pushed into the linked list. A pointer to the first element 
 * of the linked list is input to this function.
 * @param head Pointer to root of tree.
 */
void printTree(treePtr tree);

/**
 * Print out the linked list residues to a file.
 * Prints out the residues in the whole linked list in the 
 * order they are pushed into the linked list. A pointer to the first element 
 * of the linked list is input to this function.
 * @param file output file
 * @param head pointer to head of linked list
 */
void fprintTree(FILE *file, treePtr tree);

#endif
