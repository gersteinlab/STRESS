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

#ifndef LIST_H
#define LIST_H

/**
 * Linked list with pointers to both head and tail nodes to facilitate
 * both pushing and popping.
 */
struct List {
  struct Node *head;
  struct Node *tail;
};

typedef struct List *listPtr;

/**
 * Nodes for a generic doubly-linked list.
 * It has been used to point to a predecessor of the current residue "j"
 * (residue) in the shortest path between 2 given residues i and j.
 */
struct Node {
  int residue;
  struct Node *prev;
  struct Node *next;
};

typedef struct Node *nodePtr;


/**
 * Copy constructor.
 * Copy all of the member data into another List.
 * @param list List to copy.
 * @return newly created copy.
 */
listPtr copyList2(listPtr list);

/**
 * Push a new element onto a linked list.
 * So adds a new member to a linked list that already exists or creates the
 * first member of a linked list. The residue to be inserted is integer i
 * and the pointer to the first element is input to this linked list as i
 * might be the first element of the linked list.
 * @param head pointer to head of linked list
 * @param i int to be added to linked list
 */
void push(listPtr *list, int i);

/**
 * Delete the last member of a linked list and return it.
 * As the pointer may be deleted in this function, a pointer to the first
 * pointer is input to this function.
 * @param head pointer to head of linked list
 * @return value of the node that is popped
 */
int pop(listPtr *list);







/**
 * Copy constructor.
 * Copy all of the member data into another List.
 * @param list List to copy.
 * @return newly created copy.
 */
nodePtr copyList(nodePtr list);

/**
 * Push a new element onto a linked list.
 * So adds a new member to a linked list that already exists or creates the
 * first member of a linked list. The residue to be inserted is integer i
 * and the pointer to the first element is input to this linked list as i
 * might be the first element of the linked list.
 * @param head pointer to head of linked list
 * @param i int to be added to linked list
 */
void pushNode(nodePtr *head ,int i);

/**
 * Delete the last member of a linked list and return it.
 * As the pointer may be deleted in this function, a pointer to the first
 * pointer is input to this function.
 * @param head pointer to head of linked list
 * @return value of the node that is popped
 */
int popNode(nodePtr *head);

/**
 * Return the next element in the list.
 * @param head pointer to head of linked list
 * @return next node in the list
 */
nodePtr next(nodePtr head);

/**
 * Return the next element in the list.
 * @param head pointer to head of linked list
 * @return last node in the list
 */
nodePtr last(nodePtr head);

/**
 * Append the elements of one list onto another list.
 * As the head pointer may be created in this function, a pointer to the first
 * pointer is input to this function.
 * @param head pointer to head of linked list
 * @param appendList pointer to linked list that will be appended onto head
 */
void append(nodePtr *head, nodePtr appendList);

/**
 * Create a new list everse the elements of one list.
 * @param head pointer to head of linked list
 * @return pointer to the head of the newly created reverse list
 */
nodePtr reverse(nodePtr head);

/**
 * Delete the last member of a linked list.
 * As the pointer may be deleted in this function, a pointer to the first
 * pointer is input to this function.
 * @param head pointer to head of linked list
 */
void delete(nodePtr *head);

/**
 * Delete an entire linked list.
 * As the pointer may be deleted in this function, a pointer to the first
 * pointer is input to this function.
 * @param head pointer to head of linked list
 */
void deleteList(nodePtr *head);

/**
 * Determine if there are cycles in this path.
 * A cycle occurs when the same integer appears more than once in the list.
 * @param head pointer to head of linked list
 * @return 0 if there are no cycles, 1 if at least one cycle exists
 */
int hasCycle(nodePtr head);

/**
 * Determine if a list contains a second list.
 * The first list is searched for an exact copy of the second list
 * @param list pointer to head of linked list
 * @param matchList pointer to head of linked list to search for
 * @return 0 if there are no matches, 1 if at least one match exists
 */
int match(nodePtr list, nodePtr matchList);

/**
 * Print out the linked list residues.
 * Prints out the residues in the whole linked list in the 
 * order they are pushed into the linked list. A pointer to the first element 
 * of the linked list is input to this function.
 * @param head pointer to head of linked list
 */
void printList(nodePtr head);

/**
 * Print out the linked list residues to a file.
 * Prints out the residues in the whole linked list in the 
 * order they are pushed into the linked list. A pointer to the first element 
 * of the linked list is input to this function.
 * @param file output file
 * @param head pointer to head of linked list
 */
void fprintList(FILE *file, nodePtr head);

#endif
