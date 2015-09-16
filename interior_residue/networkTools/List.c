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
#include "List.h"


listPtr copyList2(listPtr list)
{
  listPtr newList = (listPtr) malloc(sizeof(struct List));
  return newList;
}

void push(listPtr *list, int i)
{
  
  
  return;
}

int pop(listPtr *list)
{
  return 0;
}

nodePtr copyList(nodePtr list)
{
  //printf("->copyList\n");
  if (list == NULL) {
    return NULL;
  }

  nodePtr currentListNode = list;
  nodePtr newList = (nodePtr) malloc(sizeof(struct Node));
  nodePtr currentNewListNode = newList;
  nodePtr tempNode;

  while(currentListNode != NULL) {
    currentNewListNode->residue = currentListNode->residue;
    if(currentListNode->next != NULL) {
      tempNode = (nodePtr) malloc(sizeof(struct Node));
      currentNewListNode->next = tempNode;
      currentNewListNode = tempNode;
    } else {
      currentNewListNode->next = NULL;
    }
    currentListNode = currentListNode->next;
  }

  //printf("<-copyList\n");
  return newList;
}


void pushNode(nodePtr *head, int i)
{
  //nodePtr newPtr, currPtr, prevPtr;
  //printf("->pushNode\n");
  nodePtr newPtr, lastPtr;
  
  newPtr = (nodePtr)malloc(sizeof(struct Node));  
  newPtr->residue = i;
  newPtr->next = NULL;
  
  lastPtr = last(*head);
  if (lastPtr == NULL) {
    *head = newPtr;
  } else {
    lastPtr->next = newPtr;
  }
  
  //printf("<-pushNode\n");
  return;
}


int popNode(nodePtr *head)
{
  //printf("->popNode\n");
  int value = -1;
  nodePtr prevPtr, currPtr;
  
  prevPtr = *head;
  currPtr = prevPtr->next;
  if (currPtr != NULL) {
    while (currPtr->next != NULL) {
      prevPtr = currPtr;
      currPtr = currPtr->next;
    }
    value = currPtr->residue;
    free(prevPtr->next);
    prevPtr->next = NULL;
  } else {
    value = (*head)->residue;
    free(*head);
    *head = NULL;
  }

  //printf("<-popNode\n");
  return value;
}


nodePtr next(nodePtr head)
{
  return head->next;
}


nodePtr last(nodePtr head)
{
  //printf("->last\n");
  nodePtr currPtr = head;

  if (currPtr != NULL) {
    while (currPtr->next != NULL) {
      currPtr = currPtr->next;
    }
  }

  //printf("<-last\n");
  return currPtr;
}


void append(nodePtr *head, nodePtr appendList)
{
  //printf("->append\n");
  nodePtr currPtr = appendList;

  while (currPtr != NULL) {
    pushNode(&*head,currPtr->residue);
    currPtr = next(currPtr);
  }

  //printf("<-append\n");
  return;
}


nodePtr reverse(nodePtr head)
{
  nodePtr reverseList = NULL;
  nodePtr copy = copyList(head);

  while (copy != NULL) {
    pushNode(&reverseList,popNode(&copy));
  }

  return reverseList;
}


void delete(nodePtr *head)
{
  //printf("->delete\n");
  nodePtr prevPtr, currPtr;
  
  prevPtr = *head;
  currPtr = prevPtr->next;
  if (currPtr != NULL) {
    while (currPtr->next != NULL) {
      prevPtr = currPtr;
      currPtr = currPtr->next;
    }
    free(prevPtr->next);
    prevPtr->next = NULL;
  } else {
    free(*head);
    *head = NULL;
  }

  //printf("<-delete\n");
  return;
}


void deleteList(nodePtr *head)
{
  nodePtr prevPtr, currPtr;

  if (*head != NULL) {
    prevPtr = *head;
    currPtr = prevPtr->next;
    if (currPtr != NULL) {
      deleteList(&currPtr);
    }
    
    free(*head);
    //free(prevPtr);
    *head = NULL;
  }

  return;
}


int hasCycle(nodePtr head)
{
  int cycle = 0;

  nodePtr checkPtr = head;
  nodePtr currPtr = head;

  while (currPtr != NULL) {
    checkPtr = currPtr->next;
    while (checkPtr != NULL) {
      if (checkPtr->residue == currPtr->residue) {
	cycle = 1;
      }
      checkPtr = checkPtr->next;
    }
    currPtr=currPtr->next;
  }

  return cycle;
}


int match(nodePtr list, nodePtr matchList)
{
  int match = 0;
  nodePtr listNode = list;
  nodePtr matchListNode = matchList;

  // Loop through list looking for match to matchList
  while (listNode != NULL && match == 0) {
    while (listNode != NULL &&
	   matchListNode != NULL &&
	   listNode->residue == matchListNode->residue) {
      //printf("(%d,%d) ", listNode->residue, matchListNode->residue);
      listNode = listNode->next;
      matchListNode = matchListNode->next;
      if (matchListNode == NULL) {
	match = 1;
      }
    }
    listNode = listNode->next;
    matchListNode = matchList;
  }
  printf("match: %d\n",match);

  return match;
}


void printList(nodePtr head)
{
  nodePtr currPtr;
	
  currPtr = head;
  while (currPtr!=NULL) {
    printf("%d, ", currPtr->residue);
    currPtr = currPtr->next;
  }
  return;
}

void fprintList(FILE *file, nodePtr head)
{
  nodePtr currPtr;
	
  currPtr = head;
  while (currPtr!=NULL) {
    fprintf(file, "%d, ", currPtr->residue);
    currPtr = currPtr->next;
  }
  return;
}
