#ifndef COMMUNITY_H
#define COMMUNITY_H

#include "List.h"


/**
 * Stores the community structure of the network.
 * Leaves is an array of pointers into a tree that displays how the different
 * nodes in the network are split into different communities at various
 * iterations.  The
 * community variable is an array that contains the community information of
 * each node at the optimal community distribution numComm.  Q is an array with
 * the Q value for different number of communities and optQ is the optimal Q
 * value.
 */
typedef struct {
  int communityNumber;       /**< Optimal number of communities */
  int *community;    /**< Community information of each node at the optimal community distribution */
  float optQ;        /**< Optimal Q */
  float *Q;          /**< Array of Q values for different numbers of communities */
  nodePtr *leaf;     /**< Tree that displays how different nodes in the network are split into different communities */
} commStr;

typedef commStr *commPtr;

#endif
