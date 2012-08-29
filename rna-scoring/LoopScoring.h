/*
 A class to score an RNA structure reading in from a file.
 */

#ifndef RNASCROING_H
#define RNASCORING_H
#include "TreeScoring.h"
#include "data.h"
#include "options.h"

#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))
#define auPen(i, j) (( ((i)==3 || (j)==3) && ( (i)==0 || (i)==2 || (j)==0 || (j)==2 )) ? 1 : 0)
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

int eMUnpairedRegion(int i1, int j1, int i2, int j2, int* RNA, nndb_constants* param);
int eH(int i, int j, int* RNA, nndb_constants* param);
int eS(int i, int j, int* RNA, nndb_constants* param);
int eL(int i, int j, int ip, int jp, int* RNA, nndb_constants* param);
int eM(TreeNode* node, int* pairedChildren, int numPairedChildren, int* RNA, nndb_constants* param);
int eE(TreeNode* node, int* pairedChildren, int numPairedChildren, int* RNA, nndb_constants* param, int length);

#endif
