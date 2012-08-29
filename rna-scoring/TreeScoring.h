/*
	A class to score a RNA sructure stored in a tree (node) data structure.
*/

#ifndef TREESCORING_H
#define TREESCORING_H
#include "StructureReader.h" 
#include "data.h"



int ScoreNode(TreeNode* node, int* RNA, nndb_constants* param, int length);

#endif
