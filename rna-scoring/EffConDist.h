/*
	A class to record effective contact distances of pairs in
	a sructure stored in a tree (node) data structure.
*/

#ifndef EFFCONDIST_H
#define EFFCONDIST_H
#include "StructureReader.h" 
#include "data.h"

int GetConDist(TreeNode* node, int* RNA, nndb_constants* param);
//int main(int argc, char* argv[]);

#endif
