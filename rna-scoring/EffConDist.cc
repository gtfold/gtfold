extern "C" {
#include "StructureReader.h"
}
//#include <limits.h>
#include "EffConDist.h"
#include "loader.h"
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <cstdlib>
#include <iostream>
//#include <math.h>
using namespace std;


//Shel: Function for scoring a node (recursive)
int GetConDist(TreeNode* node, int* RNA, nndb_constants* param){
	int numPairedChildren;
	numPairedChildren = 0;
	int i;
	
	for (i = 0 ; i < node->numChildren ; i++){ 
       // find location and number of paired children 
       //and add scores of associated loops
		if ((node->children[i])->isPair) {
			GetConDist(node->children[i], RNA, param);
			numPairedChildren += 1;
		}
	}
	printf("%d\t%d\t%d\n", node->lowBase.index, node->highBase.index, node->numChildren + numPairedChildren +1);
	
	return 0;
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
		fprintf(stderr, "USAGE: EffConDist <filename>\n");
		return 1;
    }
	
	char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U',
		'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    ResultBundle* resultBundle = CreateFromFile(argv[1]);
    TreeNode* tree = resultBundle->treenode;
    int* RNA = resultBundle->RNA_seq;
    
	// PrintTree(tree, 0);
	nndb_constants* param = populate("data/Turner99", 1);
	
	int one; int two; int three; int four; int score; 
	
	GetConDist(tree, resultBundle->RNA_seq, param);
	
    return 0;
}



