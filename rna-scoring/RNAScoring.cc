extern "C" {
    #include "StructureReader.h"
    #include "RNAScoring.h"
    #include "TreeScoring.h"    
}
//#include <ctype.h>
//#include <limits.h>
#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include "loader.h"
//#include "data.h"
//#include "constants.h"
#include <cstdlib>
#include <iostream>
#include<string.h>
#include "LoopScoring.h"
#include "options.h"
#include <math.h>
#include "PFtest.h"

using namespace std;
extern char paramDir[200];

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
		//fprintf(stderr, "USAGE: RNAScoring <filename>\n");
		help();
		return 1;
    }
	if(strcmp(argv[1],"--pf_test")==0){
		PFtest(argc,argv);
		return 0;
	}

//printf("prg is not in pftest\n");
	parse_options(argc, argv);
	
	//printf("hi");
	char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U',
                      'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    //ResultBundle* resultBundle = CreateFromFile(argv[1]);
    ResultBundle* resultBundle = CreateFromFile(seqfile);
	int length = resultBundle->length;
    TreeNode* tree = resultBundle->treenode;
    int* RNA = resultBundle->RNA_seq;
    
   // PrintTree(tree, 0);
	 nndb_constants* param = populate(strcat(paramDir, "data/Turner99"), 1);
	
	int one; int two; int three; int four; int score; 
	
	int tree_score = ScoreNode(tree, resultBundle->RNA_seq, param, length);
	printf("Tree score is ");
	printf("%.2f\n", (double)tree_score/100);
	

    return 0;
}
