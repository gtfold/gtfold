#include "TreeScoring.h"
#include "StructureReader.h"
#include "LoopScoring.h"
//#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int printOn2=1;
//Shel: Function for scoring a node (recursive)
int ScoreNode(TreeNode* node, int* RNA, nndb_constants* param, int length){
	int result;
	result = 0;
	int *pairedChildren;
	pairedChildren = NULL;
	int numPairedChildren;
	numPairedChildren = 0;
	int i;
	
	for (i = 0 ; i < node->numChildren ; i++){ 
       // find location and number of paired children 
       //and add scores of associated loops
		if ((node->children[i])->isPair) {
			//Manoj: Changed the code to generate warnings in case of pairing in structure which are not valid
			 /*Base lb = (node->children[i])->lowBase.base;
                                char lbChar;
                                if(lb==1)lbChar='A';else if(lb==2)lbChar='C';else if(lb==4)lbChar='G';else if(lb==8)lbChar='U';else lbChar="X";
                                Base hb = (node->children[i])->highBase.base;
                                char hbChar;
                                if(hb==1)hbChar='A';else if(hb==2)hbChar='C';else if(hb==4)hbChar='G';else if(hb==8)hbChar='U';else hbChar="X";
                                printf("CHECKING: bases %d and %d (%c%c) if they can pair!\n",(node->children[i])->lowBase.index, (node->children[i])->highBase.index,lbChar, hbChar);*/
			if(!canPair((node->children[i])->lowBase.base, (node->children[i])->highBase.base)){
				Base lb = (node->children[i])->lowBase.base;
				char lbChar;
				if(lb==1)lbChar='A';else if(lb==2)lbChar='C';else if(lb==4)lbChar='G';else if(lb==8)lbChar='U';else lbChar="X";
	                        Base hb = (node->children[i])->highBase.base;
				char hbChar;
                                if(hb==1)hbChar='A';else if(hb==2)hbChar='C';else if(hb==4)hbChar='G';else if(hb==8)hbChar='U';else hbChar="X";
				printf("WARNING: bases %d and %d (%c%c) can't pair; structure cannot be scored. Exiting.\n",(node->children[i])->lowBase.index, (node->children[i])->highBase.index,lbChar, hbChar);
				exit(-1);
				//continue;
			}
			//printf("Yes they can pair\n");
			result += ScoreNode(node->children[i], RNA, param, length);
			numPairedChildren += 1;
			pairedChildren = realloc(pairedChildren, sizeof(int) * numPairedChildren);
			pairedChildren[numPairedChildren - 1] = i;
		}
	}
	
	if (node->lowBase.index != 0) 
	{
		
		if (numPairedChildren == 0)  // must be a hairpin
		{
			
			 int energy = eH(node->lowBase.index,node->highBase.index, RNA, param);
          result += energy;
			if(printOn2)printf("%d \t %d: Hairpin Loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
		}
		else if (numPairedChildren == 1)  // must be stack, bulge, or internal
		{
			if (node->numChildren == 1)  // must be stack 
			{
          	int energy = eS(node->lowBase.index, node->highBase.index, RNA, param);
	       	result += energy; 
	     		if(printOn2)printf("%d \t %d: Stacked pair with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
			}
			else 
			{  // must be bulge or internal 
			   
				int energy = eL(node->lowBase.index, node->highBase.index, 
				              node->children[pairedChildren[0]]->lowBase.index,
				              node->children[pairedChildren[0]]->highBase.index,
								  RNA, param);
            result += energy;
				if(printOn2)printf("%d \t %d: Bulge or Inernal Loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
			}
		}
		else  // must be a multi-loop
		{	
			int energy = eM(node, pairedChildren, numPairedChildren, RNA, param);
			result += energy;
			if(printOn2)printf("%d \t %d: Multi-loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
		}
	}
	else { // must be external
      int energy = eE(node, pairedChildren, numPairedChildren, RNA, param, length);
		result += energy; 
		if(printOn2)printf("%d \t %d: External loop with energy %.2f\n",  node->lowBase.index, node->highBase.index, (double)energy/100);
	}

	return result;
}


