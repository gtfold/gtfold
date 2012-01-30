#include "StructureReader.h"
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

inline int canPair(int a, int b) {
	int result;
	//The possible base pairs are (A,U), (U,A), (C,G), (G,C), (G,U) and (U,G).
	//AU,UA
	if((a==1 && b==8) || (b==1 && a==8)) result=1;//return 1;
	//CG,GC
	else if((a==2 && b==4) || (b==2 && a==4)) result=1;//return 1;
	//GU,UG
	else if((a==4 && b==8) || (b==4 && a==8)) result=1;//return 1;
	else result=0;//return 0;
/* Base lb = a;
                                char lbChar;
                                if(lb==1)lbChar='A';else if(lb==2)lbChar='C';else if(lb==4)lbChar='G';else if(lb==8)lbChar='U';else lbChar="X";
                                Base hb = b;
                                char hbChar;
                                if(hb==1)hbChar='A';else if(hb==2)hbChar='C';else if(hb==4)hbChar='G';else if(hb==8)hbChar='U';else hbChar="X";
                                	
	printf("Inside canPair, a=%c, b=%c, result=%d\n",lbChar,hbChar,result);*/
	return result;
        //return ((1 << (((a)<<2) + (b))));
}


/*ZS: Converts to GT fold base code*/
		
int convertToGTfoldBase(int value){
	 if(value==1){
       return 0; 
	 }
	 else if(value==2){
		 return 1;
	 }
	 else if(value==4){
		 return 2;
	 }
	 else if(value==8){
		 return 3;
	 }
	 else{
       char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
		 fprintf(stderr, "Ambiguous base: %c\n", bases[value]);
	 	 return -1;
    }
}

/*ZS: Reads the length of the sequence from the .ct file. 
Assumes the first line of the file is 

sometabsandspaces thenumber sometext

If it isn't like that, it will break. 

Someone more competent than me should fix this someday, for now it works. 
*/
int readLength(FILE* filePtr){
   int junk = 0;
	
	char* junktext; 
	
	/* Read the first number. */
	int numRead = fscanf(filePtr, "%d%s\n", &junk, &junktext);
	/* Close file and Return the number */
	printf("The number read was %i \n", junk);
	fclose(filePtr);
	return junk;
}

/* Read data on a single base, filling baseData. Returns false if nothing could be read. Puts any pair
** index in pairIndex, or -1 if no pair.
*/
unsigned char ReadBase(FILE* filePtr, int index, BaseData* baseData, int* pairIndex, unsigned char isBPSEQ)
{
    while (1) /* Read as much as we have to to get the next base */
    {
	int junk = index - 1;
	int nextChar = 0;

	/* Check for a number. If not, ignore the line, or maybe the file is done. */
	int numRead = fscanf(filePtr, "%d", &junk);
	if (numRead != 1)
	{
	    /* Ignore this line, if it exists */
	    nextChar = fgetc(filePtr);
	    while (1)
	    {
		if (nextChar == EOF || nextChar == '\n')
		    break;
	    	nextChar = fgetc(filePtr);
	    }

	    /* If we are out of input, just return NULL */
	    if (nextChar == EOF)
	    {
		return 0;
	    }

	    /* Try the next line. */
	    continue;
	}

	/* Check for the next ID. If not, ignore the line and try the next one. */
	if (junk != index)
	{
	    while (1)
	    {
		nextChar = fgetc(filePtr);
		if (nextChar == EOF || nextChar == '\n')
		    break;
	    }
	    continue;
	}

	/* We're here if we have a valid index. */
	baseData->index = index;

	nextChar = fgetc(filePtr);
	while (isspace(nextChar))
	    nextChar = fgetc(filePtr);

	switch (nextChar)
	{
	    case 'a':
	    case 'A':
		baseData->base = A;
		break;

	    case 'c':
	    case 'C':
		baseData->base = C;
		break;

	    case 'g':
	    case 'G':
		baseData->base = G;
		break;

	    case 't':
	    case 'T':
	    case 'u':
	    case 'U':
		baseData->base = U;
		break;

/* NO SUPPORT FOR AMBIGUOUS BASES CURRENTLY */
/*
	    case 'm':
	    case 'M':
		baseData->base = M;
		break;
		
        case 'r':
	    case 'R':
		baseData->base = R;
		break;
		
        case 's':
	    case 'S':
		baseData->base = S;
		break;
		
		case 'v':
	    case 'V':
		baseData->base = V;
		break;
		
        case 'w':
	    case 'W':
		baseData->base = W;
		break;
		
        case 'y':
	    case 'Y':
		baseData->base = Y;
		break;
		
	    case 'h':
	    case 'H':
		baseData->base = H;
		break;
		
        case 'k':
	    case 'K':
		baseData->base = K;
		break;
		
        case 'd':
	    case 'D':
		baseData->base = D;
		break;
		
        case 'b':
	    case 'B':
		baseData->base = B;
		break;
		
        case 'n':
	    case 'N':
		baseData->base = N;
		break;

*/
	    default:
		fprintf(stderr, "Bad base: index %d\n", index);
		return 0;
	}

	nextChar = fgetc(filePtr);
	while (isspace(nextChar))
	    nextChar = fgetc(filePtr);
	ungetc(nextChar, filePtr);

	if (!isBPSEQ)
	{
	    if (!fscanf(filePtr, "%d", &junk))
	    {
		fprintf(stderr, "Bad prev id: id %d\n", index);
		return 0;
	    }

	    nextChar = fgetc(filePtr);
	    while (isspace(nextChar))
		nextChar = fgetc(filePtr);
	    ungetc(nextChar, filePtr);

	    if (!fscanf(filePtr, "%d", &junk))
	    {
		fprintf(stderr, "Bad next id: id %d\n", index);
		return 0;
	    }

	    nextChar = fgetc(filePtr);
	    while (isspace(nextChar))
		nextChar = fgetc(filePtr);
	    ungetc(nextChar, filePtr);
	}

	if (!fscanf(filePtr, "%d", pairIndex))
	{
	    fprintf(stderr, "Bad pair %d\n", index);
	    return 0;
	}

	if (!isBPSEQ)
	{
	    nextChar = fgetc(filePtr);
	    while (isspace(nextChar))
		nextChar = fgetc(filePtr);
	    ungetc(nextChar, filePtr);

	    if (!fscanf(filePtr, "%d", &junk))
	    {
		fprintf(stderr, "Bad trailing id: id %d\n", index);
		return 0;
	    }
	}

	return 1;
    }
}

TreeNode* CreateNode(
    int lowIndex,
    int* nextIndex,
    FILE* filePtr,
    unsigned char isBPSEQ,
    BaseData* returnData,
    int* returnPair,
	 int* RNA)
{
    if (!ReadBase(filePtr, lowIndex, returnData, returnPair, isBPSEQ))
    {
	return NULL;
    }

    TreeNode* result = (TreeNode*)malloc(sizeof(TreeNode));
    result->numChildren = 0;
    result->children = NULL;
    result->lowBase.index = returnData->index;
    result->lowBase.base = returnData->base;
	 RNA[returnData->index] = convertToGTfoldBase(returnData->base);
	 
    if (!(*returnPair))
    {
	  	result->isPair = 0;
		*nextIndex = lowIndex + 1;
		return result;
    }

    result->isPair = 1;
    result->highBase.index = *returnPair;
    *nextIndex = lowIndex + 1;
    TreeNode* child;
    do
    {
	child = CreateNode(*nextIndex, nextIndex, filePtr, isBPSEQ, returnData, returnPair, RNA);
	result->numChildren += 1;
	result->children = (TreeNode**)realloc(result->children, sizeof(TreeNode*) * result->numChildren);
	result->children[result->numChildren - 1] = child;
    } while (child && *nextIndex < result->highBase.index);

    if (*nextIndex == result->highBase.index)
    {
	ReadBase(filePtr, *nextIndex, returnData, returnPair, isBPSEQ);
	result->highBase.base = returnData->base;
	RNA[result->highBase.index] = convertToGTfoldBase(returnData->base);
	*nextIndex += 1;
    }

    return result;
}

ResultBundle* CreateFromFile(char* filename)
{
    // Figure out what kind of file we have and try to load it.
    char* extension = strrchr(filename, '.');
    unsigned char isBPSEQ = 0;
    if (extension && !strncmp(extension, ".bpseq", 6))
    {
	isBPSEQ = 1;
    }
    else if (extension && !strncmp(extension, ".ct", 3))
    {
	isBPSEQ = 0;
    }
    else
    {
	fprintf(stderr, "Unknown file type: %s\n", filename);
	return NULL;
    }
    FILE* filePtr = fopen(filename, "r");
    if (!filePtr)
    {fprintf(stderr, "Unable to open file: %s\n", filename);
	  return NULL;}

	 //ZS: Open a second copy briefly, which we only use to read the length 
	 //of the sequence 
	 //(I know it's an ugly solution but I don't want to mess up Steven's code
	 //and I don't really understand it) 
	 FILE* filePtr_second = fopen(filename, "r");
	 int value = readLength(filePtr_second); 
	 printf("Length of sequence: %d\n", value);

	 int *RNA = (int *)malloc(sizeof(int)*(value+1)); 
	 //ZS: In GTfold, RNA[1] is the first base, so complying with that here
 
    TreeNode* result = (TreeNode*)malloc(sizeof(TreeNode));
    result->numChildren = 0;
    result->children = NULL;
    result->lowBase.index = 0;
    result->lowBase.base = 0;
    result->highBase.index = 0;
    result->highBase.base = 0;
    result->isPair = 0;
    int nextIndex = 1;
    int pairIndex;
    BaseData bData;
    
    TreeNode* child = CreateNode(nextIndex, &nextIndex, filePtr, isBPSEQ, &bData, &pairIndex, RNA);
    
    while (child)
    {
	result->numChildren += 1;
	result->children = (TreeNode**)realloc(result->children, sizeof(TreeNode*) * result->numChildren);
	result->children[result->numChildren - 1] = child;
	child = CreateNode(nextIndex, &nextIndex, filePtr, isBPSEQ, &bData, &pairIndex, RNA);
    }

    fclose(filePtr);

    if (!result)
    {
	fprintf(stderr, "Empty or malformed file: %s\n", filename);
    }
    
   ResultBundle* resultBundle = (ResultBundle*)malloc(sizeof(ResultBundle));
   resultBundle->RNA_seq = RNA;
	resultBundle->length = value;
   resultBundle->treenode = result;
   
	return resultBundle;
}



#define PrintTabs(indent) for (i = 0; i < indent; ++i) printf("   ");

void PrintTree(TreeNode* tree, int indent)
{
    char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    int i;

    PrintTabs(indent);
    if (tree->lowBase.index == 0)
    {
	printf("Root\n");
    }
    else if (tree->isPair)
    {
	printf("(%d %c - %c %d)\n",
	       tree->lowBase.index, bases[tree->lowBase.base],
	       bases[tree->highBase.base], tree->highBase.index);
    }
    else
    {
	printf("%d %c\n", tree->lowBase.index, bases[tree->lowBase.base]);
    }

    for (i = 0 ; i < tree->numChildren; ++i)
    {
	PrintTree(tree->children[i], indent + 1);
    }
}
 
