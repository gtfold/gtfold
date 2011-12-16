/*This program is used for testing single stranded constraints. Functions written in this file can be used for other scripts as well
Author: Manoj Soni
Email: manoj6891@gmail.com
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#include<fstream.h>
#define LENGTH 10
#define MAX_CONSTRS 10

void executeCommand(char* cmd){
	printf("Executing command: %s\n", cmd);
	int returnVal = system(cmd);
	printf("return value: %d\n", returnVal);
	if(returnVal != 0){
		printf("Error: returnVal:%d, Exiting...\n", returnVal);
		exit(-1);
	}
}

char* getToken(char* str, char* sep, int index){
  char * pch;
  pch = strtok (str,sep);
  char* fName = pch;
  if(index==0) return fName;
  if(index != -1){
	printf("Error: index can be only 0 or -1. 0 means first and -1 means last token\n");
	exit(-1);
  }
  while (pch != NULL)
  {
    //printf ("%s\n",pch);
    fName=pch;
    pch = strtok (NULL, sep);
  }	
  printf("%s\n",fName);
  return fName;
}

char* getSeqName(char* seqFileName1){
	char* seqFileName = malloc(200);
	strcpy(seqFileName,seqFileName1);
	char* seq = getToken(seqFileName, "/", -1);
        char* seqName = getToken(seq, ".", 0);
	printf("seqName: %s\n",seqName);
	return seqName;

}

void parseCtFileToFindForcedConstraintRegion(char* ctFileName, char* constrs[MAX_CONSTRS]){//printf("Entering parse method\n");
//ct file contents are
//dG = -5.8
//1       C       0       2       12      1
//2       C       1       3       11      2
//... ... ...
//then 1 is pairing with 12, 2 is pairing with 11 .. ?
//in other words, first column is i and 5th column is j for a pair (i,j)
  //char* constrs[MAX_CONSTRS];
  int index=0;
  FILE* ctFile = fopen(ctFileName, "r");
  char str[100];
  int arr[5];
  char c;
  fscanf(ctFile, "%[^\n]", str);printf("first line: %s\n",str);
  int start=-1, end=-1;
 while(feof(ctFile)==0){
   fscanf(ctFile, "%d %c %d %d %d %d\n", &arr[0],&c,&arr[1],&arr[2],&arr[3],&arr[4]);
   if(arr[3]!=0){
	if(start==-1){ start=arr[0]; end=arr[0];}
	else end=arr[0];
	if (start!=-1 && end-start+1 >= LENGTH){
		//char* ss_constr;
		char* ss_constr = (char*)malloc(30);//new char[30];
		sprintf(ss_constr, "%s%d%s%d", "P ", start, " 0 ", (end-start+1));
		constrs[index++] = ss_constr;
		printf("constraint is: %s\n",ss_constr);
		if(index>=MAX_CONSTRS) break;
		start=-1;
	}   
   }
  else{
	if(start==-1) continue;
	//	  char* ss_constr;
                char* ss_constr = malloc(30);//new char[30];
               sprintf(ss_constr, "%s%d%s%d", "P ", start, " 0 ", (end-start+1));
                constrs[index++] = ss_constr;
                printf("constraint is: %s\n",ss_constr);
                if(index>=MAX_CONSTRS) break;
                start=-1;

  }
  // printf("line is: %d %c %d %d %d %d\n", arr[0],c,arr[1],arr[2],arr[3],arr[4]);
  }
  if(index>=MAX_CONSTRS) return ;//constrs;
  //char* ss_constr;
  char* ss_constr = malloc(30);//new char[30];
  if(start!=-1){
	sprintf(ss_constr, "%s%d%s%d", "P ", start, " 0 ", (end-start+1));
	printf("constraint is: %s\n",ss_constr);
	constrs[index++] = ss_constr;
  }
  fclose(ctFile);
  return ;//ss_constr;
}

char* createConstraintFile(char* constrs[MAX_CONSTRS]){
	char* constrFileName = "constraints.txt";
	FILE* constrFile = fopen(constrFileName, "w+");
	int i;
	for(i=0; i<MAX_CONSTRS; ++i){
		if(constrs[i]==NULL)break;
		fprintf(constrFile, "%s\n", constrs[i]);
		//rewind(constrFile);
		printf("writing constraint: %s\n", constrs[i]);
	}
	fclose(constrFile);

}

int validateCtFileForSSconstraint(char* ctFileName, char* constrs[MAX_CONSTRS]){
	return 1;
}
int main(int argc, char* argv[]){
	//ar
	if(argc < 2 ){
		printf("Error: argc < 2, Exiting...\n");
		exit(-1);
	}
	char* seqFileName = argv[1];
	char cmd[200];
	printf("seqFileName: %s\n", seqFileName);
	printf("running the program gtfold on this seq without constraints\n");
	strcat(cmd, "./gtfold ");strcat(cmd, seqFileName);
	executeCommand(cmd);
	
 	char* seqName = getSeqName(seqFileName);	
	
	char* ctFileName = strcat(seqName, ".ct");
	printf("ctFileName: %s\n", ctFileName);
	char* constrs[MAX_CONSTRS];
	parseCtFileToFindForcedConstraintRegion(ctFileName, constrs);

	createConstraintFile(constrs);

	//now run same program with this constraint
	cmd[0] = '\0';
	strcat(cmd, "./gtfold -v -c constraints.txt ");strcat(cmd, seqFileName);
        executeCommand(cmd);	
	
	int passed = validateCtFileForSSconstraint(ctFileName, constrs);
	if(passed==1){
		printf("Single constraint test passed\n");
	}
}
