/*This program is used for testing single stranded constraints. Functions written in this file can be used for other scripts as well
Author: Manoj Soni
Email: manoj6891@gmail.com
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
//#include<fstream.h>
int LENGTH = 10;//default value
int CHANCE = 50;
char* seqFileName;
int VERBOSE=0;
#define MAX_CONSTRS 750

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
				if(VERBOSE==1)printf("constraint is: %s\n",ss_constr);
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
			if(VERBOSE==1)printf("constraint is: %s\n",ss_constr);
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
		if(VERBOSE==1)printf("constraint is: %s\n",ss_constr);
		constrs[index++] = ss_constr;
	}
	if(index<MAX_CONSTRS) constrs[index] = NULL;

	fclose(ctFile);
	return ;//ss_constr;
}

char* createConstraintFile(char* constrs[MAX_CONSTRS]){
	char* constrFileName = "constraints.txt";
	FILE* constrFile = fopen(constrFileName, "w+");
	int i;

	for(i=0; i<MAX_CONSTRS; ++i){
		if(constrs[i]==NULL)break;
		int chance = rand()%100;
		if(chance>CHANCE) continue;	
		if(VERBOSE==1)printf("writing constraint: %s\n", constrs[i]);
		fprintf(constrFile, "%s\n", constrs[i]);
		//rewind(constrFile);
	}
	fclose(constrFile);

}

int validateCtFileForSSconstraint(char* ctFileName, char* constrs[MAX_CONSTRS]){
	return 1;
}

void help(){
	printf("\nUsage: test_constraints [OPTION]... FILE\n\n");
	printf("   FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");
	printf("OPTIONS\n");
	printf("   --chance INT	int value between 1 and 100, to be considered as chance or percentage probability for considering any generated constraint (default=50).\n");
	printf("   --ssmaxlen INT	Restricts maximum length for any single stranded constraint region (default=10).\n");
	printf("   -h, --help	Output help (this message) and exit.\n");
	printf("\n");
	exit(-1);
}

void parse_options(int argc, char* argv[]){
	int i;
	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-') {
			if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
				help();
			} else if(strcmp(argv[i], "--chance") == 0 ) {
				if(i < argc)
					CHANCE = atoi(argv[++i]);
				else
					help();	
			} else if(strcmp(argv[i], "--ssmaxlen") == 0 ) {
				if(i < argc)
					LENGTH = atoi(argv[++i]);
				else
					help();
			} else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
				VERBOSE = 1;
			}	
		}
		else {
			seqFileName = argv[i];
		}

	}
	if(seqFileName==NULL) {
		help();
		printf("Missing input file.\n");
	}
}

int main(int argc, char* argv[]){
	parse_options(argc, argv);

	srand(time(NULL));
	char cmd[200];
	printf("seqFileName: %s\n", seqFileName);
	printf("running the program gtmfe on this seq without constraints\n");
	strcat(cmd, "./gtmfe ");strcat(cmd, seqFileName);
	executeCommand(cmd);

	char* seqName = getSeqName(seqFileName);	

	char* ctFileName = strcat(seqName, ".ct");
	printf("ctFileName: %s\n", ctFileName);
	char* constrs[MAX_CONSTRS];
	parseCtFileToFindForcedConstraintRegion(ctFileName, constrs);
	printf("completed creating constraints and now creating constraints file\n");
	createConstraintFile(constrs);
	printf("completed creating constraints file and now executing gtmfe command with constraints\n");

	//now run same program with this constraint
	cmd[0] = '\0';
	strcat(cmd, "./gtmfe -v -c constraints.txt ");strcat(cmd, seqFileName);
	executeCommand(cmd);	

	int passed = validateCtFileForSSconstraint(ctFileName, constrs);
	if(passed==1){
		printf("Single constraint test passed\n");
	}
}
