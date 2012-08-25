#include<iostream>
#include<fstream>
#include<sstream>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <stack>
using namespace std;

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
        char* seqFileName = (char*)malloc(200);
        strcpy(seqFileName,seqFileName1);
        char* seq = getToken(seqFileName, "/", -1);
        char* seqName = getToken(seq, ".", 0);
        printf("seqName: %s\n",seqName);
        return seqName;

}

int* getStructure(const char* ensemble, int length){
        std::stack<int> openBracketStack;
        int* structure = new int[length+1];
        for(int j=0; j<=length; ++j) structure[j]=-1;
        int i=0;
        while(i<length){
                if(ensemble[i]=='('){
                        openBracketStack.push(i);
                }
                else if(ensemble[i]==')'){
                        if(openBracketStack.empty()){
                                printf("%s structure is not a valid structure in dot-bracket notation, returning empty structures\n\n", ensemble);
				for(int i=1; i<=length; ++i) structure[i]=-1;
		                return structure;
                        }
                        int openingBracketIndexForI = openBracketStack.top();
                        openBracketStack.pop();
                        structure[i+1]=openingBracketIndexForI+1;
                        structure[openingBracketIndexForI+1]=i+1;
                }
                i++;
        }
        if(!openBracketStack.empty()){
                printf("%s structure is not a valid structure in dot-bracket notation, returning empty structure\n\n", ensemble);
                for(int i=1; i<=length; ++i) structure[i]=-1;
		return structure;
        }

        return structure;
}

bool validateLCDforStruc(int* struc, int length, int lcd){
	bool isValid = true;
	for(int i=1; i<=length; ++i){
		int j=struc[i];
		if(j<0) continue;
		if(i>j) continue;
		if(j-i>lcd){ isValid=false; break;}
	}
	return isValid;
}

int main(int argc, char** argv){
	char cmd[200];
	cmd[0]='\0';
	if(argc<4){
		cout<<"Usage: ./a.out numsamples lcd <seqFile>\n";
		exit(-1);
	}
	int numsamples = atoi(argv[1]);
	int lcd = atoi(argv[2]);
        char* seqfile = argv[3];
	strcat(cmd, "./gtboltzmann -t 1 -l ");strcat(cmd, argv[2]); strcat(cmd, " -s ");strcat(cmd, argv[1]); strcat(cmd, " ");strcat(cmd, seqfile);
        executeCommand(cmd);

	char* seqName = getSeqName(seqfile);
        char* samplesfile = strcat(seqName, ".samples");
	bool pass = true;
	ifstream fin(samplesfile);
	for(int i=0; i<numsamples; ++i){
		//((((....))))    -6.4    1 12 4,
		string line;
		getline(fin, line, '\n');
		stringstream ss(line);
		string dotBracketStructure;
		ss>>dotBracketStructure;
		//cout<<"dotBracketStructure: "<<dotBracketStructure<<endl;
		int length = dotBracketStructure.length();
		int* struc = getStructure(dotBracketStructure.c_str(), length);
		bool isValid = validateLCDforStruc(struc, length, lcd);
		if(!isValid){
			cout<<"sampled structure "<<dotBracketStructure<<" is not valid for lcd "<<lcd<<endl;
			pass=false;
		}
	}
	fin.close();
	if(pass) cout<<"Test completed successfully.\n";
	else cout<<"Test did not complete successfully\n";
	return 0;
}
