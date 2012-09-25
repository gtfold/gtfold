#include<iostream>
#include<fstream>
#include <stack>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stack>
#include <map>
#include<sstream>
static int verbose=0;
using namespace std;

int* getStructureArray(const char* ensemble, int length){
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
                                printf("%s structure is not a valid structure in dot-bracket notation, exiting...\n\n", ensemble);
                                exit(-1);
                        }
                        int openingBracketIndexForI = openBracketStack.top();
                        openBracketStack.pop();
                        structure[i+1]=openingBracketIndexForI+1;
                        structure[openingBracketIndexForI+1]=i+1;
                }
                i++;
        }
        if(!openBracketStack.empty()){
                printf("%s structure is not a valid structure in dot-bracket notation, exiting...\n\n", ensemble);
                exit(-1);
        }

        return structure;
}


void updatePnum(int * strucArr, int** pnumArr, int* pnumLen, int length){
	for(int i=1; i<=length; ++i){
		int j1=strucArr[i];
		if(j1==-1) continue;
		int isAreadyExists = 0;
		for(int j=0; j<pnumLen[i]; ++j){
			if(pnumArr[i][j]==j1){
				isAreadyExists = 1;
				break;
			}
		}
		if(isAreadyExists==0){
			pnumArr[i][ pnumLen[i] ] = j1;
			pnumLen[i]++;
		}
	}
}	

int main(int argc, char** argv){
	if(argc<2){
		cout<<"Usage: ./calc_pnum file_ss.txt\n";
		exit(-1);
	}
	char* filename = argv[1];
	ifstream fin(filename);
	string seq;
	string struc;
	double energy;
	fin>>seq;fin>>energy;
	if(verbose==1)cout<<"seq="<<seq<<endl;
	int length = seq.length();
	if(verbose==1)cout<<"length="<<length<<endl;
	int** pnumArr = new int*[length+1];
	for(int i=1; i<=length; ++i) pnumArr[i] = new int[length];
	int* pnumLen = new int[length+1];
	for(int i=1; i<=length; ++i) pnumLen[i]=0;
	int numStrucs = 0;
	while(fin>>struc){
		if(verbose==1)cout<<"struc="<<struc<<endl;
		int* strucArr = getStructureArray(struc.c_str(), struc.length());
		if(verbose==1){
			for(int i=1; i<=struc.length(); ++i)cout<<i<<" "<<strucArr[i]<<"\n";
		}
		updatePnum(strucArr, pnumArr, pnumLen, length);	
		if(verbose==1)cout<<"\n\n";
		fin>>energy;
		numStrucs++;
	}
	if(numStrucs==0){
		cout<<"zero structures\n";
		return 0;
	}
	if(verbose==1)cout<<"Success\n";
	for(int i=1; i<=length; ++i){
		cout<<i<<" "<<(double)pnumLen[i]/numStrucs<<endl;
	}
	
	fin.close();
	return 0;
}
