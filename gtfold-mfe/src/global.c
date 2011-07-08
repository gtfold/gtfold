#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "constraints.h"

unsigned char *RNA; 
int *structure; 
unsigned int chPairKey;

int g_nthreads;
int g_dangles;
int g_unamode;
int g_mismatch;
int g_verbose;
int g_prefilter_mode;
int g_prefilter1;
int g_prefilter2;

void init_global_params(int len) {
	RNA = (unsigned char *) malloc((len+1)* sizeof(unsigned char));
	if (RNA == NULL) {
		perror("Cannot allocate variable 'RNA'");
		exit(-1);
	}
	structure = (int *) malloc((len+1) * sizeof(int));
	if (structure == NULL) {
		perror("Cannot allocate variable 'structure'");
		exit(-1);
	}

	init_checkPair();
}

void free_global_params() {
	free(structure);
	free(RNA);
}

void print_sequence(int len) {
	int i;
	for (i = 1; i <= len; i++) {
		if (RNA[i] == BASE_A)
			printf("A");
		else if (RNA[i] == BASE_C)
			printf("C");
		else if (RNA[i] == BASE_G)
			printf("G");
		else if (RNA[i] == BASE_U)
			printf("U");
		else
			printf("N");
	}
	printf("\n");
}

void print_structure(int len) {
	int i = 1;
	for (i = 1; i <= len; i++) 
	{
		if (structure[i] > 0 && structure[i] > i)
			printf("(");
		else if (structure[i] > 0 && structure[i] < i)
			printf(")");
		else
			printf(".");
	}
	printf("\n");
}

void init_checkPair() {
	int i, j;
	chPairKey = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			chPairKey += checkPair(i, j) << ((i << 2) + j);
}

int update_checkPair(int i, int j) {
	int r = 0;
	if (!((i >= 0 && i <=3 )&&(j >=0 && j <=3)))
		return r;
	if (!(chPairKey & (1 << ((i << 2) + j)))) {
		chPairKey += 1 << ((i << 2) + j);	
		r = 1;
	}
	return r;
}

inline int canPair(int a, int b) {
	return (chPairKey & (1 << (((a)<<2) + (b))));
}
