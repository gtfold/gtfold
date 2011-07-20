#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "constants.h"

extern unsigned char *RNA; 
extern int *structure; 
extern int* constraints;

extern int g_nthreads;
extern int g_unamode;
extern int g_dangles;
extern int g_mismatch;
extern int g_verbose;
extern int g_prefilter_mode;
extern int g_prefilter1;
extern int g_prefilter2;

extern unsigned int chPairKey;

// The possible base pairs are (A,U), (U,A), (C,G), (G,C), (G,U) and (U,G). 
#define checkPair(i, j) (((((i)-(j)) % 2) == 1 || (((i)-(j)) % 2)== -1) && (!( ((i)==BASE_A && (j)==BASE_C) || ((i)==BASE_C && (j)==BASE_A) )))


#ifdef __cplusplus
extern "C" {
#endif
int canPair(int a, int b);
void init_global_params(int len);
void free_global_params();
void print_sequence(int len); 
void print_structure(int len); 
#ifdef __cplusplus
}
#endif

void init_checkPair(); 
int  update_checkPair(int i, int j);

#endif
