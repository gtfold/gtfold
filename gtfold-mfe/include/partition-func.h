#ifndef _PARTITION_DANGLE_H
#define _PARTITION_DANGLE_H


#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

extern	double ** u;
extern	double ** up;
extern	double ** upm;
extern	double ** ud;
extern	double ** u1d;

extern	double ** s1;
extern	double ** s2;
extern	double ** s3;
extern	double ** u1;

extern int part_len;
double ED3_new(int i, int j, int k);
double ED5_new(int i, int j, int k);
double EA_new();
double EB_new();
double EC_new();
double eS_new(int i, int j);
double eL_new(int i, int j, int p, int q);
double eH_new(int i, int j);
double auPenalty_new(int i, int j);
double f(int j, int h, int l);
double calculate_partition(int len, int pf_count_mode, int no_dangle_mode);
void free_partition();


#ifdef __cplusplus
}
#endif

#endif
