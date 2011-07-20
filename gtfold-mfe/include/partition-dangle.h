#ifndef _PARTITION_DANGLE_H
#define _PARTITION_DANGLE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct partition_d{
	double ** u;
	double ** up;
	double ** upm;
	double ** s1;
	double ** s2;
	double ** s3;
	double ** u1;
	int length;
}dangle_struct;

double cond_dangle(int j, int h, int l);
dangle_struct malloc_partition_arrays_d(int len);
void fill_partition_arrays_d(dangle_struct part_struct);

#ifdef __cplusplus
}
#endif

#endif
