#ifndef _PARTITION_DANGLE_H
#define _PARTITION_DANGLE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct partition_struct {
	double ** u;
	double ** up;
	double ** upm;
	double ** ud;
	double ** u1d;

	double ** s1;
	double ** s2;
	double ** s3;
	double ** u1;

	int length;
} partition_t;

extern partition_t partition;

void calculate_partition();

#ifdef __cplusplus
}
#endif

#endif
