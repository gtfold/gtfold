#ifndef _RANDOM_SAMPLE_H
#define _RANDOM_SAMPLE_H

#include <list>

using namespace std;

typedef struct sub_seq_t{
	int start;
	int end;
	int paired;
}sub_seq;

void multi_loop_strand(int i, int j, dangle_struct d_struct, int * structure, list<sub_seq> * stack);
#endif
