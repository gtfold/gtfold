
#ifndef _ALGORITHMS_PARTITION_H
#define _ALGORITHMS_PARTITION_H

#ifdef __cplusplus
extern "C" {
#endif
/*
typedef struct _pFuncData {
    int len;
    double** QB;
    double** Q;
    double** QM;
} p_func_data;
*/
void fill_partition_fn_arrays(int len, double** QB, double** Q, double** QM);
void fillBasePairProbabilities(int length, double **Q, double **QB, double **QM, double**P);
void printBasePairProbabilities(int n, int *structure, double **P, const char* bppfile);
void printBasePairProbabilitiesDetail(int n, int *structure, double **P1, const char* outfile);
void printBasePairProbabilitiesComparison(int n, int *structure, double **P1, double** P2, const char* outfile);
double probabilityUnpaired(int length, int i, double **P);

double **mallocTwoD(int r, int c);
void freeTwoD(double** arr, int r, int c);



#ifdef __cplusplus
}
#endif

#endif
