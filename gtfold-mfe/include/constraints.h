#ifndef _CONSTRAINTS_H_
#define _CONSTRAINTS_H_

extern int* BP;
extern int** PBP;
extern int** FBP;
extern int* ind;

extern int nPBP;
extern int nFBP;

//static int load_constraints(const char* constr_file, int verbose=0);

#define BP(i,j) BP[ind[j]+i]

int init_constraints(const char* constr_file, int length) ;
void free_constraints(int length) ;
void print_constraints(int length) ;


#ifdef __cplusplus
extern "C" {
#endif
//int is_ss(int i, int j); 
//int prohibit_base(int i) ;
int check_ssregion(int i, int j);
//int check_base(int i) ;
//int force_pair(int i, int j) ;
int force_pair1(int i, int j) ;
int force_ss1(int i);
int force_ssregion1(int i, int j);
int check_iloop(int i, int j, int p, int q) ;
int check_pair(int i, int j) ;
int check_stack(int i, int j) ;
int check_hairpin(int i, int j) ;
int can_dangle(int i);
int withinCD(int i, int j);
#ifdef __cplusplus
}
#endif

#endif
