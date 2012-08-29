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
int canStack(int i, int j);
int canSS(int i);
int canSSregion(int i, int j);
int canHairpin(int i, int j);
int canILoop(int i, int j, int p, int q);

int forceSS(int i);
int forceSSregion(int i, int j);
int forcePair(int i, int j);
int forcePaired(int i);

int withinCD(int i, int j);
int verify_structure();

void enable_constraints(int b);
void enable_limit_distance(int b);
void set_contact_distance(int dist);

#ifdef __cplusplus
}
#endif

#endif
