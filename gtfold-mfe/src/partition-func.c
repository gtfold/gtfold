#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "partition-func.h"
#include "energy.h"
#include "algorithms-partition.h"
#include "global.h"
#include "utils.h"
#include<omp.h>

double ** u;
double ** up;
double ** upm;
double ** ud;
double ** u1d;
double ** s1;
double ** s2;
double ** s3;
double ** u1;
int part_len;
int PF_COUNT_MODE_;
static int NO_DANGLE_MODE_=0;
static void create_partition_arrays();
static void init_partition_arrays();
static void fill_partition_arrays();
static void free_partition_arrays();
static void calc_u(int i, int j);
static void calc_ud(int i, int j);
static void calc_up(int i, int j);
void calc_up_parallel(int i, int j);
static void calc_upm(int i, int j);
static void calc_u1(int i, int j);
static void calc_u1d(int i, int j);
static void calc_s1(int i, int j);
static void calc_s2(int i, int j);
static void calc_s3(int i, int j);
static double get_u(int i, int j);
static double get_ud(int i, int j);
static double get_up(int i, int j);
static double get_upm(int i, int j);
static double get_u1(int i, int j);
static double get_u1d(int i, int j);
static double get_s1(int i, int j);
static double get_s2(int i, int j);
static double get_s3(int i, int j);
static void set_u(int i, int j, double val);
static void set_ud(int i, int j, double val);
static void set_up(int i, int j, double val);
static void set_upm(int i, int j, double val);
static void set_u1(int i, int j, double val);
static void set_u1d(int i, int j, double val);
static void set_s1(int i, int j, double val);
static void set_s2(int i, int j, double val);
static void set_s3(int i, int j, double val);

void errorAndExit(char* msg, int i, int j, double oldVal, double newVal){
	printf("%s\n", msg);
	printf("i=%d,j=%d,oldVal=%0.1f,newVal=%0.1f\n",i,j,oldVal,newVal);
	printf("%s","\nprogram is exiting now due to above error\n");
	exit(-1);
}
//OPTIMIZED CODE STARTS
/*
double get_u(int i, int j) {if(u[i][j]==-1) errorAndExit("get_u entry is -1.",i,j,-1,0); return u[i][j];}
double get_ud(int i, int j) {if(ud[i][j]==-1) errorAndExit("get_ud entry is -1.\n",i,j,-1,0); return ud[i][j];}
double get_up(int i, int j) {if(up[i][j]==-1) errorAndExit("get_up entry is -1.\n",i,j,-1,0); return up[i][j];}
double get_upm(int i, int j) {if(upm[i][j]==-1) errorAndExit("get_upm entry is -1.\n",i,j,-1,0); return upm[i][j];}
double get_u1(int i, int j) {if(u1[i][j]==-1) errorAndExit("get_u1 entry is -1.\n",i,j,-1,0); return u1[i][j];}
double get_u1d(int i, int j) {if(u1d[i][j]==-1) errorAndExit("get_u1d entry is -1.\n",i,j,-1,0); return u1d[i][j];}
double get_s1(int i, int j) {if(s1[i][j]==-1) errorAndExit("get_s1 entry is -1.\n",i,j,-1,0); return s1[i][j];}
double get_s2(int i, int j) {if(s2[i][j]==-1) errorAndExit("get_s2 entry is -1.\n",i,j,-1,0); return s2[i][j];}
double get_s3(int i, int j) {if(s3[i][j]==-1) errorAndExit("get_s3 entry is -1.\n",i,j,-1,0); return s3[i][j];}
void set_u(int i, int j, double val) {if(u[i][j]!=-1 && u[i][j]!=val) errorAndExit("set_u entry is not -1.\n",i,j,u[i][j],val); u[i][j]=val;}
void set_ud(int i, int j, double val) {if(ud[i][j]!=-1 && ud[i][j]!=val) errorAndExit("set_ud entry is not -1.\n",i,j,ud[i][j],val); ud[i][j]=val;}
void set_up(int i, int j, double val) {if(up[i][j]!=-1 && up[i][j]!=val) errorAndExit("set_up entry is not -1.\n",i,j,up[i][j],val); up[i][j]=val;}
void set_upm(int i, int j, double val) {if(upm[i][j]!=-1 && upm[i][j]!=val) errorAndExit("set_upm entry is not -1.\n",i,j,upm[i][j],val); upm[i][j]=val;}
void set_u1(int i, int j, double val) {if(u1[i][j]!=-1 && u1[i][j]!=val) errorAndExit("set_u1 entry is not -1.\n",i,j,u1[i][j],val); u1[i][j]=val;}
void set_u1d(int i, int j, double val) {if(u1d[i][j]!=-1 && u1d[i][j]!=val) errorAndExit("set_u1d entry is not -1.\n",i,j,u1d[i][j],val); u1d[i][j]=val;}
void set_s1(int i, int j, double val) {if(s1[i][j]!=-1 && s1[i][j]!=val) errorAndExit("set_s1 entry is not -1.\n",i,j,s1[i][j],val); s1[i][j]=val;}
void set_s2(int i, int j, double val) {if(s2[i][j]!=-1 && s2[i][j]!=val) errorAndExit("set_s2 entry is not -1.\n",i,j,s2[i][j],val); s2[i][j]=val;}
void set_s3(int i, int j, double val) {if(s3[i][j]!=-1 && s3[i][j]!=val) errorAndExit("set_s3 entry is not -1.\n",i,j,s3[i][j],val); s3[i][j]=val;}
*/
/*int u_g_printOn=1;
int ud_g_printOn=1;
int up_g_printOn=1;
int upm_g_printOn=1;
int u1_g_printOn=1;
int u1d_g_printOn=1;
int s1_g_printOn=1;
int s2_g_printOn=1;
int s3_g_printOn=1;*/
inline double get_u(int i, int j) {return u[i][j];}
inline double get_ud(int i, int j) {return ud[i][j];}
inline double get_up(int i, int j) {return up[i][j];}
inline double get_upm(int i, int j) {return upm[i][j];}
inline double get_u1(int i, int j) {return u1[i][j];}
inline double get_u1d(int i, int j) {return u1d[i][j];}
inline double get_s1(int i, int j) {return s1[i][j];}
inline double get_s2(int i, int j) {return s2[i][j];}
inline double get_s3(int i, int j) {return s3[i][j];}
inline void set_u(int i, int j, double val) {u[i][j]=val;}
inline void set_ud(int i, int j, double val) {ud[i][j]=val;}
inline void set_up(int i, int j, double val) {up[i][j]=val;}
inline void set_upm(int i, int j, double val) {upm[i][j]=val;}
inline void set_u1(int i, int j, double val) {u1[i][j]=val;}
inline void set_u1d(int i, int j, double val) {u1d[i][j]=val;}
inline void set_s1(int i, int j, double val) {s1[i][j]=val;}
inline void set_s2(int i, int j, double val) {s2[i][j]=val;}
inline void set_s3(int i, int j, double val) {s3[i][j]=val;}
/*
inline double get_u(int i, int j) {if(u_g_printOn) printf("get_u entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return u[i][j];}
inline double get_ud(int i, int j) {if(ud_g_printOn) printf("get_ud entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return ud[i][j];}
inline double get_up(int i, int j) {if(up_g_printOn) printf("get_up entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return up[i][j];}
inline double get_upm(int i, int j) {if(upm_g_printOn) printf("get_upm entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return upm[i][j];}
inline double get_u1(int i, int j) {if(u1_g_printOn) printf("get_u1 entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return u1[i][j];}
inline double get_u1d(int i, int j) {if(u1d_g_printOn) printf("get_u1d entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return u1d[i][j];}
inline double get_s1(int i, int j) {if(s1_g_printOn) printf("get_s1 entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return s1[i][j];}
inline double get_s2(int i, int j) {if(s2_g_printOn) printf("get_s2 entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return s2[i][j];}
inline double get_s3(int i, int j) {if(s3_g_printOn) printf("get_s3 entry is accessed,i=%d,j=%d,index=%d\n",i,j,i*part_len+j); return s3[i][j];}
inline void set_u(int i, int j, double val) {u[i][j]=val;}
inline void set_ud(int i, int j, double val) {ud[i][j]=val;}
inline void set_up(int i, int j, double val) {up[i][j]=val;}
inline void set_upm(int i, int j, double val) {upm[i][j]=val;}
inline void set_u1(int i, int j, double val) {u1[i][j]=val;}
inline void set_u1d(int i, int j, double val) {u1d[i][j]=val;}
inline void set_s1(int i, int j, double val) {s1[i][j]=val;}
inline void set_s2(int i, int j, double val) {s2[i][j]=val;}
inline void set_s3(int i, int j, double val) {s3[i][j]=val;}
*/
//OPTIMIZED CODE ENDS
inline double myExp(double arg){
	double posArgMultRT = (-1)*arg*RT;
	if(posArgMultRT>=INFINITY_){
		return 0;
	}
	return exp(arg);
}
inline double eS_new(int i, int j){
	if(PF_COUNT_MODE_) return 0;
	return eS(i,j);
	//return eS(i,j)/100;
}
inline double eH_new(int i, int j){
	if(PF_COUNT_MODE_) return 0;
	return eH(i,j);
	//return eH(i,j)/100;
}
inline double eL_new(int i, int j, int p, int q){
	if(PF_COUNT_MODE_) return 0;
	return eL(i,j,p,q);
	//return eL(i,j,p,q)/100;
}
inline double ED3_new(int i, int j, int k){
	if(NO_DANGLE_MODE_) return 0;
	if(PF_COUNT_MODE_) return 0;
	return Ed5(j,i,k);
	//return Ed5(j,i,k)/100;
}
inline double ED5_new(int i, int j, int k){
	if(NO_DANGLE_MODE_) return 0;
	if(PF_COUNT_MODE_) return 0;
	if (k<1) return 0;
	return Ed3(j,i,k);
	//return Ed3(j,i,k)/100;
}
inline double EA_new(){
	if(PF_COUNT_MODE_) return 0;
	return Ea;
	//return Ea/100;
}
inline double EB_new(){
	if(PF_COUNT_MODE_) return 0;
	return Ec;
	//return Ec/100;
}
inline double EC_new(){
	if(PF_COUNT_MODE_) return 0;
	return Eb;
	//return Eb/100;
}
inline double auPenalty_new(int i, int j){
	if(PF_COUNT_MODE_) return 0;
	return auPenalty(i,j);
	//return auPenalty(i,j)/100;
}
inline 	double f(int j, int h, int l){
	if(j - 1 == l)
		return 1;
	else
		return myExp(-ED3_new(h,l,l+1)/RT);
}
void printMatrix(double** u, int part_len){
	int i,j;
	for (i = 0; i <= part_len+1; ++i)
	{
		for (j = 0; j <= part_len+1; ++j)
			printf("%0.1f ",u[i][j]);
		printf("\n");
	}
}
void printAllMatrixes(){
	printf("\n\nAfter calculation, u matrix:\n\n");
	printMatrix(u,part_len);
	printf("\n\nAfter calculation, ud matrix:\n\n");
	printMatrix(ud,part_len);
	printf("\n\nAfter calculation, up matrix:\n\n");
	printMatrix(up,part_len);
	printf("\n\nAfter calculation, upm matrix:\n\n");
	printMatrix(upm,part_len);
	printf("\n\nAfter calculation, u1 matrix:\n\n");
	printMatrix(u1,part_len);
	printf("\n\nAfter calculation, u1d matrix:\n\n");
	printMatrix(u1d,part_len);
	printf("\n\nAfter calculation, s1 matrix:\n\n");
	printMatrix(s1,part_len);
	printf("\n\nAfter calculation, s2 matrix:\n\n");
	printMatrix(s2,part_len);
	printf("\n\nAfter calculation, s3 matrix:\n\n");
	printMatrix(s3,part_len);
}
double calculate_partition(int len, int pf_count_mode, int no_dangle_mode)
{
	PF_COUNT_MODE_ = pf_count_mode;
	NO_DANGLE_MODE_ = no_dangle_mode;
	part_len = len;
	//OPTIMIZED CODE STARTS
        #ifdef _OPENMP
        if (g_nthreads > 0) omp_set_num_threads(g_nthreads);
        #endif

        #ifdef _OPENMP
        #pragma omp parallel
        #pragma omp master
        fprintf(stdout,"Thread count: %3d \n",omp_get_num_threads());
	#endif
	//OPTIMIZED CODE ENDSS

	create_partition_arrays();
	init_partition_arrays();
	fill_partition_arrays();
	//printAllMatrixes();//TODO uncomment it
	printf("Partition Function value is: ");printf("%4.4f\n",u[1][part_len]);printf("\n");
	return u[1][part_len];

}
void free_partition()
{
	free_partition_arrays();
}
void init_part_arrays_zeros(){
	int i,j,n;
	n = part_len+1;
	for(i=0; i<=n; ++i){
		for(j=0; j<=n; ++j){
			u[i][j]=0;
			up[i][j]=0;
			upm[i][j]=0;
			ud[i][j]=0;
			u1d[i][j]=0;
			s1[i][j]=0;
			s2[i][j]=0;
			s3[i][j]=0;
			u1[i][j]=0;
		}
	}
}
void init_part_arrays_ones(){
	int i,j,n;
	n = part_len+1;
	for(i=0; i<=n; ++i){
		for(j=0; j<=n; ++j){
			u[i][j]=1;
			up[i][j]=1;
			upm[i][j]=1;
			ud[i][j]=1;
			u1d[i][j]=1;
			s1[i][j]=1;
			s2[i][j]=1;
			s3[i][j]=1;
			u1[i][j]=1;
		}
	}
}
void init_part_arrays_negatives(){
	int i,j,n;
	n = part_len+1;
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i,j) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	for(i=0; i<=n; ++i){
		for(j=0; j<=n; ++j){
			u[i][j]=-1;
			up[i][j]=-1;
			upm[i][j]=-1;
			ud[i][j]=-1;
			u1d[i][j]=-1;
			s1[i][j]=-1;
			s2[i][j]=-1;
			s3[i][j]=-1;
			u1[i][j]=-1;
		}
	}
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i,j) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	for(i=0; i<=n+1; ++i){
		for(j=0; j<=n+1; ++j){
			u1[i][j]=-1;
		}
	}
}
void init_partition_arrays()
{  init_part_arrays_negatives();
	int i, j;
	int n = part_len;
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i,j) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	for(i=1; i<=n; ++i){
		for(j=i; j<=i+TURN && j<=n; ++j){
			u[i][j] = 1;
			up[i][j] = 0;
			ud[i][j] = 0;
			u1[i][j] = 0;
			u1d[i][j] =0 ;
			s1[i][j] = 0;
			s2[i][j] = 0;
			s3[i][j] = 0;
		}
	}
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	for(i=1; i<=n; ++i){
		u[i+1][i] = 1;
		u1[i+1][i] = 0;
		u1d[i+1][i] = 0;
	}
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	for(i=1; i<=n; i++){
		u1[i+2][i] = 0;
	}
}
void fill_partition_arrays()
{
	int b,i,j;
	int n=part_len;
	int numThds = omp_get_num_threads();
	int* i_canPair = (int*)malloc((n+1)*sizeof(int));
	int* i_cannotPair = (int*)malloc((n+1)*sizeof(int));
	int len_i_canPair=0, len_i_cannotPair=0;
	int index1=0;
	int b_threshold = n-numThds/2;//TURN+1;//n-numThds/2;//n;
	for(b=TURN+1; b<b_threshold; ++b){
		//OPTIMIZED CODE STARTS
		len_i_canPair=0;
	       	len_i_cannotPair=0;
		for(i=1; i<=n-b; ++i){
			j=i+b;
			if(canPair(RNA[i],RNA[j])){i_canPair[len_i_canPair++]=i;}
			else {i_cannotPair[len_i_cannotPair++]=i;}
		}
		#ifdef _OPENMP
		#pragma omp parallel for private (index1,i,j) schedule(guided)
		#endif
		for(index1=0; index1<len_i_canPair; ++index1){
			i=i_canPair[index1];
			j=i+b;
			calc_s1(i,j);
			calc_s2(i,j);
			calc_upm(i,j);
			calc_up(i,j);
			calc_s3(i,j);
			calc_u1d(i,j);
			calc_u1(i,j);
			calc_ud(i,j);
			calc_u(i,j);
		}
		#ifdef _OPENMP
		#pragma omp parallel for private (index1,i,j) schedule(guided)
		#endif
		for(index1=0; index1<len_i_cannotPair; ++index1){
			i=i_cannotPair[index1];
			j=i+b;
			calc_s1(i,j);
			calc_s2(i,j);
			calc_upm(i,j);
			calc_up(i,j);
			calc_s3(i,j);
			calc_u1d(i,j);
			calc_u1(i,j);
			calc_ud(i,j);
			calc_u(i,j);
		}
		//OPTIMIZED CODE ENDS
	}
	for(b=b_threshold; b<n; ++b){
		for(i=1; i<=n-b; ++i){
			j=i+b;
			calc_s1(i,j);
			calc_s2(i,j);
			calc_upm(i,j);
			calc_up_parallel(i,j);
			calc_s3(i,j);
			calc_u1d(i,j);
			calc_u1(i,j);
			calc_ud(i,j);
			calc_u(i,j);
		}
	}
	/*//OLD Implementation
	for(b=TURN+1; b<n; ++b){
		#ifdef _OPENMP
		#pragma omp parallel for private (i,j) schedule(guided)
		#endif
		for(i=1; i<=n-b; ++i){
			j=i+b;
			calc_s1(i,j);
			calc_s2(i,j);
			calc_upm(i,j);
			calc_up(i,j);
			calc_s3(i,j);
			calc_u1d(i,j);
			calc_u1(i,j);
			calc_ud(i,j);
			calc_u(i,j);
		}
	}
	*/

}
void calc_s1(int h, int j)
{
	int l;
	double s1_val = 0.0;
	for (l = h+1; l < j; ++l)
	{
		double v1 = (get_up(h,l)*(myExp(-(ED5_new(h,l,h-1)+auPenalty_new(h,l))/RT)));
		double v2 = (myExp(-ED3_new(h,l,l+1)/RT)*get_u(l+2,j));
		double v3 = get_ud(l+1,j);
		double v4 = (get_up(l+1,j)*myExp(-(auPenalty_new(l+1,j)/RT)));
		double val = v1*(v2+v3+v4);
		s1_val += val;
	}
	set_s1(h,j,s1_val);
}
void calc_s2(int h, int j)
{
	int l;
	double s2_val = 0.0;
	for (l = h+1; l < j; ++l)
	{
		double v1 = (get_up(h,l)*(myExp(-(ED5_new(h,l,h-1)+auPenalty_new(h,l))/RT)));
		double v2 = (myExp(-(ED3_new(h,l,l+1)+EB_new())/RT)*get_u1(l+2,j-1));
		double v3 = get_u1d(l+1,j-1);
		double val = v1*(v2+v3);
		s2_val += val;
	}
	set_s2(h, j, s2_val);
}
void calc_s3(int h, int j)
{int l;
	double s3_val = 0.0;
	for (l = h+1; l <= j && l+2<=part_len; ++l){
		double v1 = (get_up(h,l)*(myExp(-(ED5_new(h,l,h-1)+auPenalty_new(h,l))/RT)));
		double v2 = (f(j+1,h,l)*myExp(-((j-l)*EB_new())/RT));
		double v3 = (myExp(-(ED3_new(h,l,l+1)+EB_new())/RT)*get_u1(l+2,j));
		double v4 = get_u1d(l+1,j);
		double val = v1*(v2+v3+v4);
		s3_val += val;
	}
	set_s3(h, j, s3_val);
}
void create_partition_arrays()
{
	int len = part_len + 2;
	u = mallocTwoD(len,len);
	up = mallocTwoD(len,len);
	upm = mallocTwoD(len,len);
	ud = mallocTwoD(len,len);
	u1d = mallocTwoD(len,len);
	s1 = mallocTwoD(len,len);
	s2 = mallocTwoD(len,len);
	s3 = mallocTwoD(len,len);
	u1 = mallocTwoD(len+1,len+1);
}
void free_partition_arrays()
{
	int len = part_len + 2;
	freeTwoD(u,len,len);
	freeTwoD(up,len,len);
	freeTwoD(upm,len,len);
	freeTwoD(ud,len,len);
	freeTwoD(u1d,len,len);
	freeTwoD(s1,len,len);
	freeTwoD(s2,len,len);
	freeTwoD(s3,len,len);
	freeTwoD(u1,len+1,len+1);
}
void calc_upm(int i, int j){
	double a = EA_new();
	double b = EB_new();
	double c = EC_new();
	double p_val = 0;
	int l,h;
	double quadraticSum = 0;
	if (canPair(RNA[i],RNA[j]))
	{
		for(l=i+2; l<j; ++l){
			double v1 = (get_up(i+1,l) * myExp((-1)*(a+2*c+auPenalty_new(i+1,l))/RT));
			double v2 = (myExp((-1)*(ED3_new(i+1,l,l+1)+b)/RT) * get_u1(l+2,j-1));
			double v3 = get_u1d(l+1,j-1);
			p_val = p_val + (v1*(v2+v3));
		}
		for(l=i+3; l<j; ++l){
			double v1 = (get_up(i+2,l)*myExp((-1)*(a+2*c+b+ED3_new(j,i,i+1)+auPenalty_new(i+2,l))/RT));
			double v2 = (myExp((-1)*(ED3_new(i+2,l,l+1)+b)/RT)*get_u1(l+2,j-1));
			double v3 = get_u1d(l+1,j-1);
			p_val = p_val + (v1*(v2+v3));
		}
		for(h=i+3; h<j-1; ++h){
			quadraticSum += (get_s2(h,j) * myExp((-1)*(a+2*c+(h-i-1)*b)/RT));
		}
		quadraticSum *= (myExp((-1)*ED3_new(j,i,i+1)/RT));
		p_val += quadraticSum;
		set_upm(i, j, p_val);  }
	else {
		set_upm(i, j, 0.0);  }
}
void calc_u1(int i, int j){
	double b = EB_new();
	double c = EC_new();
	double p_val = get_u1d(i,j);
	int h;
	double quadraticSum = 0;
	for(h=i+1; h<j; ++h){
		quadraticSum += (get_s3(h,j) * myExp((-1)*(c+(h-i)*b)/RT));
	}
	p_val += quadraticSum;
	set_u1(i, j, p_val);
}
void calc_u1d(int i, int j){
	double b = EB_new();
	double c = EC_new();
	double p_val = 0;
	int l;
	for(l=i+1; l<=j; ++l){
		double v1 = (get_up(i,l)*myExp((-1)*(c+auPenalty_new(i,l))/RT));
		double v2 = (f(j+1,i,l)*myExp((-1)*(j-l)*b/RT));
		double v3 = (myExp((-1)*(ED3_new(i,l,l+1)+b)/RT)*get_u1(l+2,j));
		double v4 = get_u1d(l+1,j);
		p_val += (v1*(v2+v3+v4));
	}
	set_u1d(i, j, p_val);
}
void calc_u(int i, int j)
{
	double uval = 1 + get_up(i,j)*myExp(-auPenalty_new(i,j)/RT);
	int h;
	int ctr;
	uval +=  get_ud(i,j);
	for (h = i+1; h < j; ++h) {
		uval += (get_up(h,j) * myExp( -(ED5_new(h,j,h-1) + auPenalty_new(h,j)) / RT ));
	}
	for (ctr = i+1; ctr < j-1; ++ctr) {
		uval += get_s1(ctr,j);
	}
	set_u(i, j, uval);
}
void calc_ud(int i, int j)
{
	int l;
	double udij = 0;
	for (l = i+1; l < j; ++l)
	{
		double val1, val2, val3;
		val1 = get_up(i,l);
		val1 = val1 * myExp(-auPenalty_new(i,l) / RT);
		val2 = get_u(l+2,j);
		val2 = val2 * myExp(-ED3_new(i,l,l+1)/RT);
		val3 = get_ud(l+1,j);
		val3 = val3 + get_up(l+1,j) * myExp( -auPenalty_new(l+1,j) / RT );
		udij += (val1 * (val2 + val3));
	}
	set_ud(i, j, udij);
}
void printUPprobabilities(int i, int j){
	int h,l;
	double maxIntLoopProb = 0.0;
	int h_max=-1, l_max=-1;
	double sumIntLoopProb=0.0;
	for (h = i+1; h < j ; h++) {
		for (l = h+1; l < j; l++) {
			if (canPair(RNA[h],RNA[l])==0) continue;
			if(h==(i+1) && l==(j-1)) continue;
			double intLoopProb = (get_up(h,l) * myExp(-((double)eL_new(i,j,h,l))/RT))/up[i][j];
			sumIntLoopProb+=intLoopProb;
			if(intLoopProb > maxIntLoopProb){ maxIntLoopProb = intLoopProb; h_max=h; l_max=l;}
		}
	}
	double hpProb = myExp(-((double)eH_new(i,j))/RT )/up[i][j];
	double stackProb = (myExp(-((double)eS_new(i,j))/RT ) * get_up(i+1,j-1))/up[i][j];
	double upmProb = get_upm(i,j)/up[i][j];
	if(sumIntLoopProb>=hpProb && sumIntLoopProb>=stackProb && sumIntLoopProb >=upmProb) printf("INT ");
	else if(hpProb>=sumIntLoopProb && hpProb>=stackProb && hpProb>=upmProb) printf("HPL ");
	else if(stackProb>=hpProb && stackProb>=sumIntLoopProb && stackProb>=upmProb) printf("STK ");
	else if(upmProb>=hpProb && upmProb>=stackProb && upmProb>=sumIntLoopProb) printf("UPM ");
	printf("printing probabilities: i=%d, j=%d, upmProb=%.6f, stackProb=%.6f, hpProb=%.6f, maxIntLoopProb=%.6f,  sumIntLoopProbs=%.6f, h_max=%d, l_max=%d\n",i,j, upmProb, stackProb, hpProb, maxIntLoopProb,sumIntLoopProb,h_max,l_max);
}
void calc_up(int i, int j)
{
	double up_val = 0.0;
	if (canPair(RNA[i],RNA[j]))
	{
		if (g_LIMIT_DISTANCE && j-i > g_contactDistance){
			set_up(i,j,0.0);
		}
		else {
			int h,l;
			for (h = i+1; h < j ; h++) {
				for (l = h+1; l < j; l++) {
					if (canPair(RNA[h],RNA[l])==0) continue;
					if(h==(i+1) && l==(j-1)) continue;
					up_val += (get_up(h,l) * myExp(-((double)eL_new(i,j,h,l))/RT));
				}
			}
			up_val = up_val + myExp(-((double)eH_new(i,j))/RT );
			up_val = up_val + (myExp(-((double)eS_new(i,j))/RT ) * get_up(i+1,j-1));
			up_val = up_val + get_upm(i,j);
			set_up(i, j, up_val);
			//printUPprobabilities(i,j);
		}
	}
	else  {
		set_up(i, j, 0.0);
	}
}
void calc_up_parallel(int i, int j)
{
	double up_val = 0.0;
	if (canPair(RNA[i],RNA[j]))
	{
		if (g_LIMIT_DISTANCE && j-i > g_contactDistance){
			set_up(i,j,0.0);
		}
		else {
			int h,l;
			#ifdef _OPENMP
			#pragma omp parallel for private (h,l) schedule(guided) reduction(+ : up_val)
			#endif
			for (h = i+1; h < j ; h++) {
				double my_up_val=0.0;
				for (l = h+1; l < j; l++) {
					if (canPair(RNA[h],RNA[l])==0) continue;
					if(h==(i+1) && l==(j-1)) continue;
					my_up_val += (get_up(h,l) * myExp(-((double)eL_new(i,j,h,l))/RT));
				}
				up_val += my_up_val;
			}
			up_val = up_val + myExp(-((double)eH_new(i,j))/RT );
			up_val = up_val + (myExp(-((double)eS_new(i,j))/RT ) * get_up(i+1,j-1));
			up_val = up_val + get_upm(i,j);
			set_up(i, j, up_val);
			//printUPprobabilities(i,j);
		}
	}
	else  {
		set_up(i, j, 0.0);
	}
}
