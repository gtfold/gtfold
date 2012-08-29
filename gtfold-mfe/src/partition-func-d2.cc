#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "partition-func-d2.h"
#include "energy.h"
#include "algorithms-partition.h"
#include "global.h"
#include "utils.h"
#include<omp.h>
#include <assert.h>

#define PAIRABLE_POINTS_GATHER_OPTIMIZATION_DISABLED true

static void errorAndExit(char* msg, int i, int j, MyDouble oldVal, MyDouble newVal){
	printf("%s\n", msg);
	printf("i=%d,j=%d,oldVal=",i,j);oldVal.print();printf(",newVal=");newVal.print();printf("\n");
	printf("%s","\nprogram is exiting now due to above error\n");
	exit(-1);
}
MyDouble PartitionFunctionD2::myExp(double arg){
	double posArgMultRT = (-1)*arg*RT;
	if(posArgMultRT>=INFINITY_){
	//if(posArgMultRT>=999999999999999){
		return MyDouble(0.0);
	}
	//Varifed it, working correctly
	return MyDouble(exp(arg));
}

/*void PartitionFunctionD2::printMatrix(MyDouble** u, int part_len){
	int i,j;
	for (i = 0; i <= part_len+1; ++i)
	{
		for (j = 0; j <= part_len+1; ++j){
			//printf("%0.1f ",u[i][j]);
			u[i][j].print();
			printf(",");
		}
		printf("\n");
	}
}*/

void PartitionFunctionD2::printMatrix(MyDouble** u, int part_len, FILE* pfarraysoutfile){
	int i,j;
	for (i = 0; i <= part_len+1; ++i)
	{
		for (j = 0; j <= part_len+1; ++j){
			//printf("%0.1f ",u[i][j]);
			u[i][j].print(pfarraysoutfile);
			fprintf(pfarraysoutfile, ",");
		}
		fprintf(pfarraysoutfile, "\n");
	}
}

//Functions to retrieve partition function array entries
//operator == varified with MyDouble
MyDouble PartitionFunctionD2::get_u(int i, int j) {if(u[i][j]==-1.0) errorAndExit("get_u entry is -1.",i,j,-1.0,0.0); return u[i][j];}
MyDouble PartitionFunctionD2::get_up(int i, int j) {if(up[i][j]==-1.0) errorAndExit("get_up entry is -1.\n",i,j,-1.0,0.0); return up[i][j];}
MyDouble PartitionFunctionD2::get_upm(int i, int j) {if(upm[i][j]==-1.0) errorAndExit("get_upm entry is -1.\n",i,j,-1.0,0.0); return upm[i][j];}
MyDouble PartitionFunctionD2::get_u1(int i, int j) {if(u1[i][j]==-1.0) errorAndExit("get_u1 entry is -1.\n",i,j,-1.0,0.0); return u1[i][j];}
MyDouble PartitionFunctionD2::get_s1(int i, int j) {if(s1[i][j]==-1.0) errorAndExit("get_s1 entry is -1.\n",i,j,-1.0,0.0); return s1[i][j];}
MyDouble PartitionFunctionD2::get_s2(int i, int j) {if(s2[i][j]==-1.0) errorAndExit("get_s2 entry is -1.\n",i,j,-1.0,0.0); return s2[i][j];}
MyDouble PartitionFunctionD2::get_s3(int i, int j) {if(s3[i][j]==-1.0) errorAndExit("get_s3 entry is -1.\n",i,j,-1.0,0.0); return s3[i][j];}

//Functions to set partition function array entries
//operator != varified with MyDouble
void PartitionFunctionD2::set_u(int i, int j, MyDouble val) {if(u[i][j]!=-1 && u[i][j]!=val) errorAndExit("set_u entry is not -1.\n",i,j,u[i][j],val); u[i][j]=val;}
void PartitionFunctionD2::set_up(int i, int j, MyDouble val) {if(up[i][j]!=-1 && up[i][j]!=val) errorAndExit("set_up entry is not -1.\n",i,j,up[i][j],val); up[i][j]=val;}
void PartitionFunctionD2::set_upm(int i, int j, MyDouble val) {if(upm[i][j]!=-1 && upm[i][j]!=val) errorAndExit("set_upm entry is not -1.\n",i,j,upm[i][j],val); upm[i][j]=val;}
void PartitionFunctionD2::set_u1(int i, int j, MyDouble val) {if(u1[i][j]!=-1 && u1[i][j]!=val) errorAndExit("set_u1 entry is not -1.\n",i,j,u1[i][j],val); u1[i][j]=val;}
void PartitionFunctionD2::set_s1(int i, int j, MyDouble val) {if(s1[i][j]!=-1 && s1[i][j]!=val) errorAndExit("set_s1 entry is not -1.\n",i,j,s1[i][j],val); s1[i][j]=val;}
void PartitionFunctionD2::set_s2(int i, int j, MyDouble val) {if(s2[i][j]!=-1 && s2[i][j]!=val) errorAndExit("set_s2 entry is not -1.\n",i,j,s2[i][j],val); s2[i][j]=val;}
void PartitionFunctionD2::set_s3(int i, int j, MyDouble val) {if(s3[i][j]!=-1 && s3[i][j]!=val) errorAndExit("set_s3 entry is not -1.\n",i,j,s3[i][j],val); s3[i][j]=val;}

/*//Following functions are to be used once testing completes as they are quicker then their test counterparts
inline double get_u(int i, int j) {return u[i][j];}
inline double get_up(int i, int j) {return up[i][j];}
inline double get_upm(int i, int j) {return upm[i][j];}
inline double get_u1(int i, int j) {return u1[i][j];}
inline double get_s1(int i, int j) {return s1[i][j];}
inline double get_s2(int i, int j) {return s2[i][j];}
inline double get_s3(int i, int j) {return s3[i][j];}
inline void set_u(int i, int j, double val) {u[i][j]=val;}
inline void set_up(int i, int j, double val) {up[i][j]=val;}
inline void set_upm(int i, int j, double val) {upm[i][j]=val;}
inline void set_u1(int i, int j, double val) {u1[i][j]=val;}
inline void set_s1(int i, int j, double val) {s1[i][j]=val;}
inline void set_s2(int i, int j, double val) {s2[i][j]=val;}
inline void set_s3(int i, int j, double val) {s3[i][j]=val;}

inline void set_u(int i, int j, double val) {u[i][j]=val;}
inline void set_up(int i, int j, double val) {up[i][j]=val;}
inline void set_upm(int i, int j, double val) {upm[i][j]=val;}
inline void set_u1(int i, int j, double val) {u1[i][j]=val;}
inline void set_s1(int i, int j, double val) {s1[i][j]=val;}
inline void set_s2(int i, int j, double val) {s2[i][j]=val;}
inline void set_s3(int i, int j, double val) {s3[i][j]=val;}
*/

//Functions providing general utilities related to energy
double PartitionFunctionD2::eS_new(int i, int j){
	if(PF_COUNT_MODE_) return 0;
	return eS(i,j);
	//return eS(i,j)/100;
}
double PartitionFunctionD2::eH_new(int i, int j){
	if(PF_COUNT_MODE_) return 0;
	return eH(i,j);
	//return eH(i,j)/100;
}
double PartitionFunctionD2::eL_new(int i, int j, int p, int q){
	if(PF_COUNT_MODE_) return 0;
	return eL(i,j,p,q);
	//return eL(i,j,p,q)/100;
}
double PartitionFunctionD2::ED3_new(int i, int j, int k){
	if(NO_DANGLE_MODE_) return 0;
	if(PF_COUNT_MODE_) return 0;
	//if(k > part_len) return 0;//This is to take care of round robin way of d2, this is shel's suggestion and rnafold also seems to follow this
	if(k > part_len) k=1;//This is to take care of round robin way of d2, rnascoring code follows this means, you need to do this in case you want to pass scoring test
	return Ed5(j,i,k);
	//return Ed5(j,i,k)/100;
}
double PartitionFunctionD2::ED5_new(int i, int j, int k){
	if(NO_DANGLE_MODE_) return 0;
	if(PF_COUNT_MODE_) return 0;
	//if (k<1) return 0;//This is to take care of round robin way of d2, this is shel's suggestion and rnafold also seems to follow this
	if (k<1) k=part_len;//This is to take care of round robin way of d2, rnascoring code follows this means, you need to do this in case you want to pass scoring test
	return Ed3(j,i,k);
	//return Ed3(j,i,k)/100;
}
double PartitionFunctionD2::EA_new(){
	if(PF_COUNT_MODE_) return 0;
	return Ea;
	//return Ea/100;
}
double PartitionFunctionD2::EB_new(){
	if(PF_COUNT_MODE_) return 0;
	return Ec;
	//return Ec/100;
}
double PartitionFunctionD2::EC_new(){
	if(PF_COUNT_MODE_) return 0;
	return Eb;
	//return Eb/100;
}
double PartitionFunctionD2::auPenalty_new(int i, int j){
	if(PF_COUNT_MODE_) return 0;
	return auPenalty(i,j);
	//return auPenalty(i,j)/100;
}
MyDouble PartitionFunctionD2::f(int j, int h, int l){
	//if(j - 1 == l || PF_COUNT_MODE_ || NO_DANGLE_MODE_)//TODO: if(j - 1 == l)
	if(PF_COUNT_MODE_ || NO_DANGLE_MODE_)//New: please confirm it
	return MyDouble(1.0);
	
	if(j - 1 == l)//TODO: if(j - 1 == l)
		return MyDouble(1.0);
	else{
		//return myExp(-ED3_new(h,l,l+1)/RT);
		return MyDouble(1.0);
	}
}

//Functions to calculate partition, and other partition function related utilities exposed to outside world
void PartitionFunctionD2::printAllMatrixes(){
	printf("\n\nAfter calculation, u matrix:\n\n");
	printMatrix(u,part_len,stdout);
	printf("\n\nAfter calculation, up matrix:\n\n");
	printMatrix(up,part_len,stdout);
	printf("\n\nAfter calculation, upm matrix:\n\n");
	printMatrix(upm,part_len,stdout);
	printf("\n\nAfter calculation, u1 matrix:\n\n");
	printMatrix(u1,part_len,stdout);
	printf("\n\nAfter calculation, s1 matrix:\n\n");
	printMatrix(s1,part_len,stdout);
	printf("\n\nAfter calculation, s2 matrix:\n\n");
	printMatrix(s2,part_len,stdout);
	printf("\n\nAfter calculation, s3 matrix:\n\n");
	printMatrix(s3,part_len,stdout);
}

void PartitionFunctionD2::printAllMatrixesToFile(string pfArraysOutputFile){
	FILE* pfarraysoutfile = fopen(pfArraysOutputFile.c_str(), "w");
        if(pfarraysoutfile==NULL){
        	cerr<<"Error in opening file: "<<pfarraysoutfile<<endl;
                exit(-1);
        }
	
	fprintf(pfarraysoutfile, "u matrix:\n\n");
        printMatrix(u,part_len,pfarraysoutfile);
        fprintf(pfarraysoutfile, "\n\nup matrix:\n\n");
        printMatrix(up,part_len,pfarraysoutfile);
        fprintf(pfarraysoutfile, "\n\nupm matrix:\n\n");
        printMatrix(upm,part_len,pfarraysoutfile);
        fprintf(pfarraysoutfile, "\n\nu1 matrix:\n\n");
        printMatrix(u1,part_len,pfarraysoutfile);
        fprintf(pfarraysoutfile, "\n\ns1 matrix:\n\n");
        printMatrix(s1,part_len,pfarraysoutfile);
        fprintf(pfarraysoutfile, "\n\ns2 matrix:\n\n");
        printMatrix(s2,part_len,pfarraysoutfile);
        fprintf(pfarraysoutfile, "\n\ns3 matrix:\n\n");
        printMatrix(s3,part_len,pfarraysoutfile);

	fclose(pfarraysoutfile);
        printf("\nPartition function arrays are printed to file %s\n\n", pfArraysOutputFile.c_str());


}
MyDouble PartitionFunctionD2::calculate_partition(int len, int pf_count_mode, int no_dangle_mode, bool PF_D2_UP_APPROX_ENABLED1)
{
	PF_COUNT_MODE_ = pf_count_mode;
	NO_DANGLE_MODE_ = no_dangle_mode;
	part_len = len;
	PF_D2_UP_APPROX_ENABLED = PF_D2_UP_APPROX_ENABLED1;
		
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
	/*if(g_verbose==1){
		printf("Printing partition function table...\n");
		printAllMatrixes();
	}*/
	//printf("%4.4f\n",u[1][part_len]);
	//if(pf_count_mode==1){ printf("Possible structure count: ");(u[1][part_len]).print();}
	if(pf_count_mode==1){ printf("Possible structure count: ");(u[1][part_len]).printInt();}
	else{ printf("Partition Function Value: ");(u[1][part_len]).print();}

	printf("\n");
	return u[1][part_len];

}
//partition arrays management related functions
void PartitionFunctionD2::free_partition()
{
	free_partition_arrays();
}
void PartitionFunctionD2::init_part_arrays_negatives(){
	int i,j,n;
	n = part_len+1;
	
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i,j) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	for(i=0; i<=n; ++i){
		for(j=0; j<=n; ++j){
			u[i][j]=MyDouble(-1.0);
			up[i][j]=MyDouble(-1.0);
			upm[i][j]=MyDouble(-1.0);
			s1[i][j]=MyDouble(-1.0);
			s2[i][j]=MyDouble(-1.0);
			s3[i][j]=MyDouble(-1.0);
			u1[i][j]=MyDouble(-1.0);
		}
	}
	
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i,j) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	
	for(i=0; i<=n+1; ++i){
		for(j=0; j<=n+1; ++j){
			u1[i][j]=MyDouble(-1.0);
		}
	}
}
void PartitionFunctionD2::init_partition_arrays()
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
			u[i][j] = MyDouble(1.0);
			up[i][j] = MyDouble(0.0);
			u1[i][j] = MyDouble(0.0);
			s1[i][j] = MyDouble(0.0);
			s2[i][j] = MyDouble(0.0);
			s3[i][j] = MyDouble(0.0);
		}
	}
	for(i=1; i<=n-4; ++i){
                s1[i][i+4] = MyDouble(0.0);
                s2[i][i+4] = MyDouble(0.0);
        }

	
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	
	for(i=1; i<=n; ++i){
		u[i+1][i] = MyDouble(1.0);
		u1[i+1][i] = MyDouble(0.0);
	}
	
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	#pragma omp parallel for private (i) schedule(guided)
	#endif
	//OPTIMIZED CODE ENDS
	//for(i=1; i<=n; i++){//OLD
	for(i=1; i<=n-1; i++){//NEW
		u1[i+2][i] = MyDouble(0.0);//TODO Uncomment it, as of now i have tested it that commenting it does not impact correctness
	}
}
void PartitionFunctionD2::create_partition_arrays()
{
	int len = part_len + 2;
	u = mallocTwoD(len,len);
	up = mallocTwoD(len,len);
	upm = mallocTwoD(len,len);
	s1 = mallocTwoD(len,len);
	s2 = mallocTwoD(len,len);
	s3 = mallocTwoD(len,len);
	u1 = mallocTwoD(len+1,len+1);
}
void PartitionFunctionD2::free_partition_arrays()
{
	int len = part_len + 2;
	freeTwoD(u,len,len);
	freeTwoD(up,len,len);
	freeTwoD(upm,len,len);
	freeTwoD(s1,len,len);
	freeTwoD(s2,len,len);
	freeTwoD(s3,len,len);
	freeTwoD(u1,len+1,len+1);
}

//general utility functions
MyDouble** PartitionFunctionD2::mallocTwoD(int r, int c) {
    MyDouble** arr = (MyDouble **)malloc(r*sizeof(MyDouble));
    int i,j;
    for(i=0; i<r; i++) {
        arr[i] = (MyDouble *)malloc(c*sizeof(MyDouble));

        // failed allocating a row, so free all previous rows, free the main
        // array, and return NULL
        if(arr[i] == NULL) {
            for(j=0; j<i; j++)
                free(arr[j]);
            free(arr);
            return NULL;
        }
    }
/*
    for(i=0; i<r; ++i){
	for(j=0; j<c; ++j){
	  //arr[i][j].init();
	}
    }
*/
    return arr;
}
void PartitionFunctionD2::freeTwoD(MyDouble** arr, int r, int c) {
    int i,j;
     for(i=0; i<r; ++i){
        for(j=0; j<c; ++j){
          arr[i][j].deallocate();
        }
    }

    for(i=0; i<r; i++)
        free(arr[i]);
    free(arr);
}

void PartitionFunctionD2::fill_partition_arrays()
{
	int b,i,j;
	int n=part_len;
	
	//Optimized Parallel Implementation
	int numThds = omp_get_num_threads();
	int* i_canPair = (int*)malloc((n+1)*sizeof(int));
	int* i_cannotPair = (int*)malloc((n+1)*sizeof(int));
	int len_i_canPair=0, len_i_cannotPair=0;
	int index1=0;
	int b_threshold = n;//n-numThds/2;//TURN+1;//n-numThds/2;//n;//TODO it should be equal to n only as parallel up is not defined completely
	for(b=TURN+1; b<b_threshold; ++b){
		//OPTIMIZED CODE STARTS
		len_i_canPair=0;
	       	len_i_cannotPair=0;
		for(i=1; i<=n-b; ++i){
			j=i+b;
			if(PAIRABLE_POINTS_GATHER_OPTIMIZATION_DISABLED || canPair(RNA[i],RNA[j])){i_canPair[len_i_canPair++]=i;}
			else {i_cannotPair[len_i_cannotPair++]=i;}
		}

		if(PAIRABLE_POINTS_GATHER_OPTIMIZATION_DISABLED) assert(len_i_cannotPair==0);

		#ifdef _OPENMP
		//#pragma omp parallel for private (index1,i,j) schedule(guided)
		#pragma omp parallel for private (index1,i,j) schedule(dynamic)
		#endif
		for(index1=0; index1<len_i_canPair; ++index1){
			i=i_canPair[index1];
			j=i+b;
			calc_s1(i,j);
			calc_s2(i,j);
			calc_upm(i,j);
			if(PF_D2_UP_APPROX_ENABLED){ calc_up_serial_and_approximate(i, j);}
			else calc_up(i,j);
			calc_s3(i,j);
			calc_u1(i,j);
			calc_u(i,j);
		}
		#ifdef _OPENMP
		//#pragma omp parallel for private (index1,i,j) schedule(guided)
		#pragma omp parallel for private (index1,i,j) schedule(dynamic)
		#endif
		for(index1=0; index1<len_i_cannotPair; ++index1){
			i=i_cannotPair[index1];
			j=i+b;
			calc_s1(i,j);
			calc_s2(i,j);
			calc_upm(i,j);
			if(PF_D2_UP_APPROX_ENABLED) calc_up_serial_and_approximate(i, j);
			else calc_up(i,j);
			calc_s3(i,j);
			calc_u1(i,j);
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
			if(PF_D2_UP_APPROX_ENABLED) calc_up_parallel_and_approximate(i, j);
			else calc_up_parallel(i,j);
			calc_s3(i,j);
			calc_u1(i,j);
			calc_u(i,j);
		}
	}
	

	/*
	//parallel implementation, sequential Implementation
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
			calc_u1(i,j);
			calc_u(i,j);
		}
	}*/

}

//Functions to calculate partition function array entries
void PartitionFunctionD2::calc_s1(int h, int j)
{
	int l;
	MyDouble s1_val(0.0);
	for (l = h+1; l < j; ++l)
	{
		MyDouble val = (get_up(h,l)*(myExp(-(ED5_new(h,l,h-1)+ED3_new(h,l,l+1)+auPenalty_new(h,l))/RT))*get_u(l+1,j));
		s1_val = s1_val + val;
	}
	set_s1(h,j,s1_val);
}
void PartitionFunctionD2::calc_s2(int h, int j)
{
	int l;
	MyDouble s2_val(0.0);
	for (l = h+1; l < j; ++l)
	{
		MyDouble val = (get_up(h,l)*(myExp(-(auPenalty_new(h,l)+ED5_new(h,l,h-1)+ED3_new(h,l,l+1))/RT))*get_u1(l+1,j-1));
		s2_val = s2_val + val;
	}
	set_s2(h, j, s2_val);
}
void PartitionFunctionD2::calc_s3(int h, int j)
{
	int l;
	MyDouble s3_val(0.0);
	//for (l = h+1; l <= j && l+2<=part_len; ++l){//TODO: old
	for (l = h+1; l <= j && l+1<=part_len; ++l){//TODO: new, comment it
		MyDouble v1 = (get_up(h,l)*(myExp(-(auPenalty_new(h,l)+ED5_new(h,l,h-1)+ED3_new(h,l,l+1))/RT)));
		MyDouble v2 = (f(j+1,h,l)*myExp(-((j-l)*EB_new())/RT));
		MyDouble val = v1*(v2 + get_u1(l+1,j));//TODO verify it as it is different from dS, in dS it is get_u1(l+2,j) and when we use it, we get error of re-using entry without initialization
		s3_val = s3_val + val;
	}
	set_s3(h, j, s3_val);
}
void PartitionFunctionD2::calc_upm(int i, int j){
	double a = EA_new();
	double b = EB_new();
	double c = EC_new();
	int l,h;
	MyDouble quadraticSum(0.0);//Default constructor of MyDouble will be called, which creates a double with value zero.
	if (canPair(RNA[i],RNA[j]))
	{
		//for(h=i+3; h<j-1; ++h){//TODO According to Shel's document
		for(h=i+1; h<j-1; ++h){//Manoj has changed it
			quadraticSum = quadraticSum + (get_s2(h,j) * myExp((-1)*((h-i-1)*b)/RT));
		}
		//quadraticSum = quadraticSum * (myExp((-1)*(a+ auPenalty_new(i,j) + ED5_new(j,i,j-1)/RT + ED3_new(j,i,i+1)/RT)));//TODO: make sure which one out of ed3(i,j,j-1) or ed3(i,j,j-1) is correct, similarly for ed5
		quadraticSum = quadraticSum * (myExp((-1)*(a+auPenalty_new(i,j) + ED5_new(j,i,j-1) + ED3_new(j,i,i+1) +2*c)/RT));//TODO Old impl using ed3(j,i) instead of ed3(i,j)
		//quadraticSum = quadraticSum * (myExp((-1)*(a+auPenalty_new(i,j) + ED5_new(i,j,j-1) + ED3_new(i,j,i+1) +2*c)/RT));//TODO New impl using ed5(i,j) instead of ed3(j,i)
		
		set_upm(i, j, quadraticSum);  
	}
	else {
		set_upm(i, j, 0.0);  
	}
}
void PartitionFunctionD2::calc_u1(int i, int j){
	double b = EB_new();
	double c = EC_new();
	int h;
	MyDouble quadraticSum(0.0);
	//for(h=i+1; h<j; ++h){//OLD
	for(h=i; h<j; ++h){//NEW, suggested by Shel
		quadraticSum = quadraticSum + (get_s3(h,j) * myExp((-1)*(c+(h-i)*b)/RT));
	}
	set_u1(i, j, quadraticSum);
}
void PartitionFunctionD2::calc_u(int i, int j)
{
	MyDouble uval(1.0);//uval will be initialized to double with value of 1.0
	int h;
	int ctr;
	//for (h = i+1; h < j; ++h) {//OLD: comment it
	for (h = i; h < j; ++h) {//TODO: New, check if OLD one is correct
		uval = uval + (get_up(h,j) * myExp( -(ED5_new(h,j,h-1) + ED3_new(h,j,j+1) + auPenalty_new(h,j)) / RT ));
	}
	//for (ctr = i+1; ctr < j-1; ++ctr) {//Shel's doc
	for (ctr = i; ctr < j-1; ++ctr) {//TODO Manoj corrected it
		uval = uval + get_s1(ctr,j);
	}
	set_u(i, j, uval);
}
void PartitionFunctionD2::calc_up(int i, int j)
{
	MyDouble up_val(0.0);
	if (canPair(RNA[i],RNA[j]))
	{
		if (g_LIMIT_DISTANCE && j-i > g_contactDistance){
			set_up(i,j,0.0);
		}
		else {
			int h,l;
			//for (h = i+1; h < j-1 ; h++) {
			for (h = i+1; h <= j-2-TURN ; h++) {
				for (l = h+1+TURN; l < j; l++) {
				//for (l = h+1; l < j; l++) {
					if (canPair(RNA[h],RNA[l])==0) continue;
					if(h==(i+1) && l==(j-1)) continue;
					up_val = up_val + (get_up(h,l) * myExp(-((double)eL_new(i,j,h,l))/RT));
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

void PartitionFunctionD2::calc_up_serial_and_approximate(int i, int j)
{
        MyDouble up_val(0.0);
        if (canPair(RNA[i],RNA[j]))
        {
        	if (g_LIMIT_DISTANCE && j-i > g_contactDistance){
				set_up(i,j,0.0);
			}
			else {
                int p,q;
                //for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
                for (p = i+1; p <= j-2-TURN ; p++) {
                        MyDouble my_up_val(0.0);
                        int minq = j-i+p-MAXLOOP-2;
                        if (minq < p+1+TURN) minq = p+1+TURN;
                        int maxq = (p==(i+1))?(j-2):(j-1);
                        for (q = minq; q <= maxq; q++) {
                                if (canPair(RNA[p],RNA[q])==0) continue;
                                my_up_val = my_up_val + (get_up(p,q) * myExp(-((double)eL_new(i,j,p,q))/RT));
                        }
                        up_val = up_val + my_up_val;
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
//TODO complete it
void PartitionFunctionD2::calc_up_parallel(int i, int j)
{
	MyDouble up_val(0.0);
	if (canPair(RNA[i],RNA[j]))
	{
		if (g_LIMIT_DISTANCE && j-i > g_contactDistance){
			set_up(i,j,0.0);
		}
		else {
			int h;
			/*#ifdef _OPENMP
        	        #pragma omp parallel for private (h) schedule(guided) reduction(+ : up_val)
            	    #endif*/
			for (h = i+1; h < j-1 ; h++) {
				int l;
				MyDouble my_up_val(0.0);
				for (l = h+1; l < j; l++) {
					if (canPair(RNA[h],RNA[l])==0) continue;
					if(h==(i+1) && l==(j-1)) continue;
					my_up_val = my_up_val + (get_up(h,l) * myExp(-((double)eL_new(i,j,h,l))/RT));
				}
				up_val = up_val + my_up_val;
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
//TODO complete it
void PartitionFunctionD2::calc_up_parallel_and_approximate(int i, int j)
{
        MyDouble up_val(0.0);
        if (canPair(RNA[i],RNA[j]))
        {
        	if (g_LIMIT_DISTANCE && j-i > g_contactDistance){
				set_up(i,j,0.0);
			}
			else {
                int p;
		/*#ifdef _OPENMP
                #pragma omp parallel for private (p) schedule(guided) reduction(+ : up_val)
                #endif*/
                for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
			int q;
                        MyDouble my_up_val(0.0);
                        int minq = j-i+p-MAXLOOP-2;
                        if (minq < p+1+TURN) minq = p+1+TURN;
                        int maxq = (p==(i+1))?(j-2):(j-1);
                        for (q = minq; q <= maxq; q++) {
                                if (canPair(p,q)==0) continue;
                                my_up_val = my_up_val + (get_up(p,q) * myExp(-((double)eL_new(i,j,p,q))/RT));
                        }
                        up_val = up_val + my_up_val;
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

/*
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
}*/

