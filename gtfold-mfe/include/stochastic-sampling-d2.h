#ifndef _STOCHASTIC_SAMPLING_D2_H
#define _STOCHASTIC_SAMPLING_D2_H

#include "pf-shel-check.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stack>
#include <stdlib.h>
#include "partition-func-d2.h"
#include "energy.h"
#include <math.h>

using namespace std;
//#include "MyDouble.cc"
/*
#ifdef __cplusplus
extern "C" {
#endif
 */
template <class MyDouble>
class StochasticTracebackD2{
public:
		struct base_pair 
		{
			int i;
			int j;
			int t;

			base_pair(int i_, int j_, int t_) : i(i_), j(j_), t(t_) {}
			base_pair(const base_pair& bp) :i(bp.i), j(bp.j), t(bp.t) { }
			base_pair& operator = (const base_pair& bp)  
			{
				if (this != &bp) 
				{
					i = bp.i;
					j = bp.j;
					t = bp.t;
				}
				return *this;
			}

			int type() const { return t ;} 

			bool isPaired() const 
			{
				return t == UP;
			}

			friend std::ostream& operator << (std::ostream& out, const base_pair& bp)
			{
				out << '(' << bp.i << '-' << bp.j << ')' << ' ' << bp.t << std::endl;
				return out;
			}
		};
	private:		
		pf_shel_check fraction;
		bool checkFraction;
		bool PF_D2_UP_APPROX_ENABLED;
		enum {U=0,UP,U1};
		//int* structure;

		PartitionFunctionD2<MyDouble> pf_d2;

		int print_energy_decompose;
		FILE* energy_decompose_outfile;
		//std::stack<base_pair> g_stack;
		//double energy;
		int length;
                int PF_COUNT_MODE;
                int NO_DANGLE_MODE;

 
		MyDouble randdouble();
                bool feasible(int i, int j);
		
		MyDouble U_0(int i, int j);
		MyDouble U_ihj(int i, int h, int j);
		MyDouble U_s1_ihj(int i, int h, int j);

		MyDouble S1_ihlj(int i, int h, int l, int j);

		MyDouble Q_H_ij(int i, int j);
		MyDouble Q_S_ij(int i, int j);
		MyDouble Q_M_ij(int i, int j);
		MyDouble Q_BI_ihlj(int i, int h, int l, int j);

		MyDouble UPM_S2_ihj(int i, int h, int j);

		MyDouble S2_ihlj(int i, int h, int l, int j);

		MyDouble U1_s3_ihj(int i, int h, int j);

		MyDouble S3_ihlj(int i, int h, int l, int j);
		MyDouble S3_MB_ihlj(int i, int h, int l, int j);

		void set_single_stranded(int i, int j, int* structure);
		void set_base_pair(int i, int j, int* structure);
		void rnd_u(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_s1(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_up(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_up_approximate(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_upm(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_s2(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_u1(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_s3(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		void rnd_s3_mb(int i, int h, int l, int j, int* structure, double & energy, std::stack<base_pair>& g_stack);
		double rnd_structure(int* structure);
		double rnd_structure_parallel(int* structure, int threads_for_one_sample);
		void updateBppFreq(std::string struc_str, int struc_freq, int ** bpp_freq, int length, int& total_bpp_freq);
		void printEnergyAndStructureInDotBracketAndTripletNotation(int* structure, std::string ensemble, int length, double energy, ostream& outfile);
		std::string getStructureStringInTripletNotation(int* structure, int length);
		std::string getStructureStringInTripletNotation(const char* ensemble, int length);
	public:
		void initialize(int length1, int PF_COUNT_MODE1, int NO_DANGLE_MODE1, int print_energy_decompose, bool PF_D2_UP_APPROX_ENABLED, bool checkFraction1, std::string energy_decompose_output_file, double scaleFactor);
		void free_traceback();
		void batch_sample(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION, bool ST_D2_ENABLE_UNIFORM_SAMPLE, double ST_D2_UNIFORM_SAMPLE_ENERGY, bool ST_D2_ENABLE_BPP_PROBABILITY, std::string sampleOutFile, std::string estimateBppOutputFile, std::string scatterPlotOutputFile);
		void batch_sample_parallel(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION, bool ST_D2_ENABLE_BPP_PROBABILITY, std::string sampleOutFile, std::string estimateBppOutputFile, std::string scatterPlotOutputFile);
		void batch_sample_and_dump(int num_rnd, std::string ctFileDumpDir, std::string stochastic_summery_file_name, std::string seq, std::string seqfile);
		void printPfMatrixesToFile(std::string pfArraysOutputFile);
};


//Now the definitions of all these method starts, earlier they were in separate C++ source file and then I needed to move them to header file only because of this
//http://www.parashift.com/c%2B%2B-faq-lite/templates-defn-vs-decl.html


#include "global.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stack>
#include <map>
#include<sstream>
#include<fstream>
#include "utils.h"
#include<omp.h>
#include<time.h>
#include<unistd.h>

//Basic utility functions
template <class MyDouble>
void StochasticTracebackD2<MyDouble>::initialize(int length1, int PF_COUNT_MODE1, int NO_DANGLE_MODE1, int print_energy_decompose1, bool PF_D2_UP_APPROX_ENABLED1, bool checkFraction1, std::string energy_decompose_output_file, double scaleFactor){
	checkFraction = checkFraction1;
	length = length1;
	//if(checkFraction) fraction = pf_shel_check(length);
	print_energy_decompose = print_energy_decompose1; 
	if(print_energy_decompose==1){
		//open file handler with file energy_decompose_output_file
		//FILE* energy_decompose_outfile;
        	energy_decompose_outfile = fopen(energy_decompose_output_file.c_str(), "w");
        	if(energy_decompose_outfile==NULL){
                	cerr<<"Error in opening file: "<<energy_decompose_output_file<<endl;
                	exit(-1);
        	}
		printf("\nEnergy decomposition for Stochastically sampled structure will be saved to %s\n\n", energy_decompose_output_file.c_str());
	}
	//energy = 0.0;
	//structure = new int[length+1];
	//std::stack<base_pair> g_stack;
	PF_COUNT_MODE = PF_COUNT_MODE1;
	NO_DANGLE_MODE = NO_DANGLE_MODE1;
	PF_D2_UP_APPROX_ENABLED = PF_D2_UP_APPROX_ENABLED1;
	pf_d2.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED, scaleFactor);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::free_traceback(){
	pf_d2.free_partition();
	if(print_energy_decompose==1){
		fclose(energy_decompose_outfile);
	}
	//delete[] structure;
}

template <class MyDouble>
inline MyDouble StochasticTracebackD2<MyDouble>::randdouble()
{
	return MyDouble( rand()/(double(RAND_MAX)+1) );
}

template <class MyDouble>
inline bool StochasticTracebackD2<MyDouble>::feasible(int i, int j)
{
	return j-i > TURN && canPair(RNA[i],RNA[j]);
}

//Probability calculation functions
template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::U_0(int i, int j){
	//return (MyDouble(1.0))/pf_d2.get_u(i,j);
	return pf_d2.myExp( pf_d2.get_M_RT()*(j-i+1) / RT) /pf_d2.get_u(i,j);
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::U_ihj(int i, int h, int j){
	//return (feasible(h,j) == true) ? (pf_d2.get_up(h,j)) * (pf_d2.myExp(-((pf_d2.ED5_new(h,j,h-1))+(pf_d2.ED3_new(h,j,j+1))+(pf_d2.auPenalty_new(h,j)))/RT)) / (pf_d2.get_u(i,j)) : MyDouble(0.0);
	return (feasible(h,j) == true) ? (pf_d2.get_up(h,j)) * (pf_d2.myExp(-((pf_d2.ED5_new(h,j,h-1))+(pf_d2.ED3_new(h,j,j+1))+(pf_d2.auPenalty_new(h,j))-pf_d2.get_M_RT()*(h-i))/RT)) / (pf_d2.get_u(i,j)) : MyDouble(0.0);
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::U_s1_ihj(int i, int h, int j){
	//return (pf_d2.get_s1(h,j)) / (pf_d2.get_u(i,j));
	return (pf_d2.get_s1(h,j)) * (pf_d2.myExp(pf_d2.get_M_RT()*(h-i)/RT)) / (pf_d2.get_u(i,j));
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::S1_ihlj(int i, int h, int l, int j){
	//return (feasible(h,l) == true) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.ED5_new(h,l,h-1))+(pf_d2.ED3_new(h,l,l+1))+(pf_d2.auPenalty_new(h,l)))/RT)) * (pf_d2.get_u(l+1,j)) / (pf_d2.get_s1(h,j)) : MyDouble(0.0);
	return (feasible(h,l) == true) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.ED5_new(h,l,h-1))+(pf_d2.ED3_new(h,l,l+1))+(pf_d2.auPenalty_new(h,l)))/RT)) * (pf_d2.get_u(l+1,j)) / (pf_d2.get_s1(h,j)) : MyDouble(0.0);
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::Q_H_ij(int i, int j){
	//return (pf_d2.myExp(-(pf_d2.eH_new(i,j))/RT)) / (pf_d2.get_up(i,j));
	return (pf_d2.myExp(-(pf_d2.eH_new(i,j)-pf_d2.get_M_RT()*(j-i+1))/RT)) / (pf_d2.get_up(i,j));
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::Q_S_ij(int i, int j){
	//return (pf_d2.myExp(-(pf_d2.eS_new(i,j))/RT)) * (pf_d2.get_up(i+1,j-1)) / (pf_d2.get_up(i,j));
	return (pf_d2.myExp(-(pf_d2.eS_new(i,j)-pf_d2.get_M_RT()*2)/RT)) * (pf_d2.get_up(i+1,j-1)) / (pf_d2.get_up(i,j));
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::Q_M_ij(int i, int j){
	//return (pf_d2.get_upm(i,j)) / (pf_d2.get_up(i,j));
	return (pf_d2.get_upm(i,j)) / (pf_d2.get_up(i,j));
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::Q_BI_ihlj(int i, int h, int l, int j){
	//return feasible(h,l) ? (pf_d2.myExp(-1*(pf_d2.eL_new(i,j,h,l))/RT)) * (pf_d2.get_up(h,l)) / (pf_d2.get_up(i,j)) : MyDouble(0.0);
	return feasible(h,l) ? (pf_d2.myExp(-1*(pf_d2.eL_new(i,j,h,l)-pf_d2.get_M_RT()*(h-i+j-l))/RT)) * (pf_d2.get_up(h,l)) / (pf_d2.get_up(i,j)) : MyDouble(0.0);
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::UPM_S2_ihj(int i, int h, int j){
	//return (pf_d2.get_s2(h,j)) * (pf_d2.myExp(-((pf_d2.EA_new())+ 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(j,i,j-1)) + (pf_d2.ED3_new(j,i,i+1)))/RT)) / (pf_d2.get_upm(i,j)); //TODO: Old impl, using ed3(j,i) instead of ed3(i,j), similarly in ed5
	return (pf_d2.get_s2(h,j)) * (pf_d2.myExp(-((pf_d2.EA_new())+ 2*(pf_d2.EB_new()) + (h-i-1)*(pf_d2.EC_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(j,i,j-1)) + (pf_d2.ED3_new(j,i,i+1))-pf_d2.get_M_RT()*(h-i))/RT)) / (pf_d2.get_upm(i,j)); //TODO: Old impl, using ed3(j,i) instead of ed3(i,j), similarly in ed5
	//return (pf_d2.get_s2(h,j)) * (pf_d2.myExp(-((pf_d2.EA_new())+ 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(i,j,j-1)) + (pf_d2.ED3_new(i,j,i+1)))/RT)) / (pf_d2.get_upm(i,j)); //TODO: New impl, using ed3(i,j) instead of ed3(j,i), similarly in ed5
	//return exp((-1)*ED3_new(j,i,i+1)/RT)* (s2[h][j] * exp((-1)*(EA_new()+2*EC_new()+(h-i-1)*EB_new())/RT))/upm[i][j];//TODO: New impl
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::S2_ihlj(int i, int h, int l, int j){
	//return feasible(h,l) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)))/RT)) * (pf_d2.get_u1(l+1,j-1)) / pf_d2.get_s2(h,j) : MyDouble(0.0);
	return feasible(h,l) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1))-pf_d2.get_M_RT())/RT)) * (pf_d2.get_u1(l+1,j-1)) / pf_d2.get_s2(h,j) : MyDouble(0.0);
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::U1_s3_ihj(int i, int h, int j){
	//return (pf_d2.get_s3(h,j)) * (pf_d2.myExp((-1)*((pf_d2.EC_new())+(h-i)*(pf_d2.EB_new()))/RT)) / (pf_d2.get_u1(i,j));
	//return (pf_d2.get_s3(h,j)) * (pf_d2.myExp((-1)*((pf_d2.EC_new())+(h-i)*(pf_d2.EB_new())-pf_d2.get_M_RT()*(h-i))/RT)) / (pf_d2.get_u1(i,j));
	return (pf_d2.get_s3(h,j)) * (pf_d2.myExp((-1)*((pf_d2.EB_new())+(h-i)*(pf_d2.EC_new())-pf_d2.get_M_RT()*(h-i))/RT)) / (pf_d2.get_u1(i,j));//Manoj111
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::S3_ihlj(int i, int h, int l, int j){
	//return feasible(h,l) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)))/RT)) * ( (pf_d2.myExp(-(j-l)*(pf_d2.EB_new())/RT)) * (pf_d2.f(j+1,h,l)) + (pf_d2.get_u1(l+1,j)) ) / (pf_d2.get_s3(h,j)) : MyDouble(0.0);
	return feasible(h,l) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)))/RT)) * ( (pf_d2.myExp(-(j-l)*(pf_d2.EC_new()-pf_d2.get_M_RT())/RT)) * (pf_d2.f(j+1,h,l)) + (pf_d2.get_u1(l+1,j)) ) / (pf_d2.get_s3(h,j)) : MyDouble(0.0);
}

template <class MyDouble>
MyDouble StochasticTracebackD2<MyDouble>::S3_MB_ihlj(int i, int h, int l, int j){
	MyDouble term1 =  (pf_d2.myExp(-(j-l)*(pf_d2.EB_new())/RT)) * (pf_d2.f(j+1,h,l));
	//MyDouble term2 = (pf_d2.get_u1(l+1,j);
	MyDouble term2 = (pf_d2.get_u1(l+1,j)*pf_d2.myExp(pf_d2.get_M_RT()*(l-j)/RT));
	return  term1 / (term1+term2);
}

/*
   double Q_ijBI(int i, int j)
   {
   double sum = 0;
   for (int h = i+1; h < j-1; ++h)
   for (int l = h+1; l < j; ++l)
   {
   if (h == i+1 && l == j-1) continue;
   sum += feasible(h,l)?(exp(-1*eL_new(i,j,h,l)/RT)*up[h][l]):0; 
   }
   return sum/up[i][j];
   }*/

//Functions related to sampling
template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_u(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	cum_prob = cum_prob + U_0(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction){
			//printf("U_0(i, j)");
			//printf("\n");
			fraction.add(0, i, j, false);
		}
		return;
	}

	for (int h = i; h < j; ++h)
	{
		cum_prob = cum_prob + U_ihj(i,h,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("U_ihj(i, h, j)");
				//printf("\n");
				fraction.add(1, h, j, true);
				fraction.add(0, i, j, false);
			}
			double e2 = ( (pf_d2.ED5_new(h,j,h-1)) + (pf_d2.ED3_new(h,j,j+1)) + (pf_d2.auPenalty_new(h,j)) );
			if (print_energy_decompose == 1) {
				fprintf(energy_decompose_outfile, " (pf_d2.ED5_new(h,j,h-1))=%f, (pf_d2.ED3_new(h,j,j+1))=%f, (pf_d2.auPenalty_new(h,j))=%f\n", (pf_d2.ED5_new(h,j,h-1))/100.0, (pf_d2.ED3_new(h,j,j+1))/100.0, (pf_d2.auPenalty_new(h,j))/100.0);
				fprintf(energy_decompose_outfile, " U_ihj(i=%d,h=%d,j=%d)= %lf\n",i,h,j, e2/100.0);
				
			}
			energy += e2;
			base_pair bp(h,j,UP);
			//set_single_stranded(i,h-1,structure);
			g_stack.push(bp);
			
			return;
		}
	}

	int h1 = -1;
	for (int h = i;  h < j-1; ++h)
	{
		cum_prob = cum_prob + U_s1_ihj(i,h,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("U_s1_ihj(i, h, j)");
				//printf("\n");
				fraction.add(2, h, j, true);
				fraction.add(0, i, j, false);
			}
			h1 = h;
			rnd_s1(i,h1,j, structure, energy, g_stack);
			return;
		}
	}
	//printf("rnd=");rnd.print();printf(",cum_prob=");cum_prob.print();
	assert (0) ;
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_s1(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	for (int l = h+1; l < j; ++l)
	{
		cum_prob = cum_prob + S1_ihlj(i,h,l,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("S1_ihlj(i,h,l,j)");
				//printf("\n");
				fraction.add(1, h, l, true);
				fraction.add(0, l+1, j, true);
				fraction.add(2, h, j, false);
			}
			double e2 = (pf_d2.ED5_new(h,l,h-1))+ (pf_d2.auPenalty_new(h,l)) + (pf_d2.ED3_new(h,l,l+1));
			if (print_energy_decompose == 1) {
				fprintf(energy_decompose_outfile, "(pf_d2.ED5_new(h,l,h-1))=%f, (pf_d2.auPenalty_new(h,l))=%f, (pf_d2.ED3_new(h,l,l+1))=%f\n",(pf_d2.ED5_new(h,l,h-1))/100.0, (pf_d2.auPenalty_new(h,l))/100.0, (pf_d2.ED3_new(h,l,l+1))/100.0);
				fprintf(energy_decompose_outfile, "(%d %d) %lf\n",i,j, e2/100.0);
			}
			energy += e2;
			base_pair bp1(h,l,UP);
			base_pair bp2(l+1,j,U);
			g_stack.push(bp1);
			g_stack.push(bp2);
			return ;
		}
	}
	assert(0);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_up(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	assert(structure[i] == 0);
	assert(structure[j] == 0);

	set_base_pair(i,j, structure);

	//for (int h = i+1; h < j-1; ++h)
	for (int h = i+1; h <= j-2-TURN ; h++)
                for (int l = h+1+TURN; l < j; l++)	
		//for (int l = h+1; l < j; ++l)
		{
			if (h == i+1 && l == j-1) continue;
			cum_prob = cum_prob + Q_BI_ihlj(i,h,l,j);
			if (rnd < cum_prob)
			{
				if(checkFraction) {
					//printf("Q_BI_ihlj(i, h, l, j)");
					//printf("\n");
					fraction.add(1, h, l, true);
					fraction.add(1, i, j, false);
				}
				double e2 = (pf_d2.eL_new(i,j,h,l));
				if (print_energy_decompose == 1) 
					fprintf(energy_decompose_outfile, "IntLoop(%d %d) %lf\n",i,j, e2/100.0);
				energy += e2;
				base_pair bp(h,l,UP);
				g_stack.push(bp);
				
				return;
			}
		}

	cum_prob = cum_prob + Q_H_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction){
			//printf("Q_H_ij(i, j)");
			//printf("\n");
			fraction.add(1, i, j, false);
		}
		double e2 = (pf_d2.eH_new(i,j));
		if (print_energy_decompose == 1) 
			fprintf(energy_decompose_outfile, "Hairpin(%d %d) %lf\n",i,j, e2/100.0);
		energy += e2;
		//set_single_stranded(i+1,j-1,structure);
		
		return ;
	}

	cum_prob = cum_prob + Q_S_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("Q_S_ij(i,j)");
			//printf("\n");
			fraction.add(1, i+1, j-1, true);
			fraction.add(1, i, j, false);
		}
		double e2 = (pf_d2.eS_new(i,j));
		if (print_energy_decompose == 1) 
			fprintf(energy_decompose_outfile, "Stack(%d %d) %lf\n",i,j, e2/100.0);
		energy+=e2;
		base_pair bp(i+1,j-1,UP);
		g_stack.push(bp);
		
		return ;
	}

	cum_prob = cum_prob + Q_M_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("Q_M_ij(i,j)");
			//printf("\n");
			fraction.add(3, i, j,true);
			fraction.add(1, i, j, false);
		}
		rnd_upm(i,j, structure, energy, g_stack);
		return;
	}
/*
	for (int h = i+1; h < j-1; ++h)
		for (int l = h+1; l < j; ++l)
		{
			if (h == i+1 && l == j-1) continue;
			cum_prob = cum_prob + Q_BI_ihlj(i,h,l,j);
			if (rnd < cum_prob)
			{
				if(checkFraction) {
					//printf("Q_BI_ihlj(i,h,l,j)");
					//printf("\n");
					fraction.add(1, h, l, true);
					fraction.add(1, i, j, false);
				}
				double e2 = (pf_d2.eL_new(i,j,h,l));
				if (print_energy_decompose == 1) 
					fprintf(energy_decompose_outfile, "IntLoop(%d %d) %lf\n",i,j, e2/100.0);
				energy += e2;
				base_pair bp(h,l,UP);
				g_stack.push(bp);
				return;
			}
		}
*/
	assert(0);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_up_approximate(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	assert(structure[i] == 0);
	assert(structure[j] == 0);

	set_base_pair(i,j, structure);

	//for (int p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1); ++p){
	for (int p = i+1; p <= j-2-TURN; ++p){
		int minq = j-i+p-MAXLOOP-2;
		if (minq < p+1+TURN) minq = p+1+TURN;
		int maxq = (p==(i+1))?(j-2):(j-1);
		for (int q = minq; q <=maxq ; ++q)
		{
			cum_prob = cum_prob + Q_BI_ihlj(i,p,q,j);
			if (rnd < cum_prob)
			{
				if(checkFraction) {
                                        //printf("Q_BI_ihlj(i, h, l, j)");
                                        //printf("\n");
                                        fraction.add(1, p, q, true);
                                        fraction.add(1, i, j, false);
                                }
				double e2 = (pf_d2.eL_new(i,j,p,q));
				if (print_energy_decompose == 1) 
					fprintf(energy_decompose_outfile, "IntLoop(%d %d) %lf\n",i,j, e2/100.0);
			
				energy += e2;
				base_pair bp(p,q,UP);
				g_stack.push(bp);
				return;
			}
		}
	}

	cum_prob = cum_prob + Q_H_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction){
                        //printf("Q_H_ij(i, j)");
                        //printf("\n");
                        fraction.add(1, i, j, false);
                }
		double e2 = (pf_d2.eH_new(i,j));
		if (print_energy_decompose == 1) 
			fprintf(energy_decompose_outfile, "Hairpin(%d %d) %lf\n",i,j, e2/100.0);
		energy += e2;
		//set_single_stranded(i+1,j-1,structure);
		return ;
	}

	cum_prob = cum_prob + Q_S_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("Q_S_ij(i,j)");
			//printf("\n");
			fraction.add(1, i+1, j-1, true);
			fraction.add(1, i, j, false);
		}
		double e2 = (pf_d2.eS_new(i,j));
		if (print_energy_decompose == 1) 
			fprintf(energy_decompose_outfile, "Stack(%d %d) %lf\n",i,j, e2/100.0);
		energy+=e2;
		base_pair bp(i+1,j-1,UP);
		g_stack.push(bp);
		return ;
	}

	cum_prob = cum_prob + Q_M_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("Q_M_ij(i,j)");
			//printf("\n");
			fraction.add(3, i, j,true);
			fraction.add(1, i, j, false);
		}
		rnd_upm(i,j, structure, energy, g_stack);
		return;
	}

	assert(0);
}


template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_u1(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);

	int h1 = -1;
	//for (int h = i+1; h < j-1; ++h)//TODO OLD 
	for (int h = i; h < j; ++h)//TODO NEW 
	{
		cum_prob = cum_prob + U1_s3_ihj(i,h,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("U1_s3_ihj(i,h,j)");
				//printf("\n");
				fraction.add(5, h, j, true);
				fraction.add(6, i, j, false);
			}
			double e2 = (pf_d2.EB_new()) + (h-i)*(pf_d2.EC_new());
			if (print_energy_decompose == 1){ 
				fprintf(energy_decompose_outfile, "(pf_d2.EB_new())=%f (h-i)*(pf_d2.EC_new())=%f\n",(pf_d2.EB_new())/100.0, (h-i)*(pf_d2.EC_new())/100.0);
				fprintf(energy_decompose_outfile, "U1_s3_ihj(%d %d %d) %lf\n",i,h,j, e2/100.0);
			}
			energy += e2;
			h1 = h;
			rnd_s3(i, h1, j, structure, energy, g_stack);
			return;
		}
	}
	assert(0);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_s3(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){
	// sample l given h1 
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	for (int l = h+1; l <= j &&  l+1<=length ; ++l)
	{
		cum_prob = cum_prob +  S3_ihlj(i,h,l,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("S3_ihlj(i,h,l,j)");
				//printf("\n");
				fraction.add(1, h, l, true);
				fraction.add(6, l+1, j, true);
				fraction.add(5, h, j, false);
			}
			double e2 = ((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)));
			if (print_energy_decompose == 1) {
				fprintf(energy_decompose_outfile, "(pf_d2.auPenalty_new(h,l))=%f, (pf_d2.ED5_new(h,l,h-1))=%f, (pf_d2.ED3_new(h,l,l+1))=%f\n",(pf_d2.auPenalty_new(h,l))/100.0, (pf_d2.ED5_new(h,l,h-1))/100.0, (pf_d2.ED3_new(h,l,l+1))/100.0);
				fprintf(energy_decompose_outfile, "S3_ihlj(%d %d %d %d) %lf\n",i,h,l,j,e2/100.0);
			}
			energy += e2;
			base_pair bp(h,l,UP);
			g_stack.push(bp);
			rnd_s3_mb(i,h,l,j, structure, energy, g_stack);
			return;
		}
	}
	assert(0);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_s3_mb(int i, int h, int l, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){//shel's document call this method with arguments i,h,l,j+1 therefore one will see difference of 1 in this code and shel's document
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	cum_prob = cum_prob +  S3_MB_ihlj(i,h,l,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("S3_MB_ihlj(i,h,l,j)");
			//printf("\n");
			fraction.add(6, l+1, j, false);
		}
		double tt =  0;//(j == l)? 0 : (pf_d2.ED3_new(h,l,l+1));//this term is corresponding to f(j+1,h,l)
		double e2 = tt + (j-l)*(pf_d2.EC_new());
		if (print_energy_decompose == 1){
			fprintf(energy_decompose_outfile, "j=%d,l=%d,tt=(j == l)?0:(pf_d2.ED3_new(h,l,l+1))=%f,(j-l)*(pf_d2.EC_new())=%f\n",j,l,tt/100.0,(j-l)*(pf_d2.EC_new())/100.0);
			fprintf(energy_decompose_outfile, "S3_MB_ihlj(%d %d %d %d) %lf\n",i,h,l,j, e2/100.0);
		}
		energy += e2;
		
		return;
	}
	else{
		base_pair bp1(l+1,j,U1);
		g_stack.push(bp1);
		return;
	}
	assert(0);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_upm(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	if (print_energy_decompose == 1)
		fprintf(energy_decompose_outfile, "Multiloop (%d %d)\n",i,j);

	int h1 = -1;
	for (int h = i+1; h < j-1; ++h)
	{
		cum_prob = cum_prob + UPM_S2_ihj(i,h,j);
		if (rnd < cum_prob )
		{
			if(checkFraction) {
				//printf("UPM_S2_ihj(i,h,j)");
				//printf("\n");
				fraction.add(4, h, j, true);
				fraction.add(3, i, j, false);
			}
			double e2 = (pf_d2.EA_new()) + 2*(pf_d2.EB_new()) + (h-i-1)*(pf_d2.EC_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(j,i,j-1)) + (pf_d2.ED3_new(j,i,i+1));//TODO Old impl using ed3(j,i) instead of ed3(i,j)
			//double e2 = (pf_d2.EA_new()) + 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(i,j,j-1)) + (pf_d2.ED3_new(i,j,i+1));//TODO New impl using ed3(i,j( instead of ed3(j,i)
			energy += e2;
			h1 = h;
			if (print_energy_decompose == 1) {
				fprintf(energy_decompose_outfile, "(pf_d2.EA_new())=%f, (pf_d2.EC_new())=%f, (pf_d2.EB_new())=%f\n",(pf_d2.EA_new())/100.0, (pf_d2.EC_new())/100.0, (pf_d2.EB_new())/100.0);
				fprintf(energy_decompose_outfile, "(pf_d2.EA_new()) + 2*(pf_d2.EB_new()) + (h-i-1)*(pf_d2.EC_new())=%f, (pf_d2.auPenalty_new(i,j))=%f, (pf_d2.ED5_new(j,i,j-1))=%f, (pf_d2.ED3_new(j,i,i+1))=%f\n",((pf_d2.EA_new()) + 2*(pf_d2.EB_new()) + (h-i-1)*(pf_d2.EC_new()))/100.0, (pf_d2.auPenalty_new(i,j))/100.0, (pf_d2.ED5_new(j,i,j-1))/100.0, (pf_d2.ED3_new(j,i,i+1))/100.0);
				fprintf(energy_decompose_outfile, "%s(%d %d %d) %lf\n", "UPM_S2_ihj",i,h1,j,e2/100.0);
			}
			rnd_s2(i,h1,j, structure, energy, g_stack);
			return;
		}
	}
	//assert(h1!=-1);
	assert(0);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::rnd_s2(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	for (int l = h+1; l < j; ++l)
	{
		cum_prob = cum_prob + S2_ihlj(i,h,l,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("S2_ihlj(i,h,l,j)");
				//printf("\n");
				fraction.add(1, h, l, true);
				fraction.add(6, l+1, j-1, true);
				fraction.add(4, h, j, false);
			}
			double e2 = (pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1));       
			energy += e2;
			if (print_energy_decompose == 1){
				fprintf(energy_decompose_outfile, "(pf_d2.auPenalty_new(h,l))=%f, (pf_d2.ED5_new(h,l,h-1))=%f, (pf_d2.ED3_new(h,l,l+1))=%f\n",(pf_d2.auPenalty_new(h,l))/100.0, (pf_d2.ED5_new(h,l,h-1))/100.0, (pf_d2.ED3_new(h,l,l+1))/100.0);
				fprintf(energy_decompose_outfile, "%s(%d %d %d %d) %lf\n"," S2_ihlj",i,h,l,j, e2/100.0);
			}
			base_pair bp1(h,l,UP);
			base_pair bp2(l+1,j-1,U1);
			g_stack.push(bp1);
			g_stack.push(bp2);
			
			return;
		}
	}
	assert(0);
}

template <class MyDouble>
double StochasticTracebackD2<MyDouble>::rnd_structure(int* structure)
{
	//printf("%lf %lf %lf\n", EA_new(), EB_new(), EC_new());
	srand(rand());
	//MyDouble U = pf_d2.get_u(1,len);
	base_pair first(1,length,U);
	std::stack<base_pair> g_stack;
	g_stack.push(first);
	double energy = 0.0;
	if(checkFraction) fraction = pf_shel_check(length);
	while (!g_stack.empty())
	{
		base_pair bp = g_stack.top();
		//   std::cout << bp;
		g_stack.pop();

		if (bp.type() == U)
			rnd_u(bp.i,bp.j, structure, energy, g_stack);
		else if (bp.type() == UP){
			if(pf_d2.PF_D2_UP_APPROX_ENABLED) {rnd_up_approximate(bp.i,bp.j, structure, energy, g_stack);}
			else rnd_up(bp.i,bp.j, structure, energy, g_stack);
		}
		else if (bp.type() == U1)
			rnd_u1(bp.i,bp.j, structure, energy, g_stack);
	}
	if(checkFraction){
		int remains = fraction.count();
		cout<<"Fraction Check test: remains:"<<remains<<endl;
		if(remains>0) cout<<"Fraction Check Error: remains="<<remains<<endl;
		fraction.clear();
	}
	return (double)energy/100.0;
}

template <class MyDouble>
double StochasticTracebackD2<MyDouble>::rnd_structure_parallel(int* structure, int threads_for_one_sample)
{
	//printf("%lf %lf %lf\n", EA_new(), EB_new(), EC_new());
	srand(rand());
	//MyDouble U = pf_d2.get_u(1,len);
	base_pair first(1,length,U);
	//std::stack<base_pair> g_stack;
	std::stack<base_pair> g_stack;
	g_stack.push(first);
	double energy = 0.0;
	std::stack<base_pair> g_stack_threads[threads_for_one_sample];
	double* energy_threads = new double[threads_for_one_sample];
	for(int index=0; index<threads_for_one_sample; ++index)energy_threads[index]=0.0;

	while (!g_stack.empty())
	{
		if(g_stack.size()%threads_for_one_sample != 0){
			base_pair bp = g_stack.top();
			//   std::cout << bp;
			g_stack.pop();

			if (bp.type() == U)
				rnd_u(bp.i,bp.j, structure, energy, g_stack);
			else if (bp.type() == UP){
				if(pf_d2.PF_D2_UP_APPROX_ENABLED) { rnd_up_approximate(bp.i,bp.j, structure, energy, g_stack);}
				else rnd_up(bp.i,bp.j, structure, energy, g_stack);
			}
			else if (bp.type() == U1)
				rnd_u1(bp.i,bp.j, structure, energy, g_stack);
		}
		else{
			std::deque<base_pair> g_deque;
			while(!g_stack.empty()){
				base_pair bp = g_stack.top();
				g_stack.pop();
				g_deque.push_back(bp);
			}
			int index;
			
			#ifdef _OPENMP
			//#pragma omp parallel for private(index) shared(energy_threads, g_stack_threads, structure) schedule(guided) num_threads(threads_for_one_sample)
			#pragma omp parallel for private(index) shared(energy_threads, g_stack_threads, structure) schedule(dynamic) num_threads(threads_for_one_sample)
			//#pragma omp parallel for private(index) shared(energy_threads, g_stack_threads, structure) schedule(guided)
			#endif
			for (index = 0; index < (int)g_deque.size(); ++index) {
				int thdId = omp_get_thread_num();
				base_pair bp = g_deque[index];
				if (bp.type() == U)
					rnd_u(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
				else if (bp.type() == UP){
					if(pf_d2.PF_D2_UP_APPROX_ENABLED) rnd_up_approximate(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
					else rnd_up(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
				}
				else if (bp.type() == U1)
					rnd_u1(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
			}

			for(index=0; index<threads_for_one_sample; ++index){
				energy += energy_threads[index];
				energy_threads[index] = 0.0;
				while(!g_stack_threads[index].empty()){
					base_pair bp = g_stack_threads[index].top();
					g_stack_threads[index].pop();
					g_stack.push(bp);
				}
			}
		}
	}
	delete[] energy_threads;
	return (double)energy/100.0;
}


/*
   void batch_sample(int num_rnd, int length, double U)
   {
   int* structure = new int[length+1];
   srand(time(NULL));
   std::map<std::string,std::pair<int,double> >  uniq_structs;

   if (num_rnd > 0 ) {
   printf("\nSampling structures...\n");
   int count; //nsamples =0;
   for (count = 1; count <= num_rnd; ++count) 
   {
   memset(structure, 0, (length+1)*sizeof(int));
   double energy = rnd_structure(structure, length);

   std::string ensemble(length+1,'.');
   for (int i = 1; i <= (int)length; ++ i) {
   if (structure[i] > 0 && ensemble[i] == '.')
   {
   ensemble[i] = '(';
   ensemble[structure[i]] = ')';
   }
   }
///double myEnegry = -88.4;
//++nsamples;
//if (fabs(energy-myEnegry)>0.0001) continue; //TODO: debug
//++count;

std::map<std::string,std::pair<int,double> >::iterator iter ;
if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
{
std::pair<int,double>& pp = iter->second;
pp.first++;
}
else {
uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
}

// std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
}
//std::cout << nsamples << std::endl;
int pcount = 0;
int maxCount = 0; std::string bestStruct;
double bestE = INFINITY;

std::map<std::string,std::pair<int,double> >::iterator iter ;
for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
{
const std::string& ss = iter->first;
const std::pair<int,double>& pp = iter->second;
const double& estimated_p =  (double)pp.first/(double)num_rnd;
const double& energy = pp.second;
double actual_p = pow(2.718281,-1.0*energy/RT_)/U;

printf("%s %lf %lf %lf %d\n",ss.c_str(),energy,actual_p,estimated_p,pp.first);
pcount += pp.first;
if (pp.first > maxCount)
{
maxCount = pp.first;
bestStruct  = ss;
bestE = pp.second;
}
}
assert(num_rnd == pcount);
printf("\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);

}

delete [] structure;
}
 */

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::batch_sample(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION ,bool ST_D2_ENABLE_UNIFORM_SAMPLE, double ST_D2_UNIFORM_SAMPLE_ENERGY, bool ST_D2_ENABLE_BPP_PROBABILITY, std::string samplesOutputFile, std::string estimateBppOutputFile, std::string scatterPlotOutputFile)
{cout<<"ST_D2_ENABLE_UNIFORM_SAMPLE="<<ST_D2_ENABLE_UNIFORM_SAMPLE<<",ST_D2_UNIFORM_SAMPLE_ENERGY="<<ST_D2_UNIFORM_SAMPLE_ENERGY<<endl;
	MyDouble U;
	
	ofstream outfile;
        outfile.open(samplesOutputFile.c_str());
	if(!outfile.good()){
                cerr<<"Error in opening file: "<<samplesOutputFile<<endl;
                exit(-1);
        }

	/*if(PF_D2_UP_APPROX_ENABLED){
	  double t1 = get_seconds();
	  PartitionFunctionD2 pf_d2_exact_up;
	  bool PF_D2_UP_APPROX_ENABLED2 = false;
	  pf_d2_exact_up.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED2);
	  U = pf_d2_exact_up.get_u(1,length);
	  pf_d2_exact_up.free_partition();
	  t1 = get_seconds() - t1;
	  printf("D2 Exact UP partition function computation running time: %9.6f seconds\n", t1);
	  }
	  else U = pf_d2.get_u(1,length);*/
	//U = pf_d2.get_u(1,length);
	U = pf_d2.unscale(1,length,pf_d2.get_u(1,length));
	srand(time(NULL));

	int threads_for_one_sample = 1;
	#ifdef _OPENMP
	#pragma omp parallel
	//#pragma omp master
	{
		int thdId1 = omp_get_thread_num();
		if(thdId1==0){
			if(g_nthreads < 0) threads_for_one_sample = omp_get_num_threads();//TODO move this line to mfe_main.cc
			else threads_for_one_sample = g_nthreads;
		}
	}
	#endif
	if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION) fprintf(stdout,"Stochastic Traceback: Thread count for one sample parallelization: %3d \n",threads_for_one_sample);

	std::map<std::string,std::pair<int,double> >  uniq_structs;
	int* structure = new int[length+1];

	if (num_rnd > 0 ) {
		printf("\nSampling structures...\n");
		int count, nsamples =0;
		for (count = 1; count <= num_rnd; ++count) 
		{
			nsamples++;
			memset(structure, 0, (length+1)*sizeof(int));
			double energy;
			if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
				energy = rnd_structure_parallel(structure, threads_for_one_sample);
			}
			else{
				energy = rnd_structure(structure);
			}

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			
			//Below line of codes is for finding samples with particular energy
			if(ST_D2_ENABLE_UNIFORM_SAMPLE){
				//double myEnegry = -92.1;//-91.3;//-94.8;//dS=-88.4;//d2=-93.1
				if (fabs(energy-ST_D2_UNIFORM_SAMPLE_ENERGY)>0.0001){ count--;continue;} //TODO: debug
			}

			std::map<std::string,std::pair<int,double> >::iterator iter ;
			if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
			{
				std::pair<int,double>& pp = iter->second;
				pp.first++;
				assert(energy==pp.second);
			}
			else {
				if(ST_D2_ENABLE_SCATTER_PLOT) uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
			}

			//if(!ST_D2_ENABLE_SCATTER_PLOT){
				//std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
				printEnergyAndStructureInDotBracketAndTripletNotation(structure, ensemble, (int)length, energy, outfile);
			//}
		}
		//std::cout << nsamples << std::endl;
		if(ST_D2_ENABLE_SCATTER_PLOT && !ST_D2_ENABLE_BPP_PROBABILITY){
			FILE* scatterPlotoutfile;
			scatterPlotoutfile = fopen(scatterPlotOutputFile.c_str(), "w");
                	if(scatterPlotoutfile==NULL){
                        	cerr<<"Error in opening file: "<<scatterPlotOutputFile<<endl;
                        	exit(-1);
                	}

			int pcount = 0;
			int maxCount = 0; std::string bestStruct;
			double bestE = INFINITY;
			fprintf(scatterPlotoutfile, "nsamples=%d\n",nsamples);
			fprintf(scatterPlotoutfile, "%s,%s,%s","structure","energy","boltzman_probability");
                        fprintf(scatterPlotoutfile, ",%s,%s\t%s\n","estimated_probability","frequency","structure in triplet notation");
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			int index=0;
			for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
			{
				index++;
				const std::string& ss = iter->first;
				const std::pair<int,double>& pp = iter->second;
				const double& estimated_p =  (double)pp.first/(double)num_rnd;
				const double& energy = pp.second;
				//MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy/RT_)))/U;
				//MyDouble actual_p = (pf_d2.myExp(-(energy)/(RT_)))/U;
				MyDouble actual_p = (pf_d2.myExp(-(energy*100)/(RT)))/U;
				//MyDouble actual_p(-(energy)/(RT_));///U;
				//fprintf(scatterPlotoutfile, "%s,%f,",ss.c_str(),energy);actual_p.print(scatterPlotoutfile);
				fprintf(scatterPlotoutfile, "S%d,%f,",index,energy);actual_p.print(scatterPlotoutfile);
				//fprintf(scatterPlotoutfile, ",%f,%d,\t",estimated_p,pp.first);
				fprintf(scatterPlotoutfile, ",%.20f,%d,\t",estimated_p,pp.first);
				std::string tripletNotationStructureString = getStructureStringInTripletNotation(ss.c_str(), length);
                                fprintf(scatterPlotoutfile, "%s\n", tripletNotationStructureString.c_str());

				//printf("%s %lf\n",ss.c_str(),energy);actual_p.print();
				//printf("%lf %d\n",estimated_p,pp.first);
				pcount += pp.first;
				if (pp.first > maxCount)
				{
					maxCount = pp.first;
					bestStruct  = ss;
					bestE = pp.second;
				}
			}
			assert(num_rnd == pcount);
			fprintf(scatterPlotoutfile, "\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);
			printf("\nScatter plot frequency data for stochastic samples, saved to %s\n", scatterPlotOutputFile.c_str());
			fclose(scatterPlotoutfile);
		}
		else{
			//printf("nsamples=%d\n",nsamples);
		}
		if(ST_D2_ENABLE_BPP_PROBABILITY){
			ofstream estimateBppoutfile;
                        estimateBppoutfile.open(estimateBppOutputFile.c_str());
                        if(!estimateBppoutfile.good()){
                                cerr<<"Error in opening file: "<<estimateBppOutputFile<<endl;
                                exit(-1);
                        }
			
			int** bpp_freq = new int*[length+1];
			for(int p=1; p<=length; ++p) bpp_freq[p] = new int[length+1];
			for(int p=1; p<=length; ++p) for(int q=p+1; q<=length; ++q) bpp_freq[p][q]=0;
			//for(int p=1; p<=length; ++p) for(int q=1; q<=length; ++q) bpp_freq[p][q]=0;
			int total_bpp_freq=0;
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
			{
				const std::string& struc_str = iter->first;
				const std::pair<int,double>& pp = iter->second;
				const int& struc_freq =  pp.first;
				updateBppFreq(struc_str, struc_freq, bpp_freq, length, total_bpp_freq);
			}
			//cout<<"\nBPP Probabilities are\ni,j,bppFreq,totalBppFreq\n";
			estimateBppoutfile<<"BPP Probabilities are\ni,j,bppFreq,totalSamples\n";
			for(int p=1; p<=length; ++p) for(int q=p+1; q<=length; ++q){
				//if(bpp_freq[p][q]>0) cout<<p<<","<<q<<","<<bpp_freq[p][q]<<","<<total_bpp_freq<<endl;
				if(bpp_freq[p][q]>0) estimateBppoutfile<<p<<","<<q<<","<<bpp_freq[p][q]<<","<<num_rnd<<endl;
			}
			for(int p=1; p<=length; ++p) delete[] bpp_freq[p];
			delete[] bpp_freq;
			printf("\nEstimated base pair probabilities for stochastic samples, saved to %s\n", estimateBppOutputFile.c_str());
			estimateBppoutfile.close();
		}

	}
	delete[] structure;
	printf("\nStochastic samples saved to %s\n", samplesOutputFile.c_str());
        outfile.close();
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::updateBppFreq(std::string struc_str, int struc_freq, int** bpp_freq, int length, int& total_bpp_freq){
	//cout<<"Entering updateBppFreq\n";
	std::stack<int> pos_stack;
	int i=0;
	while(i<length){
		while(i<length && struc_str[i]=='('){ pos_stack.push(i);i++;}
		while(i<length && struc_str[i]=='.'){ i++;}
		while(i<length && struc_str[i]==')'){
			int j=pos_stack.top();
			pos_stack.pop();
			//cout<<"pair: "<<j+1<<" "<<i+1<<endl;
			bpp_freq[j+1][i+1] += struc_freq;
			total_bpp_freq += struc_freq;
			i++;
		}
	}
	assert(pos_stack.size()==0);
	//cout<<"Exiting updateBppFreq\n";
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::batch_sample_parallel(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION, bool ST_D2_ENABLE_BPP_PROBABILITY, std::string samplesOutputFile, std::string estimateBppOutputFile, std::string scatterPlotOutputFile)
{
	//MyDouble U = pf_d2.get_u(1,length);
	MyDouble U;
	ofstream outfile;
        outfile.open(samplesOutputFile.c_str());
	if(!outfile.good()){
		cerr<<"Error in opening file: "<<samplesOutputFile<<endl;
		exit(-1);
	}
	/*if(PF_D2_UP_APPROX_ENABLED){
	  double t1 = get_seconds();
	  PartitionFunctionD2 pf_d2_exact_up;
	  bool PF_D2_UP_APPROX_ENABLED2 = false;
	  pf_d2_exact_up.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED2);
	  U = pf_d2_exact_up.get_u(1,length);
	  pf_d2_exact_up.free_partition();
	  t1 = get_seconds() - t1;
	  printf("D2 Exact UP partition function computation running time: %9.6f seconds\n", t1);
	  }
	  else U = pf_d2.get_u(1,length);
	 */
	//U = pf_d2.get_u(1,length);
	U = pf_d2.unscale(1,length,pf_d2.get_u(1,length));

	srand(time(NULL));
	/*
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	if (g_nthreads > 0) omp_set_num_threads(g_nthreads);
	#endif

	#ifdef _OPENMP
	#pragma omp parallel
	//#pragma omp master
	{
	int thdId1 = omp_get_thread_num();
	if(thdId1==0){
	fprintf(stdout,"Stochastic Traceback: Thread count: %3d \n",omp_get_num_threads());
	if(g_nthreads < 0) g_nthreads = omp_get_num_threads();//TODO move this line to mfe_main.cc
	}
	}
	#endif
	//OPTIMIZED CODE ENDSS
	 */
	int total_used_threads = 1;
	#ifdef _OPENMP
	#pragma omp parallel
	//#pragma omp master
	{
		int thdId1 = omp_get_thread_num();
		if(thdId1==0){
			if(g_nthreads < 0) total_used_threads = omp_get_num_threads();//TODO move this line to mfe_main.cc
			else total_used_threads = g_nthreads;
		}
	}
	#endif
	int threads_for_counts = total_used_threads;
	int threads_for_one_sample = 1;

	if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
		int sample_count_per_thread = num_rnd/total_used_threads;
		if(total_used_threads%2 !=0 || total_used_threads < 4 || length < 200 || sample_count_per_thread > 40){
			threads_for_counts = total_used_threads;
			threads_for_one_sample = 1;
		}
		else if(total_used_threads%4 !=0 || total_used_threads < 8 || length < 600 || sample_count_per_thread > 20 ){
			threads_for_counts = total_used_threads/2;
			threads_for_one_sample = 2;
		}
		else if(total_used_threads%8 !=0 || total_used_threads < 16 || length < 1500 || sample_count_per_thread > 10){
			threads_for_counts = total_used_threads/4;
			threads_for_one_sample = 4;
		}
		else if(total_used_threads%16 !=0 || total_used_threads < 32 || length < 3000 || sample_count_per_thread > 5 ){
			threads_for_counts = total_used_threads/8;
			threads_for_one_sample = 8;
		}
		else if(total_used_threads%32 !=0 || total_used_threads < 64 || length < 6000 || sample_count_per_thread > 2 ){
			threads_for_counts = total_used_threads/16;
			threads_for_one_sample = 16;
		}
		else{
			threads_for_counts = total_used_threads/32;
			threads_for_one_sample = 32;
		}

		fprintf(stdout,"Stochastic Traceback: Thread count for one sample parallelization: %3d \n",threads_for_one_sample);
	}
	fprintf(stdout,"Stochastic Traceback: Thread count for counts parallelization: %3d \n",threads_for_counts);



	std::map<std::string,std::pair<int,double> >  uniq_structs;
	//std::map<std::string,std::pair<int,double> >  uniq_structs_thread[g_nthreads];
	//g_nthreads=4;//TODO remove this line
	//cout<<"Manoj after: g_nthreads="<<g_nthreads<<endl;
	std::map<std::string,std::pair<int,double> > *  uniq_structs_thread = new std::map<std::string,std::pair<int,double> >[threads_for_counts];
	int* structures_thread = new int[threads_for_counts*(length+1)];

	if (num_rnd > 0 ) {
		printf("\nSampling structures...\n");
		int count;// nsamples =0;
		int* countArr = new int [threads_for_counts];
		for(int ind=0; ind<threads_for_counts; ind++) countArr[ind]=0;
		#ifdef _OPENMP
		//#pragma omp parallel for private (count) shared(structures_thread) schedule(guided) num_threads(threads_for_counts)
		#pragma omp parallel for private (count) shared(structures_thread, countArr, uniq_structs_thread, outfile) schedule(guided) num_threads(threads_for_counts)
		#endif
		for (count = 1; count <= num_rnd; ++count) 
		{
			//nsamples++;
			int thdId = omp_get_thread_num();
			countArr[thdId]++;
			//cout<<"thdId="<<thdId<<endl;
			int* structure = structures_thread + thdId*(length+1);
			memset(structure, 0, (length+1)*sizeof(int));
			//double energy = rnd_structure(structure);
			double energy;
			if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
				energy = rnd_structure_parallel(structure, threads_for_one_sample);
			}
			else{
				energy = rnd_structure(structure);
			}

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			/*
			//Below line of codes is for finding samples with particular energy
			double myEnegry = -92.1;//-91.3;//-94.8;//dS=-88.4;//d2=-93.1
			if (fabs(energy-myEnegry)>0.0001){ count--;continue;} //TODO: debug
			 */

			if(ST_D2_ENABLE_SCATTER_PLOT){
				std::map<std::string,std::pair<int,double> >::iterator iter ;
				if ((iter =uniq_structs_thread[thdId].find(ensemble.substr(1))) != uniq_structs_thread[thdId].end())
				{
					std::pair<int,double>& pp = iter->second;
					pp.first++;
					//cout<<"energy="<<energy<<",pp.second="<<pp.second<<endl;
					assert(energy==pp.second);
				}
				else{
					std::pair< std::string, std::pair<int,double> > new_pp = make_pair(ensemble.substr(1),std::pair<int,double>(1,energy));
					uniq_structs_thread[thdId].insert(new_pp); 
				}
			}
			//uniq_structs_thread[thdId].insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 

			//if(!ST_D2_ENABLE_SCATTER_PLOT){
				//std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
				//printEnergyAndStructureInDotBracketAndTripletNotation(structure, ensemble, (int)length, energy, std::cout);
				printEnergyAndStructureInDotBracketAndTripletNotation(structure, ensemble, (int)length, energy, outfile);
			//}
		}

		int finalCount=0;
		for(int ind=0; ind<threads_for_counts; ind++) finalCount+=countArr[ind];
		for (count = finalCount+1; count <= num_rnd; ++count) 
		{
			//nsamples++;
			int thdId = 0;//omp_get_thread_num();
			//countArr[thdId]++;
			//cout<<"thdId="<<thdId<<endl;
			int* structure = structures_thread + thdId*(length+1);
			memset(structure, 0, (length+1)*sizeof(int));
			//double energy = rnd_structure(structure);
			double energy;
			if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
				energy = rnd_structure_parallel(structure, threads_for_one_sample);
			}
			else{
				energy = rnd_structure(structure);
			}

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			/*
			//Below line of codes is for finding samples with particular energy
			double myEnegry = -92.1;//-91.3;//-94.8;//dS=-88.4;//d2=-93.1
			if (fabs(energy-myEnegry)>0.0001){ count--;continue;} //TODO: debug
			 */

			std::map<std::string,std::pair<int,double> >::iterator iter ;
			if ((iter =uniq_structs_thread[thdId].find(ensemble.substr(1))) != uniq_structs_thread[thdId].end())
			{
				std::pair<int,double>& pp = iter->second;
				pp.first++;
				//cout<<"energy="<<energy<<",pp.second="<<pp.second<<endl;
				assert(energy==pp.second);
			}
			else {
				if(ST_D2_ENABLE_SCATTER_PLOT){
					std::pair< std::string, std::pair<int,double> > new_pp = make_pair(ensemble.substr(1),std::pair<int,double>(1,energy));
					uniq_structs_thread[thdId].insert(new_pp); 
				}
				//uniq_structs_thread[thdId].insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
			}

			//if(!ST_D2_ENABLE_SCATTER_PLOT){
				//std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
				//printEnergyAndStructureInDotBracketAndTripletNotation(structure, ensemble, (int)length, energy, std::cout);
				printEnergyAndStructureInDotBracketAndTripletNotation(structure, ensemble, (int)length, energy, outfile);
			//}
		}



		if(ST_D2_ENABLE_SCATTER_PLOT){
			for(int thd_id=0; thd_id<threads_for_counts; thd_id++){
				std::map<std::string,std::pair<int,double> >::iterator thd_iter ;
				for (thd_iter = uniq_structs_thread[thd_id].begin(); thd_iter != uniq_structs_thread[thd_id].end();  ++thd_iter)
				{
					const std::string& thd_ss = thd_iter->first;
					const std::pair<int,double>& thd_pp = thd_iter->second;
					const int thd_freq = thd_pp.first;
					const double thd_e = thd_pp.second;

					std::map<std::string,std::pair<int,double> >::iterator iter ;
					if ((iter =uniq_structs.find(thd_ss)) != uniq_structs.end())
					{
						std::pair<int,double>& pp = iter->second;
						(pp.first)+=thd_freq;
						const double e = pp.second;
						//cout<<"e="<<e<<",thd_e="<<thd_e<<endl;
						assert(e==thd_e);
					}
					else {
						uniq_structs.insert(make_pair(thd_ss,std::pair<int,double>(thd_freq,thd_e)));
					}
				}
			}
		}
		if(ST_D2_ENABLE_SCATTER_PLOT && !ST_D2_ENABLE_BPP_PROBABILITY){
			FILE* scatterPlotoutfile;
                        scatterPlotoutfile = fopen(scatterPlotOutputFile.c_str(), "w");
                        if(scatterPlotoutfile==NULL){
                                cerr<<"Error in opening file: "<<scatterPlotOutputFile<<endl;
                                exit(-1);
                        }
			//std::cout << nsamples << std::endl;
			int pcount = 0;
			int maxCount = 0; std::string bestStruct;
			double bestE = INFINITY;
			fprintf(scatterPlotoutfile, "nsamples=%d\n",num_rnd);
			fprintf(scatterPlotoutfile, "%s,%s,%s","structure","energy","boltzman_probability");
			fprintf(scatterPlotoutfile, ",%s,%s\t%s\n","estimated_probability","frequency","structure in triplet notation");
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			int index=0;
			for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
			{
				index++;
				const std::string& ss = iter->first;
				const std::pair<int,double>& pp = iter->second;
				const double& estimated_p =  (double)pp.first/(double)num_rnd;
				const double& energy = pp.second;
				//MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy/RT_)))/U;
				//MyDouble actual_p = (pf_d2.myExp(-(energy)/(RT_)))/U;
				MyDouble actual_p;
				actual_p = (pf_d2.myExp(-(energy*100)/(RT)))/U;
				//MyDouble actual_p(-(energy)/(RT_));///U;
				//fprintf(scatterPlotoutfile, "%s,%f,",ss.c_str(),energy);actual_p.print(scatterPlotoutfile);
				fprintf(scatterPlotoutfile, "S%d,%f,",index,energy);actual_p.print(scatterPlotoutfile);
				//fprintf(scatterPlotoutfile, ",%f,%d,\t",estimated_p,pp.first);
				fprintf(scatterPlotoutfile, ",%.20f,%d,\t",estimated_p,pp.first);
				std::string tripletNotationStructureString = getStructureStringInTripletNotation(ss.c_str(), length);
				fprintf(scatterPlotoutfile, "%s\n", tripletNotationStructureString.c_str());
				//printf("%s %lf\n",ss.c_str(),energy);actual_p.print();
				//printf("%lf %d\n",estimated_p,pp.first);
				pcount += pp.first;
				if (pp.first > maxCount)
				{
					maxCount = pp.first;
					bestStruct  = ss;
					bestE = pp.second;
				}
			}
			assert(num_rnd == pcount);
			fprintf(scatterPlotoutfile, "\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);
			printf("\nScatter plot frequency data for stochastic samples, saved to %s\n", scatterPlotOutputFile.c_str());
                        fclose(scatterPlotoutfile);	
		}
		else{
			//printf("nsamples=%d\n",num_rnd);
		}
		if(ST_D2_ENABLE_BPP_PROBABILITY){
			ofstream estimateBppoutfile;
        		estimateBppoutfile.open(estimateBppOutputFile.c_str());
        		if(!estimateBppoutfile.good()){
                		cerr<<"Error in opening file: "<<estimateBppOutputFile<<endl;
        	        	exit(-1);
        		}

                        int** bpp_freq = new int*[length+1];
                        for(int p=1; p<=length; ++p) bpp_freq[p] = new int[length+1];
                        for(int p=1; p<=length; ++p) for(int q=p+1; q<=length; ++q) bpp_freq[p][q]=0;
                        //for(int p=1; p<=length; ++p) for(int q=1; q<=length; ++q) bpp_freq[p][q]=0;
                        int total_bpp_freq=0;
                        std::map<std::string,std::pair<int,double> >::iterator iter ;
                        for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
                        {
                                const std::string& struc_str = iter->first;
                                const std::pair<int,double>& pp = iter->second;
                                const int& struc_freq =  pp.first;
                                updateBppFreq(struc_str, struc_freq, bpp_freq, length, total_bpp_freq);
                        }
                        //cout<<"\nBPP Probabilities are\ni,j,bppFreq,totalBppFreq\n";
                        estimateBppoutfile<<"BPP Probabilities are\ni,j,bppFreq,totalSamples\n";
                        for(int p=1; p<=length; ++p) for(int q=p+1; q<=length; ++q){
                                //if(bpp_freq[p][q]>0) cout<<p<<","<<q<<","<<bpp_freq[p][q]<<","<<total_bpp_freq<<endl;
                                if(bpp_freq[p][q]>0) estimateBppoutfile<<p<<","<<q<<","<<bpp_freq[p][q]<<","<<num_rnd<<endl;
                        }
                        for(int p=1; p<=length; ++p) delete[] bpp_freq[p];
                        delete[] bpp_freq;
			printf("\nEstimated base pair probabilities for stochastic samples, saved to %s\n", estimateBppOutputFile.c_str());
			estimateBppoutfile.close();
                }
	}
	delete [] structures_thread;
	delete [] uniq_structs_thread;
	printf("\nStochastic samples saved to %s\n", samplesOutputFile.c_str());
        outfile.close();
}


template <class MyDouble>
void StochasticTracebackD2<MyDouble>::batch_sample_and_dump(int num_rnd, std::string ctFileDumpDir, std::string stochastic_summery_file_name, std::string seq, std::string seqfile)
{
	//MyDouble U = pf_d2.get_u(1,length);
	 MyDouble U;
        /*if(PF_D2_UP_APPROX_ENABLED){
		double t1 = get_seconds();
                PartitionFunctionD2 pf_d2_exact_up;
                bool PF_D2_UP_APPROX_ENABLED2 = false;
                pf_d2_exact_up.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED2);
                U = pf_d2_exact_up.get_u(1,length);
                pf_d2_exact_up.free_partition();
                t1 = get_seconds() - t1;
                printf("D2 Exact UP partition function computation running time: %9.6f seconds\n", t1);
       }
        else U = pf_d2.get_u(1,length);
	*/
         U = pf_d2.get_u(1,length);
	//data dump preparation code starts here
	if(ctFileDumpDir.compare("")==0){
		char abspath[1000];
		char* tmp = getcwd(abspath, 1000);
		if(tmp!=abspath){//TODO debug
			cout<<"Error in getcwd, exiting...\n";
			exit(-1);
		}
		ctFileDumpDir = abspath;
	}
	cout<<"Using ctFileDumpDir = "<<ctFileDumpDir<<endl;
	std::stringstream ss;
	ss<<ctFileDumpDir<<"/"<<stochastic_summery_file_name;
	stochastic_summery_file_name = ss.str();
	cout<<"Using stochastic_summary_file_name = "<<stochastic_summery_file_name<<endl;
	std::ofstream summaryoutfile;
	summaryoutfile.open(stochastic_summery_file_name.c_str());
	std::string seqname = seqfile.substr(seqfile.find_last_of("/\\") + 1, seqfile.length() - 1);
	cout<<"Sequence Name = "<<seqname<<endl;
	//data dump preparation code ends here

	srand(time(NULL));
	std::map<std::string,std::pair<int,double> >  uniq_structs;
	int* structure = new int[length+1];
	if (num_rnd > 0 ) {
		printf("\nSampling structures...\n");
		int count;
		//int nsamples =0;
		for (count = 1; count <= num_rnd; ++count) 
		{
			memset(structure, 0, (length+1)*sizeof(int));
			double energy = rnd_structure(structure);

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			/*
			//below code is for uniform sampling test
			double myEnegry = -88.4;
			++nsamples;
			if(nsamples>1000)break;
			if (fabs(energy-myEnegry)>0.0001) continue; //TODO: debug
			//++count;
			*/
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
			{
				std::pair<int,double>& pp = iter->second;
				pp.first++;
			}
			else {
				uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
			}

			// std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
			//data dump code starts here
			std::stringstream ss;
			ss<<ctFileDumpDir<<"/"<<seqname<<"_"<<count<<".ct";
			save_ct_file(ss.str(), seq, energy, structure);
			summaryoutfile<<ss.str()<<" "<<ensemble.substr(1)<<" "<<energy<< std::endl;
			//data dump code ends here
		}
		//std::cout << nsamples << std::endl;
		int pcount = 0;
		int maxCount = 0; std::string bestStruct;
		double bestE = INFINITY;

		printf("%s,%s,%s","structure","energy","boltzman_probability");
		printf(",%s,%s\n","estimated_probability","frequency");
		std::map<std::string,std::pair<int,double> >::iterator iter ;
		for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
		{
			const std::string& ss = iter->first;
			const std::pair<int,double>& pp = iter->second;
			const double& estimated_p =  (double)pp.first/(double)num_rnd;
			const double& energy = pp.second;
	
			//MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy/RT_)))/U;
			MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy*100/RT)))/U;
			printf("%s,%lf,",ss.c_str(),energy);actual_p.print();
			printf(",%lf,%d\n",estimated_p,pp.first);
	
			pcount += pp.first;
			if (pp.first > maxCount)
			{
				maxCount = pp.first;
				bestStruct  = ss;
				bestE = pp.second;
			}
		}
		assert(num_rnd == pcount);
		printf("\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);
	}
	summaryoutfile.close();
	delete[] structure;

}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::set_single_stranded(int i, int j, int* structure)
{
	for(;i<=j;++i) 
		structure[i] = 0;
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::set_base_pair(int i, int j, int* structure)
{
	bool cond = j-i > TURN && canPair(RNA[i],RNA[j]);
	assert(cond);
	structure[i] = j;
	structure[j] = i;
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::printEnergyAndStructureInDotBracketAndTripletNotation(int* structure, std::string ensemble, int length, double energy, ostream & outfile){
	//std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
	//ssobj << energy << "\t";	
	std::string tripletNotationStructureString = getStructureStringInTripletNotation(structure, length);
	stringstream printline;
	//outfile << ensemble.substr(1) << "\t" << energy << "\t" << tripletNotationStructureString<<endl;
	printline << ensemble.substr(1) << "\t" << energy << "\t" << tripletNotationStructureString<<endl;
	#pragma omp critical
	{
	outfile << printline.str();
	if(print_energy_decompose==1){
		fprintf(energy_decompose_outfile, "%s\t%f\t%s\n\n\n", ensemble.substr(1).c_str(), energy, tripletNotationStructureString.c_str());
	}
	}
}

template <class MyDouble>
std::string StochasticTracebackD2<MyDouble>::getStructureStringInTripletNotation(int* structure, int length){
	stringstream ssobj;
	int i=1;
        while( i <= (int)length ) {
                int myI = i;
                int myJ = structure[i];
                int myCount = 1;
                if (myJ > 0 && myI<myJ)
                {
                        for(int tempI=myI+1; tempI <= (int)length; ++tempI){
                                int tempJ = structure[tempI];
                                if(tempJ>0  && tempI<tempJ  && (tempI-myI)==(myJ-tempJ)){
                                        myCount++;
                                        i = tempI+1;
                                }
                                else{
                                        //std::cout<<myI<<" "<<myJ<<" "<<myCount<<", "; 
                                        ssobj<<myI<<" "<<myJ<<" "<<myCount<<", ";
                                        i = tempI;
                                        break;
                                }
                        }
                }
                else{
                        i++;
                }
        }
        return ssobj.str();
}

template <class MyDouble>
std::string StochasticTracebackD2<MyDouble>::getStructureStringInTripletNotation(const char* ensemble, int length){
	std::stack<int> openBracketStack;
	int* structure = new int[length+1];
	for(int j=0; j<=length; ++j) structure[j]=-1;
	int i=0;
	while(i<length){
		if(ensemble[i]=='('){
			openBracketStack.push(i);
		}
		else if(ensemble[i]==')'){
			if(openBracketStack.empty()){
				printf("%s structure is not a valid structure in dot-bracket notation, exiting...\n\n", ensemble);
				exit(-1);
			}
			int openingBracketIndexForI = openBracketStack.top();
			openBracketStack.pop();
			structure[i+1]=openingBracketIndexForI+1;
			structure[openingBracketIndexForI+1]=i+1;
		}
		i++;
	}
	if(!openBracketStack.empty()){
        	printf("%s structure is not a valid structure in dot-bracket notation, exiting...\n\n", ensemble);
        	exit(-1);
        }

	return getStructureStringInTripletNotation(structure, length);
}

template <class MyDouble>
void StochasticTracebackD2<MyDouble>::printPfMatrixesToFile(std::string pfArraysOutputFile){
	pf_d2.printAllMatrixesToFile(pfArraysOutputFile);
}


#endif
