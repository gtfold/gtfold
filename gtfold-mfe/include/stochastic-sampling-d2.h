#ifndef _STOCHASTIC_SAMPLING_D2_H
#define _STOCHASTIC_SAMPLING_D2_H

#include "pf-shel-check.h"
#include <iostream>
#include <fstream>
#include<stdio.h>
#include <stack>
#include <stdlib.h>
#include "partition-func-d2.h"
#include "energy.h"
#include <math.h>
#include "MyDouble.cc"
/*
#ifdef __cplusplus
extern "C" {
#endif
 */
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
		PartitionFunctionD2 pf_d2;
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
		string getStructureStringInTripletNotation(int* structure, int length);
		string getStructureStringInTripletNotation(const char* ensemble, int length);
	public:
		void initialize(int length1, int PF_COUNT_MODE1, int NO_DANGLE_MODE1, int print_energy_decompose, bool PF_D2_UP_APPROX_ENABLED, bool checkFraction1, string energy_decompose_output_file);
		void free_traceback();
		void batch_sample(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION, bool ST_D2_ENABLE_UNIFORM_SAMPLE, double ST_D2_UNIFORM_SAMPLE_ENERGY, bool ST_D2_ENABLE_BPP_PROBABILITY, string sampleOutFile, string estimateBppOutputFile, string scatterPlotOutputFile);
		void batch_sample_parallel(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION, bool ST_D2_ENABLE_BPP_PROBABILITY, string sampleOutFile, string estimateBppOutputFile, string scatterPlotOutputFile);
		void batch_sample_and_dump(int num_rnd, std::string ctFileDumpDir, std::string stochastic_summery_file_name, std::string seq, std::string seqfile);
		void printPfMatrixesToFile(string pfArraysOutputFile);
};
#endif
