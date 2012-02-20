#ifndef _STOCHASTIC_SAMPLING_D2_H
#define _STOCHASTIC_SAMPLING_D2_H

#include <iostream>
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
		enum {U=0,UP,U1};
		int* structure;
		PartitionFunctionD2 pf_d2;
		int ss_verbose;
		std::stack<base_pair> g_stack;
		double energy;
		int length;
                int PF_COUNT_MODE;
                int NO_DANGLE_MODE;
 
		MyDouble randdouble();
                bool feasible(int i, int j);
		
		MyDouble U_0(int i, int j);
		MyDouble U_hj(int i, int j);
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

		void set_single_stranded(int i, int j);
		void set_base_pair(int i, int j);
		
		void rnd_u(int i, int j);
		void rnd_s1(int i, int h, int j);
		void rnd_up(int i, int j);
		void rnd_upm(int i, int j);
		void rnd_s2(int i, int h, int j);
		void rnd_u1(int i, int j);
		void rnd_s3(int i, int h, int j);
		void rnd_s3_mb(int i, int h, int l, int j);

		double rnd_structure();
	public:
		void initialize(int length1, int PF_COUNT_MODE1, int NO_DANGLE_MODE1, int ss_verbose1);
		void free_traceback();
		void batch_sample(int num_rnd);
		void batch_sample_and_dump(int num_rnd, std::string ctFileDumpDir, std::string stochastic_summery_file_name, std::string seq, std::string seqfile);


};
#endif
