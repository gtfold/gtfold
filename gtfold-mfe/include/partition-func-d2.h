#ifndef _PARTITION_FUNCTION_D2_H
#define _PARTITION_FUNCTION_D2_H

#include <math.h>
#include "MyDouble.cc"
/*
#ifdef __cplusplus
extern "C" {
#endif
 */
class PartitionFunctionD2{
	public:
		bool PF_D2_UP_APPROX_ENABLED;
	private:
		//Arrays to store partition function values
		MyDouble ** u;
		MyDouble ** up;
		MyDouble ** upm;
		MyDouble ** s1;
		MyDouble ** s2;
		MyDouble ** s3;
		MyDouble ** u1;
		//Different Modes and other variables
		int part_len;
		int PF_COUNT_MODE_;
		int NO_DANGLE_MODE_;
		//partition arrays management related functions
		void create_partition_arrays();
		void init_part_arrays_negatives();
		void init_partition_arrays();
		void fill_partition_arrays();
		void free_partition_arrays();
		//Functions to set partition function array entries
		void set_u(int i, int j, MyDouble val);
		void set_up(int i, int j, MyDouble val);
		void set_upm(int i, int j, MyDouble val);
		void set_u1(int i, int j, MyDouble val);
		void set_s1(int i, int j, MyDouble val);
		void set_s2(int i, int j, MyDouble val);
		void set_s3(int i, int j, MyDouble val);
		//Functions to calculate partition function array entries
		void calc_u(int i, int j);
		void calc_up(int i, int j);
		void calc_upm(int i, int j);
		void calc_u1(int i, int j);
		void calc_s1(int i, int j);
		void calc_s2(int i, int j);
		void calc_s3(int i, int j);
		void calc_up_serial_and_approximate(int i, int j);
		void calc_up_parallel_and_approximate(int i, int j);
		void calc_up_parallel(int i, int j);
		//general utility functions
		MyDouble **mallocTwoD(int r, int c);
		void freeTwoD(MyDouble** arr, int r, int c);
		//void printMatrix(MyDouble** u, int part_len);
		void printMatrix(MyDouble** u, int part_len, FILE* pfarraysoutputfile);//pfarraysoutputfile can be stdin in order to make it to print to standard output
	public:
		//Functions providing general utilities related to energy
		MyDouble myExp(double arg);
		double ED3_new(int i, int j, int k);
		double ED5_new(int i, int j, int k);
		double EA_new();
		double EB_new();
		double EC_new();
		double eS_new(int i, int j);
		double eL_new(int i, int j, int p, int q);
		double eH_new(int i, int j);
		double auPenalty_new(int i, int j);
		MyDouble f(int j, int h, int l);
		//Functions to retrieve partition function array entries
		MyDouble get_u(int i, int j);
		MyDouble get_up(int i, int j);
		MyDouble get_upm(int i, int j);
		MyDouble get_u1(int i, int j);
		MyDouble get_s1(int i, int j);
		MyDouble get_s2(int i, int j);
		MyDouble get_s3(int i, int j);
		//Functions to calculate partition, and other partition function related utilities exposed to outside world
		MyDouble calculate_partition(int len, int pf_count_mode, int no_dangle_mode, bool PF_D2_UP_APPROX_ENABLED);
		void free_partition();
		void printAllMatrixes();
		void printAllMatrixesToFile(string pfArraysOutputFile);
};
/*
#ifdef __cplusplus
}
#endif
 */
#endif
