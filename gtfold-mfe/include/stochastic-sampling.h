#ifndef __STOCHASTIC_SAMPLING__
#define __STOCHASTIC_SAMPLING__

#include <iostream>
#include <stdlib.h>
#include "partition-func.h"
#include "energy.h"

enum {U=0,UP,UD,U1,U1D};

double randdouble();
double U_0(int i, int j);
double U_ij(int i, int j);
double U_hj(int i, int h, int j);
double U_il(int i, int j);
double U_s1h(int i, int h, int j);
double U_ihlj_case1(int i, int h, int l, int j);
double U_ihlj_case2(int i, int h, int l, int j);  
double U_ihlj_case3(int i, int h, int l, int j);
double UD_il_case1(int i, int l, int j);
double UD_il_case2(int i, int l, int j);
double UD_il_case3(int i, int l, int j);
double Q_ijH(int i, int j);
double Q_ijS(int i, int j);
double Q_ijM(int i, int j);
double UPM_ip1l_case1(int i, int l, int j);
double UPM_ip1l_case2(int i, int l, int j);
double UPM_ip2l_case1(int i, int l , int j);
double UPM_ip2l_case2(int i, int l , int j);
double UPM_ijs2h(int i, int h , int j);
double UPM_ijhl_case1(int i, int h, int l, int j);
double UPM_ijhl_case2(int i, int h, int l, int j);
double U1D_ij_il_case1(int i, int l, int j);
double U1D_ij_il_case2(int i, int l, int j);
double U1D_ij_il_case3(int i, int l, int j);
double U1_ij(int i, int j);
double U1_ij_s3h(int i, int h, int j);
double U1_j_hl_case1(int h, int l, int j);
double U1_j_hl_case2(int h, int l, int j);
double U1_j_hl_case3(int h, int l, int j);

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

void set_single_stranded(int i, int j, int* structure);
void set_base_pair(int i, int j, int* structure);

void rnd_upm(int i, int j, int* structure);
void rnd_u1d(int i, int j, int* structure) ;
void  rnd_u1(int i, int j, int* structure);
void  rnd_up(int i, int j, int* structure);
void  rnd_ud(int i, int j, int* structure);
void   rnd_u(int i, int j, int* structure);

double rnd_structure(int* structure, int len);
void batch_sample(int num_rnd, int length, double U);
void batch_sample_and_dump(int num_rnd, int length, double U, std::string ctFileDumpDir, std::string stochastic_summery_file_name, std::string seq, std::string seqfile);

#endif
