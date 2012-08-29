/**
GTfold: compute minimum free energy of RNA secondary structure
Copyright (C) 2008  David A. Bader
http://www.cc.gatech.edu/~bader

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@author prashant {pgaurav@gatech.edu}

 */

#include <cstdio>

#include "constants.h"
#include "energy.h"
#include "utils.h"
#include "global.h"
#include "subopt_traceback.h"

#include <iostream>
#include <iterator>
#include <cstdlib>

//#define DEBUG 1

using std::pair;
using std::cout;
using std::endl;

const char* lstr[] = {"W", "V", "VBI", "VM", "WM", "WMPrime", "fm", "fm1"};

void (*trace_func[8]) (int i, int j, ps_t& ps, ps_stack_t& gs);
static int delta = 0;
static int mfe = INFINITY_;
static int length = -1;
static int gflag = 0;

#ifdef UNIQUE_MULTILOOP_DECOMPOSITION

static int FM1[1500][1500] = {{0}};
static int FM[1500][1500] = {{0}};

static inline int Ed5_new(int i, int j, int k) {
  return (k!=0) ? Ed3(j, i, k) : Ed3(j,i,length);
}

static inline int Ed3_new(int i, int j, int k) {
  return (k!=length+1)?Ed5(j, i, k) : Ed5(j, i, 1);
}

void calculate_fm1() {

  for (int i = 1; i <= length; ++i) {
    for (int j = i+1; j <= length; ++j) {
      int min = INFINITY_;
      for (int l = i+TURN+1; l <= j; ++l) {
        int d5 = Ed5_new(i,l,i-1);
        int d3 = Ed3_new(i,l,l+1);

        int fm1 = V(i,l) + auPenalty(i,l) + d5 + d3 + Ec*(j-l) + Eb;

        min = MIN(min, fm1);
      }
      FM1[i][j] = min;
    }
  }

}

void calculate_fm() {

  for (int i = 1; i <= length; ++i) {
    for (int j = i+1; j <= length;++j) {
      int min1 = INFINITY_;
      for (int k = i+TURN+1; k <= j-TURN-1; ++k) {
        int x = FM[i][k-1] + FM1[k][j];
        min1 = MIN(min1, x);
      }
      int min2 = INFINITY_;
      for (int k = i; k <= j-TURN-1; ++k) {
        int x = FM1[k][j] + Ec*(k-i);
        min2 = MIN(min2, x);
      }
      FM[i][j] = MIN(min1, min2);
    }
  }
}

#endif

void process(ss_map_t& subopt_data, int len) {
        int count = 0;
        length = len;

#ifdef UNIQUE_MULTILOOP_DECOMPOSITION
        calculate_fm1();
        calculate_fm();
#endif

        ps_stack_t gstack;


        // initialize the partial structure, segment stack = {[1,n]}, label = W, list_bp = {} 
        ps_t first(0, len);
        first.push(segment(1, len, lW, W[len]));	
        gstack.push(first); // initialize the partial structure stacka

        while (1) {
                if (gstack.empty()) break; // exit
                ps_t ps = gstack.top();
                gstack.pop();

                if (ps.empty()) {
                        count++;
                        pair<ss_map_t::iterator, bool> ins_result;
                        ins_result = subopt_data.insert(std::make_pair<std::string,int>(ps.str,ps.ae_));
                        if (ins_result.second == false) {
                          printf("Duplicate Structure!!!");
                          exit(1);
                        }
                        //cout << ps.str << endl;
                        continue;
                }	
                else {
                       segment smt = ps.top();
                        ps.pop();

                        gflag = 0;
                        if (smt.j_ - smt.i_ > TURN) {
                                (*trace_func[smt.label_])(smt.i_, smt.j_, ps, gstack);
                        }

                        // discarded current segment, using remaining ones
                        if (!gflag) {
                                ps_t ps1(ps);
                                gstack.push(ps1);
                        }
                }
        }

#ifdef DEBUG 
        //printf("# SS = %d\n", count);
#endif
}

ss_map_t subopt_traceback(int len, int _delta) {
        trace_func[0] = traceW;
        trace_func[1] = traceV;
        trace_func[2] = traceVBI;
        trace_func[3] = traceVM;
        trace_func[4] = traceWM;
        trace_func[5] = traceWMPrime;
#ifdef UNIQUE_MULTILOOP_DECOMPOSITION
        trace_func[6] = traceM;
        trace_func[7] = traceM1;
#endif

        mfe = W[len];
        delta = _delta;
        length = len;

        ss_map_t subopt_data;
        process(subopt_data, len);

        return subopt_data;
}

void traceV(int i, int j, ps_t& ps, ps_stack_t& gstack) {

        // Hairpin Loop
        if (eH(i,j) + ps.total()  <= mfe + delta) {
                ps_t ps1(ps); 
                ps1.accumulate(eH(i,j));
                ps1.update(i, j, '(', ')');
                push_to_gstack(gstack, ps1);
        }

        // Stack
        if (eS(i, j) + V(i+1, j-1) + ps.total() <= mfe + delta) {
                ps_t ps1(ps);
                ps1.push(segment(i+1, j-1, lV, V(i+1, j-1)));
                ps1.accumulate(eS(i,j));
                ps1.update(i, j , '(', ')');
                push_to_gstack(gstack, ps1);
        }

        // Internal Loop
        if (VBI(i,j) + ps.total() <= mfe + delta) {
                traceVBI(i,j,ps,gstack);
        }

#ifdef UNIQUE_MULTILOOP_DECOMPOSITION
        int k;

        for (k = i+2; k <= j-TURN-1; ++k) {

          int kenergy1 = FM[i+1][k] + FM1[k+1][j-1];
          int d5 = Ed5(i, j, i+1);
          int d3 = Ed3(i, j, j-1);
          int aup = auPenalty(i,j);
          int kenergy2 = d5 + d3 + aup + Ea + Eb;
          int kenergy_total = kenergy1 + kenergy2;
          if (kenergy_total + ps.total() <= mfe + delta) {
            ps_t ps1(ps);
            ps1.push(segment(i+1,k, lM, FM[i+1][k]));
            ps1.push(segment(k+1,j-1, lM1, FM1[k+1][j-1]));
            ps1.accumulate(kenergy2);
            ps1.update(i,j,'(',')');
            push_to_gstack(gstack, ps1);
          }
        }
#else

        // Multiloop
        if (VM(i,j) + ps.total() <= mfe + delta) {
                ps_t ps1(ps);
                ps1.push(segment(i, j, lVM, VM(i,j)));
                ps1.update(i, j, '(', ')');
                push_to_gstack(gstack, ps1);
        }

#endif

}

void traceVBI(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        int p,q;

        for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
                int minq = j-i+p-MAXLOOP-2;
                if (minq < p+1+TURN) minq = p+1+TURN;
                int maxq = (p==(i+1))?(j-2):(j-1);
                for (q = minq; q <= maxq; q++) {
                        if (V(p, q) + eL(i, j, p, q) + ps.total() <= mfe + delta) {
                                ps_t ps1(ps);
                                ps1.push(segment(p, q, lV, V(p, q)));
                                ps1.update(i, j , '(', ')');
                                ps1.accumulate(eL(i, j, p, q));
                                push_to_gstack(gstack, ps1);
                        }
                }
        }
}

void traceW(int i, int j, ps_t& ps, ps_stack_t& gstack) {

        for (int l = i; l < j-TURN; ++l) {
                int wim1 =  MIN(0, W[l-1]);
                int d3 = (l>i)?Ed3(j,l,l-1):0;
                int d5 = (j<length)?Ed5(j,l,j+1):0;

                int Wij = V(l,j) + auPenalty(l, j) + d3 + d5 + wim1;
                if (Wij + ps.total() <= mfe + delta ) {
                        ps_t ps1(ps);
                        ps1.push(segment(l, j, lV, V(l,j)));
                        if (wim1 == W[l-1]) ps1.push(segment(i, l-1, lW, W[l-1]));
                        ps1.accumulate(auPenalty(l, j) + d3 + d5);
                        push_to_gstack(gstack, ps1);
                }
        }

        if (W[j-1] + ps.total() <= mfe + delta) {
                ps_t ps1(ps);
                ps1.push(segment(i, j-1, lW, W[j-1]));
                push_to_gstack(gstack, ps1);
        }
}

#ifdef UNIQUE_MULTILOOP_DECOMPOSITION
void traceM1(int i, int j, ps_t& ps, ps_stack_t& gstack) {

  if (FM1[i][j-1] + Ec + ps.total() <= mfe + delta) {
    ps_t ps1(ps);
    ps1.push(segment(i, j-1, lM1, FM1[i][j-1]));
    ps1.accumulate(Ec);
    push_to_gstack(gstack, ps1);
  }

  int d5;
  int d3;
  int aup;
  //d5 = (i != 1) ? Ed5(i,j,i-1) : Ed5(i,j,length);
  //d3 = (j != length) ? Ed3(i,j,j+1) : Ed3(i,j,1);
  d5 = Ed5_new(i, j, i-1);
  d3 = Ed3_new(i, j, j+1);
  aup = auPenalty(i,j);

  if (V(i,j) + d5 + d3 + aup + Eb + ps.total() <= mfe + delta) {
    ps_t ps1(ps);
    ps1.push(segment(i, j, lV, V(i,j)));
    ps1.accumulate(d5 + d3 + aup + Eb);
    ps1.update(i,j,'(',')');
    push_to_gstack(gstack, ps1);
  }
}

void traceM(int i, int j, ps_t& ps, ps_stack_t& gstack) {

  int d5, d3;
  int aup;

  if (FM[i][j-1] + Ec + ps.total() <= mfe + delta) {
    ps_t ps1(ps);
    ps1.push(segment(i, j-1, lM, FM[i][j-1]));
    ps1.accumulate(Ec);
    push_to_gstack(gstack, ps1);
  }

  //d5 = (i != 1) ? Ed3(i,j,i-1) : Ed3(i,j,length);
  //d3 = (j != length) ? Ed5(i,j,j+1) : Ed5(i,j,1);
  d5 = Ed5_new(i, j, i-1);
  d3 = Ed3_new(i, j, j+1);
  aup = auPenalty(i,j);
  if (V(i,j) + d5 + d3 + Eb + aup + ps.total() <= mfe + delta) {
    ps_t ps1(ps);
    ps1.push(segment(i, j, lV, V(i,j)));
    ps1.accumulate(d5 + d3 + Eb + aup);
    ps1.update(i,j,'(',')');
    push_to_gstack(gstack, ps1);
  }


  for (int k = i+TURN+1; k <= j-TURN-1; ++k) {

    d5 = Ed5_new(k+1, j, k);
    d3 = Ed3_new(k+1, j, j+1);
    aup = auPenalty(k+1, j);

    if (FM[i][k] + V(k+1,j) + d5 + d3
         + Eb + aup + ps.total() <= mfe + delta) {
      ps_t ps1(ps);
      ps1.push(segment(i, k, lM, FM[i][k]));
      ps1.push(segment(k+1, j, lV, V(k+1,j)));
      ps1.accumulate(d5 + d3 + Eb + aup);
      ps1.update(k+1,j,'(',')');
      push_to_gstack(gstack, ps1);
    }
  }

  for (int k = i; k <= j-TURN-1; ++k) {

    d5 = Ed5_new(k+1, j, k);
    d3 = Ed3_new(k+1, j, j+1);
    aup = auPenalty(k+1, j);

    if (V(k+1,j) + d5 + d3
         + Eb + Ec*(k-i+1) + aup + ps.total() <= mfe+delta) {
      ps_t ps1(ps);
      ps1.push(segment(k+1, j, lV, V(k+1,j)));
      ps1.accumulate(d5 + d3 + Eb + Ec*(k-i+1) + aup);
      ps1.update(k+1, j, '(', ')');
      push_to_gstack(gstack, ps1);
    }
  }
}

#endif

void traceWM(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        int d3 = (i==1)?Ed3(j,i,length):Ed3(j,i,i-1);
        int d5 = Ed5(j,i,j+1);

        if (V(i,j) + auPenalty(i, j) + Eb + d3 + d5 + ps.total() <= mfe + delta) {
                ps_t ps_new(ps);
                ps_new.accumulate(auPenalty(i, j) + Eb + d3 + d5);
                ps_new.push(segment(i,j, lV, V(i,j)));
                push_to_gstack(gstack, ps_new);
        }

        if (WMPrime[i][j] + ps.total() <= mfe + delta) {
                ps_t ps_new(ps);
                ps_new.push(segment(i,j, lWMPrime, WMPrime[i][j]));
                push_to_gstack(gstack, ps_new);
        }

        if (WM(i+1,j) + Ec + ps.total() <= mfe + delta) {
                ps_t ps_new(ps);
                ps_new.accumulate(Ec);
                ps_new.push(segment(i+1,j, lWM, WM(i+1,j)));
                push_to_gstack(gstack, ps_new);
        }

        if (WM(i,j-1) + Ec + ps.total() <= mfe + delta) {
                ps_t ps_new(ps);
                ps_new.accumulate(Ec);
                ps_new.push(segment(i,j-1, lWM, WM(i,j-1)));
                push_to_gstack(gstack, ps_new);
        }
}

void traceWMPrime(int i, int j, ps_t& ps, ps_stack_t& gstack) {
        for (int h = i+TURN+1 ; h <= j-TURN-2; h++) {
                if (WM(i,h-1) + WM(h,j) + ps.total() <= mfe + delta) {
                        ps_t ps_new(ps);
                        ps_new.push(segment(i,h-1, lWM, WM(i,h-1)));
                        ps_new.push(segment(h,j, lWM, WM(h,j)));
                        push_to_gstack(gstack, ps_new);
                }
        }
}

void traceVM(int i, int j, ps_t& ps, ps_stack_t& gstack) {

        int d3 = Ed3(i,j,j-1);
        int d5 = Ed5(i,j,i+1);

        if (WMPrime[i+1][j-1] + Ea + Eb + auPenalty(i, j) + d3 + d5 + ps.total() <= mfe + delta) {
                ps_t ps_new(ps);
                ps_new.accumulate(Ea + Eb + auPenalty(i, j) + d3 + d5); 
                ps_new.push(segment(i+1,j-1, lWMPrime,WMPrime[i+1][j-1] ));
                push_to_gstack(gstack, ps_new);
        }
}

void push_to_gstack(ps_stack_t& gstack, const ps_t& v) {
        gflag = 1;
        gstack.push(v);
}
