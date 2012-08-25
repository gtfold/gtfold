/*
GTfold: compute minimum free energy of RNA secondary structure
Copyright (C) 2008 David A. Bader
http://www.cc.gatech.edu/~bader
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
 */ 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include "constants.h"
#include "utils.h"
#include "energy.h"
#include "global.h"
#include "algorithms.h"
#include "constraints.h"
#include "shapereader.h"
#ifdef _OPENMP 
#include "omp.h"
#endif

//#define DEBUG 1

void initializeMatrix(int len) {
  int i, j;

  for (i = 1; i <= len; ++i) 
    for (j = len; j >= i; --j) 
      if (canPair(RNA[i],RNA[j]) && j-i > TURN) 
        PP[i][j]  = 1;
}

void prefilter(int len, int prefilter1, int prefilter2) {
  char** in;
  int i, j, k, count;

  in = (char**)malloc(len*sizeof(char*));
  for (i = 1; i <= len; ++i) in[i - 1] = (char*)malloc(len*sizeof(char));

  for (i = 1; i <= len - prefilter2 + 1; ++i)
    for (j = len; j >= prefilter2 && j >= i; --j) {
      count = 0;
      for (k = 0; k < prefilter2 && k <= (j - i) / 2; ++k)
        if (PP[i + k][j - k] == 1) ++count;
      if (count >= prefilter1)
        for (k = 0; k < prefilter2 && k <= (j - i) / 2; ++k)
          ++in[i + k - 1][j - k - 1];
    }

  for (i = 1; i <= len; ++i) {
    for (j = len; j >= i; --j)
      if (!in[i - 1][j - 1]) PP[i][j] = 0;
    free(in[i - 1]);
  }

  free(in);
}

int calcVBI(int i, int j) {
  int p=0, q=0;
  int VBIij = INFINITY_;

  for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
    int minq = j-i+p-MAXLOOP-2;
    if (minq < p+1+TURN) minq = p+1+TURN;
    int maxq = (p==(i+1))?(j-2):(j-1);

    for (q = minq; q <= maxq; q++) {
      if (PP[p][q]==0) continue;
      if (!canILoop(i,j,p,q)) continue;
      VBIij = MIN(eL(i, j, p, q) + V(p,q), VBIij);
    }
  }

  return VBIij;
}

int calcVBI1(int i, int j) {
  int p=0, q=0;
  int VBIij = INFINITY_;

  for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
    int minq = j-i+p-MAXLOOP-2;
    if (minq < p+1+TURN) minq = p+1+TURN;
    int maxq = (p==(i+1))?(j-2):(j-1);

    for (q = minq; q <= maxq; q++) {
      if (PP[p][q]==0) continue;
      if (!canILoop(i,j,p,q)) continue;
      VBIij = MIN(eL1(i, j, p, q) + V(p,q), VBIij);
    }
  }

  return VBIij;
}

int calcVBI2(int i, int j, int  len) {
  int d, ii, jj; 
  int energy = INFINITY_;

  for (d = j-i-3; d >= TURN+1 && d >= j-i-2-MAXLOOP; --d) 
    for (ii = i + 1; ii < j - d && ii <= len; ++ii)
    {    
      jj = d + ii;
      if (PP[ii][jj]==1)
        energy = MIN(energy, eL1(i, j, ii, jj) + V(ii, jj));
    }    

  return energy;
}

int calculate(int len) { 
  int b, i, j;
#ifdef _OPENMP
  if (g_nthreads > 0) omp_set_num_threads(g_nthreads);
#endif

#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
  //if (g_verbose) fprintf(stdout,"Thread count: %3d \n",omp_get_num_threads());
  fprintf(stdout,"Thread count: %3d \n",omp_get_num_threads());
#endif

  initializeMatrix(len);
  if (g_unamode || g_prefilter_mode) {
    prefilter(len,g_prefilter1,g_prefilter2);
  }

  for (b = TURN+1; b <= len-1; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
    for (i = 1; i <= len - b; i++) {
      j = i + b;

      if (PP[i][j] == 1) {
        int eh = canHairpin(i,j)?eH(i,j):INFINITY_; //hair pin
        int es = canStack(i,j)?eS(i,j)+V(i+1,j-1):INFINITY_; // stack

        // Internal Loop BEGIN
        if (g_unamode) 
          VBI(i,j) = calcVBI1(i,j);
        else
          VBI(i,j) = calcVBI(i,j);
        // Internal Loop END

        // Multi Loop BEGIN
        int d3 = canSS(j-1)?Ed3(i,j,j-1):INFINITY_;
        int d5 = canSS(i+1)?Ed5(i,j,i+1):INFINITY_;

        if (g_unamode || g_mismatch) { // unamode, terminal mismatch
          VM(i,j) = MIN(VM(i,j), WMPrime[i+1][j-1] + auPenalty(i,j) + Ea + Eb);
          VM(i,j) = MIN(VM(i,j), WMPrime[i+2][j-1] + d5 + auPenalty(i,j) + Ea + Eb + Ec);
          VM(i,j) = MIN(VM(i,j), WMPrime[i+1][j-2] + d3 + auPenalty(i,j) + Ea + Eb + Ec);
          VM(i,j) = MIN(VM(i,j), WMPrime[i+2][j-2] + Estackm(i,j) + auPenalty(i,j) + Ea + Eb + 2*Ec);
        } else if (g_dangles == 2) { // -d2
          VM(i,j) = MIN(VM(i,j), WMPrime[i+1][j-1] + d3 + d5 + auPenalty(i,j) + Ea + Eb);
        } else if (g_dangles == 0) { // -d0
          VM(i,j) = MIN(VM(i,j), WMPrime[i+1][j-1] + auPenalty(i,j) + Ea + Eb);
        }	else { // default 
          VM(i,j) = MIN(VM(i,j), WMPrime[i+1][j-1] + auPenalty(i,j) + Ea + Eb);
          VM(i,j) = MIN(VM(i,j), WMPrime[i+2][j-1] + d5 + auPenalty(i,j) + Ea + Eb + Ec);
          VM(i,j) = MIN(VM(i,j), WMPrime[i+1][j-2] + d3 + auPenalty(i,j) + Ea + Eb + Ec);
          VM(i,j) = MIN(VM(i,j), WMPrime[i+2][j-2] + d3 + d5 + auPenalty(i,j) + Ea + Eb + 2*Ec);
        }
        VM(i,j) = canStack(i,j)?VM(i,j):INFINITY_;
        // Multi Loop END

        V(i,j) = MIN4(  eh,
            es,
            VBI(i,j),
            VM(i,j));
      }
      else {
        V(i,j) = INFINITY_;
      }

      // Added auxillary storage WMPrime to speedup multiloop calculations
      int h;
      for (h = i+TURN+1 ; h <= j-TURN-2; h++) {
        WMPrime[i][j] = MIN(WMPrime[i][j], WMU(i,h-1) + WML(h,j)); 
      }

      // WM begin
      int newWM = INFINITY_; 

      //ZS: This sum corresponds to when i,j are NOT paired with each other.
      //So we need to make sure only terms where i,j aren't pairing are considered. 
      newWM = (!forcePair(i,j))?MIN(newWM, WMPrime[i][j]):newWM;

      if (g_unamode || g_mismatch) { // unamode
        newWM = MIN(V(i,j) + auPenalty(i,j) + Eb, newWM); 
        newWM = canSS(i)?MIN(V(i+1,j) + Ed3(j,i+1,i) + auPenalty(i+1,j) + Eb + Ec, newWM):newWM; //i dangle
        newWM = canSS(j)?MIN(V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) + Eb + Ec, newWM):newWM;  //j dangle
        if (i<j-TURN-2)
          newWM = (canSS(i)&&canSS(j))?MIN(V(i+1,j-1) + Estackm(j-1,i+1) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM):newWM; 
      } else if (g_dangles == 2) {
        int energy = V(i,j) + auPenalty(i,j) + Eb;
        energy += (i==1)?Ed3(j,i,len):Ed3(j,i,i-1);
        /*if (j<len)*/ energy += Ed5(j,i,j+1);
        newWM = (canSS(i)&&canSS(j))?MIN(energy, newWM):newWM; //i,j dangle
      } else if (g_dangles == 0) {
        newWM = MIN(V(i,j) + auPenalty(i,j) + Eb, newWM); 
      } else { // default
        newWM = MIN(V(i,j) + auPenalty(i,j) + Eb, newWM); 
        newWM = canSS(i)?MIN(V(i+1,j) + Ed3(j,i+1,i) + auPenalty(i+1,j) + Eb + Ec, newWM):newWM; //i dangle
        newWM = canSS(j)?MIN(V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) + Eb + Ec, newWM):newWM;  //j dangle
        newWM = (canSS(i)&&canSS(j))?MIN(V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM):newWM; //i,j dangle
      }
      newWM = canSS(i)?MIN(WMU(i+1,j) + Ec, newWM):newWM; //i dangle
      newWM = canSS(j)?MIN(WML(i,j-1) + Ec, newWM):newWM; //j dangle

      WMU(i,j) = WML(i,j) = newWM;
      // WM end
    }
  }

  W[0] = 0;
  for (j = 1; j <= len; j++) {
    int i, Wj, Widjd, Wijd, Widj, Wij, Wim1;
    Wj = 0;
    for (i = 1; i < j-TURN; i++) {
      Wij = Widjd = Wijd = Widj = INFINITY_;
      Wim1 = MIN(0, W[i-1]); 

      if (g_unamode || g_mismatch) { // unafold option
        Wij = V(i, j) + auPenalty(i, j) + Wim1;
        Widj = canSS(i)?V(i+1, j) + auPenalty(i+1,j) + Ed3(j,i + 1,i) + Wim1:Widj;
        Wijd = canSS(j)?V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + Wim1:Wijd;
        Widjd = (canSS(i)&&canSS(j))?V(i+1,j-1) + auPenalty(i+1,j-1) + Estacke(j-1,i+1) + Wim1:Widjd;
        Wij = MIN4(Wij, Widjd, Wijd, Widj);
      } else if (g_dangles == 2) { // -d2 option
        int energy = V(i,j) +	 auPenalty(i,j) + Wim1;
        if (i>1) energy +=  Ed3(j,i,i-1);
        if (j<len) energy += Ed5(j,i,j+1);
        Widjd = (canSS(i)&&canSS(j))? energy:Widjd;
        Wij = MIN(Wij, Widjd);
      }	else if (g_dangles == 0) { // -d0 option
        Wij = V(i, j) + auPenalty(i, j) + Wim1;
      } else { // default
        Wij = V(i, j) + auPenalty(i, j) + Wim1;
        Widj = canSS(i)?V(i+1, j) + auPenalty(i+1,j) + Ed3(j,i + 1,i) + Wim1:INFINITY_;
        Wijd = canSS(j)?V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + Wim1:INFINITY_;
        Widjd = (canSS(i)&&canSS(j))?V(i+1,j-1) + auPenalty(i+1,j-1) + Ed3(j-1,i + 1,i) + Ed5(j-1,i+1,j) + Wim1:INFINITY_;
        Wij = MIN4(Wij, Widjd, Wijd, Widj);
      }

      Wj = MIN(Wj,Wij); 
    }
    W[j] = canSS(j)?MIN(Wj, W[j-1]):Wj;
  }

#ifdef DEBUG
  FILE* file = fopen("VBI.txt", "w");
  int ii, jj;
  for (ii = 1; ii <= len; ++ii) {    
    for (jj = len; jj > ii; --jj) {
      fprintf(file, "%d %d %d\n",ii,jj,VBI(ii,jj));
    }
  }    
  fclose(file);

  file = fopen("Eh.txt", "w");
  for (ii = 1; ii <= len; ++ii) {    
    for (jj = len; jj > ii; --jj) {
      int eh = INFINITY_;
      if (PP[ii][jj])	eh = eH(ii,jj);
      fprintf(file, "%d %d %d\n",ii,jj,eh>=INFINITY_?INFINITY_:eh);
    }
  }    
  fclose(file);

  file = fopen("Es.txt", "w");
  for (ii = 1; ii <= len; ++ii) {    
    for (jj = len; jj > ii; --jj) {
      int es = INFINITY_;
      if (PP[ii][jj] && PP[ii+1][jj-1]) es = eS(ii,jj);
      fprintf(file, "%d %d %d\n",ii,jj,es>=INFINITY_?INFINITY_:es);
    }
  }    
  fclose(file);

  file = fopen("BP.txt", "w");
  for (ii = 1; ii <= len; ++ii) {    
    for (jj = len; jj > ii; --jj) {
      fprintf(file, "%d %d %d\n",ii,jj,PP[ii][jj]);
    }
  }    
  fclose(file);

  file = fopen("VM.txt", "w");
  for (ii = 1; ii <= len; ++ii) {    
    for (jj = len; jj > ii; --jj) {
      fprintf(file, "%d %d %d\n",ii,jj,VM(ii,jj));
    }
  }    
  fclose(file);

  file = fopen("WM.txt", "w");
  for (ii = 1; ii <= len; ++ii) {    
    for (jj = len; jj > ii; --jj) {
      fprintf(file, "%d %d %d\n",ii,jj,WM(ii,jj));
    }
  }    
  fclose(file);
#endif

  return W[len];
}
