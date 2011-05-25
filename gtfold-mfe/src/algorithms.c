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
int calculate(int len, int nThreads) { 
	int b, i, j;
#ifdef _OPENMP
	if (nThreads>0) omp_set_num_threads(nThreads);
#endif
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
	fprintf(stdout,"Thread count: %3d \n",omp_get_num_threads());
#endif


	for (b = TURN+1; b <= len-1; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
		for (i = 1; i <= len - b; i++) {
			j = i + b;
			int flag = 0, newWM = INFINITY_; 
			if (canPair(RNA[i], RNA[j])) {
				flag = 1;
				int eh = canHairpin(i,j)?eH(i,j):INFINITY_; //hair pin
				int es = canStack(i,j)?eS(i,j)+getShapeEnergy(i)+getShapeEnergy(j)+V(i+1,j-1):INFINITY_; // stack
				if (j-i > 6) { // Internal Loop BEGIN
					int p=0, q=0;
					int VBIij = INFINITY_;
					for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
						int minq = j-i+p-MAXLOOP-2;
						if (minq < p+1+TURN) minq = p+1+TURN;
						int maxq = (p==i+1)?(j-2):(j-1);
						for (q = minq; q <= maxq; q++) {
							if (!canPair(RNA[p], RNA[q])) continue;
							if (!canILoop(i,j,p,q)) continue;
							VBIij = MIN(eL(i, j, p, q) + V(p,q), VBIij);
						}
					}
					VBI(i,j) = VBIij;
					V(i,j) = V(i,j) + getShapeEnergy(i) + getShapeEnergy(j);

				} // Internal Loop END
				if (j-i > 10) { // Multi Loop BEGIN
					int h;
					int VMij, VMijd, VMidj, VMidjd;
					VMij = VMijd = VMidj = VMidjd = INFINITY_;
					for (h = i+TURN+1; h <= j-1-TURN; h++) { 
						VMij = MIN(VMij, WMU(i+1,h-1) + WML(h,j-1)); 
						VMidj = MIN(VMidj, WMU(i+2,h-1) + WML(h,j-1)); 
						VMijd = MIN(VMijd, WMU(i+1,h-1) + WML(h,j-2)); 
						VMidjd = MIN(VMidjd, WMU(i+2,h-1) + WML(h,j-2)); 
					}
					int d3 = canSS(j-1)?Ed3(i,j,j-1):INFINITY_;
					int d5 = canSS(i+1)?Ed5(i,j,i+1):INFINITY_;
					VMij = MIN(VMij, (VMidj + d5 +Ec)) ;
					VMij = MIN(VMij, (VMijd + d3 +Ec));
					VMij = MIN(VMij, (VMidjd + d5 + d3+ 2*Ec));
					VMij = VMij + Ea + Eb + auPenalty(i,j);
					VM(i,j) = canStack(i,j)?VMij:INFINITY_;
				} // Multi Loop END
				V(i,j) = MIN4(eh,es,VBI(i,j),VM(i,j));
			}
			else V(i,j) = INFINITY_;
			if (j-i > 4) { // WM BEGIN
				int h; 
				for (h = i+TURN+1 ; h <= j-TURN-1; h++) {
					//ZS: This sum corresponds to when i,j are NOT paired with each other.
					//So we need to make sure only terms where i,j aren't pairing are considered. 
					newWM = (!forcePair(i,j))?MIN(newWM, WMU(i,h-1) + WML(h,j)):newWM;
				}
				newWM = MIN(V(i,j) + auPenalty(i,j) + Eb, newWM); 
				newWM = canSS(i)?MIN(V(i+1,j) + Ed3(j,i+1,i) + auPenalty(i+1,j) + Eb + Ec, newWM):newWM; //i dangle
				newWM = canSS(j)?MIN(V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) + Eb + Ec, newWM):newWM;  //j dangle
				newWM = (canSS(i)&&canSS(j))?MIN(V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM):newWM; //i,j dangle
				newWM = canSS(i)?MIN(WMU(i+1,j) + Ec, newWM):newWM; //i dangle
				newWM = canSS(j)?MIN(WML(i,j-1) + Ec, newWM):newWM; //j dangle
				WMU(i,j) = WML(i,j) = newWM;
			} // WM END
		}
	}
	for (j = TURN+2; j <= len; j++) {
		int i, Wj, Widjd, Wijd, Widj, Wij, Wim1;
		Wj = INFINITY_;
		for (i = 1; i < j-TURN; i++) {
			Wij = Widjd = Wijd = Widj = INFINITY_;
			Wim1 = MIN(0, W[i-1]); 
			Wij = V(i, j) + auPenalty(i, j) + Wim1;
			Widjd = (canSS(i)&&canSS(j))?V(i+1,j-1) + auPenalty(i+1,j-1) + Ed3(j-1,i + 1,i) + Ed5(j-1,i+1,j) + Wim1:Widjd;
			Wijd = canSS(j)?V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + Wim1:Wijd;
			Widj = canSS(i)?V(i+1, j) + auPenalty(i+1,j) + Ed3(j,i + 1,i) + Wim1:Widj;
			Wj = MIN(MIN4(Wij, Widjd, Wijd, Widj), Wj); 
		}
		W[j] = canSS(j)?MIN(Wj, W[j-1]):Wj;
	}
	return W[len];
}

