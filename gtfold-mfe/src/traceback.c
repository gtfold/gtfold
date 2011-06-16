/*
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "data.h"
#include "constants.h"
#include "constraints.h"
#include "energy.h"
#include "global.h"
#include "traceback.h"
#include "utils.h"
#include "shapereader.h"

int verbose = -1;
int total_en = 0;
int total_ex = 0;
int t_mismatch = 0;

void trace(int len, int vbose, int t_mismatch_) {
	int i;
	verbose = vbose;
	t_mismatch = t_mismatch_;
	for (i = 0; i < len+1; i++)
		structure[i] = 0;

	if (W[len] >= MAXENG) {
		printf("\n No Structure ");
		return;
	}

	printf("\n");
	
	traceW(len);
	printf("- sum of energy of Loops:   	  %12.2f kcal/mol\n", total_en/100.0);
	printf("- sum of energy of External Loop: %12.2f kcal/mol\n", total_ex/100.0);
	return;
}

void traceW(int j) {
	int done, i, Wj;
	int wim1, flag, Widjd, Wijd, Widj, Wij;
	Wj = INFINITY_;
	flag = 1;
	done = 0; 
	
	if (j == 0 || j == 1) return;

	for (i = 1; i < j && !done; i++) {
		if (j-i < TURN) continue;

		wim1 = MIN(0, W[i-1]);
		flag = 1;
		if ( wim1 != W[i-1] && canSSregion(0,i)) flag = 0;

		Widjd = Wijd =  Widj = INFINITY_;
		Wij = V(i,j) + auPenalty(i, j) + wim1;
		
		if (t_mismatch) {
			Widjd = V(i+1,j-1) + auPenalty(i+1, j-1) + Estacke(j-1,i+1) + wim1;
		}
		else {
			Widjd = V(i+1,j-1) + auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + wim1;
		}

		Wijd = V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + wim1;
		Widj = V(i+1,j) + auPenalty(i+1,j) + Ed3(j,i+1,i) + wim1;

		if ((W[j] == Wij && canStack(i,j)) || forcePair(i,j)) { 
			done = 1;
			if (verbose == 1) 
				printf("i %5d j %5d ExtLoop   %12.2f\n", i, j, auPenalty(i, j)/100.00);
			total_ex += auPenalty(i, j);
			structure[i] = j;
			structure[j] = i;
			traceV(i, j);
			if (flag ) traceW(i - 1);
			break;
		} else if ((W[j] == Widjd && t_mismatch && canSS(i) && canSS(j) && canStack(i+1,j-1)) || forcePair(i+1,j-1)) { 
			done = 1;
			if (verbose == 1) 
				printf("i %5d j %5d ExtLoop   %12.2f\n", i+1, j-1, (auPenalty(i+1, j-1) + Estacke(j-1,i+1))/100.00);
			total_ex += (auPenalty(i+1, j-1) + Estacke(j-1,i+1));
			structure[i + 1] = j - 1;
			structure[j - 1] = i + 1;
			traceV(i + 1, j - 1);
			if (flag ) traceW(i - 1);
			break;
		}
	   	else if ((W[j] == Widjd && canSS(i) && canSS(j) && canStack(i+1,j-1)) || forcePair(i+1,j-1)) { 
			done = 1;
			if (verbose == 1) 
				printf("i %5d j %5d ExtLoop   %12.2f\n", i+1, j-1, (auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j))/100.00);
			total_ex += (auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j));
			structure[i + 1] = j - 1;
			structure[j - 1] = i + 1;
			traceV(i + 1, j - 1);
			if (flag ) traceW(i - 1);
			break;
		} else if ((W[j] == Wijd && canSS(j) && canStack(i,j-1)) || forcePair(i, j-1)) { 	
			done = 1;
			if (verbose == 1) 
				printf("i %5d j %5d ExtLoop   %12.2f\n", i, j-1, (auPenalty(i,j-1) + Ed5(j-1,i,j))/100.00);
			total_ex += (auPenalty(i,j-1) + Ed5(j-1,i,j));
			structure[i] = j - 1;
			structure[j - 1] = i;
			traceV(i, j - 1);
			if (flag ) traceW(i - 1);
			break;
		} else if ((W[j] == Widj && canSS(i) && canStack(i+1,j)) || forcePair(i+1,j)){
			done = 1;
			if (verbose == 1) 
				printf("i %5d j %5d ExtLoop   %12.2f\n", i+1, j, (auPenalty(i+1,j) + Ed3(j,i+1,i))/100.00);
			total_ex += (auPenalty(i+1,j) + Ed3(j,i+1,i));
			structure[i + 1] = j;
			structure[j] = i + 1;
			traceV(i + 1, j);
			if (flag ) traceW(i - 1);
			break;
		}
	}
		
	if (W[j] == W[j - 1] && !done) traceW(j-1);

	return;
}

int traceV(int i, int j) {
	int a, b, c, d, Vij;
	if (j-i < TURN)  return INFINITY_;

	a = canHairpin(i,j)?eH(i, j):INFINITY_;
	// TODO 
	// if (eS(i, j) == 0) b = INFINITY_;
	b = canStack(i,j)?eS(i, j) + V(i + 1, j - 1):INFINITY_;
	c = canStack(i,j)?VBI(i,j):INFINITY_;
	d = canStack(i,j)?VM(i,j):INFINITY_;
	
	Vij = V(i,j);

	if (Vij == a ) { 
		if (verbose == 1) 
			printf("i %5d j %5d Hairpin   %12.2f\n", i, j, eH(i, j)/100.00);
		total_en += eH(i,j);
		return Vij;
	} else if (Vij == b) { 
		if (verbose == 1) 
			printf("i %5d j %5d Stack     %12.2f\n", i, j, eS(i, j)/100.00);
		total_en += eS(i,j);
		structure[i + 1] = j - 1;
		structure[j - 1] = i + 1;
		traceV(i + 1, j - 1);
		return Vij;
	} else if (Vij == c) { 
		if (verbose == 1) 
			printf("i %5d j %5d IntLoop  ", i, j);
		traceVBI(i, j);
		return Vij;
	} else if (Vij == d) { 
		int eVM = traceVM(i, j);
		if (verbose ==1) 
			printf("i %5d j %5d MultiLoop %12.2f\n", i, j, (Vij-eVM)/100.0);
		total_en += (Vij-eVM);
		return Vij;
	}

	return 0;
}

int traceVBI(int i, int j) {
	
	int VBIij_temp;
	int ip, jp, el, v;
	int ifinal, jfinal;

	ifinal = 0;
	jfinal = 0;

	for (ip = i + 1; ip < j - 1; ip++) {
		for (jp = ip + 1; jp < j; jp++) {
			el = eL(i, j, ip, jp);
			v = V(ip, jp);
			VBIij_temp = el + v;
			if (VBIij_temp == VBI(i,j) || forcePair(ip,jp)){
				ifinal = ip;
				jfinal = jp;
				break;
			}
		}
		if (jp != j)
			break;
	}

	structure[ifinal] = jfinal;
	structure[jfinal] = ifinal;
	if (verbose==1) 
		printf(" %12.2f\n", eL(i, j, ifinal, jfinal)/100.00);
	total_en += eL(i, j, ifinal, jfinal);

	int eVI = traceV(ifinal, jfinal);
	return eVI ;
}

int traceVM(int i, int j) {
	int done;
	int h;
	int A_temp;
	int eVM = 0;

	done = 0;
	int VMij = VM(i,j);

	if (t_mismatch && i<j-TURN-2) {
		for (h = i + 3; h <= j - 2 && !done; h++) { 
			A_temp = WM(i + 2,h - 1) + WM(h,j - 2) + Ea + Eb + auPenalty(i,j) + Estackm(i,j);
			if (A_temp == VMij && canSS(i+1) && canSS(j-1)) {
				done = 1;
				eVM += traceWM(i + 2, h - 1);
				eVM += traceWM(h, j - 2);
				break;
			}
		}
	}
	else {
		for (h = i + 3; h <= j - 2 && !done; h++) { 
			A_temp = WM(i + 2,h - 1) + WM(h,j - 2) + Ea + Eb + auPenalty(i,j) + Ed5(i,j,i + 1) + Ed3(i,j,j - 1);
			if (A_temp == VMij && canSS(i+1) && canSS(j-1)) {
				done = 1;
				eVM += traceWM(i + 2, h - 1);
				eVM += traceWM(h, j - 2);
				break;
			}
		}
	}

	for (h = i + 2; h <= j - 1 && !done; h++) {
		A_temp = WM(i+1,h-1) + WM(h,j - 1) + Ea + Eb + auPenalty(i, j);
		if (A_temp == VMij) { 
			done = 1;
			eVM += traceWM(i + 1, h - 1);
			eVM += traceWM(h, j - 1);
			break;
		}
	}

	for (h = i + 3; h <= j - 1 && !done; h++) {
		A_temp = WM(i + 2,h - 1) + WM(h,j - 1) + Ea + Eb + auPenalty(i,j) + Ed5(i,j,i + 1); 
		if (A_temp == VMij && canSS(i+1) ) {
			done = 1;
			eVM += traceWM(i + 2, h - 1);
			eVM += traceWM(h, j - 1);
			break;
		}
	}

	for (h = i + 2; h <= j - 2 && !done; h++) { 
		A_temp = WM(i + 1,h - 1) + WM(h,j - 2) + Ea + Eb + auPenalty(i, j) + Ed3(i,j,j - 1);
		if (A_temp == VMij && canSS(j-1)) {
			done = 1;
			eVM += traceWM(i + 1, h - 1);
			eVM += traceWM(h, j - 2);
			break;
		}
	}

		return eVM;
}

int traceWM(int i, int j) {

	int done;
	int h1, h;
	int eWM = 0; 

	done = 0;
	h1 = 0;

	if (i >= j)
		return 0;
	else {
		for (h = i; h < j && !done; h++) {
			int aa = WM(i,h) + WM(h + 1,j); 
			if (aa == WM(i,j)) {
				done = 1;
				h1 = h;
				break;
			}
		}
		if (h1 != 0) {
			eWM += traceWM(i, h);
			eWM += traceWM(h + 1, j);
			done = 1;
		} else {
			if (WM(i,j) == V(i,j) + auPenalty(i, j) + Eb && canStack(i,j)) { 
				structure[i] = j;
				structure[j] = i;
				eWM += traceV(i, j);
				done = 1;
			} else if (WM(i,j) == V(i+1, j) + Ed3(j,i + 1,i) + auPenalty(i+1, j) + Eb + Ec && canSS(i) &&  canStack(i+1,j)) { 
				eWM += traceV(i + 1, j);
				structure[i + 1] = j;
				structure[j] = i + 1;
				done = 1;
			} else if (WM(i,j) == V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) +  Eb + Ec && canSS(j) && canStack(i,j-1)) { 
				done = 1;
				eWM += traceV(i, j - 1);
				structure[i] = j - 1;
				structure[j - 1] = i;
			} else if (WM(i,j) == V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1, j-1) + Eb + 2*Ec && canSS(i) && canSS(j) && canStack(i+1,j-1)) { 
				done = 1;
				eWM += traceV(i + 1, j - 1);
				structure[i + 1] = j - 1;
				structure[j - 1] = i + 1;
			}
		   	else if (WM(i,j) == V(i+1,j-1) + Estackm(j-1,i+1) + auPenalty(i+1, j-1) + Eb + 2*Ec && canSS(i) && canSS(j) && canStack(i+1,j-1) && t_mismatch) {
				done = 1;
				eWM += traceV(i + 1, j - 1);
				structure[i + 1] = j - 1;
				structure[j - 1] = i + 1;
			}
			else if (WM(i,j) == WM(i + 1,j) + Ec && canSS(i)) { 
				done = 1;
				eWM += traceWM(i + 1, j);
			} else if (WM(i,j) == WM(i,j - 1) + Ec && canSS(j)) { 
				done = 1;
				eWM += traceWM(i, j - 1);
			}
		}
	}
	if(!done) { 
		fprintf(stderr, "ERROR: WM(%d,%d) could not be traced!\n", i,j);
	}
	return eWM;
}
