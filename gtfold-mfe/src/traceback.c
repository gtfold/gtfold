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

int total_en = 0;
int total_ex = 0;
int length = 0;
int print_energy_decompose = 0;
FILE* energy_decompose_outfile;

void trace(int len, int print_energy_decompose1, const char* energy_decompose_output_file) { 
	print_energy_decompose = print_energy_decompose1;
        if(print_energy_decompose==1){
		energy_decompose_outfile = fopen(energy_decompose_output_file, "w");
                if(energy_decompose_outfile==NULL){
                        printf("Error in opening file: %s\n", energy_decompose_output_file);
                        exit(-1);
                }
                printf("\nEnergy decomposition for MFE structure will be saved to %s\n", energy_decompose_output_file);
	}
	
	int i;
	for (i = 0; i <= len; i++) structure[i] = 0;

	length = len;
	if (W[len] >= MAXENG) {
		printf("- No Structure \n");
		return;
	}

	printf("\n");
	
	traceW(len);
	if (print_energy_decompose == 1) {
		fprintf(energy_decompose_outfile, "- sum of energy of Loops:   	  %12.2f kcal/mol\n", total_en/100.0);
		fprintf(energy_decompose_outfile, "- sum of energy of External Loop: %12.2f kcal/mol\n", total_ex/100.0);
  	}
	if(print_energy_decompose==1){
		fclose(energy_decompose_outfile);
	}
  	return;
}

void traceW(int j) {
	int done = 0, i;
	int wim1, flag = 1 ;
	
  if (j == 0 || j == 1)
    return;
	
	for (i = 1; i < j && !done; i++) {
		if (j-i < TURN) continue;

		wim1 = MIN(0, W[i-1]);
		flag = 1;
		if ( wim1 != W[i-1] && canSSregion(0,i)) flag = 0;

		if (g_unamode||g_mismatch) {
			if ((W[j] == V(i,j) + auPenalty(i, j) + wim1 && canStack(i,j)) || forcePair(i,j)) { 
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i, j, auPenalty(i, j)/100.00);
							total_ex += auPenalty(i, j);
							traceV(i, j);
							if (flag ) traceW(i - 1);
							break;
			} else if ((W[j] ==  V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + wim1 && canSS(j) && canStack(i,j-1)) || forcePair(i, j-1)) { 	
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i, j-1, (auPenalty(i,j-1) + Ed5(j-1,i,j))/100.00);
							total_ex += (auPenalty(i,j-1) + Ed5(j-1,i,j));
							traceV(i, j - 1);
							if (flag ) traceW(i - 1);
							break;
			} else if ((W[j] == V(i+1,j) + auPenalty(i+1,j) + Ed3(j,i+1,i) + wim1 && canSS(i) && canStack(i+1,j)) || forcePair(i+1,j)){
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i+1, j, (auPenalty(i+1,j) + Ed3(j,i+1,i))/100.00);
							total_ex += (auPenalty(i+1,j) + Ed3(j,i+1,i));
							traceV(i + 1, j);
							if (flag ) traceW(i - 1);
							break;
			} else if ((W[j] == V(i+1,j-1) + auPenalty(i+1, j-1) + Estacke(j-1,i+1) + wim1 && canSS(i) && canSS(j) && canStack(i+1,j-1)) || forcePair(i+1,j-1)) { 
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i+1, j-1, (auPenalty(i+1, j-1) + Estacke(j-1,i+1))/100.00);
							total_ex += (auPenalty(i+1, j-1) + Estacke(j-1,i+1));
							traceV(i + 1, j - 1);
							if (flag ) traceW(i - 1);
							break;
			}
		} else if (g_dangles == 2) {
				int e_dangles = 0;
				if (i>1) e_dangles +=  Ed3(j,i,i-1);
				if (j<length) e_dangles += Ed5(j,i,j+1);
				if ((W[j] == V(i,j) + auPenalty(i, j) + e_dangles + wim1 && canSS(i) && canSS(j) && canStack(i+1,j-1)) || forcePair(i+1,j-1)) { 
												done = 1;
												if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i+1, j-1, (auPenalty(i, j) + e_dangles)/100.00);
												total_ex += (auPenalty(i, j) + e_dangles);
												traceV(i, j);
												if (flag ) traceW(i - 1);
												break;
				} 
		}	else if (g_dangles == 0) {
			if ((W[j] == V(i,j) + auPenalty(i, j) + wim1 && canStack(i,j)) || forcePair(i,j)) { 
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i, j, auPenalty(i, j)/100.00);
							total_ex += auPenalty(i, j);
							traceV(i, j);
							if (flag ) traceW(i - 1);
							break;
			}
		} else { // default
			if ((W[j] == V(i,j) + auPenalty(i, j) + wim1 && canStack(i,j)) || forcePair(i,j)) { 
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i, j, auPenalty(i, j)/100.00);
							total_ex += auPenalty(i, j);
							traceV(i, j);
							if (flag ) traceW(i - 1);
							break;
			} else if ((W[j] ==  V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + wim1 && canSS(j) && canStack(i,j-1)) || forcePair(i, j-1)) { 	
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i, j-1, (auPenalty(i,j-1) + Ed5(j-1,i,j))/100.00);
							total_ex += (auPenalty(i,j-1) + Ed5(j-1,i,j));
							traceV(i, j - 1);
							if (flag ) traceW(i - 1);
							break;
			} else if ((W[j] == V(i+1,j) + auPenalty(i+1,j) + Ed3(j,i+1,i) + wim1 && canSS(i) && canStack(i+1,j)) || forcePair(i+1,j)){
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i+1, j, (auPenalty(i+1,j) + Ed3(j,i+1,i))/100.00);
							total_ex += (auPenalty(i+1,j) + Ed3(j,i+1,i));
							traceV(i + 1, j);
							if (flag ) traceW(i - 1);
							break;
			} else if ((W[j] == V(i+1,j-1) + auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + wim1 && canSS(i) && canSS(j) && canStack(i+1,j-1)) || forcePair(i+1,j-1)) { 
							done = 1;
							if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d ExtLoop   %12.2f\n", i+1, j-1, (auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j))/100.00);
							total_ex += (auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j));
							traceV(i + 1, j - 1);
							if (flag ) traceW(i - 1);
							break;
			} 
		}
	}
		
	if (W[j] == W[j - 1] && !done) traceW(j-1);
 
	//if (!done) {
	//	fprintf(stderr, "ERROR: W(%d) could not be traced\n", j);
	//}

	return;
}

int traceV(int i, int j) {
	int a, b, c, d, Vij;
	if (j-i < TURN)  return INFINITY_;

	a = canHairpin(i,j)?eH(i, j):INFINITY_;
	b = canStack(i,j)?eS(i, j) + V(i + 1, j - 1):INFINITY_;
	c = canStack(i,j)?VBI(i,j):INFINITY_;
	d = canStack(i,j)?VM(i,j):INFINITY_;
	
	Vij = V(i,j);
	structure[i] = j;
	structure[j] = i;

	if (Vij == a ) { 
		if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d Hairpin   %12.2f\n", i, j, eH(i, j)/100.00);
		total_en += eH(i,j);
		return Vij;
	} else if (Vij == b) { 
		if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d Stack     %12.2f\n", i, j, eS(i, j)/100.00);
		total_en += eS(i,j);
		traceV(i + 1, j - 1);
		return Vij;
	} else if (Vij == c) { 
		if (print_energy_decompose == 1) fprintf(energy_decompose_outfile, "i %5d j %5d IntLoop  ", i, j);
		traceVBI(i, j);
		return Vij;
	} else if (Vij == d) { 
		int eVM = traceVM(i, j);
		if (print_energy_decompose ==1) fprintf(energy_decompose_outfile, "i %5d j %5d MultiLoop %12.2f\n", i, j, (Vij-eVM)/100.0);
		total_en += (Vij-eVM);
		return Vij;
	}

	return 0;
}

int traceVBI(int i, int j) {
	int VBIij;
	int ip, jp;
	int ifinal, jfinal;

	ifinal = 0;
	jfinal = 0;

	for (ip = i + 1; ip < j - 1; ip++) {
		for (jp = ip + 1; jp < j; jp++) {
			VBIij = eL(i, j, ip, jp)+ V(ip,jp);
			if (VBIij == VBI(i,j) || forcePair(ip,jp)){
				ifinal = ip;
				jfinal = jp;
				break;
			}
		}
		if (jp != j) break;
	}

	if (print_energy_decompose==1) fprintf(energy_decompose_outfile, " %12.2f\n", eL(i, j, ifinal, jfinal)/100.00);
	total_en += eL(i, j, ifinal, jfinal);

	return traceV(ifinal, jfinal);
}

int traceVM(int i, int j) {
	int done = 0;
	int eVM = 0;

	if (g_unamode||g_mismatch) {
		if (VM(i,j) == WMPrime[i+1][j - 1] + Ea + Eb + auPenalty(i, j) ) {
			done = 1;
			eVM += traceWMPrime(i + 1, j - 1);
		} else if (VM(i,j) == WMPrime[i + 2][j - 1] + Ea + Eb + auPenalty(i,j) + Ed5(i,j,i + 1) + Ec && canSS(i+1) ) {
			done = 1;
			eVM += traceWMPrime(i + 2, j - 1);
		}
		else if ( VM(i,j) == WMPrime[i + 1][j - 2] + Ea + Eb + auPenalty(i, j) + Ed3(i,j,j - 1) + Ec && canSS(j-1)) {
			done = 1;
			eVM += traceWMPrime(i + 1, j - 2);
		}
		else if (V(i,j) == WMPrime[i + 2][j - 2] + Ea + Eb + auPenalty(i,j) + Estackm(i,j) + 2*Ec && canSS(i+1) && canSS(j-1) ) {
			done = 1;
			eVM += traceWMPrime(i + 2, j - 2);
		}
	} else if (g_dangles == 2) {
		if (V(i,j) ==  WMPrime[i + 1][j - 1] + Ea + Eb + auPenalty(i,j) + Ed5(i,j,i + 1) + Ed3(i,j,j - 1) && canSS(i+1) && canSS(j-1) ) {
			done = 1;
			eVM += traceWMPrime(i + 1, j - 1);
		}
	}	else if (g_dangles == 0) {
		if (VM(i,j) == WMPrime[i+1][j - 1] + Ea + Eb + auPenalty(i, j) ) {
			done = 1;
			eVM += traceWMPrime(i + 1, j - 1);
		}
	} else {
		if (VM(i,j) == WMPrime[i+1][j - 1] + Ea + Eb + auPenalty(i, j) ) {
			done = 1;
			eVM += traceWMPrime(i + 1, j - 1);
		} else if (VM(i,j) == WMPrime[i + 2][j - 1] + Ea + Eb + auPenalty(i,j) + Ed5(i,j,i + 1) + Ec && canSS(i+1) ) {
			done = 1;
			eVM += traceWMPrime(i + 2, j - 1);
		}
		else if ( VM(i,j) == WMPrime[i + 1][j - 2] + Ea + Eb + auPenalty(i, j) + Ed3(i,j,j - 1) + Ec && canSS(j-1)) {
			done = 1;
			eVM += traceWMPrime(i + 1, j - 2);
		} else if (V(i,j) ==  WMPrime[i + 2][j - 2] + Ea + Eb + auPenalty(i,j) + Ed5(i,j,i + 1) + Ed3(i,j,j - 1) + 2*Ec && canSS(i+1) && canSS(j-1) ) {
			done = 1;
			eVM += traceWMPrime(i + 2, j - 2);
		}
	}
	if(!done) { 
		fprintf(stderr, "ERROR: VM(%d,%d) could not be traced\n", i,j);
	}

	return eVM;
}

int traceWMPrime(int i, int j) {
	int done=0, h, energy=0;
	
	for (h = i; h < j && !done; h++) {
			if (WM(i,h) + WM(h + 1,j) == WMPrime(i,j)) {
				energy += traceWM(i, h);
				energy += traceWM(h + 1, j);
				done = 1;
		 	  break;
		}
	}
	return energy;
}

int traceWM(int i, int j) {
	assert(i < j);
	int done=0, eWM=0;

	if (!done && WM(i,j) == WMPrime[i][j]) {
			eWM += traceWMPrime(i,j);
			done = 1;
	}

	if (!done){
		if (g_unamode||g_mismatch) {
						if (WM(i,j) == V(i,j) + auPenalty(i, j) + Eb && canStack(i,j)) { 
										eWM += traceV(i, j);
										done = 1;
						} else if (WM(i,j) == V(i+1, j) + Ed3(j,i + 1,i) + auPenalty(i+1, j) + Eb + Ec && canSS(i) &&  canStack(i+1,j)) { 
										eWM += traceV(i + 1, j);
										done = 1;
						} else if (WM(i,j) == V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) +  Eb + Ec && canSS(j) && canStack(i,j-1)) { 
										done = 1;
										eWM += traceV(i, j - 1);
						} else if (WM(i,j) == V(i+1,j-1) + Estackm(j-1,i+1) + auPenalty(i+1, j-1) + Eb + 2*Ec && canSS(i) && canSS(j) && canStack(i+1,j-1)) {
										done = 1;
										eWM += traceV(i + 1, j - 1);
						}
		} else if (g_dangles == 2) {
						int energy = V(i,j) + auPenalty(i, j) + Eb;				
						energy += (i==1)?Ed3(j,i,length):Ed3(j,i,i-1);
						/*if (j<len)*/ energy += Ed5(j,i,j+1);
						if (WM(i,j) ==  energy && canSS(i) && canSS(j) && canStack(i+1,j-1)) { 
										eWM += traceV(i, j);
										done = 1;
						}
		} else if (g_dangles == 0) {
						if (WM(i,j) == V(i,j) + auPenalty(i, j) + Eb && canStack(i,j)) { 
										eWM += traceV(i, j);
										done = 1;
						}
		} else  {
						if (WM(i,j) == V(i,j) + auPenalty(i, j) + Eb && canStack(i,j)) { 
										eWM += traceV(i, j);
										done = 1;
						} else if (WM(i,j) == V(i+1, j) + Ed3(j,i + 1,i) + auPenalty(i+1, j) + Eb + Ec && canSS(i) &&  canStack(i+1,j)) { 
										eWM += traceV(i + 1, j);
										done = 1;
						} else if (WM(i,j) == V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) +  Eb + Ec && canSS(j) && canStack(i,j-1)) { 
										eWM += traceV(i, j - 1);
										done = 1;
						} else if (WM(i,j) == V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1, j-1) + Eb + 2*Ec && canSS(i) && canSS(j) && canStack(i+1,j-1)) { 
										eWM += traceV(i + 1, j - 1);
										done = 1;
						}
		}
	}
	
	if (!done){
		if (WM(i,j) == WM(i + 1,j) + Ec && canSS(i)) { 
							done = 1;
							eWM += traceWM(i + 1, j);
		} else if (WM(i,j) == WM(i,j - 1) + Ec && canSS(j)) { 
							done = 1;
							eWM += traceWM(i, j - 1);
		}
	}

	if(!done) { 
		fprintf(stderr, "ERROR: WM(%d,%d) could not be traced\n", i,j);
	}
	return eWM;
}
