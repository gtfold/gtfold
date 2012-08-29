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

/* Amrita: some of the variables defined in this file are not uesed. */

#ifndef _DATA_H
#define _DATA_H

#include "constants.h"

typedef struct{
	int poppen[5];
	int maxpen;
	int eparam[11];
	int multConst[3]; /* for multiloop penalties. */
	int dangle[4][4][4][2]; /* Contain dangling energy values */
	int inter[31]; /* Contains size penalty for internal loops */
	int bulge[31]; /* Contain the size penalty for bulges */
	int hairpin[31]; /* Contains the size penalty for hairpin loops */
	int stack[256]; /* Stacking energy used to calculate energy of stack loops */
	int tstkh[256]; /* Terminal mismatch energy used in the calculations of hairpin loops */
	int tstki[256]; /* Terminal mismatch energy used in the calculations of internal loops */
	int tloop[maxtloop + 1][2];
	int numoftloops;
	int iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
	int iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
	int iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */
	int coax[6][6][6][6];
	int tstackcoax[6][6][6][6];
	int coaxstack[6][6][6][6];
	int tstack[6][6][6][6];
	int tstkm[6][6][6][6];
	int auend;
	int gubonus;
	int cint; /* cint, cslope, c3 are used for poly C hairpin loops */
	int cslope;
	int c3;
	int efn2a;
	int efn2b;
	int efn2c;
	int triloop[maxtloop + 1][2];
	int numoftriloops;
	int init;
	int gail; /* It is either 0 or 1. It is used for grosely asymmetric internal loops */
	float prelog;
} nndb_constants;

#endif

