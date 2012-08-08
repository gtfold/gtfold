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

/* Modified by Prashant Gaurav <pgaurav3@gatech.edu>, 09/07/2010 */
/* Fixed the incorrect reporting of multiloop energy */

#ifndef _TRACEBACK_H
#define _TRACEBACK_H

#ifdef __cplusplus
extern "C" {
#endif
	void trace(int len, int print_energy_decompose, const char* energy_decompose_output_file); //, int vv, int mode, int mismatch);
#ifdef __cplusplus
}
#endif

void traceW(int i);
int traceV(int i, int j);
int traceVM(int i, int j);
int traceVBI(int i, int j);
int traceWM(int i, int j);
int traceWMPrime(int i, int j);

#endif
