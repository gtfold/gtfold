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

#ifndef _ALGORITHMS_H
#define _ALGORITHMS_H

#ifdef __cplusplus
extern "C" {
#endif
	int calculate(int len);//, int nThreads, int unamode ,int t_mismatch);
#ifdef __cplusplus
}
#endif

void calcWM(int i, int j);
void calcW(int j);
void calcVBIVMVWM(int i, int j);
#endif
