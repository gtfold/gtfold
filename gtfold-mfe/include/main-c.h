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

#ifndef _MAIN_C_H
#define _MAIN_C_H

extern int LENGTH;
extern unsigned char *RNA; /* [LENGTH] */
extern int *structure; /* [LENGTH] */
extern int *V; /* [(LENGTH-1)*(LENGTH)/2 + 1] */
extern int *W; /* [LENGTH] */
extern int **VBI; /* [LENGTH][LENGTH] */
extern int **VM; /* [LENGTH][LENGTH] */
extern int **WM; /* [LENGTH][LENGTH] */
extern int *indx; /* [LENGTH] */
extern int *constraints;
extern int num_threads;
extern int contact_dist;

enum BOOL {
	FALSE=0, TRUE
};


#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))

#endif
