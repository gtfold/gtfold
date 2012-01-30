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

#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#define maxfil 100    /* maximum length of file names */
#define INFINITY_ 9999999  /* an arbitrary value given to infinity */
#define SMALLINFTY_ 99999 /*an arbitray value to determine safe range from infinity */
#define maxtloop 100 /* maximum tetraloops allowed (info read from tloop) */
#define maxstructures 1010 /* maximum number of structures in ct file */
#define maxbases 10000   /* maximum number of bases in a structure */
#define ctheaderlength 125 /* maximum length of string containing info on sequence */
#define ga_bonus -10 /* the value of a bonus for the "almost coaxial stacking" case in efn2 */
#define amax 400 /* this is a maximum line length for void linout (below) */
#define col 80  /* this is the number of columns in an output file */
#define numlen 8  /* maximum digits in a number */
#define maxforce 600 /* maximum number of bases that can be forced single */
#define maxgu 5 /* maximum number of u's in gu pair */
#define C_ 1 /* "c" for optimized VBI. */

#define MAXLOOP 30 /* The maximum loop size. */
#define MAXENG 1000

#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_U 3

#endif
