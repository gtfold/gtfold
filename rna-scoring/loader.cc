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

/* Modified by Amrita Mathuriya August 2007 - January 2009. 1) Functions are rewritten, removing calls of boost library functions and corrected bugs.
 * This file defines various functions to read the thermodynamic parameters from files in the data directory.
 * */
/* Modified by Sonny Hernandez May 2007 - Aug 2007. All comments added marked by "SH: "*/
/* Modified by Sainath Mallidi August 2009 -  "*/

#include <string.h>
#include <stdlib.h>
#include "loader.h"

#define xstr(s) str(s)
#define str(s) #s

#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))


using namespace std;


string EN_DATADIR;

nndb_constants * populate(const char *userdatadir,bool userdatalogic) {

	nndb_constants * ret_pointer = (nndb_constants *)malloc(sizeof(nndb_constants));
	if(ret_pointer == NULL){
		fputs("Failed to allocate constants", stderr);
		exit(-1);
	}

	cout << "Loading in GTfold data files from ";

	if (!userdatalogic) {
		EN_DATADIR.assign(xstr(DATADIR));
		EN_DATADIR += "/";
		EN_DATADIR += userdatadir;
	} else {
		EN_DATADIR.assign(userdatadir);
	}

	//Handle the ending forward slash case
	if (EN_DATADIR[EN_DATADIR.length() - 1] != '/') {
		EN_DATADIR += "/";
	}

	cout << EN_DATADIR << endl;

	initMiscloopValues("miscloop.dat", ret_pointer);
	// miscloop.dat - Miscellaneous loop file
	initStackValues("stack.dat", ret_pointer);
	// stack.dat - free energies for base pair stacking
	initDangleValues("dangle.dat", ret_pointer);
	// dangle.dat - single base stacking free energies
	initLoopValues("loop.dat", ret_pointer);
	// loop.dat - entropic component for internal, bulge and hairpin loops.
	initTstkhValues("tstackh.dat", ret_pointer);
	// tstackh.dat - free energies for terminal mismatch stacking in hairpin loops
	initTstkiValues("tstacki.dat", ret_pointer);
	// tstacki.dat - free energies for terminal mismatch stacking in internal loops
	initTloopValues("tloop.dat", ret_pointer);
	// tloop.dat - free energies for distinguished tetraloops
	initInt21Values("int21.dat", ret_pointer);
	// int21.dat - free energies for 2 x 1 interior loops
	initInt22Values("int22.dat", ret_pointer);
	// int22.dat - free energies for 2 x 2 interior loops
	initInt11Values("int11.dat", ret_pointer);
	// int11.dat - free energies for 1 x 1 interior loops

	cout << "Done loading data files." << endl << endl;
	return ret_pointer;
}

char baseToDigit(std::string base) {
	if (!strcmp(base.c_str(), "A")) {
		return '1';
	}
	if (!strcmp(base.c_str(), "C")) {
		return '2';
	}
	if (!strcmp(base.c_str(), "G")) {
		return '3';
	}
	if (!strcmp(base.c_str(), "U")) {
		return '4';
	}
	if (!strcmp(base.c_str(), "N")) {
		return '5';
	}
	return (char) NULL;
}

unsigned char getBase(std::string base) {
	//cout << base;
	if (!strcmp(base.c_str(), "A") || !strcmp(base.c_str(), "a")) {
		//cout << "1";
		return BASE_A;
	}
	if (!strcmp(base.c_str(), "C") || !strcmp(base.c_str(), "c")) {
		//cout << "2";
		return BASE_C;
	}
	if (!strcmp(base.c_str(), "G") || !strcmp(base.c_str(), "g")) {
		//cout << "3";
		return BASE_G;
	}
	if (!strcmp(base.c_str(), "U") || !strcmp(base.c_str(), "u") || !strcmp(
			base.c_str(), "T") || !strcmp(base.c_str(), "t")) {
		//cout << "4";
		return BASE_U;
	}
	if (!strcmp(base.c_str(), "N") || !strcmp(base.c_str(), "n")||!strcmp(base.c_str(), "R")|| !strcmp(base.c_str(), "r")|| !strcmp(base.c_str(), "Y")|| !strcmp(base.c_str(), "y")|| !strcmp(base.c_str(), "M")|| !strcmp(base.c_str(), "m")|| !strcmp(base.c_str(), "K")|| !strcmp(base.c_str(), "k")|| !strcmp(base.c_str(), "S")|| !strcmp(base.c_str(), "s")|| !strcmp(base.c_str(), "W")|| !strcmp(base.c_str(), "w")|| !strcmp(base.c_str(), "B")|| !strcmp(base.c_str(), "b")|| !strcmp(base.c_str(), "D")|| !strcmp(base.c_str(), "d")|| !strcmp(base.c_str(), "H")|| !strcmp(base.c_str(), "h")|| !strcmp(base.c_str(), "V")|| !strcmp(base.c_str(), "v")) {
		//cout << "4";
		return BASE_A;
	}
	return 'X';
}
unsigned char getBase1(std::string base) {
	//cout << base;
	if (!strcmp(base.c_str(), "A") || !strcmp(base.c_str(), "a")) {
		//cout << "1";
		return BASE_A;
	}
	if (!strcmp(base.c_str(), "C") || !strcmp(base.c_str(), "c")) {
		//cout << "2";
		return BASE_C;
	}
	if (!strcmp(base.c_str(), "G") || !strcmp(base.c_str(), "g")) {
		//cout << "3";
		return BASE_G;
	}
	if (!strcmp(base.c_str(), "U") || !strcmp(base.c_str(), "u") || !strcmp(
			base.c_str(), "T") || !strcmp(base.c_str(), "t")) {
		//cout << "4";
		return BASE_U;
	}
	if (!strcmp(base.c_str(), "N") || !strcmp(base.c_str(), "n")|| !strcmp(base.c_str(), "R")|| !strcmp(base.c_str(), "r")|| !strcmp(base.c_str(), "Y")|| !strcmp(base.c_str(), "y")|| !strcmp(base.c_str(), "M")|| !strcmp(base.c_str(), "m")|| !strcmp(base.c_str(), "K")|| !strcmp(base.c_str(), "k")|| !strcmp(base.c_str(), "S")|| !strcmp(base.c_str(), "s")|| !strcmp(base.c_str(), "W")|| !strcmp(base.c_str(), "w")|| !strcmp(base.c_str(), "B")|| !strcmp(base.c_str(), "b")|| !strcmp(base.c_str(), "D")|| !strcmp(base.c_str(), "d")|| !strcmp(base.c_str(), "H")|| !strcmp(base.c_str(), "h")|| !strcmp(base.c_str(), "V")|| !strcmp(base.c_str(), "v")) {
		//cout << "4";
		return 'N';
	}
	return 'X';
}

int initStackValues(string fileName, nndb_constants * container) {

	ifstream cf; //cf = current file
	int i, j, k, l;
	int ii, jj, kk, ll;
	int index;
	char currentLine[256];
	string currentString;
	string s;

	// Initialize the array with INFINITY
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 4; l++) {
					container->stack[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// Read the thermodynamic parameters.
	// The 24 first lines are junk we don't need
	for (index = 1; index <= 15; index++) {
		cf.getline(currentLine, 256);
	}

	i = 0;
	kk = 0;
	ii = 0;
	jj = 0;
	ll = 0;

	while (i < 16) {

		if (i % 4 == 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);

		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;

		ll = 0;
		jj = 0;

		int z = 0;
		int r = 0;

		while (s[z] != '\0') {

			if (s[z] == ' ')
				z++;

			else if (s[z] == '.') {
				z++;
				ll++;
				if (ll == 4)
					ll = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			} else {
				char value[10];
				int x = 0;

				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}

				value[x] = '\0';

				int temp = (int) floor(100.0 * atof(value) + .5);
				container->stack[fourBaseIndex(ii,jj,kk,ll)] = temp;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				ll++;
				if (ll == 4)
					ll = 0;
			}
		}
		i++;
		if (!(i % 4))
			ii++;

		/*    jj = 1; */
		kk = (i % 4);

	}
	cf.close();
#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}

int initMiscloopValues(string fileName, nndb_constants * container) {
	/*
	 miscloop.dat - Miscellaneous loop file. Contains :
	 1. Extrapolation for large loops based on polymer theory
	 2. Asymmetric internal loop correction parameters.
	 3. the f(m) array (see Ninio for details)
	 4. Paremeters for multibranch loops
	 5. Paremeters for multibranch loops (for efn2 only)
	 6. Terminal AU or GU penalty
	 7. Bonus for GGG hairpin
	 8,9,10. C hairpin rules: a) slope  b) intercept c) value for size 3
	 11. Intermolecular initiation free energy
	 12. GAIL Rule (Grossly Asymmetric Interior Loop Rule) (on or off)
	 */

	char currentWord[256];
	string s;
	ifstream cf; //cf = current file

	fileName = EN_DATADIR + fileName;
#if 1
	cout << "Getting miscloop values from " << fileName << endl;
#endif
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}
	s = "";

	cf >> currentWord;
	for (int index = 1; index < 13; index++) { // There are total 12 values to read in.
		while (strcmp(currentWord, "-->")) {
			cf >> currentWord;
		}
		if (index == 1) {
			cf >> currentWord;
			container->prelog = 100 * atof(currentWord);
			cout << "prelog = " << container->prelog << endl;
		}
		if (index == 2) {
			cf >> currentWord;
			container->maxpen = int(atof(currentWord) * 100.0 + .5);
			cout << "maxpen = " << container->maxpen << endl;
		}
		if (index == 3) {
			for (int count = 1; count <= 4; count++) {
				cf >> currentWord;
				s = currentWord;
				container->poppen[count] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
		}
		if (index == 4) {
			container->eparam[1] = 0;
			container->eparam[2] = 0;
			container->eparam[3] = 0;
			container->eparam[4] = 0;
			container->eparam[7] = 30;
			container->eparam[8] = 30;
			container->eparam[9] = -500;
			int table[4];
			table[1] = 5;
			table[2] = 6;
			table[3] = 10;
			for (int count = 1; count <= 3; count++) {
				cf >> currentWord;
				s = currentWord;
				container->multConst[count - 1] = (int) (atof(s.c_str()) * 100 + 0.5);
				container->eparam[table[count]] = (int) (atof(s.c_str()) * 100 + 0.5);
				//cout << "\n multi " << multConst[count-1];
			}
		}
		if (index == 5) {
			int table[4];
			for (int count = 1; count <= 3; count++) {
				cf >> currentWord;
				s = currentWord;
				table[count] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
			/*
			efn2a = table[1];
			efn2b = table[2] - 1;
			efn2c = table[3] - 1;
			//Still need to iterate through, but we don't use efn contants
			*/
		}
		if (index == 6) {
			cf >> currentWord;
			container->auend = (int) (100 * atof(currentWord));
		}
		if (index == 7) {
			cf >> currentWord;
			container->gubonus = (int) (100 * atof(currentWord));
		}
		if (index == 8) {
			cf >> currentWord;
			container->cslope = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 9) {
			cf >> currentWord;
			container->cint = (int) (100 * atof(currentWord));
		}
		if (index == 10) {
			cf >> currentWord;
			container->c3 = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 11) {
			cf >> currentWord;
			container->init = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 12) {
			cf >> currentWord;
			container->gail = (int) floor(.5 + atof(currentWord));
		}
	}

	cf.close();
	return 0;
}

int initDangleValues(string fileName, nndb_constants * container) {
	ifstream cf; //cf = current file
	char currentLine[256];
	string currentString;
	string s;
	int index;
	int i, j, k, l;
	int ii, jj, kk, ll; // ii = 1st base, jj = 2nd base, kk = 3rd base, ll = 1 up or 2 low

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 2; l++) {
					container->dangle[i][j][k][l] = INFINITY_;
				}
			}
		}
	}

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// The 8 first lines are junk
	for (index = 1; index <= 8; index++) {
		cf.getline(currentLine, 256);
	}

	// 8 lines of useful data
	i = 0;
	ii = 0;
	jj = 0;
	kk = 0;
	ll = 0;

	while (i < 8) {

		if (i != 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);

		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;

		jj = 0;

		int z = 0;
		int r = 0;

		while (s[z] != '\0') {

			if (s[z] == ' ')
				z++;

			else if (s[z] == '.') {
				z++;
				kk++;
				if (kk == 4)
					kk = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			} else {
				char value[10];
				int x = 0;

				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}

				value[x] = '\0';

				int temp = (int) floor(100.0 * atof(value) + .5);
				//if ( temp == 0 )	temp = -1;
				container->dangle[ii][jj][kk][ll] = temp;
				//cout<< "\n " << temp << "  "<< ii <<" "<< jj << " "<< kk << " " << ll ;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				kk++;
				if (kk == 4)
					kk = 0;

			}

		}

		i++;
		ii++;
		if (ii == 4)
			ii = 0;

		if (i == 4)
			ll = 1;

	}
	cf.close();
	return 0;
}

int initLoopValues(string fileName, nndb_constants * container) {
	// algorithm.c, line 2996
	ifstream cf; // current file
	char currentLine[256];
	char currentWord[256];
	string s;
	int index;
	int tempValue = 0;

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}
	// The 4 first lines are junk we don't need
	for (index = 1; index <= 4; index++) {
		cf.getline(currentLine, 256);
		//s = currentLine;
		//cout << s << endl;
	}
	while (index < 30) {
		for (int j = 1; j <= 4; j++) {
			cf >> currentWord;
			if (j == 1)
				index = atoi(currentWord);
			if (j > 1) {
				if (strcmp(currentWord, ".")) {
					tempValue = (int) (100 * atof(currentWord) + 0.5);
				} else {
					tempValue = INFINITY_;
				}
			}
			switch (j) {
			case 2:
				container->inter[index] = tempValue;
				break;
			case 3:
				container->bulge[index] = tempValue;
				break;
			case 4:
				container->hairpin[index] = tempValue;
				break;
			}
		}
	}
	cf.close();

#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}

int initTstkhValues(string fileName, nndb_constants * container) {
	ifstream cf; //cf = current file
	int i, j, k, l;
	int ii, jj, kk, ll;
	int index;
	char currentLine[256];
	string currentString;
	string s;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 4; l++) {
					container->tstkh[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// The 27 first lines are junk we don't need
	for (index = 1; index <= 15; index++) {
		cf.getline(currentLine, 256);
	}

	i = 0;
	kk = 0;
	ii = 0;
	jj = 0;
	ll = 0;

	while (i < 16) {

		if (i % 4 == 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);

		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;

		ll = 0;
		jj = 0;

		int z = 0;
		int r = 0;

		while (s[z] != '\0') {

			if (s[z] == ' ')
				z++;

			else if (s[z] == '.') {
				z++;
				ll++;
				if (ll == 4)
					ll = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			}

			else {
				char value[10];
				int x = 0;

				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}

				value[x] = '\0';

				int temp = (int) floor(100.0 * atof(value) + .5);
				container->tstkh[fourBaseIndex(ii,jj,kk,ll)] = temp;

				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				ll++;
				if (ll == 4)
					ll = 0;
			}
		}
		i++;
		if (!(i % 4))
			ii++;

		jj = 0;
		kk = (i % 4);

	}

	cf.close();
#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}

int initTstkiValues(string fileName, nndb_constants * container) {
	ifstream cf; //cf = current file
	int i, j, k, l;
	int ii, jj, kk, ll;
	int index;
	char currentLine[256];
	string currentString;
	string s;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 4; l++) {
					container->tstki[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// The 27 first lines are junk we don't need
	for (index = 1; index <= 15; index++) {
		cf.getline(currentLine, 256);
	}

	i = 0;
	kk = 0;
	ii = 0;
	jj = 0;
	ll = 0;

	while (i < 16) {

		if (i % 4 == 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);

		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;

		ll = 0;
		jj = 0;

		int z = 0;
		int r = 0;

		while (s[z] != '\0') {

			if (s[z] == ' ')
				z++;

			else if (s[z] == '.') {
				z++;
				ll++;
				if (ll == 4)
					ll = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			}

			else {
				char value[10];
				int x = 0;

				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}

				value[x] = '\0';

				int temp = (int) floor(100.0 * atof(value) + .5);
				container->tstki[fourBaseIndex(ii,jj,kk,ll)] = temp;

				//cout << "\n temp " << temp << " " <<ii << " "<< jj << " "<< kk << " " << ll;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				ll++;
				if (ll == 4)
					ll = 0;
			}
		}
		i++;
		if (!(i % 4))
			ii++;

		jj = 0;
		kk = (i % 4);

	}

	cf.close();
#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}

//SH: Rewritten as an error was being generated by this function.
int initTloopValues(string fileName, nndb_constants * container) {
	ifstream cf;
	int count;
	char currentLine[256];
	char currentSeqNumbers[7];
	char currentValue[6];

	currentSeqNumbers[6] = '\0';
	currentValue[5] = '\0';

	string s, temp;

#if 0
	cout << "Getting tloop values...";
#endif

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}
	// skip 2 lines
	cf.getline(currentLine, 256);
	cf.getline(currentLine, 256);

	container->numoftloops = 0;

	while (!cf.eof() && (++(container->numoftloops) < maxtloop)) {
		/*
		 MFold use a weird system : it's sthg like Sum(base[i]*5^i) / Algorithm.c line 3134
		 that's 97655 valuess
		 */
		int clindex=0;
		cf.getline(currentLine, 256);
		while(currentLine[clindex]== ' ') clindex++;
		for (count = 0; count < 6; count++) {
			temp = currentLine[count + clindex];
			currentSeqNumbers[count] = baseToDigit(temp);
			//cout << currentSeqNumbers[count];
			//cout << currentLine[count+1];
		}
		//cout << endl;
		clindex=clindex+7;
		while(currentLine[clindex]== ' ') clindex++;
		count = 0;
		while(currentLine[clindex+count]!=' '&&currentLine[clindex+count]!='\0') {
			currentValue[count] = currentLine[count + clindex];
			count++;
			//cout << currentValue[count];
		}
		//cout << endl;

		container->tloop[container->numoftloops][0] = (int) atoi(currentSeqNumbers);
		container->tloop[container->numoftloops][1] = (int) floor(100.0 * atof(currentValue) + 0.5);
		//cout << " --  "<< tloop[numoftloops][0] << " | " << tloop[numoftloops][1] << endl;
	}
	cf.close();
#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}

int initInt22Values(string fileName, nndb_constants * container) {

	//Read the 2x2 internal loops
	//key iloop22[a][b][c][d][j][l][k][m] =
	//a j l b
	//c k m d
	int i, j, k, r, q, t, y, z;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							for (y = 0; y < 4; y++)
								for (z = 0; z < 4; z++)
									container->iloop22[i][j][k][r][q][t][y][z] = INFINITY_;

	ifstream cf;
	int index, flag;
	char currentLine[256], currentValue[6];
	string s, s1, s2, temp;
	string sre;

	int base[4];
	int l, m;

#if 0
	cout << "Getting Int22 values...";
#endif
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// Get rid of the 38 1st lines
	for (index = 1; index < 28; index++) {
		cf.getline(currentLine, 256);
	}

	sre = "Y";
	flag = 0;

	for (index = 1; index <= 36; index++) { // 36 tables 16x16
		// Go to the beginning of the tables
		while (!flag) {
			cf.getline(currentLine, 256);

			s = currentLine;

			int z = 0;
			while (s[z] != '\0') {

				if (s[z] == 'Y')
					flag = 1;

				z++;
			}
		}
		flag = 0;

		// We skip 5 lines
		for (i = 0; i < 5; ++i) {
			cf.getline(currentLine, 256);
		}
 
		// get the closing bases
		cf.getline(currentLine, 256);
		s1 = currentLine;
		cf.getline(currentLine, 256);
		s2 = currentLine;
		s = s1 + s2;

		int z = 0;
		int k = 0;
		while (s[z] != '\0') {

			if (s[z] == 'A')
				base[k++] = BASE_A;
			else if (s[z] == 'C')
				base[k++] = BASE_C;
			else if (s[z] == 'G')
				base[k++] = BASE_G;
			else if (s[z] == 'U')
				base[k++] = BASE_U;

			z++;
		}

		cf.getline(currentLine, 256);

		//key iloop22[a][b][c][d][j][l][k][m] =
		//a j l b
		//c k m d


		for (int rowIndex = 1; rowIndex <= 16; rowIndex++) {
			for (int colIndex = 1; colIndex <= 16; colIndex++) {
				cf >> currentValue;
				// rowIndex = j*4+k
				// colIndex = l*4+m

				j = ((rowIndex - 1) - (rowIndex - 1) % 4) / 4;
				k = (rowIndex - 1) % 4;

				l = ((colIndex - 1) - (colIndex - 1) % 4) / 4;
				m = (colIndex - 1) % 4;

				container->iloop22[base[0]][base[1]][base[2]][base[3]][j][l][k][m]
				                                                     = (int) floor(100.0 * atof(currentValue) + 0.5);
				//int temp = (int) floor(100.0*atof(currentValue)+0.5);
				//printf("\n temp is : %d    %d %d %d %d   %d  %d   %d %d", temp, base[0], base[2], base[1], base[3], j, k, l, m );
			}
		}
	}

	cf.close();
#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}

int initInt21Values(string fileName, nndb_constants * container) {

	// 24x6 arrays of 4x4 values
	//      c
	//   a     f
	//   b     g
	//     d e

	ifstream cf;
	char currentLine[256];
	string sre;
	string s, s1, s2;
	int a, b, c, d, e, f, g, index;
	int i, j, k, r, q, t, y;
	int z;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							for (y = 0; y < 4; y++)
								container->iloop21[i][j][k][r][q][t][y] = INFINITY_;
	a = 0;
	b = 0;
	c = 0;
	d = 0;
	e = 0;
	f = 0;
	g = 0;
	k = 0;

	int base1[7];
	int base2[7];

	base1[1] = BASE_A + 1;
	base2[1] = BASE_U + 1;
	base1[2] = BASE_C + 1;
	base2[2] = BASE_G + 1;
	base1[3] = BASE_G + 1;
	base2[3] = BASE_C + 1;
	base1[4] = BASE_U + 1;
	base2[4] = BASE_A + 1;
	base1[5] = BASE_G + 1;
	base2[5] = BASE_U + 1;
	base1[6] = BASE_U + 1;
	base2[6] = BASE_G + 1;

	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// Get rid of the 18 1st lines
	for (index = 1; index <= 17; index++) {
		cf.getline(currentLine, 256);
	}

	i = 1;

	while (i <= 6) {

		j = 1;

		while (j <= 4) {

			k = 1;

			for (index = 1; index <= 10; index++)
				cf.getline(currentLine, 256);

			s = currentLine;

			while (k <= 4) {

				cf.getline(currentLine, 256);

				s = currentLine;

				r = 0;
				z = 0;
				int jj = 1;
				d = 1;
				while (s[z] != '\0') {

					if (s[z] == ' ')
						z++;

					else if (s[z] == '.') {
						z++;
						d++;
						if (d == 5)
							d = 1;
						r++;
						if (r % 4 == 0)
							jj++;
					}

					else {
						char value[10];
						int x = 0;

						while (s[z] != ' ' && s[z] != '\0') {
							value[x++] = s[z++];
						}

						value[x] = '\0';

						int temp = (int) floor(100.0 * atof(value) + .5);
						a = base1[i];
						b = base2[i];
						f = base1[jj];
						g = base2[jj];
						c = k;
						e = j;

						container->iloop21[a - 1][b - 1][c - 1][d - 1][e - 1][f - 1][g - 1]
						                                                  = temp;
						//printf("\n temp %d, a %d , b %d, c %d, d %d, e %d, f %d, g %d", temp, a, b, c, d, e, f, g);
						r++;
						z++;
						if (r % 4 == 0)
							jj++;
						d++;
						if (d == 5)
							d = 1;
					}
				}
				k++;
			}
			j++;
		}
		i++;
	}
	cf.close();
	return 0;
}

int initInt11Values(string fileName, nndb_constants * container) {

	//Read the 1x1 internal loops
	//key iloop11[a][b][c][d][j][l][k][m] =
	//a b c
	//d e f

	int i, j, k, r, q, t;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							container->iloop11[i][j][k][r][q][t] = INFINITY_;

	ifstream cf;
	int index;
	char currentLine[256];
	string s;
	int base1[7];
	int base2[7];

	//int j, k, l, m;
	int a, b, c, d, f;

#if 0
	cout << "Getting Int11 values...";
#endif
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << endl;
		exit(-1);
	}

	// Get rid of the 19 1st lines
	for (index = 1; index <= 17; index++) {
		cf.getline(currentLine, 256);
	}

	/*
	 Structure of the file:
	 array of 6x6 arrays
	 Order:
	 AU CG GC UA GU UG
	 */

	base1[1] = BASE_A + 1;
	base2[1] = BASE_U + 1;
	base1[2] = BASE_C + 1;
	base2[2] = BASE_G + 1;
	base1[3] = BASE_G + 1;
	base2[3] = BASE_C + 1;
	base1[4] = BASE_U + 1;
	base2[4] = BASE_A + 1;
	base1[5] = BASE_G + 1;
	base2[5] = BASE_U + 1;
	base1[6] = BASE_U + 1;
	base2[6] = BASE_G + 1;

	i = 0;
	k = 0;

	while (k < 6) {
		k++;
		index = 0;
		for (index = 1; index <= 10; index++) {
			cf.getline(currentLine, 256);
		}

		i = 0;
		b = 1;

		while (i < 4) {

			int jj = 1;
			++i;
			cf.getline(currentLine, 256);
			s = currentLine;

			j = 0;

			int r = 0;
			int z = 0;
			int e = 1;
			while (s[z] != '\0') {

				if (s[z] == ' ')
					z++;

				else if (s[z] == '.') {
					z++;
					e++;
					if (e == 5)
						e = 1;
					r++;
					if (r % 4 == 0)
						jj++;
				}

				else {
					char value[10];
					int x = 0;

					while (s[z] != ' ' && s[z] != '\0') {
						value[x++] = s[z++];
					}

					value[x] = '\0';

					int temp = (int) floor(100.0 * atof(value) + .5);
					a = base1[k];
					d = base2[k];
					c = base1[jj];
					f = base2[jj];
					container->iloop11[a - 1][b - 1][c - 1][d - 1][e - 1][f - 1] = temp;
					r++;
					z++;
					if (r % 4 == 0)
						jj++;
					e++;
					if (e == 5)
						e = 1;
				}
			}
			b++;
		}
	}

	cf.close();
#if 0
	cout << " Done!" << endl;
#endif
	return 0;
}
