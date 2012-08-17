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

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <math.h>

#include <cstring>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "data.h"
#include "utils.h"
#include "constants.h"
#include "global.h"
#include "loader.h"

#define xstr(s) str(s)
#define str(s) #s

using namespace std;

std::string EN_DATADIR;

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
int iloop21[5][5][5][5][5][5][5]; /* 2*1 internal loops */
int iloop22[5][5][5][5][5][5][5][5]; /* 2*2 internal looops */
int iloop11[5][5][5][5][5][5]; /* 1*1 internal loops */

//int coax[6][6][6][6];
//int tstackcoax[6][6][6][6];
//int coaxstack[6][6][6][6];
//int tstack[6][6][6][6];
//int tstkm[6][6][6][6];

int tstackm[5][5][6][6];
int tstacke[5][5][6][6];
int tstacki23[5][5][5][5];

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

void readThermodynamicParameters(const char *userdatadir,bool userdatalogic, int unamode = 0, int rnamode = 0, int mismatch = 0) {
	struct stat buf;

	if (!userdatalogic) {
		std::string opt1, opt2, opt3;
		char cwd[1024];
		bool found = false;

		opt1.assign(xstr(DATADIR));
		if (getenv("GTFOLDDATADIR") != 0)
			opt2.assign(getenv("GTFOLDDATADIR"));
		if (getcwd(cwd, sizeof(cwd)) != 0) {
			opt3.assign(cwd);
			opt3 += "/data";
		}

		if (stat(opt2.c_str(), &buf) == 0) {
			found = true;
			EN_DATADIR.assign(opt2);
			fprintf(stdout,"Checking for environ variable 'GTFOLDDATADIR', found \n");
		} else {
			if (found == false)	fprintf(stdout,"Checking for environ variable 'GTFOLDDATADIR', not found \n");
		}

		if (!found && stat(opt1.c_str(), &buf) == 0) {
			found = true;
			EN_DATADIR.assign(opt1); 
			fprintf(stdout,"Checking for default datadir '%s', found.\n", opt1.c_str());
		} else  {
			if (found == false)	fprintf(stdout,"Checking for default datadir '%s', not found.\n", opt1.c_str());
		}

		if (!found && stat(opt3.c_str(), &buf) == 0){
			found = true;
			EN_DATADIR.assign(opt3);
			fprintf(stdout,"Checking under current dir '%s', found.\n", opt3.c_str());

		} else {
			if (found == false)	fprintf(stdout,"Checking under current dir '%s', not found.\n", opt3.c_str());
		}

		EN_DATADIR += "/";
		if (found == false) {
			exit(-1);
		}	

		if (unamode) {
			EN_DATADIR += "UNAParams";
		}else if (rnamode) {
			EN_DATADIR += "RNAParams";
		} else {
			EN_DATADIR += "Turner99";
		}
	}
	else {
		EN_DATADIR.assign(userdatadir);
		if (stat(EN_DATADIR.c_str(), &buf) == -1){
			fprintf(stdout,"Checking for parameter files in dir '%s', not found.\n\n", EN_DATADIR.c_str());
			exit(-1);		
		}
	}

	//Handle the ending forward slash case
	if (EN_DATADIR[EN_DATADIR.length() - 1] != '/') {
		EN_DATADIR += "/";
	}

	initMiscloopValues("miscloop.DAT", EN_DATADIR);
	initDangleValues("dangle.DAT", EN_DATADIR);
	initStackValues("stack.DAT", EN_DATADIR);
	initLoopValues("loop.DAT", EN_DATADIR);
	initTstkhValues("tstackh.DAT", EN_DATADIR);
	initTstkiValues("tstacki.DAT", EN_DATADIR);
	initTloopValues("tloop.DAT", EN_DATADIR);

	if (unamode) {
		initInt21Values("asint1x2.DAT", EN_DATADIR);
		initInt22Values("sint4.DAT", EN_DATADIR);
		initInt11Values("sint2.DAT", EN_DATADIR);
		initTstkmValues("tstackm.DAT", EN_DATADIR);
		initTstkeValues("tstacke.DAT", EN_DATADIR);
		initTstk23Values("tstacki23.DAT", EN_DATADIR);
	} else if (mismatch) {
		initTstkmValues("tstackm.DAT", EN_DATADIR);
		initTstkeValues("tstacke.DAT", EN_DATADIR);
		initInt21Values("int21.DAT", EN_DATADIR);
		initInt22Values("int22.DAT", EN_DATADIR);
		initInt11Values("int11.DAT", EN_DATADIR);

	}	
	else {
		initInt21Values("int21.DAT", EN_DATADIR);
		initInt22Values("int22.DAT", EN_DATADIR);
		initInt11Values("int11.DAT", EN_DATADIR);
	}
}

int initStackValues(const string& fileName, const string& dirPath)  { 
	std::string filePath;
	std::ifstream cf;

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	if (!cf.good()) {
		std::cerr << "File open failed - " << filePath << std::endl;
		cf.close();
		exit(-1);		
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 4; l++) {
					stack[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}
	
	int count = 0;

	while (cf.good()) {
		int i, j, k, l;
		
		std::string line;
		std::stringstream ss;

		getline(cf, line);
		ss << line;
		if (line.empty()) continue;

		i = count/4;
		k = count%4;
		
		for (int cols = 0; cols < 16; cols++) {
			std::string val;
			ss >> val;
			j = cols/4;
			l = cols%4; 		
			if (!(val == "inf")) {
				stack[fourBaseIndex(i,j,k,l)]  = (int) floor(100.0 * atof(val.c_str()) + .5);
			}
		}	
		++count;
	}
	
	return 0;
}

int initMiscloopValues(const string& fileName, const string& dirpath) {
	std::string filePath;
	std::ifstream cf;

	filePath = dirpath + fileName;
	cf.open(filePath.c_str(), ios::in);

	if (!cf.good()) {
		std::cerr << "File open failed - " << filePath << std::endl;
		cf.close();
		exit(-1);		
	}

	char currentWord[256];
	std::string s;

	for (int index = 1; index < 13; index++) { // There are total 12 values to read in.
		if (index == 1) {
			cf >> currentWord;
			prelog = 100 * atof(currentWord);
		}
		if (index == 2) {
			cf >> currentWord;
			maxpen = int(atof(currentWord) * 100.0 + .5);
		}
		if (index == 3) {
			for (int count = 1; count <= 4; count++) {
				cf >> currentWord;
				s = currentWord;
				poppen[count] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
		}
		if (index == 4) {
			eparam[1] = 0;
			eparam[2] = 0;
			eparam[3] = 0;
			eparam[4] = 0;
			eparam[7] = 30;
			eparam[8] = 30;
			eparam[9] = -500;
			int table[4];
			table[1] = 5;
			table[2] = 6;
			table[3] = 10;
			for (int count = 1; count <= 3; count++) {
				cf >> currentWord;
				s = currentWord;
				multConst[count - 1] = (int) (atof(s.c_str()) * 100 + 0.5);
				eparam[table[count]] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
		}
		if (index == 5) {
			int table[4];
			for (int count = 1; count <= 3; count++) {
				cf >> currentWord;
				s = currentWord;
				table[count] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
			efn2a = table[1];
			efn2b = table[2] - 1;
			efn2c = table[3] - 1;
		}
		if (index == 6) {
			cf >> currentWord;
			auend = (int) (100 * atof(currentWord));
		}
		if (index == 7) {
			cf >> currentWord;
			gubonus = (int) (100 * atof(currentWord));
		}
		if (index == 8) {
			cf >> currentWord;
			cslope = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 9) {
			cf >> currentWord;
			cint = (int) (100 * atof(currentWord));
		}
		if (index == 10) {
			cf >> currentWord;
			c3 = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 11) {
			cf >> currentWord;
			init = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 12) {
			cf >> currentWord;
			gail = (int) floor(.5 + atof(currentWord));
		}
	}

	cf.close();	 

	return 0;
}

int initDangleValues(const std::string& fileName, 
					  const std::string& dirPath) { 
	std::string filePath;
	std::ifstream cf;

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 2; l++) {
					dangle[i][j][k][l] = INFINITY_;
				}
			}
		}
	}
	
	int count = 0;

	while (cf.good()) {
		int i, j, k, l;
		
		std::string line;
		std::stringstream ss;

		getline(cf, line);
		ss << line;
		if (line.empty()) continue;

		i = count%4;
		l = count/4;
		
		for (int cols = 0; cols < 16; cols++) {
			std::string val;

			j = cols/4;
			ss >> val;
			k = cols%4; 		
			if (!(val == "inf")) {
				dangle[i][j][k][l] = (int) floor(100.0 * atof(val.c_str()) + .5);
			}
		}	
		++count;
	}
	
	return 0;
}

int initLoopValues( const string& fileName, const string& dirPath) {
	// algorithm.c, line 2996
	ifstream cf; // current file
	std::string filePath;

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);
	
	//char currentLine[256];
	char currentWord[256];
	string s;
	int index= 0;
	int tempValue = 0;

	if (cf.fail()) {
		cerr << "File open failed " << filePath << endl;
		exit(-1);
	}

	while (index < 30) {
		for (int j = 1; j <= 4; j++) {
			cf >> currentWord;
			if (j == 1)
				index = atoi(currentWord);
			if (j > 1) {
				if (strcmp(currentWord, "inf")) {
					tempValue = (int) (100 * atof(currentWord) + 0.5);
				} else {
					tempValue = INFINITY_;
				}
			}
			switch (j) {
			case 2:
				inter[index] = tempValue;
				break;
			case 3:
				bulge[index] = tempValue;
				break;
			case 4:
				hairpin[index] = tempValue;
				break;
			}
		}
	}
	cf.close();

	return 0;
}

int initTstk23Values(const std::string& fileName, const std::string& dirPath) {
	int i1, j1, i2, j2;                                                              
	std::string filePath;
	std::ifstream cf;
	std::string val;						

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	if (cf.fail()) {
		cerr << "File open failed " << filePath << endl;
		exit(-1);
	}

	for (i1 = 0; i1 < 5; ++i1)
		for (i2 = 0; i2 < 5; ++i2)
			for (j1 = 0; j1 < 5; ++j1)
				for (j2 = 0; j2 < 5; ++j2)
					if (i1 == 4 || j1 == 4)
						tstacki23[i1][j1][i2][j2] = INFINITY_;
					else if (i2 == 4 || j2 == 4)
						tstacki23[i1][j1][i2][j2] = 0;
					else {
						cf >> val;
						tstacki23[i1][j1][i2][j2] = (val == "inf")? (INFINITY_): 
							((int) floor(100.0 * atof(val.c_str()) + .5));
					}
	cf.close();

	return 0;
}




int initTstkeValues(const std::string& fileName, const std::string& dirPath) {
	int i1, j1, i2, j2;                                                              
	std::string filePath;
	std::ifstream cf;
	std::string val;						

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	if (cf.fail()) {
		cerr << "File open failed " << filePath << endl;
		exit(-1);
	}

	for (i1 = 0; i1 < 5; ++i1) {                                    
		for (i2 = 0; i2 < 6; ++i2) {                                                     
			for (j1 = 0; j1 < 5; ++j1) {                                                   
				for (j2 = 0; j2 < 6; ++j2) {                                                      
					if (i1 == 4 || j1 == 4)                                                       
						tstacke[i1][j1][i2][j2] = INFINITY_;                                         
					else if (i2 == 5 || j2 == 5)                                                  
						tstacke[i1][j1][i2][j2] = INFINITY_;                                         
					else if (i2 == 4 || j2 == 4)                                                  
						tstacke[i1][j1][i2][j2] = 0;                                                
					else { 
						cf >> val;
						tstacke[i1][j1][i2][j2] = (val == "inf")? (INFINITY_): 
							((int) floor(100.0 * atof(val.c_str()) + .5));
					}
				}
			}
		}			
	}

	cf.close();

	return 0;
}

int initTstkmValues(const std::string& fileName, const std::string& dirPath) {
	int i1, j1, i2, j2;                                                              
	std::string filePath;
	std::ifstream cf;
	std::string val;						

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);
	
	if (cf.fail()) {
		cerr << "File open failed " << filePath << endl;
		exit(-1);
	}

	for (i1 = 0; i1 < 5; ++i1) {                                    
		for (i2 = 0; i2 < 6; ++i2) {                                                     
			for (j1 = 0; j1 < 5; ++j1) {                                                   
				for (j2 = 0; j2 < 6; ++j2) {                                                      
					if (i1 == 4 || j1 == 4)                                                       
						tstackm[i1][j1][i2][j2] = INFINITY_;                                         
					else if (i2 == 5 || j2 == 5)                                                  
						tstackm[i1][j1][i2][j2] = INFINITY_;                                         
					else if (i2 == 4 || j2 == 4)                                                  
						tstackm[i1][j1][i2][j2] = 0;                                                
					else { 
						cf >> val;
						tstackm[i1][j1][i2][j2] = (val == "inf")? (INFINITY_): 
							((int) floor(100.0 * atof(val.c_str()) + .5));
					}
				}
			}
		}			
	}

	cf.close();

	return 0;
}

int initTstkhValues(const std::string& fileName, const std::string& dirPath) {
	std::string filePath;
	std::ifstream cf;

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 2; l++) {
					tstkh[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}
	
	int count = 0;

	while (cf.good()) {
		int i, j, k, l;
		
		std::string line;
		std::stringstream ss;

		getline(cf, line);
		ss << line;
		if (line.empty()) continue;

		i = count/4;
		k = count%4;
		
		for (int cols = 0; cols < 16; cols++) {
			std::string val;
			ss >> val;
			j = cols/4;
			l = cols%4; 		
			if (!(val == "inf")) {
				tstkh[fourBaseIndex(i,j,k,l)]= (int) floor(100.0 * atof(val.c_str()) + .5);
			}
		}	
		++count;
	}
	
	return 0;
}


int initTstkiValues(const std::string& fileName, const std::string& dirPath) {
	std::string filePath;
	std::ifstream cf;

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 2; l++) {
					tstki[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}

	int count = 0;

	while (cf.good()) {
		int i, j, k, l;

		std::string line;
		std::stringstream ss;

		getline(cf, line);
		ss << line;
		if (line.empty()) continue;

		i = count/4;
		k = count%4;

		for (int cols = 0; cols < 16; cols++) {
			std::string val;
			ss >> val;
			j = cols/4;
			l = cols%4; 		
			if (!(val == "inf")) {
				tstki[fourBaseIndex(i,j,k,l)]= (int) floor(100.0 * atof(val.c_str()) + .5);
			}
		}	
		++count;
	}

	return 0;
}

int initTloopValues(const std::string& fileName, const std::string& dirPath) {
	std::string filePath;
	std::ifstream cf;

	filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);

	int count;
	char currentLine[256];
	char currentSeqNumbers[7];
	char currentValue[6];

	currentSeqNumbers[6] = '\0';
	currentValue[5] = '\0';

	string s, temp;

	numoftloops = 0;

	while (!cf.eof() && (++(numoftloops) < maxtloop)) {
		/*
		 MFold use a weird system : it's sthg like Sum(base[i]*5^i) / Algorithm.c line 3134
		 that's 97655 valuess
		 */
		int clindex=0;
		cf.getline(currentLine, 256);
		while(currentLine[clindex]== ' ') clindex++;
		for (count = 0; count < 6; count++) {
			temp = currentLine[count + clindex];
			currentSeqNumbers[count] = baseToDigit(temp.c_str());
		}
		clindex=clindex+7;
		while(currentLine[clindex]== ' ') clindex++;
		count = 0;
		while(currentLine[clindex+count]!=' '&&currentLine[clindex+count]!='\0') {
			currentValue[count] = currentLine[count + clindex];
			count++;
		}
		
		tloop[numoftloops][0] = (int) atoi(currentSeqNumbers);
		
		if (!(strcmp(currentValue,"inf")==0)) 
			tloop[numoftloops][1] = (int) floor(100.0 * atof(currentValue) + 0.5);
	}
	cf.close();
	return 0;
}

int initInt22Values(const std::string& fileName, const std::string& dirPath) {
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
									iloop22[i][j][k][r][q][t][y][z] = INFINITY_;

	ifstream cf;
	char currentLine[256], currentValue[6];

	int base[4];
	int l, m;

	string filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);
	if (!cf.good()) {
		cerr << "File open failed " << filePath << endl;
		exit(-1);
	}

	const char* str[] = {"A-U","C-G","G-C","U-A","G-U","U-G"};
	

	int outer, inner;

	for (int index = 1; index <= 36; index++) { // 36 tables 16x16
		std::stringstream ss1;
		std::string s ;

		outer = (index-1)/6;
		inner = (index-1)%6;
		ss1 << str[outer][0] << '-' << str[inner][0] << '-' << str[outer][2] << '-' <<
			str[inner][2];
	    ss1 >> s;	
		
		int k = 0;
		int z = 0;
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

		for (int rowIndex = 1; rowIndex <= 16; rowIndex++) {
			cf.getline(currentLine, 256);
			std::stringstream ss;
			ss << currentLine;
			for (int colIndex = 1; colIndex <= 16; colIndex++) {
				ss >> currentValue;	

				j = ((rowIndex - 1) - (rowIndex - 1) % 4) / 4;
				k = (rowIndex - 1) % 4;

				l = ((colIndex - 1) - (colIndex - 1) % 4) / 4;
				m = (colIndex - 1) % 4;

				if (!(strcmp(currentValue,"inf")==0)) 
					iloop22[base[0]][base[1]][base[2]][base[3]][j][l][k][m] 
						= (int) floor(100.0 * atof(currentValue) + 0.5);
			}
		}
	}

	cf.close();
	return 0;
}

int initInt21Values(const std::string& fileName, const std::string& dirPath) {
	// 24x6 arrays of 4x4 values
	//      c
	//   a     f
	//   b     g
	//     d e

	ifstream cf;
	char currentLine[256];
	string sre;
	string s, s1, s2;
	int a, b, c, d, e, f, g;
	int i, j, k, r, q, t, y;
	int z;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							for (y = 0; y < 4; y++)
								iloop21[i][j][k][r][q][t][y] = INFINITY_;
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

	std::string filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed" << filePath << endl;
		exit(-1);
	}

	i = 1;
	while (i <= 6) {
		j = 1;
		while (j <= 4) {
			k = 1;
			while (k <= 4) {
				cf.getline(currentLine, 256);
				r = 0;
				int jj = 1;
				d = 1;

				std::stringstream ss;
				ss << currentLine;

				for (z = 1; z <=24 ; z++) {
					char value[32];
					ss >> value;
					int temp = (int) floor(100.0 * atof(value) + .5);
					a = base1[i];
					b = base2[i];
					f = base1[jj];
					g = base2[jj];
					c = k;
					e = j;

					if (!(strcmp(value,"inf")==0)) 
						iloop21[a - 1][b - 1][c - 1][d - 1][e - 1][f - 1][g - 1] = temp;
					r++;
					if (r % 4 == 0) jj++;
					d++;
					if (d == 5) d = 1;
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


int initInt11Values(const std::string& fileName, const std::string& dirPath) {
	int i, j, k, r, q, t;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							iloop11[i][j][k][r][q][t] = INFINITY_;

	ifstream cf;
	char currentLine[256];
	string s;
	int base1[7];
	int base2[7];

	int a, b, c, d, f;

	std::string filePath = dirPath + fileName;
	cf.open(filePath.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "File open failed " << filePath << endl;
		exit(-1);
	}

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
		i = 0;
		b = 1;

		while (i < 4) {
			int jj = 1;
			++i;
			cf.getline(currentLine, 256);
			j = 0;

			int r = 0;
			int e = 1;
			std::stringstream ss;
			ss << currentLine;

			for (int z=1; z <= 24; ++z) {
				char value[32];
				ss >> value;	
				int temp = (int) floor(100.0 * atof(value) + .5);
				a = base1[k];
				d = base2[k];
				c = base1[jj];
				f = base2[jj];
				
				if (!(strcmp(value,"inf")==0)) 
					iloop11[a - 1][b - 1][c - 1][d - 1][e - 1][f - 1] = temp;
				
				r++;
				if (r % 4 == 0) jj++;
				e++;
				if (e == 5) e = 1;
			}
			b++;
		}
	}

	cf.close();
	
	return 0;
}
