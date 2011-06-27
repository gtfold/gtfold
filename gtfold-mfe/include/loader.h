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

#ifndef _LOADER_H
#define _LOADER_H

#include <string>

#include "constants.h"
#include "data.h"

void readThermodynamicParameters(const char *userdatadir,bool userdatalogic, int unamode, int rnamode, int t_mismatch);

int initStackValues(const std::string& fileName, const std::string& dirPath);
int initMiscloopValues(const std::string& fileName, const std::string& dirPath);
int initDangleValues(const std::string& fileName, const std::string& dirPath);
int initLoopValues(const std::string& fileName, const std::string& dirPath);
int initTstkhValues(const std::string& fileName, const std::string& dirPath);
int initTstkiValues(const std::string& fileName, const std::string& dirPath);
int initTloopValues(const std::string& fileName, const std::string& dirPath);
int initInt21Values(const std::string& fileName, const std::string& dirPath);
int initInt22Values(const std::string& fileName, const std::string& dirPath);
int initInt11Values(const std::string& fileName, const std::string& dirPath);
int	initTstkmValues(const std::string& fileName, const std::string& dirPath);
int	initTstkeValues(const std::string& fileName, const std::string& dirPath);
int	initTstk23Values(const std::string& fileName, const std::string& dirPath);

extern std::string EN_DATADIR;

#endif
