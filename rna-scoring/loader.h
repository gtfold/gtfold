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

/* This file externs the prototypes of functions defined in loader.cc*/

#ifndef _LOADER_H
#define _LOADER_H

#include <iostream>
#include <fstream>
#include <istream>
#include <math.h>
#include <string>
#include "constants.h"
#include "data.h"

nndb_constants* populate(const char *userdatadir,bool userdatalogic);
unsigned char getBase(std::string base);
unsigned char getBase1(std::string base);
int initStackValues(std::string fileName, nndb_constants * container);
int initMiscloopValues(std::string fileName, nndb_constants * container);
int initDangleValues(std::string fileName, nndb_constants * container);
int initLoopValues(std::string fileName, nndb_constants * container);
int initTstkhValues(std::string fileName, nndb_constants * container);
int initTstkiValues(std::string fileName, nndb_constants * container);
int initTloopValues(std::string fileName, nndb_constants * container);
int initInt21Values(std::string fileName, nndb_constants * container);
int initInt22Values(std::string fileName, nndb_constants * container);
int initInt11Values(std::string fileName, nndb_constants * container);

#endif
