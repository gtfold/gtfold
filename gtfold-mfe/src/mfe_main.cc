/*
 GTfold: compute minimum free energy of RNA secondary structure
 Copyright (C) 2007-2011  David A. Bader, Christine E. Heitsch, 
 and Steve C. Harvey
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
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>

#include "main.h"
#include "mfe_main.h"
#include "utils.h"
#include "loader.h"
//#include "options.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "constraints.h"
#include "traceback.h"
#include "shapereader.h"

using namespace std;

static bool PARAM_DIR = false;
//static bool LIMIT_DISTANCE;
static bool CONS_ENABLED = false;
static bool VERBOSE = false;
static bool SHAPE_ENABLED = false;
static bool T_MISMATCH = false;
static bool UNAMODE = false;
static bool RNAMODE = false;
static bool b_prefilter = false;

static string seqfile = "";
static string constraintsFile = "";
static string outputPrefix = "";
static string outputFile = "";
static string outputDir = "";
static string shapeFile = "";
static string paramDir; // default value

static int dangles=-1;
static int prefilter1=2;
static int prefilter2=2;

static int nThreads = -1;
static int contactDistance = -1;

static void help();
void printRunConfiguration(string seq);
void parse_mfe_options(int argc, char** argv);

void init_fold(const char* seq) {
  assert(seq != NULL);
  int len = strlen(seq);

  init_global_params(len);

  if (!encodeSequence(seq)) {
    free_fold(len);
    exit(0);
  }

  create_tables(len);

  if (CONS_ENABLED) {
    init_constraints(constraintsFile.c_str(), len);
  }

  if (SHAPE_ENABLED) {
    readSHAPEarray(shapeFile.c_str(),len);
  }

  if (UNAMODE) {
    if (T_MISMATCH) printf("Ignoring -m option, using --unafold\n");
    if (PARAM_DIR) printf("Ignoring -p option, using --unafold\n");
    if (dangles == 0 || dangles == 2) 
      printf("Ignoring -d option, using --unafold\n");
    if (b_prefilter == 1) 
      printf("Ignoring --prefilter option, using --unafold\n");
    T_MISMATCH = false;
    PARAM_DIR = false;
    dangles = -1;
    b_prefilter = false;
  }

  if (RNAMODE) {
    if (T_MISMATCH) printf("Ignoring -m option, using --rnafold\n");
    if (PARAM_DIR) printf("Ignoring -p option, using --rnafold\n");
    if (dangles == 0 || dangles == 2) 
      printf("Ignoring -d option, using --rnafold\n");
    if (b_prefilter == 1) 
      printf("Ignoring --prefilter option, using --rnafold\n");
    T_MISMATCH = false;
    PARAM_DIR = false;
    dangles = -1;
    b_prefilter = false;
  }

  if ((dangles == 0 || dangles == 2) && !UNAMODE && !RNAMODE) {
    if (T_MISMATCH) printf("Ignoring -m option, using -d option\n");
    T_MISMATCH = false;
  } else {
    if (dangles != -1 && !UNAMODE && !RNAMODE) printf("Ignoring -d as it accept 0 or 2 only\n");	
    dangles = -1;
  }

  g_nthreads = nThreads;
  g_unamode  = UNAMODE;
  g_mismatch = T_MISMATCH;
  g_verbose  = VERBOSE;
  g_prefilter_mode  = b_prefilter;
  g_prefilter1  = prefilter1;
  g_prefilter2  = prefilter2;
  g_dangles = dangles;

#ifdef DEBUG
  printf("g_nthreads = %d\n", g_nthreads);
  printf("g_unamode = %d\n", g_unamode);
  printf("g_mismatch = %d\n", g_mismatch);
  printf("g_prefilter_mode = %d\n", g_prefilter_mode);
  printf("g_dangles = %d\n", g_dangles);

#endif

}

void free_fold(int len) {
	if (CONS_ENABLED) 
		free_constraints(len);
	if (SHAPE_ENABLED){
		free_shapeArray(len);
	}

	free_tables(len);
	free_global_params();
}


int mfe_main(int argc, char** argv) {
	std::string seq;
	int energy;
	
	parse_mfe_options(argc, argv);

	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}

	init_fold(seq.c_str());
	
	// Read in thermodynamic parameters. Always use Turner99 data (for now)
  readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, UNAMODE, RNAMODE, T_MISMATCH);
	printRunConfiguration(seq);

	printf("\nComputing minimum free energy structure...\n");
	fflush(stdout);

	double t1 = get_seconds();
	energy = calculate(seq.length()) ; 
	t1 = get_seconds() - t1;
	
	printf("Done.\n\n");
	printf("Results:\n");
	if (energy >= MAXENG)	
		printf("- Minimum Free Energy: %12.4f kcal/mol\n", 0.00);
	else
		printf("- Minimum Free Energy: %12.4f kcal/mol\n", energy/100.00);
	printf("- MFE runtime: %9.6f seconds\n", t1);

	t1 = get_seconds();
	trace(seq.length());
	t1 = get_seconds() - t1;
	
	printf("\n");
	print_sequence(seq.length());
	print_structure(seq.length());

	if (CONS_ENABLED)
		print_constraints(seq.length());

	if (SHAPE_ENABLED && VERBOSE)
		print_shapeArray(seq.length());

	save_ct_file(outputFile, seq, energy);
	printf("\nMFE structure saved in .ct format to %s\n", outputFile.c_str());

	if(CONS_ENABLED && VERBOSE){
		printf("Verifying that structure fulfills constraint criteria... ");
		if(verify_structure()){
			printf("OK\n");
		}
		else{
			printf("ERROR: NOT OK!!\n");
			fprintf(stderr, "ERROR: Structure does not fulfill constraint criteria.\n");
			fprintf(stderr, "Structure file: %s\n", outputFile.c_str());
			fprintf(stderr, "Constraint file: %s\n", constraintsFile.c_str());
		}
	}

	free_fold(seq.length());
	
  return EXIT_SUCCESS;
}


void parse_mfe_options(int argc, char** argv) {
  int i;

  for(i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
        help();
      } else if(strcmp(argv[i], "--constraints") == 0 || strcmp(argv[i], "-c") == 0) {
        if(i < argc) {
          constraintsFile = argv[++i];
          CONS_ENABLED = true;
        }
        else
          help();
      } else if(strcmp(argv[i], "--limitCD") == 0 || strcmp(argv[i], "-l") == 0) {
        if(i < argc){
          contactDistance = atoi(argv[++i]);
          stringstream ss;
          ss << contactDistance;
          if (contactDistance >= 0 && !strcmp(ss.str().c_str(),argv[i])) {
            //LIMIT_DISTANCE = true;
            enable_limit_distance(true);  
            set_contact_distance(contactDistance);
          }
          else
            help();
        }
        else
          help();
      }else if(strcmp(argv[i], "--prefix") == 0 || strcmp(argv[i], "-o") == 0) {
        if(i < argc)
          outputPrefix = argv[++i];
        else
          help();
      } else if (strcmp(argv[i], "--workdir") == 0 || strcmp(argv[i], "-w") == 0) {
        if(i < argc)
          outputDir = argv[++i];
        else
          help();
      } else if (strcmp(argv[i], "--paramdir") == 0 || strcmp(argv[i], "-p") == 0) {
        if(i < argc) {
          paramDir = argv[++i];
          PARAM_DIR = true;
        }
        else
          help();
      } else if (strcmp(argv[i], "--dangle") == 0 || strcmp(argv[i], "-d") == 0) {
        std::string cmd = argv[i];
        if(i < argc) {
          dangles = atoi(argv[++i]);
          if (!(dangles == 0 || dangles == 2)) {
            dangles = -1;
            printf("Ignoring %s option as it accepts either 0 or 2\n", cmd.c_str());
          } 
        } else
          help();
      } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mismatch") == 0) {
        T_MISMATCH = true;
      } else if (strcmp(argv[i], "--unafold") == 0) {
        UNAMODE = true;
      } else if (strcmp(argv[i], "--rnafold") == 0) {
        RNAMODE = true;
      } else if (strcmp(argv[i], "--prefilter") == 0) {
        if(i < argc) {
          prefilter1 = atoi(argv[++i]);
          if (prefilter1 <= 0 ) {
            printf("INVALID ARGUMENTS: --prefilter accepts positive integers\n\n");
            help();
          }
          b_prefilter = true;
          prefilter2 = prefilter1;
        } else 
          help();
      }
      else if(strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0) {
        if(i < argc)
          nThreads = atoi(argv[++i]);
        else
          help();	
      } else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
        VERBOSE = true;
      }
      else if (strcmp(argv[i], "--useSHAPE") == 0){
        if( i < argc){
          shapeFile = argv[++i];
          SHAPE_ENABLED = true;
        }
        else
          help();
      }				
    } else {
      seqfile = argv[i];
    }
  }

  // Must have an input file specified
  if(seqfile.empty()) {
    help();
    printf("Missing input file.\n");
  }

  // If no output file specified, create one
  if(outputPrefix.empty()) {
    // base it off the input file
    outputPrefix += seqfile;

    size_t pos;
    // extract file name from the path
    if ((pos=outputPrefix.find_last_of('/')) > 0) {
      outputPrefix = outputPrefix.substr(pos+1);
    }

    // and if an extension exists, remove it ...
    if(outputPrefix.find(".") != string::npos)
      outputPrefix.erase(outputPrefix.rfind("."));
  }

  // If output dir specified
  if (!outputDir.empty()) {
    outputFile += outputDir;
    outputFile += "/";
  }
  // ... and append the .ct
  outputFile += outputPrefix;
  outputFile += ".ct";

}


void printRunConfiguration(string seq) {
	bool standardRun = true;

	printf("Run Configuration:\n");
	if (RNAMODE == true) {
		printf("+ running in rnafold mode\n");
		standardRun = false;
	} 
	if (UNAMODE == true) {
		printf("+ running in unafold mode\n");
		standardRun = false;
	}
	if (dangles == 0 || dangles == 2) {
		printf("+ running in dangle %d mode\n", dangles);
		standardRun = false;
	} 
	if (T_MISMATCH == true) {
		printf("+ enabled terminal mismatch calculations\n");
		standardRun = false;
	}
	if (b_prefilter == true) {
		printf("+ running with prefilter value = %d\n",prefilter1);
		standardRun = false;
	}

	if(!constraintsFile.empty()) {
		printf("- using constraint file: %s\n", constraintsFile.c_str());
		standardRun = false;
	}

	if(!shapeFile.empty()){
		printf("- using SHAPE data file: %s\n", shapeFile.c_str());
	}
	if (contactDistance != -1) {
		printf("- maximum contact distance: %d\n", contactDistance);
		standardRun = false;
	}

	if(standardRun)
		printf("- standard\n");

	printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
	printf("- input file: %s\n", seqfile.c_str());
	printf("- sequence length: %d\n", (int)seq.length());
	printf("- output file: %s\n", outputFile.c_str());
}

static void help() {
    printf("Usage: gtmfe [OPTION]... FILE\n\n");

    printf("   FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");

    printf("OPTIONS\n");
    printf("   -c, --constraints FILE\n");
    printf("                        Load constraints from FILE.  See Constraint syntax below.\n");
    printf("   -d, --dangle INT     Restricts treatment of dangling energies (INT=0,2),\n"); 
    printf("                        see below for details.\n");
    printf("   -h, --help           Output help (this message) and exit.\n");
    printf("   -l, --limitCD INT    Set a maximum base pair contact distance to INT. If no\n");
    printf("                        limit is given, base pairs can be over any distance.\n");
    printf("   -m  --mismatch       Enable terminal mismatch calculations\n");
//    printf("   -n, --noisolate      Prevent isolated base pairs from forming.\n");
    printf("   -o, --output NAME    Write output files with prefix given in NAME\n");
    printf("   -p  --paramdir DIR   Path to directory from which parameters are to be read\n");
    printf("   -t, --threads INT    Limit number of threads used to INT.\n");
    printf("   -v, --verbose        Run in verbose mode (includes loop-by-loop energy decomposition\n");
    printf("                        and confirmation of constraints satisfied).\n");
    printf("   -w, --workdir DIR    Path of directory where output files will be written.\n");
    printf("   --prefilter INT      Prohibits any basepair which does not have appropriate\n");
    printf("                        neighboring nucleotides such that it could be part of\n");
    printf("                        a helix of length INT.\n");
    printf("   --rnafold            Run as RNAfold default mode (ViennaPackage version 1.8.5).\n");
    printf("   --unafold            Run as UNAfold default mode (version 3.8), subject to traceback\n");
    printf("                        implementation.\n");

    printf("   -s, --useSHAPE FILE  Use SHAPE constraints from FILE.\n");      

    printf("\nConstraint syntax:\n");
    printf("\tF i j k  # force (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs.\n");
    printf("\tP i j k  # prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs.\n");
    printf("\tP i 0 k  # make bases from i to i+k-1 single stranded bases.\n");

    printf("\nDangle:\n");
    printf("\tINT=0 ignores dangling energies (mostly for debugging).\n");
    printf("\tINT=2 dangling energies are added for nucleotides on either\n");
    printf("\tside of each branch in multi-loops and external loops.\n");
    printf("\tAll other values for INT are ignored.\n");
    exit(-1);
}

