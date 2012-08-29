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
static bool SILENT = false;
//static bool SHAPE_ENABLED = false;

//extern int SHAPE_ENABLED;

static bool T_MISMATCH = false;
static bool UNAMODE = false;
static bool RNAMODE = false;
static bool b_prefilter = false;

static string seqfile = "";
static string constraintsFile = "";
static string outputPrefix = "";
static string outputFile = "";
static string energyDecomposeOutFile = "";
static string outputDir = "";
static string shapeFile = "";
static string paramDir; // default value

static int dangles=-1;
static int prefilter1=2;
static int prefilter2=2;
static int print_energy_decompose = 0;

static int nThreads = -1;
static int contactDistance = -1;

static void help();
static void detailed_help();
static void printRunConfiguration(string seq);
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
    if (T_MISMATCH) if(!SILENT) printf("Ignoring -m option, using --unafold\n");
    if (PARAM_DIR) if(!SILENT) printf("Ignoring -p option, using --unafold\n");
    if (dangles == 0 || dangles == 1 || dangles == 2) 
      if(!SILENT) printf("Ignoring -d option, using --unafold\n");
    if (b_prefilter == 1) 
      if(!SILENT) printf("Ignoring --prefilter option, using --unafold\n");
    T_MISMATCH = false;
    PARAM_DIR = false;
    dangles = -1;
    b_prefilter = false;
  }

  if (RNAMODE) {
    if (T_MISMATCH) if(!SILENT) printf("Ignoring -m option, using --rnafold\n");
    if (PARAM_DIR) if(!SILENT) printf("Ignoring -p option, using --rnafold\n");
    if (dangles == 0 || dangles == 1 || dangles == 2) 
      if(!SILENT) printf("Ignoring -d option, using --rnafold\n");
    if (b_prefilter == 1) 
      if(!SILENT) printf("Ignoring --prefilter option, using --rnafold\n");
    T_MISMATCH = false;
    PARAM_DIR = false;
    dangles = -1;
    b_prefilter = false;
  }

  if ((dangles == 0 || dangles == 1 ||dangles == 2) && !UNAMODE && !RNAMODE) {
    if (T_MISMATCH) if(!SILENT) printf("Ignoring -m option, using -d option\n");
    T_MISMATCH = false;
  } else {
    if (dangles != -1 && !UNAMODE && !RNAMODE) if(!SILENT) printf("Ignoring -d as it accept 0 1 or 2 only\n");	
    dangles = -1;
  }
  if(dangles==1) dangles=-1;

  g_nthreads = nThreads;
  g_unamode  = UNAMODE;
  g_mismatch = T_MISMATCH;
  g_verbose  = VERBOSE;
  g_prefilter_mode  = b_prefilter;
  g_prefilter1  = prefilter1;
  g_prefilter2  = prefilter2;
  g_dangles = dangles;

#ifdef DEBUG
  if(!SILENT) printf("g_nthreads = %d\n", g_nthreads);
  if(!SILENT) printf("g_unamode = %d\n", g_unamode);
  if(!SILENT) printf("g_mismatch = %d\n", g_mismatch);
  if(!SILENT) printf("g_prefilter_mode = %d\n", g_prefilter_mode);
  if(!SILENT) printf("g_dangles = %d\n", g_dangles);

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
	
	//if(!SILENT) print_header();

	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}

	init_fold(seq.c_str());
	
	// Read in thermodynamic parameters. Always use Turner99 data (for now)
  readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, UNAMODE, RNAMODE, T_MISMATCH);
	printRunConfiguration(seq);

	if(!SILENT) printf("\nComputing minimum free energy structure...\n");
	fflush(stdout);

	double t1 = get_seconds();
	energy = calculate(seq.length()) ; 
	t1 = get_seconds() - t1;
	
	if(!SILENT) printf("Done.\n\n");
	if(!SILENT) printf("Results:\n");
	if (energy >= MAXENG)	
		printf("- Minimum Free Energy: %12.4f kcal/mol\n", 0.00);
	else
		printf("- Minimum Free Energy: %12.4f kcal/mol\n", energy/100.00);
	printf("- MFE runtime: %9.6f seconds\n", t1);

	t1 = get_seconds();
	trace(seq.length(), print_energy_decompose, energyDecomposeOutFile.c_str());
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
		if(!SILENT) printf("Verifying that structure fulfills constraint criteria... ");
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
      } else if(strcmp(argv[i], "--detailedhelp") == 0 ) {
	detailed_help();
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
      }else if(strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) {
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
          if (!(dangles == 0 || dangles == 1 || dangles == 2)) {
            dangles = -1;
            printf("Ignoring %s option as it accepts either 0 1 or 2\n", cmd.c_str());
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
      } else if (strcmp(argv[i], "--silent") == 0) {
        SILENT = true;
      } else if (strcmp(argv[i], "--energydetail") == 0 || strcmp(argv[i], "-e") == 0) {
	print_energy_decompose = 1;
      }
      else if (strcmp(argv[i], "--useSHAPE") == 0){
        if( i < argc){
          shapeFile = argv[++i];
          //SHAPE_ENABLED = true;
          SHAPE_ENABLED = 1;
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
 	energyDecomposeOutFile += outputDir;
	energyDecomposeOutFile += "/";

  }
  // ... and append the .ct
  outputFile += outputPrefix;
  outputFile += ".ct";

  energyDecomposeOutFile += outputPrefix;
  energyDecomposeOutFile += ".energy";

}


static void printRunConfiguration(string seq) {
	bool standardRun = true;

	if(!SILENT) printf("Run Configuration:\n");
	if (RNAMODE == true) {
		if(!SILENT) printf("+ running in rnafold mode\n");
		standardRun = false;
	} 
	if (UNAMODE == true) {
		if(!SILENT) printf("+ running in unafold mode\n");
		standardRun = false;
	}
	if (dangles == 0 || dangles == 1 || dangles == 2) {
		if(!SILENT) printf("+ running in dangle %d mode\n", dangles);
		standardRun = false;
	} 
	if (T_MISMATCH == true) {
		if(!SILENT) printf("+ enabled terminal mismatch calculations\n");
		standardRun = false;
	}
	if (b_prefilter == true) {
		if(!SILENT) printf("+ running with prefilter value = %d\n",prefilter1);
		standardRun = false;
	}

	if(!constraintsFile.empty()) {
		if(!SILENT) printf("- using constraint file: %s\n", constraintsFile.c_str());
		standardRun = false;
	}

	if(!shapeFile.empty()){
		if(!SILENT) printf("- using SHAPE data file: %s\n", shapeFile.c_str());
	}
	if (contactDistance != -1) {
		if(!SILENT) printf("- maximum contact distance: %d\n", contactDistance);
		standardRun = false;
	}

	if(standardRun)
		if(!SILENT) printf("- standard\n");

	if(!SILENT) printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
	if(!SILENT) printf("- input file: %s\n", seqfile.c_str());
	if(!SILENT) printf("- sequence length: %d\n", (int)seq.length());
	printf("- output file: %s\n", outputFile.c_str());
}

static void print_usage_developer_options() {
    printf("\nSetting default parameter directory:\n");
    printf("\tTo run properly, GTfold requires access to a set of parameter files. If you are using one of the prepackaged binaries, you may need (or chose) to \n");
    printf("\tset the GTFOLDDATADIR environment variable to specify the directory in whihc GTfold should look to find default parameter files. In a terminal \n");
    printf("\twindow, use either the command \n");
    printf("\t\texport GTFOLDDATADIR=DIR\n");
    printf("\t\tfor BASH shell users, or \n");
    printf("\t\tsetenv GTFOLDDATADIR=DIR\n");
    printf("\t\tfor tcsh shell users. Alternatively, you may use the --paramdir option described above. \n");
    printf("\tGTfold will by default look for parameter files in the following directories: \n");
    printf("\t\t(1)      The directory pointed to by environment variable GTFOLDDATADIR \n");
    printf("\t\t(2)      The install directory (eg. /usr/local/share/gtfold), if (1) fails. \n");
    printf("\t\t(3)      The subdirectory 'data' of the current directory, if (1) and (2) fail. \n");
    printf("\n");
}

static void print_usage() {
    printf("Usage: gtmfe [OPTION]... FILE\n\n");

    printf("   FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");

    printf("OPTIONS\n");
    printf("   -c, --constraints FILE\n");
    printf("                        Load constraints from FILE.  See Constraint syntax below.\n");
    printf("   -d, --dangle INT     Restricts treatment of dangling energies (INT=0,1,2), (with -d option, call to -m option will be ignored)\n"); 
    printf("                        see below for details.\n");
    printf("   -h, --help           Output help (this message) and exit.\n");
    printf("   --detailedhelp      Output help (this message) with detailed options and examples, and exit.\n");
    printf("   -l, --limitCD INT    Set a maximum base pair contact distance to INT. If no\n");
    printf("                        limit is given, base pairs can be over any distance.\n");
    printf("   -m  --mismatch       Enable terminal mismatch calculations\n");
//    printf("   -n, --noisolate      Prevent isolated base pairs from forming.\n");
    printf("   -o, --output NAME    Write output files with prefix given in NAME\n");
    printf("   -p  --paramdir DIR   Path to directory from which parameters are to be read\n");
    printf("   -t, --threads INT    Limit number of threads used to INT.\n");
    printf("   -v, --verbose        Run in verbose mode (includes confirmation of constraints satisfied).\n");
    //printf("   --silent		    Run in silent mode.\n");
    printf("   -w, --workdir DIR    Path of directory where output files will be written.\n");
    printf("   --prefilter INT      Prohibits any basepair which does not have appropriate\n");
    printf("                        neighboring nucleotides such that it could be part of\n");
    printf("                        a helix of length INT.\n");
    printf("   --rnafold            Run as RNAfold default mode (ViennaPackage version 1.8.5).\n");
    printf("   --unafold            Run as UNAfold default mode (version 3.8), subject to traceback\n");
    printf("                        implementation.\n");

    printf("   --useSHAPE FILE  Use SHAPE constraints from FILE.\n");      
    printf("   -e, --energydetail         prints energy decomposition for MFE structure to file output-prefix.energy.\n");
    printf("\nConstraint syntax:\n");
    printf("\tP i j k  # prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs.\n");
    printf("\tF i j k  # force (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs.\n");
    printf("\tP i 0 k  # make bases from i to i+k-1 single stranded bases.\n");
    printf("\tF i 0 k  # forces bases from i to i+k-1 single stranded bases.\n");

    printf("\nDangle:\n");
    printf("\tINT=0 ignores dangling energies (mostly for debugging).\n");
    printf("\tINT=1 default mode of dangling energies (consider energetically favourable).\n");
    printf("\tINT=2 dangling energies are added for nucleotides on either\n");
    printf("\tside of each branch in multi-loops and external loops.\n");
    printf("\tAll other values for INT are ignored.\n");
    
    printf("\nSHAPE syntax:\n");
    printf("\tSHAPE values should be given in a file with two space-delimited columns, for example\n");
    printf("\t\t1 0.1\n");
    printf("\t\t2 0.001 \n");
    printf("\t\t3 1.67 \n");
    printf("\t\tetc.,\n");
    printf("\twhere the first column is the nucleotide position (INT) and the second column is the SHAPE reactivity[1] (DOUBLE) for that position. The file \n");
    printf("\tshould have no header. Not all positions need to be included in the file, and the values do not need to be in order of increasing position. Negative\n");
    printf("\tSHAPE reactivities are ignored. \n");
}

static void print_examples(){
        printf("\n\nEXAMPLES:\n\n");
        printf("1. Calculate Minimum Free Energy Structure:\n\n");
        printf("gtmfe [-c FILE] [-d 0|1|2] [-t n] [-o outputPrefix] [-v] [-p DIR] [-w DIR] [-l] [-m] [--prefilter INT] [--useSHAPE FILE] [-e] <seq_file>\n\n");
        printf("gtmfe [--unafold] [--rnafold] [-c FILE] [-t n] [-o outputPrefix] [-v] [-p DIR] [-w DIR] [-l] [-m] [--prefilter INT] [--useSHAPE FILE] [-e] <seq_file>\n\n");
        printf("\n");
}

static void help(){
	print_usage();
        print_examples();
	exit(-1);
}

static void detailed_help(){
        print_usage();
        print_examples();
	print_usage_developer_options();
        exit(-1);
}
