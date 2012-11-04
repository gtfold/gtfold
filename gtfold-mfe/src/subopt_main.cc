#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

#include "loader.h"
#include "algorithms.h"
#include "subopt_traceback.h"
#include "global.h"
#include "utils.h"
#include "mfe_main.h"

using namespace std;

static bool SILENT = false;
static string seqfile = "";
static string suboptFile = "";
static double suboptDelta = 0.0;
static string outputPrefix = "";
static string outputFile = "";
static string outputDir = "";
static string paramDir = "";
static bool PARAM_DIR = false;
static int is_check_for_duplicates_enabled = 0;//TODO Default value need to be confirmed with Shel
static int max_structure_count = -1;//User can optionally use --maxcount max_structure_count to restrict program to generate only that many structure and quits after that 
static void help();
static void detailed_help();
static void printRunConfiguration(string seq);

void save_subopt_file(string outputFile, ss_map_t& ss_data, 
		const string& seq, int energy)
{
	ofstream outfile;
	outfile.open(outputFile.c_str());
	char buff[4096];

	sprintf(buff,"%s %6.2f", seq.c_str(), energy/100.0);
	outfile << buff << std::endl;
	for (ss_map_t::iterator it = ss_data.begin(); it!= ss_data.end(); ++it) 
	{
		sprintf(buff,"%s %6.2f", (it->first).c_str(), it->second/100.0);
		//outfile << it->first << '\t' << it->second/100.0 << std::endl;
		outfile << buff << std::endl;
	}

	outfile.close();
}

static void write_header_subopt_file(string outputFile, const string& seq, int energy)
{
	ofstream outfile;
	outfile.open(outputFile.c_str());
	char buff[4096];

	sprintf(buff,"count\t%s\t%6.2f", seq.c_str(), energy/100.0);
	outfile << buff << std::endl;

	outfile.close();
}

void parse_options(int argc, char** argv) {
  int i;
  g_dangles = 2;

  for(i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      if(strcmp(argv[i], "--help") == 0 || 
          strcmp(argv[i], "-h") == 0) {
        help(); 
      } else if(strcmp(argv[i], "--detailedhelp") == 0 ) {
	detailed_help();
      } else if(strcmp(argv[i], "--duplicatecheck") == 0 ) {
	is_check_for_duplicates_enabled = 1;
      } else if (strcmp(argv[i], "--paramdir") == 0 || 
          strcmp(argv[i], "-p") == 0) {
        if(i < argc) {
          paramDir = argv[++i];
          PARAM_DIR = true;
        }
        else {
          help();
        }
      } else if (strcmp(argv[i], "--workdir") == 0 || strcmp(argv[i], "-w") == 0) {
        if(i < argc)
          outputDir = argv[++i];
        else
          help();
      } else if(strcmp(argv[i], "--delta") == 0) {
        g_dangles = 2;
        if(i < argc)
          suboptDelta = atof(argv[++i]);
        else
          help();
      } else if(strcmp(argv[i], "-o") == 0) {
		outputPrefix.assign(argv[++i]);
      } else if (strcmp(argv[i], "--dangle") == 0 || strcmp(argv[i], "-d") == 0) {
        std::string cmd = argv[i];
	if(i < argc) {
          g_dangles = atoi(argv[++i]);
          if (g_dangles != 2) {
            g_dangles = 2;
            printf("Ignoring %s option as it accepts only 2 and program will continue with dangles value as 2\n", cmd.c_str());
          }
        } else
          help();
      } else if (strcmp(argv[i], "--maxcount") == 0) {
	if(i+1 < argc) {
          max_structure_count = atoi(argv[++i]);
          if (max_structure_count < 0) {
            max_structure_count = -1;
            printf("Ignoring --maxcount option as it accepts only positive numbers and program will continue with maxcount value as -1 which means there is no restriction on maximum count of structures\n");
          }
        } else
          help();
      }
    } else {
      seqfile = argv[i];
    }
  }

  if(seqfile.empty()) {
    printf("Missing input file.\n");
    help();
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
    suboptFile += outputDir;
    suboptFile += "/";
  }
  // ... and append the .ct
  outputFile += outputPrefix;
  outputFile += ".ct";
  suboptFile += outputPrefix;	
  suboptFile += "_ss.txt";	

	printf("Output: %s %s %s\n", outputPrefix.c_str(), outputFile.c_str(), suboptFile.c_str());
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
    printf("Usage: gtsubopt [OPTION]... FILE\n\n");

    printf("   FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");

    printf("OPTIONS\n");
    printf("   --delta DOUBLE         Calculate suboptimal structures within DOUBLE kcal/mol\n");
    printf("                        of the MFE. (Uses -d 2 treatment of dangling energies.)\n");
    printf("\n"); 
    printf("   -d, --dangle INT     Restricts treatment of dangling energies (INT=2),\n");
    printf("   -o, --output NAME    Write output files with prefix given in NAME\n");
    printf("   -p  --paramdir DIR   Path to directory from which parameters are to be read\n");
    printf("   -h, --help           Output help (this message) and exit.\n");
    printf("   --detailedhelp      Output help (this message) with detailed options and examples, and exit.\n");
    printf("   -w, --workdir DIR    Path of directory where output files will be written.\n");
    printf("   -v, --verbose        Run in verbose mode.\n");
    printf("   --duplicatecheck	    Enable the check if duplicate structure is coming or not, if duplicate is coming then warn the user and exits. By default this option will switched off. This option will slowdown program as well as consume more memory because of need to storing all structures.\n");
    printf("   --maxcount INT	    Optional option '--maxcount max_structure_count' to restrict program to generate only that many structure and quits after that. By default there is no maximum limit\n");
}

static void print_examples(){
        printf("\n\nEXAMPLES:\n\n");
        printf("1. Calculate Suboptimal Structures:\n");
        printf("gtsubopt --delta DOUBLE [-d 2] [-o outputPrefix] [--duplicatecheck] [--maxcount INT] [-v] [-w DIR] [-p DIR] <seq_file>\n\n");
        printf("\n\n");
}

static void help() {
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

void subopt_main(int argc, char** argv) {

  string seq = "";
  parse_options(argc, argv);

  if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
    printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
    exit(-1);
  }
  init_fold(seq.c_str());
  g_dangles = 2;  
  readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, 0, 1, 0);
  printRunConfiguration(seq);

  int energy = calculate(seq.length()) ; 
  write_header_subopt_file(suboptFile, seq, energy);	
  
  double t1 = get_seconds();
  ss_map_t subopt_data = subopt_traceback(seq.length(), 100.0*suboptDelta, suboptFile, is_check_for_duplicates_enabled, max_structure_count);
  t1 = get_seconds() - t1;
  
  //printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
//	printf("- input file: %s\n", seqfile.c_str());
//	printf("- sequence length: %d\n", (int)seq.length());
  //printf("\n");
  printf("Subopt traceback running time: %9.6f seconds\n", t1);

  printf("Subopt structures saved in %s\n", suboptFile.c_str());
  //save_subopt_file(suboptFile, subopt_data, seq, energy);	

  //printf("+ calculating suboptimal structures within %f kcal/mol of MFE\n", suboptDelta);
  //printf("+ suboptimal structures file: %s\n", suboptFile.c_str());

	free_fold(seq.length());
	printf("\n");
}

static void printRunConfiguration(string seq) {

        if(!SILENT) printf("\nRun Configuration:\n");
        
        if (g_dangles == 2) {
                if(!SILENT) printf("+ running in dangle %d mode\n", g_dangles);
        }
        

        if(!SILENT) printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
        if(!SILENT) printf("- input file: %s\n", seqfile.c_str());
        if(!SILENT) printf("- sequence length: %d\n", (int)seq.length());
	if(!SILENT) printf("+ calculating suboptimal structures within %f kcal/mol of MFE\n", suboptDelta);
	if(!SILENT) printf("+ suboptimal structures file: %s\n", suboptFile.c_str());
	if(!SILENT) printf("+ is_check_for_duplicates_enabled: %d\n", is_check_for_duplicates_enabled);
        //printf("- output file: %s\n", outputFile.c_str());
	printf("\n");
}



