#include "options.h"

//#include <sstream>

//using namespace std;

int PFMODE=0;//boltzman sampling or stochastic sampling and  partition function mode or dS mode, it is mode as defined and used for partition function of sfold
int NODANGLEMODE=0;//no dangling at all means d0
int D2MODE=0;//d2 mode
int DEFAULTMODE=1;//default mode

char seqfile[200];
char paramDir[200] = {'\0'};

/**
 * Print the help message and quit.
 */
void help() {
    printf("Usage: RNAScoring [OPTION]... FILE\n\n");

    printf("   FILE is an RNA sequence structure file (ct file) containing the structure in standard tabular format.\n\n");

    printf("OPTIONS\n");
    printf("   --dS            Calculate score based on partition function energy calculation rules.\n");
    printf("   --d0            Calculate score based on zero dangling energy.\n");
    printf("   --d2            Calculate score based on d2 energy calculation rules.\n");
    printf("   --pf_test <summary_file_path> <output_file_path> <error_file_path> [-dS|d2]		Here summary_file_path is the file that is summary of all samples from stochastic sampling. After program completes output_file_path will contains all different scores (dS,d0,d2) for all structures. error_file_path will be containing entries (lines) for structures, for which partition function energies are different from -dS score. if error_file_path file is empty then it has passed the test, other wise it has failed this test and you can refer end-user to err_file_path for further investigation, any option can be provided as -dS or -d2 in which partition function is calculated. default value is -dS.\n");
    printf("   --param-dir  paramDir        optional parameter to explicitly define parameter directory different than current directory.\n");
    exit(-1);
}

/**
 * Parse the options from argc and argv and save them into global state.
 */
void parse_options(int argc, char** argv) {
  int i;

  for(i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
        help();
      }  
      else if (strcmp(argv[i], "--dS") == 0) {
        PFMODE = 1;
	DEFAULTMODE=0;
      }
      else if(strcmp(argv[i], "--d0") == 0){
	NODANGLEMODE = 1;
	DEFAULTMODE=0;
      }
	else if(strcmp(argv[i], "--d2") == 0){
        D2MODE=1;
	DEFAULTMODE=0;
      }
	else if(strcmp(argv[i], "--param-dir")==0){
	 if(i<argc){
	  i++;
	  strcpy(paramDir, argv[i]);
	 }
	 else help();
	}
    }
    else {
      strcpy(seqfile, argv[i]);
    }

  }

  // Must have an input file specified
  if(strlen(seqfile)==0) {
    help();
    printf("Missing input file.\n");
  }
}
