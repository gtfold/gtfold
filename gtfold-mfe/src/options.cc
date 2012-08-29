#include "loader.h"
#include "options.h"

#include <sstream>

using namespace std;

bool ILSA;
bool NOISOLATE;
bool PARAM_DIR = false;
bool LIMIT_DISTANCE;
bool BPP_ENABLED;
bool SUBOPT_ENABLED;
bool CONS_ENABLED = false;
bool VERBOSE = false;
//bool SHAPE_ENABLED = false;
bool T_MISMATCH = false;
bool UNAMODE = false;
bool RNAMODE = false;
bool b_prefilter = false;
bool CALC_PART_FUNC = false;
bool RND_SAMPLE = false;
bool PF_COUNT_MODE = false;

string seqfile = "";
string constraintsFile = "";
string outputPrefix = "";
string outputFile = "";
string suboptFile = "";
string bppOutFile = "";
string outputDir = "";
string shapeFile = "";
string paramDir; // default value

int num_rnd = 0;
int dangles=-1;
int prefilter1=2;
int prefilter2=2;

float suboptDelta = 0.0;
int nThreads = -1;
int contactDistance = -1;

/**
* Print the help message and quit.
*/

/*
void help() {
printf("Usage: gtfold [OPTION]... FILE\n\n");

printf(" FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");

printf("OPTIONS\n");
printf(" -c, --constraints FILE\n");
printf(" Load constraints from FILE. See Constraint syntax below.\n");
printf(" -d, --dangle INT Restricts treatment of dangling energies (INT=0,2),\n");
printf(" see below for details.\n");
printf(" -h, --help Output help (this message) and exit.\n");
printf(" -l, --limitCD INT Set a maximum base pair contact distance to INT. If no\n");
printf(" limit is given, base pairs can be over any distance.\n");
printf(" -m --mismatch Enable terminal mismatch calculations\n");
// printf(" -n, --noisolate Prevent isolated base pairs from forming.\n");
printf(" -o, --output NAME Write output files with prefix given in NAME\n");
printf(" -p --paramdir DIR Path to directory from which parameters are to be read\n");
printf(" -t, --threads INT Limit number of threads used to INT.\n");
printf(" -v, --verbose Run in verbose mode (includes loop-by-loop energy decomposition\n");
printf(" and confirmation of constraints satisfied).\n");
printf(" -w, --workdir DIR Path of directory where output files will be written.\n");
printf(" --prefilter INT Prohibits any basepair which does not have appropriate\n");
printf(" neighboring nucleotides such that it could be part of\n");
printf(" a helix of length INT.\n");
printf(" --rnafold Run as RNAfold default mode (ViennaPackage version 1.8.5).\n");
printf(" --unafold Run as UNAfold default mode (version 3.8), subject to traceback\n");
printf(" implementation.\n");

printf("\nBETA OPTIONS\n");
printf(" --bpp Calculate base pair probabilities.\n");
printf(" --partition Calculate the partition function.\n");
printf(" --pf_count Calculate the structure count using partition function and zero energy value.\n");
printf(" --subopt NUM Calculate suboptimal structures within NUM kcal/mol\n");
printf(" of the MFE. (Uses -d 2 treatment of dangling energies.)\n");
printf(" -s, --useSHAPE FILE Use SHAPE constraints from FILE.\n");

printf("\nConstraint syntax:\n");
printf("\tF i j k # force (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs.\n");
printf("\tP i j k # prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs.\n");
printf("\tP i 0 k # make bases from i to i+k-1 single stranded bases.\n");

printf("\nDangle:\n");
printf("\tINT=0 ignores dangling energies (mostly for debugging).\n");
printf("\tINT=2 dangling energies are added for nucleotides on either\n");
printf("\tside of each branch in multi-loops and external loops.\n");
printf("\tAll other values for INT are ignored.\n");
exit(-1);
}
*/

/**
* Parse the options from argc and argv and save them into global state.
*/
/*
void parse_options(int argc, char** argv) {
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
if (contactDistance >= 0 && !strcmp(ss.str().c_str(),argv[i]))
LIMIT_DISTANCE = true;
else
help();
}
else
help();
} else if(strcmp(argv[i], "--noisolate") == 0 || strcmp(argv[i], "-n") == 0) {
NOISOLATE = true;
} else if(strcmp(argv[i], "--prefix") == 0 || strcmp(argv[i], "-o") == 0) {
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
else if(strcmp(argv[i], "--bpp") == 0) {
BPP_ENABLED = true;
} else if(strcmp(argv[i], "--subopt") == 0) {
SUBOPT_ENABLED = true;
dangles = 2;
if(i < argc)
suboptDelta = atof(argv[++i]);
else
help();
} else if (strcmp(argv[i],"--partition") == 0) {
CALC_PART_FUNC = true;
} else if (strcmp(argv[i],"--pf_count") == 0) {
CALC_PART_FUNC = true;
PF_COUNT_MODE = true;
} else if (strcmp(argv[i],"--sample") == 0) {
RND_SAMPLE = true;
if(i < argc)
num_rnd = atoi(argv[++i]);
else
help();
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
suboptFile += outputDir;
suboptFile += "/";
bppOutFile += outputDir;
bppOutFile += "/";
}
// ... and append the .ct
outputFile += outputPrefix;
outputFile += ".ct";

suboptFile += outputPrefix;
suboptFile += "_ss.txt";

bppOutFile += outputPrefix;
bppOutFile += "_bpp.txt";
}
*/
/**
* Prints the run configuration for this run.
*
* The lines that start with a '-' are normal options, the '+' are beta options.
*/
/*
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
if (NOISOLATE == true) {
printf("- preventing isolated base pairs\n");
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

if (BPP_ENABLED == true) {
printf("+ calculating base pair probabilities\n");
printf("+ BPP output file: %s\n", bppOutFile.c_str());
standardRun = false;
}

if (SUBOPT_ENABLED) {
printf("+ calculating suboptimal structures within %f kcal/mol of MFE\n", suboptDelta);
printf("+ suboptimal structures file: %s\n", suboptFile.c_str());
standardRun = false;
}

if(standardRun)
printf("- standard\n");

printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
printf("- input file: %s\n", seqfile.c_str());
printf("- sequence length: %d\n", (int)seq.length());
printf("- output file: %s\n", outputFile.c_str());
}
*/
