#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include<sstream>
//#include <sys/time.h>
//#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "loader.h"
#include "algorithms-partition.h"
#include "boltzmann_main.h"
#include "partition-func.h"
#include "partition-func-d2.h"
#include "mfe_main.h"
#include "stochastic-sampling.h"
#include "stochastic-sampling-d2.h"
#include "algorithms.h"
#include "traceback.h"
#include "utils.h"

using namespace std;

static bool SILENT = false;
static bool CALC_PART_FUNC = true;
static bool PF_COUNT_MODE = false;
static bool BPP_ENABLED = false;
static bool PARAM_DIR = false;
static bool RND_SAMPLE = false;
static bool DUMP_CT_FILE = false;
static bool CALC_PF_DO = false;
static bool CALC_PF_DS = false;
static bool CALC_PF_D2 = true;//making D2 mode as default option
static bool PF_D2_UP_APPROX_ENABLED = true;//making short internal loop code to run as default
static bool PF_PRINT_ARRAYS_ENABLED = false;//making short internal loop code to run as default
static bool ST_D2_ENABLE_COUNTS_PARALLELIZATION = true;//making parallelization of sample counts as default
static bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION = false;
static bool ST_D2_ENABLE_SCATTER_PLOT = false;
static bool ST_D2_ENABLE_UNIFORM_SAMPLE = false;
static double ST_D2_UNIFORM_SAMPLE_ENERGY = 0.0;
static bool ST_D2_ENABLE_CHECK_FRACTION = false;
static bool ST_D2_ENABLE_BPP_PROBABILITY = false;

static string seqfile = "";
static string outputPrefix = "";
static string outputDir = "";
static string outputFile = "";
static string paramDir; // default value
static string bppOutFile = "";
static string sampleOutFile = "";
static string energyDecomposeOutFile = "";
static string estimateBppOutputFile = "";
static string scatterPlotOutputFile = "";
static string pfArraysOutFile = "";
static string ctFileDumpDir = "";
static string stochastic_summery_file_name = "stochaSampleSummary.txt";

static int num_rnd = 0;
static int ss_verbose_global = 0;
static int print_energy_decompose = 0;
static int dangles=2;//making dangle default value as 2

static bool LIMIT_DISTANCE = false;
static int contactDistance = -1;

static void printRunConfiguration(string seq);


static void print_usage() {
	printf("Usage: gtboltzmann [OPTION]... FILE\n\n");

	printf("   FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");

	printf("OPTIONS\n");

	printf("   -d, --dangle INT	Restricts treatment of dangling energies (INT=0,2). See below for details.\n");
	//printf("   --detailedhelp      Output help (this message) with detailed options and examples, and exit.\n");
	printf("   --detailedhelp       Display detailed help message. Includes examples and additional options useful to developers.\n");
	//printf("   -e, --energydetail         prints energy decomposition for sampled structures to file with extention '.energy' (should be used with '-t 1' option, as otherwise all threads in parallel, will write to file and output will be intermixed from all threads).\n");
	printf("   -e, --energydetail   Writes loop-by-loop energy decomposition of structures to\n");
	printf("			output-prefix.energy. When using this function in combination\n");
	printf("			with --sample, number of threads must be limited to one (-t 1). .\n");
	printf("   --estimatebpp	Writed a csv file containing, for each sampled base pair, that base pair and it's frequency\n");
        printf("			to output-prefix.sbpp. This option is ignored if not using --sample.\n");
	printf("   --groupbyfreq        Write a csv file (output-prefix.frequency) containing, for each sampled structure, a line with\n");
	printf("			the structure's probability under the Boltzmann Distribution followed by the normalized frequency\n");
	printf("			of that structure, where (normalized frequency) = (structure frequency)/(number of structures sampled)\n");
	printf("			Only valid in combination with --sample.\n");
	printf("   -h, --help           Output help (this message) and exit.\n");
	printf("   -l|--limitcd  INT	Set a maximum base pair contact distance to INT. If no\n");
	printf(" 		      	limit is given, base pairs can be over any distance.\n");
	printf("   -o, --output NAME    Write output files with prefix given in NAME\n");
	printf("   -p  --paramdir DIR   Path to directory from which parameters are to be read\n");
	//printf("   --pfcount		Calculate the structure count using partition function and zero energy value.\n");
	printf("   --pfcount		Output the number of possible structures (using partition function).\n");
	//printf("   -s|--sample   INT	Sample number of structures equal to INT.\n");
	printf("   -s|--sample   INT	Sample INT structures from Boltzmann distribution. Writes structures to file output-prefix.samples.\n");
	printf("   -t|--threads INT	Limit number of threads used to INT.\n");
	printf("   -v, --verbose	Run in verbose mode (includes partition function table printing.)\n");
	//printf("   --estimatebpp	While sampling structures, Calculate base pair probabilities.\n");
        //printf("   --groupbyfreq        While sampling structures, Collect frequency of all structures and calculate estimate probability and boltzmann probability for scatter plot.\n");
	printf("   -w, --workdir DIR    Path of directory where output files will be written.\n");
	//printf("   --bpp		Calculate base pair probabilities.\n");
	printf("   --bpp		Calculate base pair probabilities for the predicted structure\n");
	printf("			and print to output-prefix.bpp. (Beta option)\n");
}

static void print_usage_developer_options() {
	printf("\n\nDeveloper OPTIONS\n");
        printf("   --partition          Calculate the partition function (default is using d2 dangling mode).\n");
        //printf("   --printarrays       Print the partition function arrays to outputPrefix.pfarrays file.\n");
        printf("   --printarrays        Writes partition function arrays to prefix.pfarrays. \n");
        //printf("   --exactintloop       Do the exact internal loop calculation while calculating partition function and traceback without any short internal loop approximation)\n");
        printf("   --exactintloop       Includes structures with abitrarily many unpaired nucleotides in internal loops.\n");
        printf("                        Note: using this option increases the running time by a factor of N,\n");
	printf("                        where N is the sequence length.\n");
	printf("   --checkfraction	While sampling structures, enable check for each structure, the probability used in stochastic\n");
        printf("   -dS                  Calculate the partition function using sfold reccurences and use them in traceback.\n");
	printf("			sampling matches the probability of that structure according to the Boltzmann Distribution.\n");
        //printf("   --sampleenergy DOUBLE      While sampling structures, Samples with Energy energy1 will only be sampled.\n");
        printf("   --sampleenergy DOUBLE      Writes only sampled structures with free energy equal to DOUBLE to file prefix.sample. Only valid in combination with --sample.\n");
        //printf("   --counts-parallel  While sampling structures, parallelize INT sample counts among available threads (this is also a default behaviour of sampling).\n");
        //printf("   --parallelsample        While sampling structures, parallelize the processing of one sample (useful when sampling large sequence with number of samples being less than available threads).\n");
        printf("   --parallelsample     Paralellizes the sampling of each individual structure.\n");
	printf("			Only valid in combination with --sample.\n");
        //printf("   -s|--sample   INT  --separatectfiles [--ctfilesdir dump_dir_path] [--summaryfile dump_summery_file_name] Sample number of structures equal to INT and dump each structure to a ct file in dump_dir_path directory (if no value provided then use current directory value for this purpose) and also create a summary file with name stochastic_summery_file_name in dump_dir_path directory (if no value provided, use stochaSampleSummary.txt value for this purpose).\n");
        printf("   --separatectfiles [--ctfilesdir DIR] [--summaryfile NAME] Writes each sampled structure to a separate .ct file \n");
	printf("			in the DIR directory. Also writes a summary of the sampled structures to NAME in DIR.\n");
	printf("                        Default directory is the working directory specified with -w, and the default summary file\n");
	 printf("                       name is stochaSampleSumary.txt Only valid in combination with --sample. \n");

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

static void print_examples(){
        printf("\n\nEXAMPLES:\n\n");
        printf("1. Sample structures stochastically:\n\n");
        printf("gtboltzmann [-s INT] [[-d 0|2]|[-dS]] [-t n] [-o outputPrefix] [-v] [--estimatebpp] [-p DIR] [-w DIR] [-l] <seq_file>\n\n");
        printf("2. Calculate base pair probabilities:\n\n");
	printf("gtboltzmann --bpp [-d 2] [-o outputPrefix] [-v] [-p DIR] [-w DIR] [-l] <seq_file>\n\n");
        printf("\n\n");
}

static void print_examples_developer_options(){
        printf("\n\nDeveloper Options EXAMPLES:\n\n");
        printf("1. Calculate Partition function:\n\n");
        printf("gtboltzmann --partition [[-d 0|2]|[-dS]] [-t n] [-o outputPrefix] [--exactintloop] [-v] [-p DIR] [-w DIR] [-l] <seq_file>\n\n");
        printf("2. Sample structures stochastically:\n\n");
        printf("gtboltzmann [-s INT] [[-d 0|2]|[-dS]] [-t n] [-o outputPrefix] [--exactintloop] [-v] [--groupbyfreq] [--sampleenergy DOUBLE] [--estimatebpp] [--parallelsample] [-p DIR] [-w DIR] [-l] <seq_file>\n\n");
        printf("gtboltzmann [-s INT] [[-d 0|2]|[-dS]] -t 1 [-o outputPrefix] [--exactintloop] [-v] [--groupbyfreq] [--sampleenergy DOUBLE] [-e] [--checkfraction] [--estimatebpp] [--parallelsample] [-p DIR] [-w DIR] [-l] <seq_file>\n\n");
        printf("gtboltzmann [-s] INT --separatectfiles [--ctfilesdir dump_dir_path] [--summaryfile dump_summery_file_name] [-d 2] [--exactintloop] [-v] [-p DIR] [-w DIR] [-l] <seq_file>\n\n");
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
	print_examples_developer_options();
	exit(-1);
}

static void printRunConfiguration(string seq) {

        if(!SILENT) printf("\nRun Configuration:\n");

	if (PF_COUNT_MODE == true) {
                if(!SILENT) printf("+ running with --pfcount option\n");
        }
        if (BPP_ENABLED == true) {
                if(!SILENT) printf("+ running with --bpp option\n");
        }
        if(dangles==0 && !PF_COUNT_MODE){//if (CALC_PF_DO && !CALC_PF_DS && !CALC_PF_D2 && !PF_COUNT_MODE) {
                if(!SILENT) printf("+ running in dangle d0 mode\n");
        }
        if(dangles==-1 && CALC_PF_DS){//if (!CALC_PF_DO && CALC_PF_DS && !CALC_PF_D2) {
                if(!SILENT) printf("+ running in dangle dS mode\n");
        }
        if(dangles==2 && !PF_COUNT_MODE){//if (!CALC_PF_DO && !CALC_PF_DS && CALC_PF_D2) {
                if(!SILENT) printf("- running in dangle d2 mode\n");
        }
        if (PARAM_DIR == true) {
                if(!SILENT) printf("+ running with customized param dir: %s\n",paramDir.c_str());
        }
        if (RND_SAMPLE == true) {
                if(!SILENT) printf("- running to calculate %d samples\n", num_rnd);
        }
        if (contactDistance != -1) {
                if(!SILENT) printf("- maximum contact distance: %d\n", contactDistance);
        }
        if(!SILENT) printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
        if(!SILENT) printf("- input sequence file: %s\n", seqfile.c_str());
        if(!SILENT) printf("- sequence length: %d\n", (int)seq.length());
        //if(!SILENT) printf("- output file: %s\n", outputFile.c_str());
        if(RND_SAMPLE) if(!SILENT) printf("- samples output file: %s\n", sampleOutFile.c_str());
        if(BPP_ENABLED) if(!SILENT) printf("- bpp output file: %s\n", bppOutFile.c_str());
        if(PF_PRINT_ARRAYS_ENABLED) if(!SILENT) printf("+ partition function array print output file: %s\n", pfArraysOutFile.c_str());
        if(print_energy_decompose==1) if(!SILENT) printf("+ energy decompose output file: %s\n", energyDecomposeOutFile.c_str());

	printf("\n");
}

static void validate_options(string seq){
        if(!SILENT) printf("\nValidating Options:\n");
	if(PF_COUNT_MODE){
		if(RND_SAMPLE){
			if(!SILENT) printf("Ignoring --sample Option with --pfcount option and Program will continue with --pfcount option only.\n\n");
			RND_SAMPLE = false;
		}
	}
	if(BPP_ENABLED){
		//do nothing
	}
	else if(CALC_PART_FUNC && !RND_SAMPLE){//partition function
		if(print_energy_decompose==1){
			if(!SILENT) printf("Ignoring the option -e or --energy, as it will be valid with --sample option.\n\n");	
		}
		if(ST_D2_ENABLE_SCATTER_PLOT && !ST_D2_ENABLE_BPP_PROBABILITY){
			if(!SILENT) printf("Ignoring the option --groupbyfreq, as it will be valid with --sample option.\n\n");
		}
		if(ST_D2_ENABLE_UNIFORM_SAMPLE){
                        if(!SILENT) printf("Ignoring the option --sampleenergy, as it will be valid with --sample option.\n\n");
                }
		if(ST_D2_ENABLE_CHECK_FRACTION){
                        if(!SILENT) printf("Ignoring the option --checkfraction, as it will be valid with --sample option.\n\n");
                }
		if(ST_D2_ENABLE_BPP_PROBABILITY){
                        if(!SILENT) printf("Ignoring the option --estimatebpp, as it will be valid with --sample option.\n\n");
                }
		//if(ST_D2_ENABLE_COUNTS_PARALLELIZATION){
                  //      printf("Ignoring the option --counts-parallel, as it will be valid with --sample option.\n\n");
                //}
		if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
                        if(!SILENT) printf("Ignoring the option --parallelsample, as it will be valid with --sample option.\n\n");
                }
	}
	else if(!CALC_PART_FUNC && RND_SAMPLE){//sample
		//nothing to check
        }
	else if(CALC_PART_FUNC && RND_SAMPLE){//both partition function and sample
        	//printf("Program proceeding with sampling as both partition function calculation and sampling calculation option are used.\n\n");
        	CALC_PART_FUNC = false;
        }
	else if(!CALC_PART_FUNC && !RND_SAMPLE){//neither partition function nor sample
        	printf("Program exiting as neither partition function calculation nor sampling calculation option is used.\n\n");
		help();
        	exit(-1);	
        }

	printf("\n");
}
/*
static bool is_int(char const* p)
{
    return strcmp(itoa(atoi(p)), p, 10) == 0;
}*/
static bool isNumeric( const char* pszInput)
{
	int nNumberBase = 10;
	string base = "0123456789ABCDEF";
	string input = pszInput;
 
	return (input.find_first_not_of(base.substr(0, nNumberBase)) == string::npos);
}

static void parse_options(int argc, char** argv) {
	int i;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-') {
			if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
				help(); 
			} else if(strcmp(argv[i], "--detailedhelp") == 0 ) {
				detailed_help();
			} else if (strcmp(argv[i], "--paramdir") == 0 || strcmp(argv[i], "-p") == 0) {
				if(i < argc) {
					paramDir = argv[++i];
					PARAM_DIR = true;
				}
				else {
					help();
				}
			} else if(strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0) {
				if(i < argc){
					g_nthreads = atoi(argv[++i]);
				}
				else help();
			} else if(strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) {
				if(i < argc)
					outputPrefix = argv[++i];
				else
					help();
			} else if (strcmp(argv[i], "--workdir") == 0 || strcmp(argv[i], "-w") == 0) {
				if(i+1 < argc)//if(i < argc)
					outputDir = argv[++i];
				else
					help();
			} else if(strcmp(argv[i], "--bpp") == 0) {
				BPP_ENABLED = true;
				CALC_PART_FUNC = false;
			} else if (strcmp(argv[i],"--partition") == 0) {
				CALC_PART_FUNC = true;
			} else if (strcmp(argv[i],"--printarrays") == 0) {
				PF_PRINT_ARRAYS_ENABLED = true;
			} else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
				g_verbose = 1;
			}
			else if (strcmp(argv[i], "--energydetail") == 0 || strcmp(argv[i], "-e") == 0) {
				print_energy_decompose = 1;
			}
			else if (strcmp(argv[i], "--dangle") == 0 || strcmp(argv[i], "-d") == 0) {
				std::string cmd = argv[i];
				if(i+1 < argc) {//if(i < argc) {
					dangles = atoi(argv[++i]);
					if (!(dangles == 0 || dangles == 2)) {
						dangles = 2;
						printf("Ignoring %s option as it accepts either 0 or 2, proceeding with taking value for dangle as 2\n", cmd.c_str());
					}
					if(dangles==0) CALC_PF_DO = true;
					else if(dangles==2) CALC_PF_D2 = true;
				} else
					help();
			}
			else if (strcmp(argv[i],"-dS") == 0) {
				CALC_PF_DS = true;  
				dangles=-1;
			} else if (strcmp(argv[i],"-d0") == 0) {
				CALC_PF_DO = true;
				dangles=0;  
			} else if (strcmp(argv[i],"-d2") == 0) {
				//help();
				CALC_PF_D2 = true;
				dangles=2;
			} else if(strcmp(argv[i],"--exactintloop") == 0){ 
				PF_D2_UP_APPROX_ENABLED = false;
			} else if(strcmp(argv[i],"--groupbyfreq") == 0){
				ST_D2_ENABLE_SCATTER_PLOT = true;
			} else if(strcmp(argv[i],"--sampleenergy") == 0){
				ST_D2_ENABLE_UNIFORM_SAMPLE = true;
				if(i < argc){ST_D2_UNIFORM_SAMPLE_ENERGY = atof(argv[++i]);}
				else help();
			} else if(strcmp(argv[i],"--checkfraction") == 0){
				ST_D2_ENABLE_CHECK_FRACTION = true;
			} else if(strcmp(argv[i],"--estimatebpp") == 0){ 
				ST_D2_ENABLE_BPP_PROBABILITY = true;
				ST_D2_ENABLE_SCATTER_PLOT = true;
			} else if(strcmp(argv[i],"--counts-parallel") == 0){
				ST_D2_ENABLE_COUNTS_PARALLELIZATION = true;
			} else if(strcmp(argv[i],"--parallelsample") == 0){ 
				ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION = true;
			} else if (strcmp(argv[i],"--pfcount") == 0) {
				CALC_PART_FUNC = true;
				PF_COUNT_MODE = true;
			} else if (strcmp(argv[i],"--sample") == 0 || strcmp(argv[i], "-s") == 0) {
				RND_SAMPLE = true;
				if(i+1 < argc){//if(i < argc){
					//if(is_int(argv[i+1])){		
					if(isNumeric(argv[i+1])){		
						num_rnd = atoi(argv[++i]);
					}
					else {
						printf("Error: %s is not a valid integer.\n\n",argv[i+1]);
						help();
					}
				}
				else
					help();
				if(i+1 < argc){//if(i < argc){//--separatectfiles [--ctfilesdir dump_dir_name] [--summaryfile dump_summery_file_name]
					//printf("i=%d, argc=%d, argv[i]=%s, \n",i,argc,argv[i]);
					if (strcmp(argv[i+1],"--separatectfiles") == 0){
						i=i+1;
						DUMP_CT_FILE = true;
						if (i+1 < argc && strcmp(argv[i+1],"--ctfilesdir") == 0){//if (i < argc && strcmp(argv[i+1],"--ctfilesdir") == 0){
							i=i+1;
							if(i+1 < argc)//if(i < argc)
								ctFileDumpDir = argv[++i];
							else
								help();
						}
						if (i+1 < argc && strcmp(argv[i+1],"--summaryfile") == 0){//if (i < argc && strcmp(argv[i+1],"--summaryfile") == 0){
							i=i+1;
							if(i+1 < argc)//if(i < argc)
								stochastic_summery_file_name = argv[++i];
							else
								help();
						}
					}
				}
			} else if (strcmp(argv[i],"--limitcd") == 0 || strcmp(argv[i], "-l") == 0) {
				if(i+1 < argc) {//if(i < argc) {
					LIMIT_DISTANCE = true;
					contactDistance = atoi(argv[++i]);
				}
				else
					help();
			}
			else{
				printf("Error: Option %s is Undefined option\n", argv[i]);
				help();
			}
		} else {
			seqfile = argv[i];
		}
	}

	if(seqfile.compare("")==0 || seqfile.empty()) {
		printf("Error: Missing input file.\n");
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
		bppOutFile += outputDir;
		bppOutFile += "/";
		sampleOutFile += outputDir;
		sampleOutFile += "/";
		energyDecomposeOutFile += outputDir;
		energyDecomposeOutFile += "/";
		estimateBppOutputFile += outputDir;
		estimateBppOutputFile += "/";
		scatterPlotOutputFile += outputDir;
		scatterPlotOutputFile += "/";
		pfArraysOutFile += outputDir;
		pfArraysOutFile += "/";

	}
	// ... and append the .ct
	outputFile += outputPrefix;
	outputFile += ".ct";

	bppOutFile += outputPrefix;	
	bppOutFile += "_bpp.txt";

	sampleOutFile += outputPrefix;	
	sampleOutFile += ".samples";

	energyDecomposeOutFile += outputPrefix;	
	energyDecomposeOutFile += ".energy";

	estimateBppOutputFile += outputPrefix;	
	estimateBppOutputFile += ".sbpp";

	scatterPlotOutputFile += outputPrefix;	
	scatterPlotOutputFile += ".frequency";

	pfArraysOutFile += outputPrefix;
	pfArraysOutFile += ".pfarrays";

}
/*
   double get_seconds() {
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
   }*/

int boltzmann_main(int argc, char** argv) {
	std::string seq;
	double t1;
	parse_options(argc, argv);
	validate_options(seqfile);
	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}
	parse_mfe_options(argc, argv);
	init_fold(seq.c_str());
	g_LIMIT_DISTANCE = LIMIT_DISTANCE;
	g_contactDistance = contactDistance;

	readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, 0, 0, 0);
	printRunConfiguration(seqfile);
	
	if (LIMIT_DISTANCE) {
		if (strlen(seq.c_str()) < contactDistance) 
			printf("\nContact distance limit is higher than the sequence length. Continuing without restraining contact distance.\n");
		else printf("\nLimiting contact distance to %d\n",contactDistance);
	}
	
	if (CALC_PART_FUNC == true && CALC_PF_DS == true) {
		int pf_count_mode = 0;
		if(PF_COUNT_MODE) pf_count_mode=1;
		int no_dangle_mode = 0;
		if(CALC_PF_DO) no_dangle_mode=1;
		//printf("\nComputing partition function in -dS mode ..., pf_count_mode=%d, no_dangle_mode=%d\n", pf_count_mode, no_dangle_mode);
		printf("\nComputing partition function...\n");
		t1 = get_seconds();
		calculate_partition(seq.length(),pf_count_mode,no_dangle_mode);
		t1 = get_seconds() - t1;
		printf("partition function computation running time: %9.6f seconds\n", t1);
		//calculate_partition(seq.length(),0,0);
		free_partition();
	}
	else if (CALC_PART_FUNC == true && CALC_PF_D2 == true) {
		int pf_count_mode = 0;
		if(PF_COUNT_MODE) pf_count_mode=1;
		int no_dangle_mode = 0;
		if(CALC_PF_DO) no_dangle_mode=1;
		//printf("\nComputing partition function in -d2 mode ..., pf_count_mode=%d, no_dangle_mode=%d, PF_D2_UP_APPROX_ENABLED=%d\n", pf_count_mode, no_dangle_mode,PF_D2_UP_APPROX_ENABLED);
		printf("\nComputing partition function...\n");
		t1 = get_seconds();
		PartitionFunctionD2 pf_d2;
		pf_d2.calculate_partition(seq.length(),pf_count_mode,no_dangle_mode,PF_D2_UP_APPROX_ENABLED);
		t1 = get_seconds() - t1;
		printf("partition function computation running time: %9.6f seconds\n", t1);
		//calculate_partition(seq.length(),0,0);
		if(PF_PRINT_ARRAYS_ENABLED) pf_d2.printAllMatrixesToFile(pfArraysOutFile);
		pf_d2.free_partition();
	}  
	else if (CALC_PART_FUNC == true && CALC_PF_DO == true) {
		//printf("\nCalculating partition function in -d0 mode ...\n");
		printf("\nCalculating partition function...\n");
		/*
		//Below method is not correct method for d0 mdoe partition function computation as discussed by Shel
		double ** Q,  **QM, **QB, **P;
		Q = mallocTwoD(seq.length() + 1, seq.length() + 1);
		QM = mallocTwoD(seq.length() + 1, seq.length() + 1);
		QB = mallocTwoD(seq.length() + 1, seq.length() + 1);
		P = mallocTwoD(seq.length() + 1, seq.length() + 1);


		fill_partition_fn_arrays(seq.length(), Q, QB, QM);
		freeTwoD(Q, seq.length() + 1, seq.length() + 1);
		freeTwoD(QM, seq.length() + 1, seq.length() + 1);
		freeTwoD(QB, seq.length() + 1, seq.length() + 1);
		freeTwoD(P, seq.length() + 1, seq.length() + 1);
		 */
		t1 = get_seconds();
		calculate_partition(seq.length(),0,1);
		t1 = get_seconds() - t1;
		printf("partition function computation running time: %9.6f seconds\n", t1);
		//calculate_partition(seq.length(),0,0);
		free_partition();

	}
	else if (CALC_PART_FUNC == true) {
		printf("\nComputing partition function...\n");
		int pf_count_mode = 0;
		if(PF_COUNT_MODE) pf_count_mode=1;
		t1 = get_seconds();
		calculate_partition(seq.length(),pf_count_mode, 0);
		t1 = get_seconds() - t1;
		printf("partition function computation running time: %9.6f seconds\n", t1);
		free_partition();
	} else if (RND_SAMPLE == true) {
		int pf_count_mode = 0;
		if(PF_COUNT_MODE) pf_count_mode=1;
		int no_dangle_mode = 0;
		if(CALC_PF_DO) no_dangle_mode=1;
		  
		if(CALC_PF_D2 == true){
			//printf("\nComputing stochastic traceback in -d2 mode ..., pf_count_mode=%d, no_dangle_mode=%d, PF_D2_UP_APPROX_ENABLED=%d\n", pf_count_mode, no_dangle_mode,PF_D2_UP_APPROX_ENABLED);
			printf("\nComputing stochastic traceback...\n");
			StochasticTracebackD2 st_d2;
			t1 = get_seconds();
                        st_d2.initialize(seq.length(), pf_count_mode, no_dangle_mode, print_energy_decompose, PF_D2_UP_APPROX_ENABLED,ST_D2_ENABLE_CHECK_FRACTION, energyDecomposeOutFile);
                        t1 = get_seconds() - t1;
                        printf("D2 Traceback initialization (partition function computation) running time: %9.6f seconds\n", t1);
			t1 = get_seconds();
			if(DUMP_CT_FILE==false){
				if(ST_D2_ENABLE_COUNTS_PARALLELIZATION && g_nthreads!=1)
					st_d2.batch_sample_parallel(num_rnd,ST_D2_ENABLE_SCATTER_PLOT,ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION,ST_D2_ENABLE_BPP_PROBABILITY, sampleOutFile, estimateBppOutputFile, scatterPlotOutputFile);
				else st_d2.batch_sample(num_rnd,ST_D2_ENABLE_SCATTER_PLOT,ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION,ST_D2_ENABLE_UNIFORM_SAMPLE,ST_D2_UNIFORM_SAMPLE_ENERGY,ST_D2_ENABLE_BPP_PROBABILITY, sampleOutFile, estimateBppOutputFile, scatterPlotOutputFile);
			}
                        else  st_d2.batch_sample_and_dump(num_rnd, ctFileDumpDir, stochastic_summery_file_name, seq, seqfile);
			t1 = get_seconds() - t1;
                        printf("D2 Traceback computation running time: %9.6f seconds\n", t1);
			if(PF_PRINT_ARRAYS_ENABLED) st_d2.printPfMatrixesToFile(pfArraysOutFile);
			st_d2.free_traceback();
		}
		else{//if(CALC_PF_DS == true){//TODO here it is making dS by default
			//printf("\nComputing stochastic traceback in -dS mode ..., pf_count_mode=%d, no_dangle_mode=%d\n", pf_count_mode, no_dangle_mode);
			printf("\nComputing stochastic traceback...\n");
			double U = calculate_partition(seq.length(),pf_count_mode,0);
			t1 = get_seconds();
			if(DUMP_CT_FILE==false) batch_sample(num_rnd, seq.length(), U);
			else batch_sample_and_dump(num_rnd, seq.length(), U, ctFileDumpDir, stochastic_summery_file_name, seq, seqfile); 
			t1 = get_seconds() - t1;
	                printf("Traceback computation running time: %9.6f seconds\n", t1);
			free_partition();
		}
	} else if(BPP_ENABLED) {
		printf("\n");
		printf("Calculating partition function\n");
		double ** Q,  **QM, **QB, **P;
		Q = mallocTwoD(seq.length() + 1, seq.length() + 1);
		QM = mallocTwoD(seq.length() + 1, seq.length() + 1);
		QB = mallocTwoD(seq.length() + 1, seq.length() + 1);
		P = mallocTwoD(seq.length() + 1, seq.length() + 1);


		fill_partition_fn_arrays(seq.length(), Q, QB, QM);
		fillBasePairProbabilities(seq.length(), Q, QB, QM, P);
		//printBasePairProbabilities(seq.length(), structure, P, bppOutFile.c_str());
		printBasePairProbabilitiesDetail(seq.length(), structure, P, bppOutFile.c_str());
		printf("Saved BPP output in %s\n",bppOutFile.c_str());

		freeTwoD(Q, seq.length() + 1, seq.length() + 1);
		freeTwoD(QM, seq.length() + 1, seq.length() + 1);
		freeTwoD(QB, seq.length() + 1, seq.length() + 1);
		freeTwoD(P, seq.length() + 1, seq.length() + 1);
	} else {
		printf("No valid option specified !\n\n");
		help();
	}

	free_fold(seq.length());
	printf("\n");
	return EXIT_SUCCESS;
}
