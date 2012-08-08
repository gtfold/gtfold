#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
//#include <sys/time.h>
//#include <time.h>

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

static bool CALC_PART_FUNC = false;
static bool PF_COUNT_MODE = false;
static bool BPP_ENABLED = false;
static bool PARAM_DIR = false;
static bool RND_SAMPLE = false;
static bool DUMP_CT_FILE = false;
static bool CALC_PF_DO = false;
static bool CALC_PF_DS = false;
static bool CALC_PF_D2 = true;//making D2 mode as default option
static bool PF_D2_UP_APPROX_ENABLED = true;//making short internal loop code to run as default
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
static string ctFileDumpDir = "";
static string stochastic_summery_file_name = "stochaSampleSummary.txt";

static int num_rnd = 0;
static int ss_verbose_global = 0;
static int print_energy_decompose = 0;
static int dangles=2;//making dangle default value as 2

static bool LIMIT_DISTANCE = false;
static int contactDistance = -1;


static void help() {
	printf("Usage: gtboltzmann [OPTION]... FILE\n\n");

	printf("   FILE is an RNA sequence file containing only the sequence or in FASTA format.\n\n");

	printf("OPTIONS\n");

	printf("   -t|--threads INT	Limit number of threads used to INT.\n");
	printf("   -v, --verbose	Run in verbose mode (includes partition function table printing.)\n");
	printf("   -d, --dangle INT	Restricts treatment of dangling energies (INT=0,2),\n");
	printf("   --partition		Calculate the partition function (default is using d2 dangling mode).\n");
	printf("   --exact-internal-loop	Do the exact internal loop calculation while calculating partition function and traceback without any short internal loop approximation)\n");
	printf("   -dS			Calculate the partition function using sfold reccurences and use them in traceback.\n");
	printf("   --sample   INT	Sample number of structures equal to INT.\n");
	printf("   --scatterPlot	While sampling structures, Collect frequency of all structures and calculate estimate probability and boltzmann probability for scatter plot.\n");
	printf("   --uniformSample energy1	While sampling structures, Samples with Energy energy1 will only be sampled.\n");
	printf("   --check-fraction	While sampling structures, enable test of check fraction.\n");
	printf("   --estimate-bpp	While sampling structures, Calculate base pair probabilities.\n");
	printf("   --counts-parallel	While sampling structures, parallelize INT sample counts among available threads (this is also a default behaviour of sampling).\n");
	printf("   --one-sample-parallel	While sampling structures, parallelize the processing of one sample (useful when sampling large sequence with number of samples being less than available threads).\n");
	printf("   --sample   INT  --dump [--dump_dir dump_dir_path] [--dump_summary dump_summery_file_name]	Sample number of structures equal to INT and dump each structure to a ct file in dump_dir_path directory (if no value provided then use current directory value for this purpose) and also create a summary file with name stochastic_summery_file_name in dump_dir_path directory (if no value provided, use stochaSampleSummary.txt value for this purpose).\n");
	printf("   --pfcount		Calculate the structure count using partition function and zero energy value.\n");
	printf("   --bpp		Calculate base pair probabilities.\n");
	printf("   -l|--limitCD  INT	Set a maximum base pair contact distance to INT. If no\n");
	printf(" 		      	limit is given, base pairs can be over any distance.\n");
	printf("   -o, --output NAME    Write output files with prefix given in NAME\n");
	printf("   -p  --paramdir DIR   Path to directory from which parameters are to be read\n");
	printf("   -h, --help           Output help (this message) and exit.\n");
	printf("   --detailed-help      Output help (this message) with detailed options and examples, and exit.\n");
	printf("   -e, --energy         prints energy decomposition for sampled structures to file with extention '.energy' (should be used with '-t 1' option, as otherwise all threads in parallel, will write to file and output will be intermixed from all threads).\n");
	printf("   -w, --workdir DIR    Path of directory where output files will be written.\n");
	exit(-1);
}
static void detailed_help(){
        printf("\n\nEXAMPLES:\n\n");
        printf("1. Calculate Partition function:\n");
        printf("gtboltzmann --partition [[-d 0|2]|[-dS]] [-t 1|2|...|N] [-o outputPrefix] [--exact-internal-loop] [-v] [-p DIR] [-p DIR] [-l] <seq_file>\n\n");
        printf("2. Sample structures stochastically:\n");
        printf("gtboltzmann --sample INT [[-d 0|2]|[-dS]] [-t 1|2|...|N] [-o outputPrefix] [--exact-internal-loop] [-v] [--scatterPlot] [--uniformSample DOUBLE] [--estimate-bpp] [--one-sample-parallel] [-p DIR] [-l] <seq_file>\n\n");
        printf("gtboltzmann --sample INT [[-d 0|2]|[-dS]] -t 1 [-o outputPrefix] [--exact-internal-loop] [-v] [--scatterPlot] [--uniformSample DOUBLE] [-e] [--check-fraction] [--estimate-bpp] [--one-sample-parallel] [-p DIR] [-l] <seq_file>\n\n");
        printf("gtboltzmann --sample INT --dump [--dump_dir dump_dir_path] [--dump_summary dump_summery_file_name] [-d 2] [--exact-internal-loop] [-v] [-p DIR] [-l] <seq_file>\n\n");
        printf("\n\n");
        help();
}

static void validate_options(){
	if(PF_COUNT_MODE){
		if(RND_SAMPLE){
			printf("Ignoring --sample Option with --pfcount option and Program will continue with --pfcount option only.\n\n");
			RND_SAMPLE = false;
		}
	}
	if(CALC_PART_FUNC && !RND_SAMPLE){//partition function
		if(print_energy_decompose==1){
			printf("Ignoring the option -e or --energy, as it will be valid with --sample option.\n\n");	
		}
		if(ST_D2_ENABLE_SCATTER_PLOT){
			printf("Ignoring the option --scatterPlot, as it will be valid with --sample option.\n\n");
		}
		if(ST_D2_ENABLE_UNIFORM_SAMPLE){
                        printf("Ignoring the option --uniformSample, as it will be valid with --sample option.\n\n");
                }
		if(ST_D2_ENABLE_CHECK_FRACTION){
                        printf("Ignoring the option --check-fraction, as it will be valid with --sample option.\n\n");
                }
		if(ST_D2_ENABLE_BPP_PROBABILITY){
                        printf("Ignoring the option --estimate-bpp, as it will be valid with --sample option.\n\n");
                }
		//if(ST_D2_ENABLE_COUNTS_PARALLELIZATION){
                  //      printf("Ignoring the option --counts-parallel, as it will be valid with --sample option.\n\n");
                //}
		if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
                        printf("Ignoring the option --one-sample-parallel, as it will be valid with --sample option.\n\n");
                }
	}
	else if(!CALC_PART_FUNC && RND_SAMPLE){//sample
		//nothing to check
        }
	else if(CALC_PART_FUNC && RND_SAMPLE){//both partition function and sample
        	printf("Program proceeding with sampling as both partition function calculation and sampling calculation option are used.\n\n");
        	CALC_PART_FUNC = false;
        }
	else if(!CALC_PART_FUNC && !RND_SAMPLE){//neither partition function nor sample
        	printf("Program exiting as neither partition function calculation nor sampling calculation option is used.\n\n");
		help();
        	exit(-1);	
        }
}

static void parse_options(int argc, char** argv) {
	int i;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-') {
			if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
				help(); 
			}
			if(strcmp(argv[i], "--detailed-help") == 0 ) {
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
				if(i < argc)
					outputDir = argv[++i];
				else
					help();
			} else if(strcmp(argv[i], "--bpp") == 0) {
				BPP_ENABLED = true;
			}
			else if (strcmp(argv[i],"--partition") == 0) {
				CALC_PART_FUNC = true;
			} else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
				g_verbose = 1;
      			}
			else if (strcmp(argv[i], "--energy") == 0 || strcmp(argv[i], "-e") == 0) {
                                print_energy_decompose = 1;
                        }
			else if (strcmp(argv[i], "--dangle") == 0 || strcmp(argv[i], "-d") == 0) {
        			std::string cmd = argv[i];
        			if(i < argc) {
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
			} else if (strcmp(argv[i],"-d0") == 0) {
				CALC_PF_DO = true;  
			} else if (strcmp(argv[i],"-d2") == 0) {
				//help();
				CALC_PF_D2 = true;
			} else if(strcmp(argv[i],"--exact-internal-loop") == 0){ 
				PF_D2_UP_APPROX_ENABLED = false;
			} else if(strcmp(argv[i],"--scatterPlot") == 0){
				ST_D2_ENABLE_SCATTER_PLOT = true;
			} else if(strcmp(argv[i],"--uniformSample") == 0){
				ST_D2_ENABLE_UNIFORM_SAMPLE = true;
				if(i < argc){ST_D2_UNIFORM_SAMPLE_ENERGY = atof(argv[++i]);}
				else help();
			} else if(strcmp(argv[i],"--check-fraction") == 0){
				ST_D2_ENABLE_CHECK_FRACTION = true;
			} else if(strcmp(argv[i],"--estimate-bpp") == 0){ 
				ST_D2_ENABLE_BPP_PROBABILITY = true;
				ST_D2_ENABLE_SCATTER_PLOT = true;
			} else if(strcmp(argv[i],"--counts-parallel") == 0){
				ST_D2_ENABLE_COUNTS_PARALLELIZATION = true;
			} else if(strcmp(argv[i],"--one-sample-parallel") == 0){ 
				ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION = true;
			} else if (strcmp(argv[i],"--pfcount") == 0) {
				CALC_PART_FUNC = true;
				PF_COUNT_MODE = true;
			} else if (strcmp(argv[i],"--sample") == 0) {
				RND_SAMPLE = true;
				if(i < argc)
					num_rnd = atoi(argv[++i]);
				else
					help();
				if(i < argc){//--dump [--dump_dir dump_dir_name] [--dump_summary dump_summery_file_name]
					if (strcmp(argv[i+1],"--dump") == 0){
						i=i+1;
						DUMP_CT_FILE = true;
						if (i < argc && strcmp(argv[i+1],"--dump_dir") == 0){
							i=i+1;
							if(i < argc)
								ctFileDumpDir = argv[++i];
							else
								help();
						}
						if (i < argc && strcmp(argv[i+1],"--dump_summary") == 0){
							i=i+1;
							if(i < argc)
								stochastic_summery_file_name = argv[++i];
							else
								help();
						}
					}
				}
			} else if (strcmp(argv[i],"--limitCD") == 0 || strcmp(argv[i], "-l") == 0) {
				if(i < argc) {
					LIMIT_DISTANCE = true;
					contactDistance = atoi(argv[++i]);
				}
				else
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
	validate_options();
	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}
	parse_mfe_options(argc, argv);
	init_fold(seq.c_str());
	g_LIMIT_DISTANCE = LIMIT_DISTANCE;
	g_contactDistance = contactDistance;

	readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, 0, 0, 0);
	
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
		printf("\nComputing partition function in -dS mode ..., pf_count_mode=%d, no_dangle_mode=%d\n", pf_count_mode, no_dangle_mode);
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
		printf("\nComputing partition function in -d2 mode ..., pf_count_mode=%d, no_dangle_mode=%d, PF_D2_UP_APPROX_ENABLED=%d\n", pf_count_mode, no_dangle_mode,PF_D2_UP_APPROX_ENABLED);
		t1 = get_seconds();
		PartitionFunctionD2 pf_d2;
		pf_d2.calculate_partition(seq.length(),pf_count_mode,no_dangle_mode,PF_D2_UP_APPROX_ENABLED);
		t1 = get_seconds() - t1;
		printf("partition function computation running time: %9.6f seconds\n", t1);
		//calculate_partition(seq.length(),0,0);
		pf_d2.free_partition();
	}  
	else if (CALC_PART_FUNC == true && CALC_PF_DO == true) {
		printf("\nCalculating partition function in -d0 mode ...\n");
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
		printf("\nComputing partition function...\n");
		int pf_count_mode = 0;
		if(PF_COUNT_MODE) pf_count_mode=1;
		int no_dangle_mode = 0;
		if(CALC_PF_DO) no_dangle_mode=1;
		  
		if(CALC_PF_D2 == true){
			printf("\nComputing stochastic traceback in -d2 mode ..., pf_count_mode=%d, no_dangle_mode=%d, PF_D2_UP_APPROX_ENABLED=%d\n", pf_count_mode, no_dangle_mode,PF_D2_UP_APPROX_ENABLED);
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
			st_d2.free_traceback();
		}
		else{//if(CALC_PF_DS == true){//TODO here it is making dS by default
			printf("\nComputing stochastic traceback in -dS mode ..., pf_count_mode=%d, no_dangle_mode=%d\n", pf_count_mode, no_dangle_mode);
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
		printBasePairProbabilities(seq.length(), structure, P, bppOutFile.c_str());
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
	return EXIT_SUCCESS;
}
