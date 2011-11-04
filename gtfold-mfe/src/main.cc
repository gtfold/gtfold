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
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main.h"
#include "utils.h"
#include "loader.h"
#include "options.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "algorithms-partition.h"
#include "partition-func.h"
#include "constraints.h"
#include "traceback.h"
#include "subopt_traceback.h"
#include "shapereader.h"
#include "stochastic-sampling.h"

//#define DEBUG 1

using namespace std;

double get_seconds() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

void init_fold(string seq) {
	int len = seq.length();

	init_global_params(len);

	if (!encodeSequence(seq)) {
		free_fold(seq.length());
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

/**
 * Read the sequence out of the given filename and store it in seq
 *
 * @param filename A c string with the file to open
 * @param seq A C++ string object (passed by reference) to save to
 * @return SUCCESS or FAILURE
 */
int read_sequence_file(const char* filename, std::string& seq) {
	seq = "";

	ifstream fs;
	fs.open(filename, ios::in);
	if (!fs.good()) return FAILURE;

	string line;
	while(fs.good()) {
		getline(fs, line);
		// exclude lines starting with FASTA comment characters
		if(line[0] != ';' && line[0] != '>' && line.length() > 0)
			seq += line;
	}

	fs.close();

    size_t loc;
    while((loc = seq.find(" ")) != string::npos)
        seq.erase(loc, 1);

	return SUCCESS;
}

/**
 * Encode the given string of letters into an RNA sequence.
 *
 * This handles IUPAC codes and places the results into the array RNA.
 *
 * @param seq The string to encode
 * @return A boolean representing success or failure of the encoding.  Failure
 *         can occur when seeing unknown IUPAC codes.
 */
bool encodeSequence(string seq) {
	unsigned int unspecifiedBaseCount=0;
	unsigned int unspecifiedBases[seq.length()];

	for(unsigned int i=1; i<=seq.length(); i++) {
		RNA[i] = encode(seq[i-1]);

		// die on non-IUPAC codes
		if (RNA[i]=='X') {
			fprintf(stderr, "ERROR: Non-IUPAC nucleotide code at position: %d (%c)\n", i, seq[i-1]);
			fprintf(stderr, "See http://www.bioinformatics.org/sms/iupac.html for valid codes.\n");
			return false;
		}

		// add non-canonical IUPAC codes to the warning list
		if(!isWatsonCrickBase(seq[i-1]))
			unspecifiedBases[unspecifiedBaseCount++]=i;
	}

	// just print a warning for non-canonical IUPAC codes
	if(unspecifiedBaseCount > 0) {
		printf("\nIncompletely-specified IUPAC codes have been detected at position%s: ", unspecifiedBaseCount == 1 ? "" : "s");

		// put [0] first for nice comma separation
		printf("%d (%c)", unspecifiedBases[0], seq.at(unspecifiedBases[0]-1));
		for(unsigned int i=1; i<unspecifiedBaseCount; i++)
			printf(", %d (%c)", unspecifiedBases[i], seq[unspecifiedBases[i]-1]);

		printf("\nPlease replace with fully-specified IUPAC codes (A,C,G,U,T) and retry.\n");
		return false;
	}

	return true;
}

void print_header() {
	printf("GTfold: A Scalable Multicore Code for RNA Secondary Structure Prediction\n");
	printf("(c) 2007-2011  D.A. Bader, C.E. Heitsch, S.C. Harvey\n");
	printf("Georgia Institute of Technology\n\n");
}

void save_subopt_file(string outputFile, ss_map_t& ss_data, 
		const string& seq, int energy)
{
	ofstream outfile;
	outfile.open(outputFile.c_str());
	char buff[1024];

	sprintf(buff,"%s\t%0.2f", seq.c_str(), energy/100.0);
	outfile << buff << std::endl;
	for (ss_map_t::iterator it = ss_data.begin(); it!= ss_data.end(); ++it) 
	{
		sprintf(buff,"%s %0.2f", (it->first).c_str(), it->second/100.0);
		//outfile << it->first << '\t' << it->second/100.0 << std::endl;
		outfile << buff << std::endl;
	}

	outfile.close();
}

/**
 * Save the output to a ct file
 *
 * @param outputFile The file to save to
 * @param energy The MFE energy (multiplied by 100)
 */
void save_ct_file(string outputFile, string seq, int energy) {

	ofstream outfile;
	outfile.open(outputFile.c_str());

	outfile << seq.length() << "\t  dG = " << energy/100.0 << endl;
	//outfile << seq.length() << "\tdG = " << energy/100.0 << "\t" << seqfile << endl;

	unsigned int i = 1;
	for(i=1; i <= seq.length(); i++)
		outfile << i << "\t" << seq[i-1] << "\t" << i-1 << "\t" << (i+1)%(seq.length()+1) << "\t" << structure[i] << "\t" << i << endl;

	outfile.close();
}

int main(int argc, char** argv) {
	std::string seq;
	int energy;
	double t1;
	
	print_header();

	parse_options(argc, argv);

	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}
	init_fold(seq);
	
	// Read in thermodynamic parameters. Always use Turner99 data (for now)
	if (SUBOPT_ENABLED) 
									readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, UNAMODE, 1, T_MISMATCH);
	else
									readThermodynamicParameters(paramDir.c_str(), PARAM_DIR, UNAMODE, RNAMODE, T_MISMATCH);

	printRunConfiguration(seq);

  if (CALC_PART_FUNC == true)
  {
    printf("\nComputing partition function...\n");
    int pf_count_mode = 0;
    if(PF_COUNT_MODE) pf_count_mode=1;
    int no_dangle_mode=0;
    if(NO_DANGLE_MODE) no_dangle_mode=1;
	t1 = get_seconds();
	calculate_partition(seq.length(),pf_count_mode,no_dangle_mode);
	t1 = get_seconds() - t1;
	printf("partition function computation running time: %9.6f seconds\n", t1);

    free_partition();
    free_fold(seq.length());
    exit(0);
  }
  if (RND_SAMPLE == true)
  { 
	//below code is for comparison of scores of different sampled structures from different methods
	int ctFileWriteOn=0;
	int summaryWriteOn=0;
        string summaryfile = "";
        ofstream summaryoutfile;      
	if(summaryWriteOn){//if(ctFileWriteOn){
		summaryfile = "/home/users/msoni/Desktop/manoj_gatech/research/gtfold/git_code/gtfold/gtfold-mfe/src/stochaSampleSummery.txt";
		summaryoutfile.open(summaryfile.c_str());
	}

	//below code is for comparison of BP probabilities 
	int calcBpProb=0;
         double** bpProb;
         if(calcBpProb){
                bpProb = new double*[seq.length()+1];
                for(int bpIndI=1; bpIndI<=seq.length(); ++bpIndI){
                        bpProb[bpIndI]=new double[seq.length()+1];
                        for(int bpIndJ=1; bpIndJ<=seq.length(); ++bpIndJ){
                                bpProb[bpIndI][bpIndJ]=0.0;
                        }
                }
          }

  
	printf("\nComputing partition function...\n");
	  int pf_count_mode = 0;
	  if(PF_COUNT_MODE) pf_count_mode=1;
	  int no_dangle_mode=0;
    	  if(NO_DANGLE_MODE) no_dangle_mode=1;
	t1 = get_seconds();
	  calculate_partition(seq.length(),pf_count_mode,no_dangle_mode);
	t1 = get_seconds() - t1;
        printf("partition function computation running time: %9.6f seconds\n", t1);
	t1 = get_seconds();
	  //int* structure = new int[seq.length()+1];
	  srand(time(NULL));
	std::map<std::string,std::pair<int,double> >  uniq_structs;

	  if (num_rnd > 0 ) {
		  printf("\nSampling structures...\n");
		  for (int count = 1; count <= num_rnd; ++count) 
		  {
			  memset(structure, 0, (seq.length()+1)*sizeof(int));
			  double energy = rnd_structure(structure, seq.length());

			  std::string ensemble(seq.length()+1,'.');
			  for (int i = 1; i <= (int)seq.length(); ++ i) {
				  //     printf("%d %d\n",i,structure[i]);
				  if (structure[i] > 0 && ensemble[i] == '.')
				  {
					  ensemble[i] = '(';
					  ensemble[structure[i]] = ')';
				  }
			  }
			if(summaryWriteOn){
				char abspath[1000];
                        	getcwd(abspath, 1000);
                       		std::stringstream ss;
                        	ss<<abspath<<"/"<<seqfile<<"_"<<count<<".ct";
                        	if(ctFileWriteOn) save_ct_file(ss.str(), seq, energy);
                          	summaryoutfile<<ss.str()<<" "<<ensemble.substr(1)<<" "<<energy<< std::endl;
			}

        		std::map<std::string,std::pair<int,double> >::iterator iter ;
        		if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
        		{
        			std::pair<int,double>& pp = iter->second;
        			pp.first++;
        		}
        		else {
        			uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
        		}

                        if(calcBpProb){
                                for(int struInd=1; struInd<=seq.length(); ++struInd){
                                        if(structure[struInd]>0){bpProb[struInd][structure[struInd]]++;}
                                }
                        }

			//std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
		  }
	}
	if(summaryWriteOn){
		if(ctFileWriteOn) printf("All sampled structure were dumped into corresponding ct files.\n");
        	printf("Saved Score Summary output in %s\n",summaryfile.c_str());
		summaryoutfile.close();
	}    
	
	int pcount = 0;
	int maxCount = 0; std::string bestStruct;
	double bestE = INFINITY;

	std::map<std::string,std::pair<int,double> >::iterator iter ;
	for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
	{
		const std::string& ss = iter->first;
		const std::pair<int,double>& pp = iter->second;
		//  printf("%s\tp=%lf, e=%lf\n",ss.c_str(),(double)pp.first/(double)num_rnd,pp.second);
		pcount += pp.first;
		if (pp.first > maxCount)
		{
			maxCount = pp.first;
			bestStruct  = ss;
			bestE = pp.second;
		}
	}
	assert(num_rnd == pcount);
	printf("Most favourable structure is : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);

	t1 = get_seconds() - t1;
        printf("Stochastic traceback (samples=%d), running time: %9.6f seconds\n", num_rnd, t1);

	if(calcBpProb){
        	printf("\nComputing Base Pair Probability Comparison Data:\n");
        	for(int bpIndI=1; bpIndI<=seq.length(); ++bpIndI){
                	for(int bpIndJ=1; bpIndJ<=seq.length(); ++bpIndJ){
                        	bpProb[bpIndI][bpIndJ]/=num_rnd;
                	}
        	}

        	double ** Q,  **QM, **QB, **P;
	        Q = mallocTwoD(seq.length() + 1, seq.length() + 1);
	        QM = mallocTwoD(seq.length() + 1, seq.length() + 1);
	        QB = mallocTwoD(seq.length() + 1, seq.length() + 1);
	        P = mallocTwoD(seq.length() + 1, seq.length() + 1);
	        fill_partition_fn_arrays(seq.length(), Q, QB, QM);
	        fillBasePairProbabilities(seq.length(), Q, QB, QM, P);

        	char abspath[1000];
	        getcwd(abspath, 1000);
	        std::stringstream bppss;
		printf("no dangle mode is %d\n", NO_DANGLE_MODE);
	        if(NO_DANGLE_MODE) bppss<<abspath<<"/"<<seqfile<<"_"<<"bppComparisonNoDangle"<<".ct";
		else bppss<<abspath<<"/"<<seqfile<<"_"<<"bppComparison"<<".ct";
	        printBasePairProbabilitiesComparison(seq.length(), structure, P, bpProb, bppss.str().c_str());
	        printf("Saved BPP comparison output in %s\n",bppss.str().c_str());
	        freeTwoD(Q, seq.length() + 1, seq.length() + 1);
	        freeTwoD(QM, seq.length() + 1, seq.length() + 1);
	        freeTwoD(QB, seq.length() + 1, seq.length() + 1);
	        freeTwoD(P, seq.length() + 1, seq.length() + 1);
	        for(int bpIndI=1; bpIndI<=seq.length(); ++bpIndI)delete[] bpProb[bpIndI];
	        delete[] bpProb;
	}

	free_partition();
	free_fold(seq.length());
	//delete [] structure;
	//t1 = get_seconds() - t1;
	//printf("Stochastic traceback (samples=%d) with partition function computation, running time: %9.6f seconds\n", num_rnd, t1);
	exit(0);
  }

	printf("\nComputing minimum free energy structure...\n");
	fflush(stdout);

	t1 = get_seconds();
	energy = calculate(seq.length()) ; 
	t1 = get_seconds() - t1;
	
	printf("Done.\n\n");
	printf("Results:\n");
	if (energy >= MAXENG)	
		printf("- Minimum Free Energy: %12.4f kcal/mol\n", 0.00);
	else
		printf("- Minimum Free Energy: %12.4f kcal/mol\n", energy/100.00);
	printf("- MFE runtime: %9.6f seconds\n", t1);


	if (SUBOPT_ENABLED) {	
		t1 = get_seconds();
		ss_map_t subopt_data = subopt_traceback(seq.length(), 100.0*suboptDelta);
		t1 = get_seconds() - t1;
		printf("\n");
		printf("Subopt traceback running time: %9.6f seconds\n", t1);
		
		printf("Subopt structures saved in %s\n", suboptFile.c_str());
		save_subopt_file(suboptFile, subopt_data, seq, energy);	
		free_fold(seq.length());
		exit(0);
	}
	
	t1 = get_seconds();
	trace(seq.length()); //, VERBOSE, UNAMODE, T_MISMATCH);
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

	if(BPP_ENABLED){
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
	}

	// release the malloc'd arrays
	free_fold(seq.length());
/*
	dangle_struct partition;
	partition = malloc_partition_arrays_d(seq.length());
	fill_partition_arrays_d(partition);
	printf("Done with the partition functioni.\n");
*/
	return EXIT_SUCCESS;
}
