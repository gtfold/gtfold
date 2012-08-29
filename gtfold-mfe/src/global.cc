#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "utils.h"
#include "global.h"
#include "constraints.h"

unsigned char *RNA; 
int *structure; 
unsigned int chPairKey;

int g_nthreads;
int g_dangles;
int g_unamode;
int g_mismatch;
int g_verbose = 0;
int g_prefilter_mode;
int g_prefilter1;
int g_prefilter2;

int SHAPE_ENABLED = 0;
int g_LIMIT_DISTANCE;
int g_contactDistance;

void init_global_params(int len) {
	RNA = (unsigned char *) malloc((len+1)* sizeof(unsigned char));
	if (RNA == NULL) {
		perror("Cannot allocate variable 'RNA'");
		exit(-1);
	}
	structure = (int *) malloc((len+1) * sizeof(int));
	if (structure == NULL) {
		perror("Cannot allocate variable 'structure'");
		exit(-1);
	}

	init_checkPair();
}

void free_global_params() {
	free(structure);
	free(RNA);
}

void print_sequence(int len) {
	int i;
	for (i = 1; i <= len; i++) {
		if (RNA[i] == BASE_A)
			printf("A");
		else if (RNA[i] == BASE_C)
			printf("C");
		else if (RNA[i] == BASE_G)
			printf("G");
		else if (RNA[i] == BASE_U)
			printf("U");
		else
			printf("N");
	}
	printf("\n");
}

void print_structure(int len) {
	int i = 1;
	for (i = 1; i <= len; i++) 
	{
		if (structure[i] > 0 && structure[i] > i)
			printf("(");
		else if (structure[i] > 0 && structure[i] < i)
			printf(")");
		else
			printf(".");
	}
	printf("\n");
}

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

void print_gtfold_usage_help() {
	//printf("There are three programs provided by GTfold. \n 1. gtmfe: Program implementing and supporting various options for calculating minimum free energy and MFE structure for a given sequence. Run 'gtmfe --help' for more help on this.\n 2. gtboltzmann: Program implementing and supporting various options and features for sampling structures for a sequence stochastically with help of partition function. Run 'gtboltzmann --help' for more help on this.\n 3. gtsubopt: Program implementing and supporting vcarious options and features for getting sub-optimal structures for a given sequence. Run 'gtsubopt --help' for more help on this.\n\nPlease visit our sourceforge page for more details,\nhttp://gtfold.sourceforge.net/\n\n");
	printf("The GTfold package consists of three programs.\n 1. gtmfe: Computes MFE structure for a given sequence. Use the command ‘gtmfe --help’ for more details.\n 2. gtboltzmann: Supports various features such as stochastic sampling and computation of base pair probabilities under the Boltzmann Distribution. Use the command ‘gtboltzmann --help’ for more details\n 3. gtsubopt: Computes all structures with in a given energy range of the MFE for a given sequence. Run ‘gtsubopt --help’ for more details.\n\n\n Please visit our sourceforge page for full details on all three programs\n http://gtfold.sourceforge.net/\n");
}


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

void save_ct_file(string outputFile, string seq, int energy, int *structure1) {
        ofstream outfile;
        outfile.open(outputFile.c_str());

        outfile << seq.length() << "\t  dG = " << energy/100.0 << endl;
        //outfile << seq.length() << "\tdG = " << energy/100.0 << "\t" << seqfile << endl;

        unsigned int i = 1;
        for(i=1; i <= seq.length(); i++)
                outfile << i << "\t" << seq[i-1] << "\t" << i-1 << "\t" << (i+1)%(seq.length()+1) << "\t" << structure1[i] << "\t" << i << endl;

        outfile.close();
}

void init_checkPair() {
	int i, j;
	chPairKey = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			chPairKey += checkPair(i, j) << ((i << 2) + j);
}

int update_checkPair(int i, int j) {
	int r = 0;
	if (!((i >= 0 && i <=3 )&&(j >=0 && j <=3)))
		return r;
	if (!(chPairKey & (1 << ((i << 2) + j)))) {
		chPairKey += 1 << ((i << 2) + j);	
		r = 1;
	}
	return r;
}

int canPair(int a, int b) {
	return (chPairKey & (1 << (((a)<<2) + (b))));
}

/*
void help() {
    printf("Usage: gtfold [OPTION]... FILE\n\n");

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

    printf("\nBETA OPTIONS\n");
    printf("   --bpp                Calculate base pair probabilities.\n");
    printf("   --partition          Calculate the partition function.\n");
    printf("   --pf_count          Calculate the structure count using partition function and zero energy value.\n");
    printf("   --subopt NUM         Calculate suboptimal structures within NUM kcal/mol\n");
    printf("                        of the MFE. (Uses -d 2 treatment of dangling energies.)\n");
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
*/
