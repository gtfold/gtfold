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
#include "constraints.h"
#include "traceback.h"
#include "subopt_traceback.h"
#include "shapereader.h"

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
	printf("(c) 2007-2011  D.A. Bader, S. Mallidi, A. Mathuriya, C.E. Heitsch, S.C. Harvey\n");
	printf("Georgia Institute of Technology\n\n");
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
	
	// Read in thermodynamic parameters. Always use Turner99 data (for now)
	readThermodynamicParameters(paramDir.c_str(), PARAM_DIR);

	printRunConfiguration(seq);
	
	init_fold(seq);

	printf("\nComputing minimum free energy structure...\n");
	fflush(stdout);

	t1 = get_seconds();
	energy = calculate(seq.length(), nThreads);
	t1 = get_seconds() - t1;
	
	printf("Done.\n\n");
	printf("Results:\n");
	printf("- Minimum Free Energy: %12.4f kcal/mol\n", energy/100.00);
	printf("- MFE runtime: %9.6f seconds\n", t1);
	
	
	if (SUBOPT_ENABLED) {	
		t1 = get_seconds();
		subopt_traceback(seq.length(), suboptDelta);
		t1 = get_seconds() - t1;
		printf("Subopt traceback running time: %9.6f seconds\n\n", t1);

		free_fold(seq.length());
		exit(0);
	}
	
	t1 = get_seconds();
	trace(seq.length(), VERBOSE);
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
	// release the malloc'd arrays
	free_fold(seq.length());

	return EXIT_SUCCESS;
}
