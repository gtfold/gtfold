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
#include <limits>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main.h"
#include "mfe_main.h"
#include "boltzmann_main.h"
#include "subopt_main.h"
#include "utils.h"
#include "global.h"


int main(int argc, char** argv) {
  print_header();
  
  std::string cmd = argv[0];

  if (cmd.find("gtmfe") != std::string::npos) {
    mfe_main(argc, argv);
  } else if (cmd.find("gtsubopt") != std::string::npos){
    subopt_main(argc,argv);
  } else if (cmd.find("gtboltzmann") != std::string::npos) {
    boltzmann_main(argc,argv);
  } else{
	print_gtfold_usage_help();
  }

	return EXIT_SUCCESS;
}

	/*
	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}*/

  /*
	//init_fold(seq);
	
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
    //free_fold(seq.length());
    exit(0);
  }
  
  if (RND_SAMPLE == true)
  {
    printf("\nComputing partition function...\n");
	  int pf_count_mode = 0;
	  if(PF_COUNT_MODE) pf_count_mode=1;
	  double U = calculate_partition(seq.length(),pf_count_mode);
 
    batch_sample(num_rnd, seq.length(), U); 

    free_partition();
	  //free_fold(seq.length());
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
		//free_fold(seq.length());
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
  */
/*
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
*/
	// release the malloc'd arrays
	//free_fold(seq.length());
/*
	dangle_struct partition;
	partition = malloc_partition_arrays_d(seq.length());
	fill_partition_arrays_d(partition);
	printf("Done with the partition functioni.\n");
*/
	//return EXIT_SUCCESS;
//}
