#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <iterator>

#include "global.h"
#include "options.h"
#include "constraints.h"

int* BP;
int* ind;

typedef pair<int,int> basepair_t;

/*
ZS: Explanation to BP array. 

BP(i,j) for i!=j can have one of the following values: 
0: Nothing is required about a pair i,j in constraints
1: Force pairing between i and j. This implies also that
   any non-nested pairings will be prohibited and any 
   other pairs i and j would be involved in are also prohibited. 
2: Prohibit pairing between i and j. (Nothing else is done)

BP(i,i) can have one of the following values: 
0: Nothing is known about position i at all
3: Position i is prohibited from pairing (forced to be single-stranded). This 
   implies also that any pairs with i will be prohibited. 
4: Position i is forced to be in a particular pair (i is not single-stranded).
5: Position i is prohibited from pairing with at least one nucleotide.
6: Position i is forced to be in some pair (i is not single-stranded), but no pairing partner is specified.
*/

int** PBP;
int** FBP;

int nPBP;
int nFBP;

static bool CONS_ENABLED = false;
static bool LIMIT_DISTANCE = false;
static int contactDistance = -1;

void enable_constraints(int b) {CONS_ENABLED = b;}
void enable_limit_distance(int b) {LIMIT_DISTANCE = b;}
void set_contact_distance(int dist) {contactDistance = dist;}

static int load_constraints(const char* constr_file, int seq_length, int verbose=0) {
	
	fprintf(stdout, "- Running with constraints\n");

	std::ifstream cfcons;
	cfcons.open(constr_file, std::ios::in);
    if (cfcons == 0) {
        fprintf(stderr, "Error opening constraint file\n\n");
        cfcons.close();
        return -1;
    }

    char cons[100];

    while (!cfcons.eof()) {
        cfcons.getline(cons, 100);
        if (cons[0] == 'F' || cons[0] == 'f') nFBP++;
        if (cons[0] == 'P' || cons[0] == 'p') nPBP++;
    }
    cfcons.close();

    //fprintf(stdout, "Number of Constraints given: %d\n\n", nFBP + nPBP);
    if (nFBP + nPBP == 0) {
        fprintf(stderr, "No Constraints found.\n\n");
        return -1;
    }

    FBP = (int**) malloc(nFBP*sizeof(int*));
    PBP = (int**) malloc(nPBP*sizeof(int*));

    int fit = 0, pit = 0, it = 0;

    for (it = 0; it < nFBP; it++) {
        FBP[it] = (int*) malloc(2* sizeof (int));
    }
    for(it=0; it < nPBP; it++) {
        PBP[it] = (int*)malloc(2*sizeof(int));
    }
    cfcons.open(constr_file, std::ios::in);

    while(!cfcons.eof()) {
        cfcons.getline(cons,100);
        char *p=strtok(cons, " ");
        p = strtok(0, " ");
        if(cons[0]=='F' || cons[0]=='f') {
            int fit1=0;
            while(p!=0) {
                FBP[fit][fit1++] = atoi(p);
                p = strtok(0, " ");
            }
            fit++;
        }
        if(cons[0]=='P' || cons[0]=='p') {
            int pit1=0;
            while(p!=0) {
                PBP[pit][pit1++] = atoi(p);
                p = strtok(0, " ");
            }
            pit++;
        }
    }

	for(it=0; it< nFBP; it++) {
		if (FBP[it][0] < 1 || FBP[it][0] >= seq_length - TURN) {
			fprintf(stderr, "\nBase %d from constraint 'F %d %d %d' is out of bounds: For constraint 'F i j k' value i must be between 1 and the sequence length - %d, where the sequence length is %d \n",
				FBP[it][0], FBP[it][0], FBP[it][1], FBP[it][2], TURN, seq_length);
			exit(1);
		}
		if (FBP[it][1] != 0 && (FBP[it][1] <= FBP[it][0] + TURN 
								|| FBP[it][1] > seq_length)) {
			fprintf(stderr, "\nBase %d from constraint 'F %d %d %d' has an illegal value: For constraint 'F i j k' value j must be either 0 or between i + %d and the sequence length of %d \n",
				FBP[it][1], FBP[it][0], FBP[it][1], FBP[it][2], TURN, seq_length);
			exit(1);
		}
		if (FBP[it][2] < 1) {
			fprintf(stderr, "\nValue %d from constraint 'F %d %d %d' is too small: For constraint 'F i j k' value k must be at least 1\n",
				FBP[it][1], FBP[it][0], FBP[it][1], seq_length);
			exit(1);
		}
		if (FBP[it][1] != 0 && (FBP[it][1] - FBP[it][0] - 2*FBP[it][2] + 2 <= TURN)) {
			fprintf(stderr, "\nValue %d from constraint 'F %d %d %d' is too large: For constraint 'F i j k' the values must satisfy the inequality (j-k+1)-(i+k-1) > %d \n",
				FBP[it][2], FBP[it][0], FBP[it][1], FBP[it][2], TURN);
			exit(1);
		}
		if (FBP[it][1] == 0 && (FBP[it][0] + FBP[it][2] - 1 > seq_length)) {
			fprintf(stderr, "\nValue %d from constraint 'F %d %d %d' is too large: For constraint 'F i 0 k' the value i + k  - 1 can be at most the sequence length, in this case %d \n",
				FBP[it][2], FBP[it][0],FBP[it][1],FBP[it][2], seq_length);
			exit(1);
		}
	}
	
	for(it=0; it< nPBP; it++) {
		if (PBP[it][0] < 1 || (PBP[it][1] != 0 && PBP[it][0] >= seq_length - TURN)) {
			fprintf(stderr, "\nBase %d from constraint 'P %d %d %d' is out of bounds: For constraint 'P i j k' value i must be between 1 and the sequence length - %d, where the sequence length is %d \n",
				PBP[it][0], PBP[it][0], PBP[it][1], PBP[it][2], TURN, seq_length);
			exit(1);
		}
		if (PBP[it][1] != 0 && (PBP[it][1] <= PBP[it][0] + TURN 
								|| PBP[it][1] > seq_length)) {
			fprintf(stderr, "\nBase %d from constraint 'P %d %d %d' has an illegal value: For constraint 'P i j k' value j must be either 0 or between i + %d and the sequence length of %d \n",
				PBP[it][1], PBP[it][0], PBP[it][1], PBP[it][2], TURN, seq_length);
			exit(1);
		}
		if (PBP[it][2] < 1) {
			fprintf(stderr, "\nValue %d from constraint 'P %d %d %d' is too small: For constraint 'P i j k' value k must be at least 1\n",
				PBP[it][1], PBP[it][0], PBP[it][1], seq_length);
			exit(1);
		}
		if (PBP[it][1] != 0 && (PBP[it][1] - PBP[it][0] - 2*PBP[it][2] + 2 <= TURN)) {
			fprintf(stderr, "\nValue %d from constraint 'P %d %d %d' is too large: For constraint 'P i j k' the values must satisfy the inequality (j-k+1)-(i+k-1) > %d \n",
				PBP[it][2], PBP[it][0], PBP[it][1], PBP[it][2], TURN);
			exit(1);
		}
		if (PBP[it][1] == 0 && (PBP[it][0] + PBP[it][2] - 1 > seq_length)) {
			fprintf(stderr, "\nValue %d from constraint 'P %d %d %d' is too large: For constraint 'P i 0 k' the value i + k  - 1 can be at most the sequence length, in this case %d \n",
				PBP[it][2], PBP[it][0],PBP[it][1],PBP[it][2], seq_length);
			exit(1);
		}
	}
	
	return 0;
}

int init_constraints(const char* constr_file,int length) {
  enable_constraints(true);
  load_constraints(constr_file, length);
  

	int i,j,it,k,a,b;
	ind = (int*) malloc((length+1) * sizeof(int));
	if (ind == NULL) {
		perror("Cannot allocate variable 'ind'");
		exit(-1);
	}
	for(i = 1; i<=length; i++){
		ind[i] = (i*(i-1)) >>  1; //n(n-1)/2
	}


	BP = (int*) malloc((((length+1)*(length))/2+1)*sizeof(int));
	if (BP == NULL) {
       	perror("Cannot allocate variable 'constraints'");
	    exit(-1);
	}
	
	int LLL = length*(length+1)/2 + 1;


	//ZS: initialize all basepairing constraints to 0 (default is nothing known)
	for(i = 0; i<LLL; i++){
		BP[i] = 0;
	}

	
	//ZS: for ambiguous bases (right now this only handles 'N'), force single-stranded. 
	for(i = 1; i <= length; i++){
		if(RNA[i]=='N'){
			//force single-stranded
			BP(i,i) = 3;
			//Prohibit pairing with anything else
			for(j = i+1; j<=length; j++){
				BP(i,j) = 2;
			}
		}
	}
	
	//CM: set forced basepairs
	for (it = 0; it < nFBP; it++){
		
		if(FBP[it][1] == 0){	// F i 0 k
			// For bases i through i+k-1, say they must be in some pair
			for (i = FBP[it][0]; i < FBP[it][0]+FBP[it][2]; i++) {
				if(BP(i,i) == 0){
					BP(i,i) = 6;
				}
				//If BP(i,i) == 4 or 6, ignore
			}
		}
		
		else {	// F i j k, j != 0
			// Require the stack from (i,j) to (i+k-1,j-k+1)
			
			for (k = 1; k <= FBP[it][2]; k++) {
				
				i = FBP[it][0] + k - 1;
				j = FBP[it][1] - k + 1;
				
				if (BP(i,i) == 4 && BP(i,j) != 1) {	// Base i is already in a pair
					printf("Constraint 'F %d %d %d' is trying to force the pair (%d,%d), but base %d is forced to be in a different pair by a previous 'F i j k' constraint\n",
						FBP[it][0],FBP[it][1],FBP[it][2],i,j,i);
					exit(1); // FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				}
				else if (BP(j,j) == 4 && BP(i,j) != 1) { // Base j is already in a pair
					printf("Constraint 'F %d %d %d' is trying to force the pair (%d,%d), but base %d is forced to be in a different pair by a previous 'F i j k' constraint\n",
						FBP[it][0],FBP[it][1],FBP[it][2],i,j,j);
					exit(1);
				}
				else if (BP(i,j) == 2) {	// Pseudoknot!
					printf("Constraint 'F %d %d %d' is trying to force the pair (%d,%d), but this pair conflicts with a previous 'F i j k' constraint\n",
						FBP[it][0],FBP[it][1],FBP[it][2],i,j);
					printf("\n");
					exit(1);
				}
				else if (!canPair(RNA[i], RNA[j])) {
					printf("Constraint 'F %d %d %d' is trying to force the pair (%d,%d), but the pair is non-canonical\n",
						FBP[it][0],FBP[it][1],FBP[it][2],i,j);
					exit(1);
				}
				else {
					
					// Require the pair (i,j)
					BP(i,i) = 4;
					BP(j,j) = 4;
					BP(i,j) = 1;
					
					// Nothing else can pair with i or j, and no pseudoknots
					for (b = i+1; b < j; b++) {
						BP(i,b) = 2;
						BP(b,j) = 2;
					}
					// To avoid rewriting, these only need to be set for k = 1
					if (k == 1) {
						for (a = 1; a < i; a++) {
							for (b = i; b <= j; b++) {
								BP(a,b) = 2;
							}
						}
						for (a = i; a <= j; a++) {
							for (b = j+1; b <= length; b++) {
								BP(a,b) = 2;
							}
						}
					}
				}
			}
		}
	}
	
	
	//CM: set prohibited basepairs
	for (it = 0; it < nPBP; it++){
	
		if (PBP[it][1] != 0) { // P i j k, j != 0
			// Prohibit the stack (i,j) to (i+k-1,j-k+1)
			
			for (k = 1; k <= PBP[it][2]; k++) {
				
				i = PBP[it][0] + k - 1;
				j = PBP[it][1] - k + 1;
				
				// Pair (i,j) is required by some constraint
				if (BP(i,j) == 1) {
					printf("Constraint 'P %d %d %d' is trying to prohibit pair (%d,%d), but this pair is required by a 'F i j k' constraint\n",
						PBP[it][0],PBP[it][1],PBP[it][2],i,j);
					exit(1);
				}
				
				else {
					
					BP(i,j) = 2;
					BP(i,i) = 5;	// This is only for the sake of print_constraints
					BP(j,j) = 5;	// This is only for the sake of print_constraints
					
				}
				
			}
		
		
		}
	
		else { // P i 0 k
			// For bases i to i+k-1, say they must be single-stranded
			for (i = PBP[it][0]; i < PBP[it][0] + PBP[it][2]; i++) {
				
				// Base i is in a required pair
				if (BP(i,i) == 4 || BP(i,i) == 6) {
					printf("Constraint 'P %d %d %d' is trying to prohibit base %d from being paired, but this base is required to be paired by a 'F i j k' constraint\n",
						PBP[it][0],PBP[it][1],PBP[it][2],i);
					exit(1);
				}
				
				else {
					BP(i,i) = 3;
					
					for (a = 1; a < i; a++) {
						BP(a,i) = 2;
					}
					for (a = i+1; a <= length; a++) {
						BP(i,a) = 2;
					}
				}
				
			}
			
		
		}
	
	}

	return 0;
}

int verify_structure(){
	//ZS: This method returns true if the structure (global.h) is consistent with
	//the constraints, and false if it is not.
	
	if(CONS_ENABLED){
	int errorhappened = 0;
	int it, k;
	//Check prohibited constraints
	if(nPBP != 0){
		for(it = 0; it < nPBP; it++){
			if(PBP[it][2] < 1 || PBP[it][1] == 0){
				//printf("Invalid entry (P: %d %d %d)\n", PBP[it][0], PBP[it][1], PBP[it][2]);
				continue;
			}
			for(k = 1; k <= PBP[it][2]; k++){
				//correct answer: strcuture(PBP[it][0]+k-1) != structure(PBP[it][1]-(k-1))
				if(structure[PBP[it][0]+k-1] == PBP[it][1]-(k-1) || structure[PBP[it][1]-(k-1)] == PBP[it][0]+k-1){	
					errorhappened = 1; 
					fprintf(stderr, "Constraint P %d %d %d is not fulfilled.\n",PBP[it][0], PBP[it][1], PBP[it][2]);
					break;
				}
			}
		}	
	}

	//Check forced constraints
	if(nFBP != 0){
		for(it = 0; it<nFBP; it++){
			if(FBP[it][2] < 1){
				//printf("Invalid entry (F: %d %d %d)\n", FBP[it][0], FBP[it][1], FBP[it][2]);
				continue;
			}

			for(k = 1; k<=FBP[it][2]; k++){
				int i1 = FBP[it][0]+k-1;
				int j1 = FBP[it][1]-k+1;
				if(FBP[it][1]!=0&&!canPair(RNA[FBP[it][0]+k-1], RNA[FBP[it][1]-k+1])){
					//printf("Can't force (%d, %d) to pair (non-canonical) \n", FBP[it][0]+k-1, FBP[it][1]-k+1);
					continue;			
				}
				if(FBP[it][1]!=0&&(j1-i1 < TURN)){
					//printf("Can't force (%d, %d) to pair (turn too tight) \n", FBP[it][0]+k-1, FBP[it][1]-k+1);
					continue;
				}
				if(FBP[it][1] == 0){
					//force single-stranded
					if(structure[FBP[it][0]] != 0){
						fprintf(stderr, "Constraint F %d %d %d is not fulfilled.\n",FBP[it][0], FBP[it][1], FBP[it][2]);
						errorhappened = 1;
					}
				}
				else{
					if(structure[FBP[it][0]+k-1] != FBP[it][1]-(k-1) || structure[FBP[it][1]-(k-1)] != FBP[it][0]+k-1){
						fprintf(stderr, "Constraint F %d %d %d is not fulfilled.\n",FBP[it][0], FBP[it][1], FBP[it][2]);
						errorhappened = 1;
					}
					
				}
			}
		}
	}
		return errorhappened?0:1;
	}
	else{
		return 1;
	}
		 
}


void free_constraints(int len) {
	free(BP);
}

void print_constraints(int len) {

    int i = 1;
    int j = 1;
/*

    printf("Printing constraints \n");


    for(j = 1; j<=len; j++){
	for(i = 1; i<=j; i++){
		printf("BP(%d,%d)=%d\t", i,j,BP(i,j));
	}
	printf("\n\n");
    }

*/
    for (i = 1; i <= len; i++) {
	    switch(BP(i,i)){
		case 3:
			//printf("%d: .\n", i);break;
			printf("x");break;
		case 4:
			//printf("%d - ", i);
			for(j = 1;j<=len;j++){ 
				if(i<j&&BP(i,j)==1){
					//printf("%d",j);
					printf("(");
				}
				else if(j<i&&BP(j,i)==1){
					//printf("%d",j);
					printf(")");
				}
			}
			//printf("\n");
			break;
		case 5: 
			//printf("%d P \n", i); break;
			printf("P");break;
		case 0:
			//printf("%d not constrained\n",i);
			printf(".");
			break;
		case 6:
			printf("|");
			break;
		default:
			fprintf(stderr,"ERROR in constraint value, debugging info: i=%d, BP(i,i)=%d",i,BP(i,i));break;
	    }  
    }
    printf("\n");
 
}

int forcePaired(int i){
//ZS: this function returns true if position i is forced to be paired, but partner is unknown.
	if(CONS_ENABLED)
		return BP(i,i)==6; 
	else
		return 0; 
}

int forcePair(int i, int j) {
//ZS: this function returns true if the pair i-j is forced.
	if (CONS_ENABLED) 
		return BP(i,j)==1;
	else
		return 0;
}

int forceSS(int i){
//ZS: this function returns true if base i is forced to be single-stranded 
	if(CONS_ENABLED)
		return BP(i,i)==3;
	else 
		return 0;
}

int forceSSregion(int i, int j){
	//ZS: This returns 1 if all nucleotides between i and j are forced to be single-stranded, 0 if there is at least one that is not
	if(CONS_ENABLED){
		for(int p = i; p<j; p++){
			if(BP(p,p)!=3){
				return 0;
			}
		}
		return 1;
	}
	else return 0;
}

int canILoop(int i, int j, int p, int q) {
//ZS: This function returns 1 if internal loop with pairs i,j and p,q is allowed, 0 if not allowed
	if (CONS_ENABLED){
		//Need to check that i,j and p,q pairs are allowed, 
		//and single-stranded regions between i,p and q,j are allowed
		return canStack(i,j) && canStack(p,q) && canSSregion(i,p) && canSSregion(q,j);
	}
	else 
		return withinCD(i,j) && withinCD(p,q);
}

int canHairpin(int i, int j) {
	//ZS: Returns 1 if a hairpin between i,j is allowed, 0 if not allowed
	if (CONS_ENABLED){
		//Need to check if pair is allowed and if anything pairs between them
		return canStack(i,j) && canSSregion(i,j);
	}	
	else return withinCD(i,j);
}

int canSSregion(int i, int j){
//ZS: This function returns 0 if any nucleotide between i and j is forced to pair with something
//(i and j are NOT included in the check)
	if(CONS_ENABLED){
		for(int p = i+1; p<j; p++){
			if(BP(p,p)==4||BP(p,p)==6) return 0;
		}
		return 1;}
	else{ return 1; }
}

int canSS(int i){
	//ZS: This function returns 0 if i is forced to pair with something, 1 if it can be single stranded. 
	if (CONS_ENABLED){
		return 1-(BP(i,i)==4||BP(i,i)==6);
	}
	else{
		return 1;
	}
}

int canStack(int i, int j) {
// ZS: This function returns 1 if an i,j pair is allowed, 0 if not allowed
// (Unfortunately the more intuitive "canPair" name is already used to mark complementary nucleotides so 
// I couldn't use that again) 
	if (CONS_ENABLED){ 
		//can't pair if i,j is prohibited, or i or j are forced to be single-stranded
		if(BP(i,j)==2||BP(i,i)==3||BP(j,j)==3) return 0;
		//can't pair if i or j are forced to pair with something other than each other
		if((BP(i,i)==4||BP(j,j)==4)&&BP(i,j)!=1) return 0;	
		else return withinCD(i,j);
	}
	else return withinCD(i,j);
}

int withinCD(int i, int j){
	//ZS: This function returns 1 if i and j are within the required contact distance from each other. 
	if (LIMIT_DISTANCE){
		return (j-i<=contactDistance);
	}
	else return 1;
}

