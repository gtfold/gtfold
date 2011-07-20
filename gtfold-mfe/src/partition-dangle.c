//Attempting to follow Ding and Lawrence as closely as possible

/**s1 -> partial_external
*  s2 -> partial_multi
*  s3 -> partial_multi2
*  u1 -> u_multi
*  f -> cond_dangle (conditional)
*  ebi -> eL
**/

#include <math.h>
#include <stdlib.h>

#include "energy.h"
#include "utils.h"
#include "global.h"
#include "data.h"
#include "algorithms-partition.h"
#include "partition-dangle.h"

#define MIN_TURN 4
#define MAX_LOOP 30

//double RT = (0.00198721 * 310.15)/100.00; 

double cond_dangle(int j, int h, int l){
	if(j - 1 == l)
		return 1;
	else 
		return exp(-Ed3(h,l,l+1)/RT);
}

dangle_struct malloc_partition_arrays_d(int length){
	dangle_struct ret_struct;
	ret_struct.length = length;
	
	//Add 1 since we're not zero indexed
	ret_struct.u = mallocTwoD(length + 1, length + 1); 
	ret_struct.up = mallocTwoD(length + 1, length + 1); 
	ret_struct.upm = mallocTwoD(length + 1, length + 1); 
	ret_struct.u1 = mallocTwoD(length + 1, length + 1); 
	ret_struct.s1 = mallocTwoD(length + 1, length + 1); 
	ret_struct.s2 = mallocTwoD(length + 1, length + 1); 
	ret_struct.s3 = mallocTwoD(length + 1, length + 1); 

	int i,j;

	for(i = 0; i < length + 1; i++){
		for(j = 0; j < length + 1; j++){
			ret_struct.u[i][j] = 1;
			ret_struct.up[i][j] = 0;
			ret_struct.upm[i][j] = 0;
			ret_struct.u1[i][j] = 0;
			ret_struct.s1[i][j] = 0;
			ret_struct.s2[i][j] = 0;
			ret_struct.s3[i][j] = 0;
		}
	}
	return ret_struct;
}

void free_partition_arrays_d(struct partition_d part)
{
	freeTwoD(part.u, part.length, part.length);
	freeTwoD(part.up, part.length, part.length);
	freeTwoD(part.upm, part.length, part.length);
	freeTwoD(part.u1, part.length, part.length);
	freeTwoD(part.s1, part.length, part.length);
	freeTwoD(part.s2, part.length, part.length);
	freeTwoD(part.s3, part.length, part.length);
}

void fill_partition_arrays_d(dangle_struct part_struct){
	//So that we don't have the /100 in every exponent
	//Seeing as we missed some last time

	double ** partial_external = part_struct.s1;
	double ** partial_multi = part_struct.s2;
	double ** partial_multi2 = part_struct.s3;
	double ** u_multi = part_struct.u1;
	double ** u = part_struct.u;
	double ** up = part_struct.up;
	double ** upm = part_struct.upm;
	int len = part_struct.length;

	int seg_length;
	int i,j;

	for(seg_length = MIN_TURN; seg_length <= len; len++){
		//Insert parallelism here.
		for(i = 1; i < len - seg_length; i++){
			j = i + seg_length - 1;

			int l;
			if(canPair(i,j)){
				for(l = i+2; l < j; l++){
					upm[i][j] += up[i+1][l] * exp(-(Ea + 2 * Ec + auPenalty(i + 1,l)/RT)) *
								(exp(-Ed3(i + 1,l,l+1)/RT) * u_multi[l + 2][j - 1] + 
								u_multi[l + 1][j - 1] -
								u_multi[l + 2][j - 1]);
					if(l != i+2){ //goes from i + 2 < l < j
						upm[i][j] += up[i + 2][l] * 
								exp(-(Ea + 2 * Ec + Eb + Ed3(j,i,i + 1) + auPenalty(i + 2, l))/RT) *
								(exp(-Ed3(i + 2,l,l + 1)/RT)*u_multi[l + 2][j - 1] + 
								u_multi[l + 1][j - 1]  - 
								u_multi[l + 2][j - 1]);

						if(l != j - 1){
							upm[i][j] += exp(-Ed3(j,i,i + 1)/RT)*exp(-(Ea + 2 * Ec + (l - i - 1) * Eb)/RT) * 
										partial_multi[l][j]; 
							//Changed h in the paper to l it makes no difference
							//other than unifying loops
						}
					}
				}


				up[i][j] += exp(-eH(i,j)/RT) + exp(-eS(i,j)/RT) * up[i + 1][j - 1] + upm[i][j];
				for(l = j - 1; l > i + 1; l--){
					if((j - l) - 2> MAX_LOOP){
						break;
					}
					int h;
					for(h = i + 1 ; h < l; h++){
						if(!(i == h - 1 && j == l - 1)){ //If this is true we have a stack
							up[i][j] += exp(-eL(i,j,h,l)/RT);
						}
						if((j - l) + (h - i) - 2 > MAX_LOOP){
							break;
						}
					}
				}
			}//End checkpair conditional

			for(l = i + 1; l <= j ; l++){
				u_multi[i][j] += up[i][l] * exp(-(Ec + auPenalty(i,l))/RT) *
						(cond_dangle(j + 1,i,l) * exp(-(j - l) * Eb/RT) + 
						exp(-Ed3(i,l,l + 1)/RT) * u_multi[l + 2][j] + 
						u_multi[l + 1][j] -
						u_multi[l + 2][j]);
				if(l != i+1){	
					u_multi[i][j] += up[i + 1][j]*exp(-(Ec + Eb + auPenalty(i + 1,l))/RT)*
							(cond_dangle(j + 1, i + 1, l) * exp(-(j - l)*Eb/RT) + 
							exp(-Ed3(i + 1,l,l + 1)/RT) * u_multi[l + 2][j] + 
							u_multi[l + 1][j] - 
							u_multi[l + 2][j]);
					if(l != j){
						u_multi[i][j] += exp(-(Ec + (l - i) * Eb)/RT) * partial_multi2[l][j];

					}
				}
			}


			int h = i; //To stay consistant with the notation in the paper.

			for(l = h+1; l < j; l++){
				partial_external[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
							(exp(-Ed3(h,l,l + 1)/RT)*u[l + 2][j - 1] + 
							u[l + 1][j] -
							u[l + 2][j]);

				partial_multi[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
							(exp(-Ed3(h,l,l + 1)/RT)*u_multi[l + 2][j - 1] + 
							u_multi[l + 1][j - 1] -
							u_multi[l + 2][j - 1]);

				partial_multi2[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
							(cond_dangle(j + 1, h, l)*exp(-(j - l) * Eb / RT) + 
							exp(-Ed3(h,l,l + 1)/RT)*u_multi[l + 2][j] + 
							u_multi[l + 1][j] - 
							u_multi[l + 2][j]);
							
			}
			/** partial_multi2 goes up to j **/
			partial_multi2[h][j] += up[h][l] * exp( -(Ed5(h,l,h - 1) + auPenalty(h,l)) / RT) *
						(cond_dangle(j + 1, h, l)*exp(-(j-l) * Eb / RT) + 
						exp(-Ed3(h,l,l + 1)/RT) * u_multi[l + 2][j] + 
						u_multi[l + 1][j] - 
						u_multi[l + 2][j]);

			//Finally we do the u matrix
			u[i][j] = 1 + up[i][j]*exp(-auPenalty(i,j) /RT);

			for(l = i + 1; l < j; l++){
				//Replaced l with h in the next line to unify loops
				u[i][j] += up[l][j] * exp(-(Ed5(l, j , l - 1) + auPenalty(h,j))/RT);

				u[i][j] += exp(-(auPenalty(i,l))/RT)*(exp(-Ed3(i,l,l + 1)/RT) * u[l + 2][j] +
					   u[l + 1][j] - u[l + 2][j]);
				if(l != j - 1) {
					u[i][j] += partial_external[l][j];
				}
			}
		}//end of minor (paralellizeable) for loop
	}//End of major for loop
}
