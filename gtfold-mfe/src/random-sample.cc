#include <math.h>
#include <stdlib.h>
#include <list>

#include "partition-dangle.h"
#include "algorithms-partition.h"
#include "energy.h"
#include "random-sample.h"

using namespace std;

#define MIN_TURN 4

//double RT=(0.00198721 * 310.15)/100.00; 


//I don't know why not having this function
//gives me a link error it should get this function
//from partition-dangle

/*double cond_dangle(int j, int h, int l){
	if(j - 1 == l)
		return 1;
	else 
		return exp(-Ed3(h,l,l+1)/RT);
}*/

double randdouble(){
	    return rand()/(double(RAND_MAX)+1);
} 

//P_0
double single_stranded_prob(int i, int j, dangle_struct d_struct){
	return 1 / d_struct.u[i][j];
}

//P_{ij}
double ends_pair_prob(int i, int j, dangle_struct d_struct){
	return d_struct.up[i][j] * exp(-auPenalty(i,j)/RT)/d_struct.u[i][j];
}

//P_{hj}
double prob_h_pairs_with_5_end
		(int i, int h, int j, dangle_struct d_struct){

	return d_struct.up[h][j] * 
			exp(-(Ed5(h,j,h - 1) + 
			auPenalty(h,j))/RT) /
			d_struct.u[i][j];
}

//P_{il}
double prob_l_pairs_with_3_end
	(int i, int l, int j, dangle_struct d_struct){

	return d_struct.up[i][l] * exp(-auPenalty(i,l)/RT) *
			(exp(-Ed3(i,l,l+1)/RT) * d_struct.u[l + 2][j] +
			d_struct.u[l + 1][j] -
			d_struct.u[l + 2][j]) /
			d_struct.u[i][j];
}

//P_{s1h}
double prob_h_pairs_in_ij(int i, int h, int j, dangle_struct d_struct){
	return d_struct.s1[h][j] / d_struct.u[i][j];
}

//P_{hl}
//This is the probability h pairs with l 
//**Given that h pairs before j**
double prob_h_pairs_with_l(int h, int l, int j, dangle_struct d_struct){
	return d_struct.up[h][l] * exp(-(Ed5(h,l,h - 1) + auPenalty(h,l))/RT) *
			(exp(-Ed3(h,l,l + 1)/RT) * d_struct.u[l + 2][j] +
			d_struct.u[l + 1][j] -
			d_struct.u[l + 2][j]) /
			d_struct.s1[h][j];
}

//Q_{ijH}
double prob_ij_pair_is_hairpin(int i, int j, dangle_struct d_struct){
	return exp(-eH(i,j)/RT)/d_struct.up[i][j];
}

//Q_{ijS}
double prob_ij_pair_is_stack(int i, int j, dangle_struct d_struct){
	return exp(-eS(i,j)/RT) * d_struct.up[i + 1][j - 1]/d_struct.u[i][j];
}
//To find Q_{ijBl} take the rest of the probabilities and subtract
//them from 1.

//Q_{ijM}
double prob_ij_pair_is_multiloop(int i, int j, dangle_struct d_struct){
	return d_struct.upm[i][j]/d_struct.up[i][j];
}

//This the probability hl closes a multiloop
//**GIVEN IJ CLOSES A MULTILOOP**
double prob_ijhl_is_internal_loop
	(int i, int j, int h, int l, double partition_sum, dangle_struct d_struct){
	return exp(-eL(i,j,h,l)/RT) * d_struct.up[h][l] / partition_sum;
}

//P_{ij(1+1)l}
double first_multi_prob_pair_next_base
	(int i, int l, int j, dangle_struct d_struct){
	return d_struct.up[i + 1][l] * 
			(exp( - (Ea + 2 * Ec + auPenalty(i + 1, l))/RT) *
			(exp(-Ed3(i + 1,l,l +1))/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 1][j - 1]) /
			d_struct.upm[i][j];
}

//P_{ij(1+2)l}
double first_multi_prob_pair_after_base
	(int i, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[i + 2][l] * 
			exp(-(Ea + 2 * Ec + Eb + Ed3(j, i, i + 1) + auPenalty(i + 2,l))/RT) *
			(exp(-Ed3(i + 2,l,l+1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 2][j - 1]) /
			d_struct.upm[i][j];
}

//P_{ijs2h)
double first_multi_prob_pair_rest_base(int i, int h, int j, dangle_struct d_struct){
	return exp(-(Ea + 2 * Ec + (h - i - 1) * Eb + Ed3(j,i,i+1))/RT) *
			d_struct.s2[h][j] / d_struct.upm[i][j];
}
//P_{ijhl}
double first_multi_prob_base_rest
	(int i, int h, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[h][l] * exp(-(Ed5(h,l,h-1) + auPenalty(h,l))/RT) * 
			(exp( -Ed3(h,l,l+1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 2][j - 1])/d_struct.s2[h][j];
}

//Q_{(l1 + 1)(j -1)(l1 +1)}
double rest_multi_prob_next_base_pairs (int l1, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[l1 +1][l] * exp(-(Ec + auPenalty(l1 + 1, l))/RT) *
			(cond_dangle(j, l1 + 1, l) * exp( -(j - 1 - l) * Eb /RT) +
			exp(-(j - 1 - l) * Eb) / RT + 
			exp(-Ed3(l1 + 1, l, l + 1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 1][j - 1]) /
			d_struct.u1[l1 + 1][j - 1];
}

//Q_{(l1 + 1)(j - 1)(l1 + 2)l}
double rest_multi_prob_after_base_pairs
	(int l1, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[l1 + 2][l] * exp(-(Ec + Eb + auPenalty(l1 + 2,l))/RT) *
			(cond_dangle(j, l1 + 2, l) * exp( -(j - 1 - l) * Eb / RT) +
			exp(-Ed3(l1 + 1,l ,l + 1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] - 
			d_struct.u1[l + 2][j - 1]) /
			d_struct.u1[l1 + 1][j - 1];
			
}

//Q_{(l1 + 1)(j - 1)s3h}
double rest_multi_prob_rest_base_pairs
	(int l1, int h, int j, dangle_struct d_struct){
	
	return exp( -(Ec + (h - l1 - 1) * Eb)/RT) * d_struct.s3[h][j - 1] /
			d_struct.u1[l1 + 1][j - 1];
}

//Q_{(j-1)hl}
//GIVEN THAT H PAIRS BEFORE J
double rest_multi_prob_h_pairs_with_l
	(int l1, int h, int l, int j, dangle_struct d_struct){
	
	return d_struct.up[h][l] * exp(-(Ed5(h,l,h - 1) + auPenalty(h,l))/RT) *
			(cond_dangle(j,h,l) * exp(-(j - 1 -l) * Eb/RT) + 
			exp(-Ed3(h,l,l + 1)/RT) * d_struct.u1[l + 2][j - 1] +
			d_struct.u1[l + 1][j - 1] -
			d_struct.u1[l + 2][j - 1]) /
			d_struct.s3[h][j - 1];
}

//P_{Bh2l2}
int there_is_another_helix(int h2, int l2, int j, dangle_struct d_struct){
	double counter = randdouble();
	double coin_flip = cond_dangle(j,h2,l2) * exp(-(j - 1 - l2) * Eb / RT) /
					  (cond_dangle(j,h2,l2) * exp(-(j - 1 - l2) * Eb / RT) + 
					exp(Ed3(h2,l2,l2 + 1)/RT) * d_struct.u1[l2 + 2][j - 1] +
					d_struct.u1[l2 + 1][j - 1] -
					d_struct.u1[l2 + 2][j - 1]);
	
	//If I was slick and mean I would have put
		
	//return counter > coin-flip

	if(counter <= coin_flip){
		return 0;
	}
	else{
		return 1;
	}
}

void single_stranded(int i, int j, int * structure){
	for(;i <= j; i++)
	{
		structure[i] = 0;
	}
}

void unknown_strand(int i, int j, dangle_struct d_struct, int * structure, list<sub_seq> * stack){
//Unknown strand from i,j
	double counter = randdouble();
	counter -= single_stranded_prob(i, j, d_struct);
	if (counter <= 0){
		single_stranded(i,j, structure);
		return;
	}
	counter -= ends_pair_prob(i, j, d_struct);
	if(counter <= 0){
		sub_seq seq;
		seq.start = i;
		seq.end = j;
		seq.paired = 1;
		//put i j on paired and push i j on stack as paired
		
		structure[i] = j;
		structure[j] = i;
		stack->push_front(seq);
		return;
	}
	for (int h = i + MIN_TURN; h < j; h++){
		counter -= prob_h_pairs_with_5_end(i,h,j,d_struct);
		if(counter < 0){
			sub_seq seq1;
			sub_seq seq2;

			seq1.start = i;
			seq1.end = h;
			seq1.paired = 1;

			seq2.start = h + 1;
			seq2.end = j;
			seq2.paired = 0;
			
			//Push i h as paired and h + 1 to j as unknown
			structure[i] = h;
			structure[h] = i;
			stack->push_front(seq1);
			stack->push_front(seq2);
			return;
		}
	}

	for(int l = i + 1; l <= j - MIN_TURN; l++){
		counter -= prob_l_pairs_with_3_end(i,l,j,d_struct);
		if (counter < 0){
			
			sub_seq seq;
			seq.start = l;
			seq.end = j;
			seq.paired = 1;

			//put from i to l - 1 as unpaired and push l and j as paired
			single_stranded(i, l - 1, structure);
			structure[l] = j;
			structure[j] = l;
			stack->push_front(seq);
			return;
		}
	}
	int h; //The loop is just to find h we need it later
	for(h = i + 1; i <= j - MIN_TURN - 1; h++){
		counter -= prob_h_pairs_in_ij(i,h,j,d_struct);
		if(counter < 0)
		{
			break;
		}
	}

	//Now to find l
	counter = randdouble(); //reset counter, this is a seperate case
	for(int l = h + MIN_TURN; l <= j - MIN_TURN; l++){
		counter -= prob_h_pairs_with_l(h,l,j,d_struct);	
		if (counter < 0){
			sub_seq seq1;
			sub_seq seq2;
			
			seq1.start = h;
			seq1.end = l;
			seq1.paired = 1;

			seq2.start = l + 1;
			seq2.end = j;
			seq2.paired = 0;

			//set from i to h - 1 as single stranded
			//Push h l as paired and l + 1 to j
			//as unknown.

			single_stranded(i, h - 1, structure);
			structure[h] = l;
			structure[l] = h;
			stack->push_front(seq1);
			stack->push_front(seq2);
			return;
		}	
	}
	//If we hit here, that's bad
}

void paired_strand(int i, int j, dangle_struct d_struct, int * structure, list<sub_seq> * stack){
	double counter = randdouble();
	double cumulative_prob = 0;

 
	cumulative_prob = prob_ij_pair_is_hairpin(i,j,d_struct);
	if(counter - cumulative_prob < 0){
		//i + 1, j - 1 is single stranded
		single_stranded(i + 1, j - 1, structure);
		return;
	}

	cumulative_prob += prob_ij_pair_is_stack(i,j,d_struct);
	if(counter - cumulative_prob < 0){
		sub_seq seq;
		seq.start = i + 1;
		seq.end = j - 1;
		seq.paired = 1;

		//i + 1, j - 1 is pushed as stack
		stack->push_front(seq);
		structure[i + 1] = j - 1;
		structure[j - 1] = i + 1;
		return;
	}

	cumulative_prob += prob_ij_pair_is_multiloop(i,j,d_struct);
	if(counter - cumulative_prob < 0){
		multi_loop_strand(i,j,d_struct,structure,stack);
		return;
	}

	//else we know we have an internal loop.
	//The probabilities should work out so that
	//we never get more than 30 base internal loops.
	//Might need to put guards in here.
	counter = randdouble();
	for(int h = i + 1; h < j; h++){
		for(int l = h + MIN_TURN; l < j; l++){
			if(!(h == i + 1 && l == j - 1)) //Don't want a stack
				counter -= exp(-eL(i,j,h,l)/RT) * d_struct.up[h][l] / (1 - cumulative_prob); 
				if(counter < 0){
					sub_seq	seq;
					seq.start = h;
					seq.end = l;
					seq.paired = 1;

					//hl goes on the stack as a pair and the rest of the internal loop
					//gets put on a single stranded

					stack->push_front(seq);
					structure[h] = l;
					structure[l] = h;
					single_stranded(i + 1, h - 1,structure);
					single_stranded(l + 1, j - 1,structure);
				}
		}
	}
}

//It returns where it stops.
int multi_first_loop(int i, int j, dangle_struct d_struct, int * structure, list<sub_seq> * stack){
	double counter = randdouble();
	for(int l= i + MIN_TURN + 1; l < j - 1; l++){
		counter -= first_multi_prob_pair_next_base(i,l,j,d_struct);
		if(counter < 0){
			sub_seq seq;
			seq.start = i + 1;
			seq.end = l;
			seq.paired = 1;

			//Under this condition i + 1 pairs with l;
			
			stack->push_front(seq);
			structure[i + 1] = l;
			structure[l] = i + 1;
			return l + 1;
		}
	}
	for(int l = i + MIN_TURN + 2; l < j - 1; l++){
		counter -= first_multi_prob_pair_after_base(i,l,j,d_struct);
		if(counter < 0){
			sub_seq seq;
			seq.start = i + 2;
			seq.end = l;
			seq.paired = 1;

			//Under this condition i + 2 pairs with l;
			//and i + 1 is unpaired
			
			stack->push_front(seq);
			structure[i + 1] = 0;
			structure[i + 2] = l;
			structure[l] = i + 2;
			return l + 1;
		}
	}
	
	//We need h after this loop so we declare it here.
	int h;
	for(h = i + 3; h < j - 1; h ++){
		counter -= first_multi_prob_pair_rest_base(i,h,j,d_struct);
		if(counter < 0){
			break;
		}
	}
	
	//Now h is the opening base in a base pair now we need to find the closing one
	//So first we start a new counter
	counter = randdouble();
	for (int l = h + MIN_TURN; l < j - 1; l++){
		counter -= first_multi_prob_base_rest(i,h,l,j,d_struct);
		if(counter < 0){
			sub_seq seq;
			seq.start = h;
			seq.end = l;
			seq.paired = 1;
			
			//H pairs with l and i + 1 to h - 1 is single stranded
			single_stranded(i + 1, h - 1, structure);
			stack->push_front(seq);
			structure[h] = l;
			structure[l] = h;
			return l + 1;
		}
	}
	
	return -1; //We have an underflow or calculation error.
}

void multi_loop_strand(int i, int j, dangle_struct d_struct, int * structure, list<sub_seq> * stack){
	int l1 = multi_first_loop(i,j,d_struct, structure,stack);//l1 is the name of the closing base of the first loop
	int h1 = i;
	int notdone = 1;
	int breakflag = 0;
	
	while(notdone){
		double counter = randdouble();
		for(int l = l1 + MIN_TURN + 1; l < j - 1; l++){
			counter -= rest_multi_prob_next_base_pairs(l1,l,j,d_struct);
			if(counter < 0){
				sub_seq seq;
				seq.start = l1 + 1;
				seq.end = l;
				seq.paired = 1;

				//Under this condition i + 1 pairs with l;
				
				stack->push_front(seq);
				structure[l1 + 1] = l;
				structure[l] = i + 1;
				h1 = l1 + 1;
				l1 = l;
				breakflag = 1;
				break;
			}
		}
		if(breakflag){
			notdone = there_is_another_helix(h1, l1, j, d_struct);
			breakflag = 0;
			continue; //Starts while loop over
		}
		for(int l = l1 + MIN_TURN + 2; l < j - 1; l++){
			counter -= rest_multi_prob_after_base_pairs(l1,l,j,d_struct);
			if(counter < 0){
				sub_seq seq;
				seq.start = l1 + 2;
				seq.end = l;
				seq.paired = 1;

				//Under this condition i + 2 pairs with l;
				//and i + 1 is unpaired
				
				stack->push_front(seq);
				structure[l1 + 1] = 0;
				structure[l1 + 2] = l;
				structure[l] = l1 + 2;
				h1 = l1 + 2;
				l1 = l;
				breakflag = 1;
				break;
			}
		}
		if(breakflag){
			notdone = there_is_another_helix(h1, l1, j, d_struct);
			breakflag = 0;
			continue; //Starts while loop over
		}
		//We need h after this loop so we declare it here.
		int h;
		for(h = l1 + 3; h < j - 1; h ++){
			counter -= rest_multi_prob_rest_base_pairs(i,h,j,d_struct);
			if(counter < 0){
				break;
			}
		}
		
		//Now h is the opening base in a base pair now we need to find the closing one
		//So first we start a new counter
		counter = randdouble();
		for (int l = h + MIN_TURN; l < j - 1; l++){
			counter -= rest_multi_prob_h_pairs_with_l(l1,h,l,j,d_struct);
			if(counter < 0){
				sub_seq seq;
				seq.start = h;
				seq.end = l;
				seq.paired = 1;
				
				//H pairs with l and i + 1 to h - 1 is single stranded
				single_stranded(l1 + 1, h - 1, structure);
				stack->push_front(seq);
				structure[h] = l;
				structure[l] = h;
				l1 = l + 1;
				h1 = h;
			}
		}
		notdone = there_is_another_helix(h1, l1, j, d_struct);
	}//End while
}//End function


void sample_structure(int * structure, dangle_struct d_struct){
	sub_seq seq;
	seq.start = 1;
	seq.end = d_struct.length;
	seq.paired = 0;

	list<sub_seq> da_stack;
	da_stack.push_front(seq);
	while(!da_stack.empty()){
		sub_seq cur_seq = da_stack.front();
		if(cur_seq.paired)
		{
			paired_strand(cur_seq.start, cur_seq.end, d_struct, structure, &da_stack);
		}
		else
		{
			unknown_strand(cur_seq.start, cur_seq.end, d_struct, structure, &da_stack);
		}
		da_stack.pop_front();
	}
	
}
