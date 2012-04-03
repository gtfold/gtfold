#include "stochastic-sampling-d2.h"
#include "global.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stack>
#include <map>
#include<sstream>
#include<fstream>
#include "utils.h"
#include<omp.h>

//Basic utility functions
void StochasticTracebackD2::initialize(int length1, int PF_COUNT_MODE1, int NO_DANGLE_MODE1, int ss_verbose1, bool PF_D2_UP_APPROX_ENABLED1){
	checkFraction = true;
	length = length1;
	if(checkFraction) fraction = pf_shel_check(length);
	ss_verbose = ss_verbose1; 
	//energy = 0.0;
	//structure = new int[length+1];
	//std::stack<base_pair> g_stack;
	PF_COUNT_MODE = PF_COUNT_MODE1;
	NO_DANGLE_MODE = NO_DANGLE_MODE1;
	PF_D2_UP_APPROX_ENABLED = PF_D2_UP_APPROX_ENABLED1;
	pf_d2.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED);
}

void StochasticTracebackD2::free_traceback(){
	pf_d2.free_partition();
	//delete[] structure;
}

MyDouble StochasticTracebackD2::randdouble()
{
	return MyDouble( rand()/(double(RAND_MAX)+1) );
}

bool StochasticTracebackD2::feasible(int i, int j)
{
	return j-i > TURN && canPair(RNA[i],RNA[j]);
}

//Probability calculation functions
MyDouble StochasticTracebackD2::U_0(int i, int j){
	return (MyDouble(1.0))/pf_d2.get_u(i,j);
}

MyDouble StochasticTracebackD2::U_ihj(int i, int h, int j){
	return (feasible(h,j) == true) ? (pf_d2.get_up(h,j)) * (pf_d2.myExp(-((pf_d2.ED5_new(h,j,h-1))+(pf_d2.ED3_new(h,j,j+1))+(pf_d2.auPenalty_new(h,j)))/RT)) / (pf_d2.get_u(i,j)) : MyDouble(0.0);
}

MyDouble StochasticTracebackD2::U_s1_ihj(int i, int h, int j){
	return (pf_d2.get_s1(h,j)) / (pf_d2.get_u(i,j));
}

MyDouble StochasticTracebackD2::S1_ihlj(int i, int h, int l, int j){
	return (feasible(h,l) == true) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.ED5_new(h,l,h-1))+(pf_d2.ED3_new(h,l,l+1))+(pf_d2.auPenalty_new(h,l)))/RT)) * (pf_d2.get_u(l+1,j)) / (pf_d2.get_s1(h,j)) : MyDouble(0.0);
}

MyDouble StochasticTracebackD2::Q_H_ij(int i, int j){
	return (pf_d2.myExp(-(pf_d2.eH_new(i,j))/RT)) / (pf_d2.get_up(i,j));
}

MyDouble StochasticTracebackD2::Q_S_ij(int i, int j){
	return (pf_d2.myExp(-(pf_d2.eS_new(i,j))/RT)) * (pf_d2.get_up(i+1,j-1)) / (pf_d2.get_up(i,j));
}

MyDouble StochasticTracebackD2::Q_M_ij(int i, int j){
	return (pf_d2.get_upm(i,j)) / (pf_d2.get_up(i,j));
}

MyDouble StochasticTracebackD2::Q_BI_ihlj(int i, int h, int l, int j){
	return feasible(h,l) ? (pf_d2.myExp(-1*(pf_d2.eL_new(i,j,h,l))/RT)) * (pf_d2.get_up(h,l)) / (pf_d2.get_up(i,j)) : MyDouble(0.0);
}

MyDouble StochasticTracebackD2::UPM_S2_ihj(int i, int h, int j){
	return (pf_d2.get_s2(h,j)) * (pf_d2.myExp(-((pf_d2.EA_new())+ 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(j,i,j-1)) + (pf_d2.ED3_new(j,i,i+1)))/RT)) / (pf_d2.get_upm(i,j)); //TODO: Old impl, using ed3(j,i) instead of ed3(i,j), similarly in ed5
	//return (pf_d2.get_s2(h,j)) * (pf_d2.myExp(-((pf_d2.EA_new())+ 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(i,j,j-1)) + (pf_d2.ED3_new(i,j,i+1)))/RT)) / (pf_d2.get_upm(i,j)); //TODO: New impl, using ed3(i,j) instead of ed3(j,i), similarly in ed5
	//return exp((-1)*ED3_new(j,i,i+1)/RT)* (s2[h][j] * exp((-1)*(EA_new()+2*EC_new()+(h-i-1)*EB_new())/RT))/upm[i][j];//TODO: New impl
}

MyDouble StochasticTracebackD2::S2_ihlj(int i, int h, int l, int j){
	return feasible(h,l) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)))/RT)) * (pf_d2.get_u1(l+1,j-1)) / pf_d2.get_s2(h,j) : MyDouble(0.0);
}

MyDouble StochasticTracebackD2::U1_s3_ihj(int i, int h, int j){
	return (pf_d2.get_s3(h,j)) * (pf_d2.myExp((-1)*((pf_d2.EC_new())+(h-i)*(pf_d2.EB_new()))/RT)) / (pf_d2.get_u1(i,j));
}

MyDouble StochasticTracebackD2::S3_ihlj(int i, int h, int l, int j){
	return feasible(h,l) ? (pf_d2.get_up(h,l)) * (pf_d2.myExp(-((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)))/RT)) * ( (pf_d2.myExp(-(j-l)*(pf_d2.EB_new())/RT)) * (pf_d2.f(j+1,h,l)) + (pf_d2.get_u1(l+1,j)) ) / (pf_d2.get_s3(h,j)) : MyDouble(0.0);
}

MyDouble StochasticTracebackD2::S3_MB_ihlj(int i, int h, int l, int j){
	MyDouble term1 =  (pf_d2.myExp(-(j-l)*(pf_d2.EB_new())/RT)) * (pf_d2.f(j+1,h,l));
	MyDouble term2 = (pf_d2.get_u1(l+1,j));
	return  term1 / (term1+term2);
}

/*
   double Q_ijBI(int i, int j)
   {
   double sum = 0;
   for (int h = i+1; h < j-1; ++h)
   for (int l = h+1; l < j; ++l)
   {
   if (h == i+1 && l == j-1) continue;
   sum += feasible(h,l)?(exp(-1*eL_new(i,j,h,l)/RT)*up[h][l]):0; 
   }
   return sum/up[i][j];
   }*/

//Functions related to sampling
void StochasticTracebackD2::rnd_u(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	cum_prob = cum_prob + U_0(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction){
			//printf("U_0(i, j)");
			//printf("\n");
			fraction.add(0, i, j, false);
		}
		return;
	}

	for (int h = i; h < j; ++h)
	{
		cum_prob = cum_prob + U_ihj(i,h,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("U_ihj(i, h, j)");
				//printf("\n");
				fraction.add(1, h, j, true);
				fraction.add(0, i, j, false);
			}
			double e2 = ( (pf_d2.ED5_new(h,j,h-1)) + (pf_d2.ED3_new(h,j,j+1)) + (pf_d2.auPenalty_new(h,j)) );
			if (ss_verbose == 1) {
				printf(" (pf_d2.ED5_new(h,j,h-1))=%f, (pf_d2.ED3_new(h,j,j+1))=%f, (pf_d2.auPenalty_new(h,j))=%f\n", (pf_d2.ED5_new(h,j,h-1))/100.0, (pf_d2.ED3_new(h,j,j+1))/100.0, (pf_d2.auPenalty_new(h,j))/100.0);
				printf(" U_ihj(i=%d,h=%d,j=%d)= %lf\n",i,h,j, e2/100.0);
			}
			energy += e2;
			base_pair bp(h,j,UP);
			//set_single_stranded(i,h-1,structure);
			g_stack.push(bp);
			
			return;
		}
	}

	int h1 = -1;
	for (int h = i;  h < j-1; ++h)
	{
		cum_prob = cum_prob + U_s1_ihj(i,h,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("U_s1_ihj(i, h, j)");
				//printf("\n");
				fraction.add(2, h, j, true);
				fraction.add(0, i, j, false);
			}
			h1 = h;
			rnd_s1(i,h1,j, structure, energy, g_stack);
			return;
		}
	}
	//printf("rnd=");rnd.print();printf(",cum_prob=");cum_prob.print();
	assert (0) ;
}

void StochasticTracebackD2::rnd_s1(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	for (int l = h+1; l < j; ++l)
	{
		cum_prob = cum_prob + S1_ihlj(i,h,l,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("S1_ihlj(i,h,l,j)");
				//printf("\n");
				fraction.add(1, h, l, true);
				fraction.add(0, l+1, j, true);
				fraction.add(2, h, j, false);
			}
			double e2 = (pf_d2.ED5_new(h,l,h-1))+ (pf_d2.auPenalty_new(h,l)) + (pf_d2.ED3_new(h,l,l+1));
			if (ss_verbose == 1) {
				printf("(pf_d2.ED5_new(h,l,h-1))=%f, (pf_d2.auPenalty_new(h,l))=%f, (pf_d2.ED3_new(h,l,l+1))=%f\n",(pf_d2.ED5_new(h,l,h-1)), (pf_d2.auPenalty_new(h,l)), (pf_d2.ED3_new(h,l,l+1)));
				printf("(%d %d) %lf\n",i,j, e2/100.0);
			}
			energy += e2;
			base_pair bp1(h,l,UP);
			base_pair bp2(l+1,j,U);
			g_stack.push(bp1);
			g_stack.push(bp2);
			return ;
		}
	}
	assert(0);
}

void StochasticTracebackD2::rnd_up(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	assert(structure[i] == 0);
	assert(structure[j] == 0);

	set_base_pair(i,j, structure);

	for (int h = i+1; h < j-1; ++h)
		for (int l = h+1; l < j; ++l)
		{
			if (h == i+1 && l == j-1) continue;
			cum_prob = cum_prob + Q_BI_ihlj(i,h,l,j);
			if (rnd < cum_prob)
			{
				if(checkFraction) {
					//printf("Q_BI_ihlj(i, h, l, j)");
					//printf("\n");
					fraction.add(1, h, l, true);
					fraction.add(1, i, j, false);
				}
				double e2 = (pf_d2.eL_new(i,j,h,l));
				if (ss_verbose == 1) 
					printf("IntLoop(%d %d) %lf\n",i,j, e2/100.0);
				energy += e2;
				base_pair bp(h,l,UP);
				g_stack.push(bp);
				
				return;
			}
		}

	cum_prob = cum_prob + Q_H_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction){
			//printf("Q_H_ij(i, j)");
			//printf("\n");
			fraction.add(1, i, j, false);
		}
		double e2 = (pf_d2.eH_new(i,j));
		if (ss_verbose == 1) 
			printf("Hairpin(%d %d) %lf\n",i,j, e2/100.0);
		energy += e2;
		//set_single_stranded(i+1,j-1,structure);
		
		return ;
	}

	cum_prob = cum_prob + Q_S_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("Q_S_ij(i,j)");
			//printf("\n");
			fraction.add(1, i+1, j-1, true);
			fraction.add(1, i, j, false);
		}
		double e2 = (pf_d2.eS_new(i,j));
		if (ss_verbose == 1) 
			printf("Stack(%d %d) %lf\n",i,j, e2/100.0);
		energy+=e2;
		base_pair bp(i+1,j-1,UP);
		g_stack.push(bp);
		
		return ;
	}

	cum_prob = cum_prob + Q_M_ij(i,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("Q_M_ij(i,j)");
			//printf("\n");
			fraction.add(3, i, j,true);
			fraction.add(1, i, j, false);
		}
		rnd_upm(i,j, structure, energy, g_stack);
		return;
	}

	for (int h = i+1; h < j-1; ++h)
		for (int l = h+1; l < j; ++l)
		{
			if (h == i+1 && l == j-1) continue;
			cum_prob = cum_prob + Q_BI_ihlj(i,h,l,j);
			if (rnd < cum_prob)
			{
				if(checkFraction) {
					//printf("Q_BI_ihlj(i,h,l,j)");
					//printf("\n");
					fraction.add(1, h, l, true);
					fraction.add(1, i, j, false);
				}
				double e2 = (pf_d2.eL_new(i,j,h,l));
				if (ss_verbose == 1) 
					printf("IntLoop(%d %d) %lf\n",i,j, e2/100.0);
				energy += e2;
				base_pair bp(h,l,UP);
				g_stack.push(bp);
				return;
			}
		}
	assert(0);
}

void StochasticTracebackD2::rnd_up_approximate(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	assert(structure[i] == 0);
	assert(structure[j] == 0);

	set_base_pair(i,j, structure);

	for (int p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1); ++p){
		int minq = j-i+p-MAXLOOP-2;
		if (minq < p+1+TURN) minq = p+1+TURN;
		int maxq = (p==(i+1))?(j-2):(j-1);
		for (int q = minq; q <=maxq ; ++q)
		{
			cum_prob = cum_prob + Q_BI_ihlj(i,p,q,j);
			if (rnd < cum_prob)
			{
				double e2 = (pf_d2.eL_new(i,j,p,q));
				if (ss_verbose == 1) 
					printf("IntLoop(%d %d) %lf\n",i,j, e2/100.0);
				energy += e2;
				base_pair bp(p,q,UP);
				g_stack.push(bp);
				return;
			}
		}
	}

	cum_prob = cum_prob + Q_H_ij(i,j);
	if (rnd < cum_prob)
	{
		double e2 = (pf_d2.eH_new(i,j));
		if (ss_verbose == 1) 
			printf("Hairpin(%d %d) %lf\n",i,j, e2/100.0);
		energy += e2;
		//set_single_stranded(i+1,j-1,structure);
		return ;
	}

	cum_prob = cum_prob + Q_S_ij(i,j);
	if (rnd < cum_prob)
	{
		double e2 = (pf_d2.eS_new(i,j));
		if (ss_verbose == 1) 
			printf("Stack(%d %d) %lf\n",i,j, e2/100.0);
		energy+=e2;
		base_pair bp(i+1,j-1,UP);
		g_stack.push(bp);
		return ;
	}

	cum_prob = cum_prob + Q_M_ij(i,j);
	if (rnd < cum_prob)
	{
		rnd_upm(i,j, structure, energy, g_stack);
		return;
	}

	assert(0);
}


void StochasticTracebackD2::rnd_u1(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);

	int h1 = -1;
	//for (int h = i+1; h < j-1; ++h)//TODO OLD 
	for (int h = i; h < j; ++h)//TODO NEW 
	{
		cum_prob = cum_prob + U1_s3_ihj(i,h,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("U1_s3_ihj(i,h,j)");
				//printf("\n");
				fraction.add(5, h, j, true);
				fraction.add(6, i, j, false);
			}
			double e2 = (pf_d2.EC_new()) + (h-i)*(pf_d2.EB_new());
			if (ss_verbose == 1){ 
				printf("(pf_d2.EC_new())=%f (h-i)*(pf_d2.EB_new())=%f\n",(pf_d2.EC_new())/100.0, (h-i)*(pf_d2.EB_new())/100.0);
				printf("U1_s3_ihj(%d %d %d) %lf\n",i,h,j, e2/100.0);
			}
			energy += e2;
			h1 = h;
			rnd_s3(i, h1, j, structure, energy, g_stack);
			return;
		}
	}
	assert(0);
}

void StochasticTracebackD2::rnd_s3(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){
	// sample l given h1 
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	for (int l = h+1; l <= j &&  l+1<=length ; ++l)
	{
		cum_prob = cum_prob +  S3_ihlj(i,h,l,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("S3_ihlj(i,h,l,j)");
				//printf("\n");
				fraction.add(1, h, l, true);
				fraction.add(6, l+1, j, true);
				fraction.add(5, h, j, false);
			}
			double e2 = ((pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1)));
			if (ss_verbose == 1) {
				printf("(pf_d2.auPenalty_new(h,l))=%f, (pf_d2.ED5_new(h,l,h-1))=%f, (pf_d2.ED3_new(h,l,l+1))=%f\n",(pf_d2.auPenalty_new(h,l))/100.0, (pf_d2.ED5_new(h,l,h-1))/100.0, (pf_d2.ED3_new(h,l,l+1))/100.0);
				printf("S3_ihlj(%d %d %d %d) %lf\n",i,h,l,j,e2/100.0);
			}
			energy += e2;
			base_pair bp(h,l,UP);
			g_stack.push(bp);
			rnd_s3_mb(i,h,l,j, structure, energy, g_stack);
			return;
		}
	}
	assert(0);
}

void StochasticTracebackD2::rnd_s3_mb(int i, int h, int l, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){//shel's document call this method with arguments i,h,l,j+1 therefore one will see difference of 1 in this code and shel's document
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	cum_prob = cum_prob +  S3_MB_ihlj(i,h,l,j);
	if (rnd < cum_prob)
	{
		if(checkFraction) {
			//printf("S3_MB_ihlj(i,h,l,j)");
			//printf("\n");
			fraction.add(6, l+1, j, false);
		}
		double tt =  0;//(j == l)? 0 : (pf_d2.ED3_new(h,l,l+1));//this term is corresponding to f(j+1,h,l)
		double e2 = tt + (j-l)*(pf_d2.EB_new());
		if (ss_verbose == 1){
			printf("j=%d,l=%d,tt=(j == l)?0:(pf_d2.ED3_new(h,l,l+1))=%f,(j-l)*(pf_d2.EB_new())=%f\n",j,l,tt/100.0,(j-l)*(pf_d2.EB_new())/100.0);				printf("S3_MB_ihlj(%d %d %d %d) %lf\n",i,h,l,j, e2/100.0);
		}
		energy += e2;
		
		return;
	}
	else{
		base_pair bp1(l+1,j,U1);
		g_stack.push(bp1);
		return;
	}
	assert(0);
}

void StochasticTracebackD2::rnd_upm(int i, int j, int* structure, double & energy, std::stack<base_pair>& g_stack)
{
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	if (ss_verbose == 1)
		printf("Multiloop (%d %d)\n",i,j);

	int h1 = -1;
	for (int h = i+1; h < j-1; ++h)
	{
		cum_prob = cum_prob + UPM_S2_ihj(i,h,j);
		if (rnd < cum_prob )
		{
			if(checkFraction) {
				//printf("UPM_S2_ihj(i,h,j)");
				//printf("\n");
				fraction.add(4, h, j, true);
				fraction.add(3, i, j, false);
			}
			double e2 = (pf_d2.EA_new()) + 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(j,i,j-1)) + (pf_d2.ED3_new(j,i,i+1));//TODO Old impl using ed3(j,i) instead of ed3(i,j)
			//double e2 = (pf_d2.EA_new()) + 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()) + (pf_d2.auPenalty_new(i,j)) + (pf_d2.ED5_new(i,j,j-1)) + (pf_d2.ED3_new(i,j,i+1));//TODO New impl using ed3(i,j( instead of ed3(j,i)
			energy += e2;
			h1 = h;
			if (ss_verbose == 1) {
				printf("(pf_d2.EA_new())=%f, (pf_d2.EC_new())=%f, (pf_d2.EB_new())=%f\n",(pf_d2.EA_new())/100.0, (pf_d2.EC_new())/100.0, (pf_d2.EB_new())/100.0);
				printf("(pf_d2.EA_new()) + 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new())=%f, (pf_d2.auPenalty_new(i,j))=%f, (pf_d2.ED5_new(j,i,j-1))=%f, (pf_d2.ED3_new(j,i,i+1))=%f\n",((pf_d2.EA_new()) + 2*(pf_d2.EC_new()) + (h-i-1)*(pf_d2.EB_new()))/100.0, (pf_d2.auPenalty_new(i,j))/100.0, (pf_d2.ED5_new(j,i,j-1))/100.0, (pf_d2.ED3_new(j,i,i+1))/100.0);
				printf("%s(%d %d %d) %lf\n", "UPM_S2_ihj",i,h1,j,e2/100.0);
			}
			rnd_s2(i,h1,j, structure, energy, g_stack);
			return;
		}
	}
	//assert(h1!=-1);
	assert(0);
}

void StochasticTracebackD2::rnd_s2(int i, int h, int j, int* structure, double & energy, std::stack<base_pair>& g_stack){
	MyDouble rnd = randdouble();
	MyDouble cum_prob(0.0);
	for (int l = h+1; l < j; ++l)
	{
		cum_prob = cum_prob + S2_ihlj(i,h,l,j);
		if (rnd < cum_prob)
		{
			if(checkFraction) {
				//printf("S2_ihlj(i,h,l,j)");
				//printf("\n");
				fraction.add(1, h, l, true);
				fraction.add(6, l+1, j-1, true);
				fraction.add(4, h, j, false);
			}
			double e2 = (pf_d2.auPenalty_new(h,l)) + (pf_d2.ED5_new(h,l,h-1)) + (pf_d2.ED3_new(h,l,l+1));       
			energy += e2;
			if (ss_verbose == 1){
				printf("(pf_d2.auPenalty_new(h,l))=%f, (pf_d2.ED5_new(h,l,h-1))=%f, (pf_d2.ED3_new(h,l,l+1))=%f\n",(pf_d2.auPenalty_new(h,l))/100.0, (pf_d2.ED5_new(h,l,h-1))/100.0, (pf_d2.ED3_new(h,l,l+1))/100.0);
				printf("%s(%d %d %d %d) %lf\n"," S2_ihlj",i,h,l,j, e2/100.0);
			}
			base_pair bp1(h,l,UP);
			base_pair bp2(l+1,j-1,U1);
			g_stack.push(bp1);
			g_stack.push(bp2);
			
			return;
		}
	}
	assert(0);
}

double StochasticTracebackD2::rnd_structure(int* structure)
{
	//printf("%lf %lf %lf\n", EA_new(), EB_new(), EC_new());
	srand(rand());
	//MyDouble U = pf_d2.get_u(1,len);
	base_pair first(1,length,U);
	std::stack<base_pair> g_stack;
	g_stack.push(first);
	double energy = 0.0;

	while (!g_stack.empty())
	{
		base_pair bp = g_stack.top();
		//   std::cout << bp;
		g_stack.pop();

		if (bp.type() == U)
			rnd_u(bp.i,bp.j, structure, energy, g_stack);
		else if (bp.type() == UP){
			if(pf_d2.PF_D2_UP_APPROX_ENABLED) rnd_up_approximate(bp.i,bp.j, structure, energy, g_stack);
			else rnd_up(bp.i,bp.j, structure, energy, g_stack);
		}
		else if (bp.type() == U1)
			rnd_u1(bp.i,bp.j, structure, energy, g_stack);
	}
	return (double)energy/100.0;
}

double StochasticTracebackD2::rnd_structure_parallel(int* structure, int threads_for_one_sample)
{
	//printf("%lf %lf %lf\n", EA_new(), EB_new(), EC_new());
	srand(rand());
	//MyDouble U = pf_d2.get_u(1,len);
	base_pair first(1,length,U);
	//std::stack<base_pair> g_stack;
	std::stack<base_pair> g_stack;
	g_stack.push(first);
	double energy = 0.0;
	std::stack<base_pair> g_stack_threads[threads_for_one_sample];
	double* energy_threads = new double[threads_for_one_sample];
	for(int index=0; index<threads_for_one_sample; ++index)energy_threads[index]=0.0;

	while (!g_stack.empty())
	{
		if(g_stack.size()%threads_for_one_sample != 0){
			base_pair bp = g_stack.top();
			//   std::cout << bp;
			g_stack.pop();

			if (bp.type() == U)
				rnd_u(bp.i,bp.j, structure, energy, g_stack);
			else if (bp.type() == UP){
				if(pf_d2.PF_D2_UP_APPROX_ENABLED) rnd_up_approximate(bp.i,bp.j, structure, energy, g_stack);
				else rnd_up(bp.i,bp.j, structure, energy, g_stack);
			}
			else if (bp.type() == U1)
				rnd_u1(bp.i,bp.j, structure, energy, g_stack);
		}
		else{
			std::deque<base_pair> g_deque;
			while(!g_stack.empty()){
				base_pair bp = g_stack.top();
				g_stack.pop();
				g_deque.push_back(bp);
			}
			int index;
			
			#ifdef _OPENMP
			#pragma omp parallel for private(index) shared(energy_threads, g_stack_threads, structure) schedule(guided) num_threads(threads_for_one_sample)
			//#pragma omp parallel for private(index) shared(energy_threads, g_stack_threads, structure) schedule(guided)
			#endif
			for (index = 0; index < g_deque.size(); ++index) {
				int thdId = omp_get_thread_num();
				base_pair bp = g_deque[index];
				if (bp.type() == U)
					rnd_u(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
				else if (bp.type() == UP){
					if(pf_d2.PF_D2_UP_APPROX_ENABLED) rnd_up_approximate(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
					else rnd_up(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
				}
				else if (bp.type() == U1)
					rnd_u1(bp.i,bp.j, structure, energy_threads[thdId], g_stack_threads[thdId]);
			}

			for(index=0; index<threads_for_one_sample; ++index){
				energy += energy_threads[index];
				energy_threads[index] = 0.0;
				while(!g_stack_threads[index].empty()){
					base_pair bp = g_stack_threads[index].top();
					g_stack_threads[index].pop();
					g_stack.push(bp);
				}
			}
		}
	}
	delete[] energy_threads;
	return (double)energy/100.0;
}


/*
   void batch_sample(int num_rnd, int length, double U)
   {
   int* structure = new int[length+1];
   srand(time(NULL));
   std::map<std::string,std::pair<int,double> >  uniq_structs;

   if (num_rnd > 0 ) {
   printf("\nSampling structures...\n");
   int count; //nsamples =0;
   for (count = 1; count <= num_rnd; ++count) 
   {
   memset(structure, 0, (length+1)*sizeof(int));
   double energy = rnd_structure(structure, length);

   std::string ensemble(length+1,'.');
   for (int i = 1; i <= (int)length; ++ i) {
   if (structure[i] > 0 && ensemble[i] == '.')
   {
   ensemble[i] = '(';
   ensemble[structure[i]] = ')';
   }
   }
///double myEnegry = -88.4;
//++nsamples;
//if (fabs(energy-myEnegry)>0.0001) continue; //TODO: debug
//++count;

std::map<std::string,std::pair<int,double> >::iterator iter ;
if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
{
std::pair<int,double>& pp = iter->second;
pp.first++;
}
else {
uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
}

// std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
}
//std::cout << nsamples << std::endl;
int pcount = 0;
int maxCount = 0; std::string bestStruct;
double bestE = INFINITY;

std::map<std::string,std::pair<int,double> >::iterator iter ;
for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
{
const std::string& ss = iter->first;
const std::pair<int,double>& pp = iter->second;
const double& estimated_p =  (double)pp.first/(double)num_rnd;
const double& energy = pp.second;
double actual_p = pow(2.718281,-1.0*energy/RT_)/U;

printf("%s %lf %lf %lf %d\n",ss.c_str(),energy,actual_p,estimated_p,pp.first);
pcount += pp.first;
if (pp.first > maxCount)
{
maxCount = pp.first;
bestStruct  = ss;
bestE = pp.second;
}
}
assert(num_rnd == pcount);
printf("\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);

}

delete [] structure;
}
 */

void StochasticTracebackD2::batch_sample(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION ,bool ST_D2_ENABLE_UNIFORM_SAMPLE, double ST_D2_UNIFORM_SAMPLE_ENERGY)
{cout<<"ST_D2_ENABLE_UNIFORM_SAMPLE="<<ST_D2_ENABLE_UNIFORM_SAMPLE<<",ST_D2_UNIFORM_SAMPLE_ENERGY="<<ST_D2_UNIFORM_SAMPLE_ENERGY<<endl;
	MyDouble U;
	/*if(PF_D2_UP_APPROX_ENABLED){
	  double t1 = get_seconds();
	  PartitionFunctionD2 pf_d2_exact_up;
	  bool PF_D2_UP_APPROX_ENABLED2 = false;
	  pf_d2_exact_up.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED2);
	  U = pf_d2_exact_up.get_u(1,length);
	  pf_d2_exact_up.free_partition();
	  t1 = get_seconds() - t1;
	  printf("D2 Exact UP partition function computation running time: %9.6f seconds\n", t1);
	  }
	  else U = pf_d2.get_u(1,length);*/
	U = pf_d2.get_u(1,length);

	int threads_for_one_sample = 1;
	#ifdef _OPENMP
	#pragma omp parallel
	//#pragma omp master
	{
		int thdId1 = omp_get_thread_num();
		if(thdId1==0){
			if(g_nthreads < 0) threads_for_one_sample = omp_get_num_threads();//TODO move this line to mfe_main.cc
			else threads_for_one_sample = g_nthreads;
		}
	}
	#endif
	if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION) fprintf(stdout,"Stochastic Traceback: Thread count for one sample parallelization: %3d \n",threads_for_one_sample);

	std::map<std::string,std::pair<int,double> >  uniq_structs;
	int* structure = new int[length+1];

	if (num_rnd > 0 ) {
		printf("\nSampling structures...\n");
		int count, nsamples =0;
		for (count = 1; count <= num_rnd; ++count) 
		{
			nsamples++;
			memset(structure, 0, (length+1)*sizeof(int));
			double energy;
			if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
				energy = rnd_structure_parallel(structure, threads_for_one_sample);
			}
			else{
				energy = rnd_structure(structure);
			}

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			
			//Below line of codes is for finding samples with particular energy
			if(ST_D2_ENABLE_UNIFORM_SAMPLE){
				//double myEnegry = -92.1;//-91.3;//-94.8;//dS=-88.4;//d2=-93.1
				if (fabs(energy-ST_D2_UNIFORM_SAMPLE_ENERGY)>0.0001){ count--;continue;} //TODO: debug
			}

			std::map<std::string,std::pair<int,double> >::iterator iter ;
			if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
			{
				std::pair<int,double>& pp = iter->second;
				pp.first++;
				assert(energy==pp.second);
			}
			else {
				if(ST_D2_ENABLE_SCATTER_PLOT) uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
			}

			if(!ST_D2_ENABLE_SCATTER_PLOT) std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
		}
		//std::cout << nsamples << std::endl;
		if(ST_D2_ENABLE_SCATTER_PLOT){
			int pcount = 0;
			int maxCount = 0; std::string bestStruct;
			double bestE = INFINITY;
			printf("nsamples=%d\n",nsamples);
			printf("%s,%s,%s","structure","energy","boltzman_probability");
			printf(",%s,%s\n","estimated_probability","frequency");
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
			{
				const std::string& ss = iter->first;
				const std::pair<int,double>& pp = iter->second;
				const double& estimated_p =  (double)pp.first/(double)num_rnd;
				const double& energy = pp.second;
				//MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy/RT_)))/U;
				//MyDouble actual_p = (pf_d2.myExp(-(energy)/(RT_)))/U;
				MyDouble actual_p = (pf_d2.myExp(-(energy*100)/(RT)))/U;
				//MyDouble actual_p(-(energy)/(RT_));///U;
				printf("%s,%f,",ss.c_str(),energy);actual_p.print();
				printf(",%f,%d\n",estimated_p,pp.first);

				//printf("%s %lf\n",ss.c_str(),energy);actual_p.print();
				//printf("%lf %d\n",estimated_p,pp.first);
				pcount += pp.first;
				if (pp.first > maxCount)
				{
					maxCount = pp.first;
					bestStruct  = ss;
					bestE = pp.second;
				}
			}
			assert(num_rnd == pcount);
			printf("\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);
		}
		else{
			printf("nsamples=%d\n",nsamples);
		}

	}
	delete[] structure;
}

void StochasticTracebackD2::batch_sample_parallel(int num_rnd, bool ST_D2_ENABLE_SCATTER_PLOT, bool ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION)
{
	//MyDouble U = pf_d2.get_u(1,length);
	MyDouble U;
	/*if(PF_D2_UP_APPROX_ENABLED){
	  double t1 = get_seconds();
	  PartitionFunctionD2 pf_d2_exact_up;
	  bool PF_D2_UP_APPROX_ENABLED2 = false;
	  pf_d2_exact_up.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED2);
	  U = pf_d2_exact_up.get_u(1,length);
	  pf_d2_exact_up.free_partition();
	  t1 = get_seconds() - t1;
	  printf("D2 Exact UP partition function computation running time: %9.6f seconds\n", t1);
	  }
	  else U = pf_d2.get_u(1,length);
	 */
	U = pf_d2.get_u(1,length);

	srand(time(NULL));
	/*
	//OPTIMIZED CODE STARTS
	#ifdef _OPENMP
	if (g_nthreads > 0) omp_set_num_threads(g_nthreads);
	#endif

	#ifdef _OPENMP
	#pragma omp parallel
	//#pragma omp master
	{
	int thdId1 = omp_get_thread_num();
	if(thdId1==0){
	fprintf(stdout,"Stochastic Traceback: Thread count: %3d \n",omp_get_num_threads());
	if(g_nthreads < 0) g_nthreads = omp_get_num_threads();//TODO move this line to mfe_main.cc
	}
	}
	#endif
	//OPTIMIZED CODE ENDSS
	 */
	int total_used_threads = 1;
	#ifdef _OPENMP
	#pragma omp parallel
	//#pragma omp master
	{
		int thdId1 = omp_get_thread_num();
		if(thdId1==0){
			if(g_nthreads < 0) total_used_threads = omp_get_num_threads();//TODO move this line to mfe_main.cc
			else total_used_threads = g_nthreads;
		}
	}
	#endif
	int threads_for_counts = total_used_threads;
	int threads_for_one_sample = 1;

	if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
		int sample_count_per_thread = num_rnd/total_used_threads;
		if(total_used_threads%2 !=0 || total_used_threads < 4 || length < 200 || sample_count_per_thread > 40){
			threads_for_counts = total_used_threads;
			threads_for_one_sample = 1;
		}
		else if(total_used_threads%4 !=0 || total_used_threads < 8 || length < 600 || sample_count_per_thread > 20 ){
			threads_for_counts = total_used_threads/2;
			threads_for_one_sample = 2;
		}
		else if(total_used_threads%8 !=0 || total_used_threads < 16 || length < 1500 || sample_count_per_thread > 10){
			threads_for_counts = total_used_threads/4;
			threads_for_one_sample = 4;
		}
		else if(total_used_threads%16 !=0 || total_used_threads < 32 || length < 3000 || sample_count_per_thread > 5 ){
			threads_for_counts = total_used_threads/8;
			threads_for_one_sample = 8;
		}
		else if(total_used_threads%32 !=0 || total_used_threads < 64 || length < 6000 || sample_count_per_thread > 2 ){
			threads_for_counts = total_used_threads/16;
			threads_for_one_sample = 16;
		}
		else{
			threads_for_counts = total_used_threads/32;
			threads_for_one_sample = 32;
		}

		fprintf(stdout,"Stochastic Traceback: Thread count for one sample parallelization: %3d \n",threads_for_one_sample);
	}
	fprintf(stdout,"Stochastic Traceback: Thread count for counts parallelization: %3d \n",threads_for_counts);



	std::map<std::string,std::pair<int,double> >  uniq_structs;
	//std::map<std::string,std::pair<int,double> >  uniq_structs_thread[g_nthreads];
	//g_nthreads=4;//TODO remove this line
	//cout<<"Manoj after: g_nthreads="<<g_nthreads<<endl;
	std::map<std::string,std::pair<int,double> > *  uniq_structs_thread = new std::map<std::string,std::pair<int,double> >[threads_for_counts];
	int* structures_thread = new int[threads_for_counts*(length+1)];

	if (num_rnd > 0 ) {
		printf("\nSampling structures...\n");
		int count;// nsamples =0;
		#ifdef _OPENMP
		#pragma omp parallel for private (count) shared(structures_thread) schedule(guided) num_threads(threads_for_counts)
		#endif
		for (count = 1; count <= num_rnd; ++count) 
		{
			//nsamples++;
			int thdId = omp_get_thread_num();
			//cout<<"thdId="<<thdId<<endl;
			int* structure = structures_thread + thdId*(length+1);
			memset(structure, 0, (length+1)*sizeof(int));
			//double energy = rnd_structure(structure);
			double energy;
			if(ST_D2_ENABLE_ONE_SAMPLE_PARALLELIZATION){
				energy = rnd_structure_parallel(structure, threads_for_one_sample);
			}
			else{
				energy = rnd_structure(structure);
			}

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			/*
			//Below line of codes is for finding samples with particular energy
			double myEnegry = -92.1;//-91.3;//-94.8;//dS=-88.4;//d2=-93.1
			if (fabs(energy-myEnegry)>0.0001){ count--;continue;} //TODO: debug
			 */

			std::map<std::string,std::pair<int,double> >::iterator iter ;
			if ((iter =uniq_structs_thread[thdId].find(ensemble.substr(1))) != uniq_structs_thread[thdId].end())
			{
				std::pair<int,double>& pp = iter->second;
				pp.first++;
				//cout<<"energy="<<energy<<",pp.second="<<pp.second<<endl;
				assert(energy==pp.second);
			}
			else {
				if(ST_D2_ENABLE_SCATTER_PLOT){
					std::pair< std::string, std::pair<int,double> > new_pp = make_pair(ensemble.substr(1),std::pair<int,double>(1,energy));
					uniq_structs_thread[thdId].insert(new_pp); 
				}
				//uniq_structs_thread[thdId].insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
			}

			if(!ST_D2_ENABLE_SCATTER_PLOT) std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
		}

		if(ST_D2_ENABLE_SCATTER_PLOT){
			for(int thd_id=0; thd_id<threads_for_counts; thd_id++){
				std::map<std::string,std::pair<int,double> >::iterator thd_iter ;
				for (thd_iter = uniq_structs_thread[thd_id].begin(); thd_iter != uniq_structs_thread[thd_id].end();  ++thd_iter)
				{
					const std::string& thd_ss = thd_iter->first;
					const std::pair<int,double>& thd_pp = thd_iter->second;
					const int thd_freq = thd_pp.first;
					const double thd_e = thd_pp.second;

					std::map<std::string,std::pair<int,double> >::iterator iter ;
					if ((iter =uniq_structs.find(thd_ss)) != uniq_structs.end())
					{
						std::pair<int,double>& pp = iter->second;
						(pp.first)+=thd_freq;
						const double e = pp.second;
						//cout<<"e="<<e<<",thd_e="<<thd_e<<endl;
						assert(e==thd_e);
					}
					else {
						uniq_structs.insert(make_pair(thd_ss,std::pair<int,double>(thd_freq,thd_e)));
					}
				}
			}
			//std::cout << nsamples << std::endl;
			int pcount = 0;
			int maxCount = 0; std::string bestStruct;
			double bestE = INFINITY;
			printf("nsamples=%d\n",num_rnd);
			printf("%s,%s,%s","structure","energy","boltzman_probability");
			printf(",%s,%s\n","estimated_probability","frequency");
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
			{
				const std::string& ss = iter->first;
				const std::pair<int,double>& pp = iter->second;
				const double& estimated_p =  (double)pp.first/(double)num_rnd;
				const double& energy = pp.second;
				//MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy/RT_)))/U;
				//MyDouble actual_p = (pf_d2.myExp(-(energy)/(RT_)))/U;
				MyDouble actual_p;
				actual_p = (pf_d2.myExp(-(energy*100)/(RT)))/U;
				//MyDouble actual_p(-(energy)/(RT_));///U;
				printf("%s,%f,",ss.c_str(),energy);actual_p.print();
				printf(",%f,%d\n",estimated_p,pp.first);

				//printf("%s %lf\n",ss.c_str(),energy);actual_p.print();
				//printf("%lf %d\n",estimated_p,pp.first);
				pcount += pp.first;
				if (pp.first > maxCount)
				{
					maxCount = pp.first;
					bestStruct  = ss;
					bestE = pp.second;
				}
			}
			assert(num_rnd == pcount);
			printf("\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);
		}
		else{
			printf("nsamples=%d\n",num_rnd);
		}
	}
	delete [] structures_thread;
	delete [] uniq_structs_thread;
}


void StochasticTracebackD2::batch_sample_and_dump(int num_rnd, std::string ctFileDumpDir, std::string stochastic_summery_file_name, std::string seq, std::string seqfile)
{
	//MyDouble U = pf_d2.get_u(1,length);
	 MyDouble U;
        /*if(PF_D2_UP_APPROX_ENABLED){
		double t1 = get_seconds();
                PartitionFunctionD2 pf_d2_exact_up;
                bool PF_D2_UP_APPROX_ENABLED2 = false;
                pf_d2_exact_up.calculate_partition(length,PF_COUNT_MODE,NO_DANGLE_MODE, PF_D2_UP_APPROX_ENABLED2);
                U = pf_d2_exact_up.get_u(1,length);
                pf_d2_exact_up.free_partition();
                t1 = get_seconds() - t1;
                printf("D2 Exact UP partition function computation running time: %9.6f seconds\n", t1);
       }
        else U = pf_d2.get_u(1,length);
	*/
         U = pf_d2.get_u(1,length);
	//data dump preparation code starts here
	if(ctFileDumpDir.compare("")==0){
		char abspath[1000];
		getcwd(abspath, 1000);
		ctFileDumpDir = abspath;
	}
	cout<<"Using ctFileDumpDir = "<<ctFileDumpDir<<endl;
	std::stringstream ss;
	ss<<ctFileDumpDir<<"/"<<stochastic_summery_file_name;
	stochastic_summery_file_name = ss.str();
	cout<<"Using stochastic_summary_file_name = "<<stochastic_summery_file_name<<endl;
	std::ofstream summaryoutfile;
	summaryoutfile.open(stochastic_summery_file_name.c_str());
	std::string seqname = seqfile.substr(seqfile.find_last_of("/\\") + 1, seqfile.length() - 1);
	cout<<"Sequence Name = "<<seqname<<endl;
	//data dump preparation code ends here

	srand(time(NULL));
	std::map<std::string,std::pair<int,double> >  uniq_structs;
	int* structure = new int[length+1];
	if (num_rnd > 0 ) {
		printf("\nSampling structures...\n");
		int count, nsamples =0;
		for (count = 1; count <= num_rnd; ++count) 
		{
			memset(structure, 0, (length+1)*sizeof(int));
			double energy = rnd_structure(structure);

			std::string ensemble(length+1,'.');
			for (int i = 1; i <= (int)length; ++ i) {
				if (structure[i] > 0 && ensemble[i] == '.')
				{
					ensemble[i] = '(';
					ensemble[structure[i]] = ')';
				}
			}
			/*
			//below code is for uniform sampling test
			double myEnegry = -88.4;
			++nsamples;
			if(nsamples>1000)break;
			if (fabs(energy-myEnegry)>0.0001) continue; //TODO: debug
			//++count;
			*/
			std::map<std::string,std::pair<int,double> >::iterator iter ;
			if ((iter =uniq_structs.find(ensemble.substr(1))) != uniq_structs.end())
			{
				std::pair<int,double>& pp = iter->second;
				pp.first++;
			}
			else {
				uniq_structs.insert(make_pair(ensemble.substr(1),std::pair<int,double>(1,energy))); 
			}

			// std::cout << ensemble.substr(1) << ' ' << energy << std::endl;
			//data dump code starts here
			std::stringstream ss;
			ss<<ctFileDumpDir<<"/"<<seqname<<"_"<<count<<".ct";
			save_ct_file(ss.str(), seq, energy, structure);
			summaryoutfile<<ss.str()<<" "<<ensemble.substr(1)<<" "<<energy<< std::endl;
			//data dump code ends here
		}
		//std::cout << nsamples << std::endl;
		int pcount = 0;
		int maxCount = 0; std::string bestStruct;
		double bestE = INFINITY;

		printf("%s,%s,%s","structure","energy","boltzman_probability");
		printf(",%s,%s\n","estimated_probability","frequency");
		std::map<std::string,std::pair<int,double> >::iterator iter ;
		for (iter = uniq_structs.begin(); iter != uniq_structs.end();  ++iter)
		{
			const std::string& ss = iter->first;
			const std::pair<int,double>& pp = iter->second;
			const double& estimated_p =  (double)pp.first/(double)num_rnd;
			const double& energy = pp.second;
	
			//MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy/RT_)))/U;
			MyDouble actual_p = (MyDouble(pow(2.718281,-1.0*energy*100/RT)))/U;
			printf("%s,%lf,",ss.c_str(),energy);actual_p.print();
			printf(",%lf,%d\n",estimated_p,pp.first);
	
			pcount += pp.first;
			if (pp.first > maxCount)
			{
				maxCount = pp.first;
				bestStruct  = ss;
				bestE = pp.second;
			}
		}
		assert(num_rnd == pcount);
		printf("\nMax frequency structure : \n%s e=%lf freq=%d p=%lf\n",bestStruct.c_str(),bestE,maxCount,(double)maxCount/(double)num_rnd);
	}
	summaryoutfile.close();
	delete[] structure;

}

void StochasticTracebackD2::set_single_stranded(int i, int j, int* structure)
{
	for(;i<=j;++i) 
		structure[i] = 0;
}

void StochasticTracebackD2::set_base_pair(int i, int j, int* structure)
{
	bool cond = j-i > TURN && canPair(RNA[i],RNA[j]);
	assert(cond);
	structure[i] = j;
	structure[j] = i;
}
