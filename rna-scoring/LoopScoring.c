#include <stdio.h>
#include <math.h>
#include "LoopScoring.h"
#include "TreeScoring.h"

char bases[16] = {0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U',
                      'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
int printOn1=1;
int eMUnpairedRegion(int i1, int j1, int i2, int j2, int* RNA, nndb_constants* param){//printf("Entering eMUnpairedRegion\n");
	//Shel: helper function for calculating muliloop energies. 
	//Computes dangling energy for unpaired region given the pairs for the stems on either side.
	//(based on Andronescu masters thesis,
	// M.Sc., Academy of Economic Studies, Bucharest, Romania, 2000, 
	// pg 32.)
	int energy;
	/*
	int PFMODE=0;//boltzman sampling or stochastic sampling and  partition function mode or dS mode, it is mode as defined and used for partition function of sfold
	int NODANGLEMODE=0;//no dangling at all means d0
	int D2MODE=0;//d2 mode
	int DEFAULTMODE=1;//default mode
	*/

	if(NODANGLEMODE==1){
		energy = 0;
	}
	else if (j1+1 < i2-1 || D2MODE==1) {
		// if there are at least two nucleotides in the unpaired region,
		// then add the energy for both a 3' and 5' dangling end 
		// for the first and last nucleotides in the unpaired region, respectively.
		//printf("A\n");
		if(D2MODE==1){
			energy = param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0] + param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1];
			if(printOn1)printf("i1=%d,j1=%d,i2=%d,j2=%d, param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0]/100=%f,  param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1]/100=%f\n", i1,j1,i2,j2,param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0]/100.0, param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1]/100.0);
		}
	} else if (j1+1 == i2-1) { //ZS: fixed limits 
		// if there is only one nucleotide in the unpaired region,
		// then add which energy is more favorable, 
		// the unpaired region as a 3' or a 5' dangling end.
		//printf("B\n");
		//energy = MIN(param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0], param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1]);
		//int e5 = param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0];
		//int e3 = param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1];
		//printf("In eMUnpairedRegion, e5=%d, e3=%d\n",e5,e3);
		if(PFMODE==1){
			energy = param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0];
		}
		else if(DEFAULTMODE==1){
			energy = MIN(param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0], param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1]);
		}
		else energy=0;
	} else {
		// if the unpaire region is empty, there can be no dangling ends.
		//	printf("No dangling\n");
		energy = 0;
	}
	//printf("between branch %d - %d and %d - %d, \n 3dangle has energy %d, and \n 5dangle energy %d\n",  
	//	   i1, j1, i2, j2, param->dangle[RNA[j1]][RNA[i1]][RNA[j1+1]][0],  param->dangle[RNA[j2]][RNA[i2]][RNA[i2-1]][1]);
	//printf("Nucleotides were: %d-%d and %d-%d \n", RNA[i1], RNA[j1], RNA[i2], RNA[j2]);
	//printf("Returning energy = %i \n", energy);
        if(printOn1)printf("eMUnpairedRegion(i1=%d,j1=%d,i2=%d,j2=%d)=%f\n", i1, j1, i2, j2, (double)energy/100);
	return energy;
}
int eL(int i, int j, int ip, int jp, int* RNA, nndb_constants* param) {
    //ZS: internal loop calculations, borrowed from GTfold, with small 
    //modifications (eg. deleted eparam's which nobody understood and 
    //they were zero to our knowledge anyway) 
	int energy;
	int size1, size2, size;
	int loginc;
	int lopsided; /* define the asymmetry of an interior loop */

	energy = INFINITY_;
	loginc = 0;

	size1 = ip - i - 1;
	size2 = j - jp - 1;
	size = size1 + size2;

	if (size1 == 0 || size2 == 0) {
		if (size > 30) {
			loginc = (int) floor(param->prelog * log((double) size / 30.0));
			//(Note: auPen is defined in LoopScoring.h and returns 0 if no AU 
			//penalty should be applied, and 1 if it should be applied. This is 
			//different from how it is in GTfold.)	
			energy = param->bulge[30] + loginc + 
						param->auend*(auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]));	
		} else if (size <= 30 && size != 1) {
			energy = param->bulge[size];
			energy += param->auend*(auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]));
		} else if (size == 1) {
			energy = param->stack[fourBaseIndex(RNA[i], RNA[j], RNA[ip], RNA[jp])]
			               + param->bulge[size];
		}
	} else {
		//Internal loop
		lopsided = abs(size1 - size2);
		if (size > 30) {
			loginc = (int) floor(param->prelog * log((double) size / 30.0));
			//ZS: Somebody else's previous comment follows, I don't know the answer: 
			/* Please check what should be the difference in the following two options. Is it correct?*/
			if (!((size1 == 1 || size2 == 1) && param->gail)) { 
				//ZS: gail = grossly assymetric interior loop rule, ON(1)/OFF(0)
				energy = param->tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])]
						 	+ param->tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], 
							 RNA[ip - 1])] + param->inter[30] + loginc + param->eparam[3]   //ZS: I have no idea about eparam here 
							  + MIN(param->maxpen, (lopsided * param->poppen[MIN(2, MIN(size1, size2))]));
			} else {
				energy = param->tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)]
				        + param->tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A,
				        		BASE_A)] + param->inter[30] + loginc + param->eparam[3] + 
								  MIN(param->maxpen, (lopsided	* param->poppen[MIN(2, MIN(size1, size2))]));
			}
		}
		/* if size is not > 30, we have a lot of cases... */
		else if (size1 == 2 && size2 == 2) {
			/* 2x2 internal loop */
			energy = param->iloop22[RNA[i]][RNA[ip]][RNA[j]][RNA[jp]][RNA[i + 1]][RNA[i+ 2]][RNA[j - 1]][RNA[j - 2]];
		} else if (size1 == 1 && size2 == 2) {
			energy = param->iloop21[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]][RNA[j - 2]][RNA[ip]][RNA[jp]];
		} else if (size1 == 2 && size2 == 1) {
			/* 1x2 internal loop */
			energy = param->iloop21[RNA[jp]][RNA[ip]][RNA[j - 1]][RNA[i + 2]][RNA[i + 1]][RNA[j]][RNA[i]];
		} else if (size == 2) {
			/* 1*1 internal loops */
			energy = param->iloop11[RNA[i]][RNA[i + 1]][RNA[ip]][RNA[j]][RNA[j - 1]][RNA[jp]];
		} else if ((size1 == 1 || size2 == 1) && param->gail) { 
			energy = param->tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)]
			         + param->tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A, BASE_A)]
                  + param->inter[size] + loginc + param->eparam[3] + MIN(param->maxpen, 
						(lopsided * param->poppen[MIN(2, MIN(size1, size2))]));
		} else { /* General Internal loops */
			energy = param-> tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] 
					   + param->tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], 
						RNA[ip - 1])] + param->inter[size] + loginc + param->eparam[3] 
					/* AM: I don't understand this eparam value, 
					I think they do not play any role currently. 
					Please look in loader.cc file, for what value 
					have been assinged to various elements of eparam array */
					//ZS: I don't understand it either and I think we should just 
					//delete it 
               + MIN(param->maxpen, (lopsided * param->poppen[MIN(2, MIN(size1, size2))]));
		}
	}
        if(printOn1)printf("eL(i=%d, j=%d, ip=%d, jp=%d)=%f\n",i,j,ip,jp,(double)energy/100);
	return energy;
}

int eH(int i, int j, int* RNA, nndb_constants* param) {
	/*  Hairpin loop for all the bases between i and j */
	/*  size for size of the loop, energy is the result, 
	loginc is for the extrapolation for loops bigger than 30 */
	int size;
	int loginc;
	int energy = INFINITY_;
	int key, index, count, tlink, kmult;
	size = j - i - 1; /*  size is the number of bases in the loop, when the closing pair is excluded */

	/*  look in hairpin, and be careful that there is only 30 values */
	if (size > 30) {
		loginc = (int) (param->prelog * log(((double) size) / 30.0));
		energy = param->hairpin[30] + loginc + param->tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + param->eparam[4]; /* size penalty + terminal mismatch stacking energy*/
	}

	else if (size <= 30 && size > 4) {
		energy = param->hairpin[size] + param->tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + param->eparam[4]; /* size penalty + terminal mismatch stacking energy*/
	}

	else if (size == 4) {
		/*  tetraloop */
		key = 0;
		tlink = 0;
		int cnt2; 
		for (index = 0; index < 6; ++index) {
			switch (RNA[i + index]) {
			case BASE_A:
				kmult = 1;
				break;
			case BASE_C:
				kmult = 2;
				break;
			case BASE_G:
				kmult = 3;
				break;
			case BASE_U:
				kmult = 4;
				break;
			default:
				kmult = 0;
				fprintf(stderr, "ERROR: in tetraloop calculation\n");
			}
			//ZS: The math.pow function didn't work for some reason on my machine.
			//It is also silly to convert to doubles when it isn't necessary.
			//So I just made this "quick and dirty". fix A more ideal solution 
			//would use bit operations for the "key", but then it has to be 
			//changed both here and in the loader. 
		   int powval=1;
			for(cnt2=5-index; cnt2>0; cnt2--){
					powval *=10;
			}
			
			//key += kmult * (int) pow(10.0, 5 - index);  //ZS: This didn't work
			key += kmult * powval; 
			
		}
		/*  if the sequence is in tloop, we use this value */
		for (count = 1; count < param->numoftloops && tlink == 0; ++count) {
			if (key == param->tloop[count][0]) {
				tlink = param->tloop[count][1];
			}
		}
		energy = tlink + param->hairpin[size] + param->tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + param->eparam[4];
	   
	}

	else if (size == 3) {
		/*  triloop... For the moment, the file triloop.dat is empty */
		/*  else, should have a treatment like the one if size==4 */
		energy = param->hairpin[size];
		/* AM: Don't include stacking energy terms for triloopls */
		/* + tstkh[RNA[i]][RNA[j]][RNA[i+1]][RNA[j-1]]  */
		/* + eparam[4]; */
		/*  Must be another penalty for terminal AU... Not sure of this */
		energy += param->auend*auPen(RNA[i], RNA[j]);
	}

	else if (size < 3 && size != 0) {
		/*  no terminal mismatch */
		energy = param->hairpin[size] + param->eparam[4];
		if ((RNA[i] == BASE_A && RNA[j] == BASE_U) || (RNA[i] == BASE_U
				&& RNA[j] == BASE_A)) {
			energy += 6; /*  Seems to be a penalty for terminal AU.  *//* Hairpin Loops of size 3 are not allowed, the term hairpin[size] will result in a very large value.  */
		}
	} else if (size == 0)
		return INFINITY_;

	/*  GGG Bonus => GU closure preceded by GG */
	/*  i-2 = i-1 = i = G, and j = U; i < j */
	if (i > 2) {
		if (RNA[i - 2] == BASE_G && RNA[i - 1] == BASE_G && RNA[i] == BASE_G
				&& RNA[j] == BASE_U) {
			energy += param->gubonus;
			/*  printf ("\n GGG bonus for i %d j %d ", i, j); */
		}
	}

	/*  Poly-C loop => How many C are needed for being a poly-C loop */
	tlink = 1;
	for (index = 1; (index <= size) && (tlink == 1); ++index) {
		if (RNA[i + index] != BASE_C)
			tlink = 0;
	}
	if (tlink == 1) {
		if (size == 3) {
			energy += param->c3;
		} else {
			energy += param->cint + size * param->cslope;
		}
	}
	if(printOn1)printf("eH(i=%d, j=%d)=%f\n",i,j,(double)energy/100); 
	return energy;
}

int eS(int i, int j, int* RNA, nndb_constants* param) {
	//ZS: Score a stack. 
   int energy;
	energy = param->stack[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])]; 
	if(printOn1)printf("eS(i=%d, j=%d, i+1=%d, j-1=%d)=%f\n",i,j,i+1,j-1, (double)energy/100);
	return energy;
}


int _eM(TreeNode* node, int* pairedChildren, int numPairedChildren, int* RNA, nndb_constants* param) {
    //ZS: Score a multiloop 
    //Here, no dangling energies are taken into account, we are using the formulas from the 
    //Turner website: for fewer than 7 unpaired nucleotides: 
    //eM(i,j, i1, j1, ... , ik, jk) = a + b*<nr of unpaired nt> + c*<nr of branches>
	 //for more than 6 unpaired nucleotides: 
	 //eM = a + 6b + (1.1 kcal/mol)×ln([number of unpaired nucleotides]/6) 
	//+ c×[number of branching helices] 
    
	int energy;
	int a = param->multConst[0];
	int b = param->multConst[1];
	int c = param->multConst[2];
	int nr_branches = numPairedChildren + 1;
	int nr_unpaired = node->numChildren - numPairedChildren;
	
	//	printf("nr branches = %i, nr unpaired = %i, a = %i, b = %i, c = %i \n",  
	//	   nr_branches, nr_unpaired, a, b, c);
	
	if(nr_unpaired <= 6){
		energy = a + nr_unpaired*b + nr_branches*c;
	}
	else{
	   energy = a + 6*b + 110*log(nr_unpaired/6) + c*nr_branches;
   }
	if(printOn1)printf("_eM(i=%d,j=%d)=%f\n",(node->lowBase).index, node->highBase.index, (double)energy/100);
	return energy;
}

int eM(TreeNode* node, int* pairedChildren, int numPairedChildren, int* RNA, nndb_constants* param) {
    //Shel: Score a multiloop 
    //Here the dangling bases are also taken into account 
	int energy;
	int a = param->multConst[0];
	int b = param->multConst[1];
	int c = param->multConst[2];
	int nr_branches = numPairedChildren + 1;
	int nr_unpaired = node->numChildren - numPairedChildren;
	energy = a + nr_unpaired*b + nr_branches*c;
	
	if(printOn1)printf("nr branches = %i, nr unpaired = %i, a/100 = %f, b/100 = %f, c/100 = %f, ENERGY=a + nr_unpaired*b + nr_branches*c /100.0 = %f \n",  
		   nr_branches, nr_unpaired, a/100.0, b/100.0, c/100.0, energy/100.0);
	
	energy += eMUnpairedRegion(node->highBase.index, node->lowBase.index, 
							   node->children[pairedChildren[0]]->lowBase.index, node->children[pairedChildren[0]]->highBase.index, 
							   RNA, param);
	int i;
	for (i = 0; i < numPairedChildren-1; i++) {
		//Scores the dangling ends in unpaired regions between paired children
		energy += eMUnpairedRegion(node->children[pairedChildren[i]]->lowBase.index, node->children[pairedChildren[i]]->highBase.index, 
								   node->children[pairedChildren[i+1]]->lowBase.index, node->children[pairedChildren[i+1]]->highBase.index, 
								   RNA, param);	
		//printf("Energy so far: %i \n", energy);
	}
	energy += eMUnpairedRegion(node->children[pairedChildren[numPairedChildren-1]]->lowBase.index, node->children[pairedChildren[numPairedChildren-1]]->highBase.index, 
									  node->highBase.index, node->lowBase.index,
									  RNA, param);
	//printf("Energy so far: %i \n", energy);
	
	i=0;
	for(i=0; i < numPairedChildren; i++){
     //ZS: Give AU penalty for every non-GC pair. 
     //the aupen function is defined in LoopScoring.h and will do this automatically
     //so we just call it 
     energy += param->auend*auPen(RNA[node->children[pairedChildren[i]]->lowBase.index], RNA[node->children[pairedChildren[i]]->highBase.index]);
     if(auPen(RNA[node->children[pairedChildren[i]]->lowBase.index], RNA[node->children[pairedChildren[i]]->highBase.index]) != 0){
				 if(printOn1)printf("AU penalty awarded for branch nr. %i (%i, %i%): %i  \n", i, 
				 RNA[node->children[pairedChildren[i]]->lowBase.index], 
				 RNA[node->children[pairedChildren[i]]->highBase.index], 
				 param->auend*auPen(RNA[node->children[pairedChildren[i]]->lowBase.index], 
				 RNA[node->children[pairedChildren[i]]->highBase.index]) );
				 }
   }
   
   //Check the same node for dangling
    energy += param->auend*auPen(RNA[node->lowBase.index], RNA[node->highBase.index]);
	 if(auPen(RNA[node->lowBase.index],RNA[node->highBase.index])>0){
        if(printOn1)printf("AU penalty awarded for root branch with bases (%i, %i%): %i  \n", 
				 RNA[node->lowBase.index], 
				 RNA[node->highBase.index],
				 param->auend*auPen(RNA[node->lowBase.index], RNA[node->highBase.index]));
	 }
        if(printOn1)printf("eM(i=%d,j=%d)=%f\n",(node->lowBase).index, node->highBase.index, (double)energy/100);
	return energy;
}


int eE(TreeNode* node, int* pairedChildren, int numPairedChildren, int* RNA, nndb_constants* param, int length){
	int energy;
	energy = 0;
	int i;
	for(i = 0; i < numPairedChildren; i++){
		//Shel: Give AU penalty for every non-GC pair. 
		//the aupen function is defined in LoopScoring.h and will do this automatically
		//so we just call it 
		energy += param->auend*auPen(RNA[node->children[pairedChildren[i]]->lowBase.index], RNA[node->children[pairedChildren[i]]->highBase.index]);
		if(auPen(RNA[node->children[pairedChildren[i]]->lowBase.index], RNA[node->children[pairedChildren[i]]->highBase.index]) != 0){
			if(printOn1)printf("AU penalty awarded for exterior branch nr. %i (%i, %i%): %i  \n", i, 
				   RNA[node->children[pairedChildren[i]]->lowBase.index], 
				   RNA[node->children[pairedChildren[i]]->highBase.index], 
				   param->auend*auPen(RNA[node->children[pairedChildren[i]]->lowBase.index], 
									  RNA[node->children[pairedChildren[i]]->highBase.index]) );
		}
	}
    //Here the dangling bases are also taken into account 
	if ((node->children[0])->isPair == 0 || D2MODE==1) {//If there is a initial dangling end
		int j;
		i = (node->children[pairedChildren[0]])->lowBase.index;
		j = (node->children[pairedChildren[0]])->highBase.index;
		if(NODANGLEMODE==1) energy += 0;
		else {
			//if(i==1){ if(D2MODE==1) energy+=0; else energy += param->dangle[RNA[j]][RNA[i]][RNA[length]][1];}//TODO Manoj Changed it for D2
			if(i==1){ 
				energy += param->dangle[RNA[j]][RNA[i]][RNA[length]][1];
				if(printOn1)printf("i=%d,j=%d,param->dangle[RNA[j]][RNA[i]][RNA[length]][1]/100=%f\n",i,j,param->dangle[RNA[j]][RNA[i]][RNA[length]][1]/100.0);
			}
			else{
				energy += param->dangle[RNA[j]][RNA[i]][RNA[i-1]][1];	
				if(printOn1)printf("i=%d,j=%d,param->dangle[RNA[j]][RNA[i]][RNA[i-1]][1]/100=%f\n",i,j,param->dangle[RNA[j]][RNA[i]][RNA[i-1]][1]/100.0);
			}
		}
	}
	
	
	if(numPairedChildren > 1) {	//If there is more than one branch on the external loop			
		for (i = 0; i < numPairedChildren-1; i++) {
			//Scores the dangling ends in unpaired regions between paired children
			energy += eMUnpairedRegion(node->children[pairedChildren[i]]->lowBase.index, node->children[pairedChildren[i]]->highBase.index, 
									   node->children[pairedChildren[i+1]]->lowBase.index, node->children[pairedChildren[i+1]]->highBase.index, 
									   RNA, param);	
			//printf("Energy so far: %i \n", energy);
		}
	}	
	if(pairedChildren[numPairedChildren-1] < node->numChildren - 1 || D2MODE==1){ //If there is a trailing dangling end
		int j;
		i = (node->children[pairedChildren[numPairedChildren-1]])->lowBase.index;
		j = (node->children[pairedChildren[numPairedChildren-1]])->highBase.index;
		if(NODANGLEMODE==1) energy += 0;
                else{
			//if(j==length){ if(D2MODE==1) energy+=0; else energy += param->dangle[RNA[j]][RNA[i]][RNA[1]][0];}//TODO Manoj Changed it for D2
			if(j==length){
				 energy += param->dangle[RNA[j]][RNA[i]][RNA[1]][0];
				if(printOn1)printf("i=%d,j=%d,param->dangle[RNA[j]][RNA[i]][RNA[1]][0]/100=%f\n",i,j,param->dangle[RNA[j]][RNA[i]][RNA[1]][0]/100.0);
			}
			else{
			 energy += param->dangle[RNA[j]][RNA[i]][RNA[j+1]][0];	
			if(printOn1)printf("i=%d,j=%d,param->dangle[RNA[j]][RNA[i]][RNA[j+1]][0]/100=%f\n",i,j,param->dangle[RNA[j]][RNA[i]][RNA[j+1]][0]/100.0);
			}
		}
	}
        if(printOn1)printf("eE(i=%d,j=%d)=%f\n",(node->lowBase).index, node->highBase.index, (double)energy/100);

	return energy; 
}

