#ifndef _MY_DOUBLE_H_
#define _MY_DOUBLE_H_

#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "gmp.h"
using namespace std;

//#define USE_DOUBLE_ONLY 1

//#ifdef USE_DOUBLE_ONLY
//#define MyDouble double
//#else
const int PRECISION = 1024;
static char BIG_NUM_ENABLED = 'Y';//'N';//'Y';
static int BIGNUM_ONLY=0;
const int PRINT_DIGITS_AFTER_DECIMAL = 10;
static int verbose=0;
/*mpf_t getBigNum(double val2){
  mpf_t* bigValue = new mpf_t;
  mpf_init2(*bigValue);
  mpf_set_d(*bigValue, val2);
  return *bigValue;
  }*/

class MyDouble{
	private: 
		mpf_t* bigValue;
		double* smallValue;
		char isBig;//'y' means it is already BigNum, 'n' means it is native double
	public:
		MyDouble(){if(verbose==1)printf("Default constructor called\n");
			//this((double)0.0);
			bigValue=0;smallValue=0;isBig='n';
			createDouble();			
		}
	/*	void init(){
			bigValue=0;smallValue=0;isBig='n';
			createDouble();
		}*/
		MyDouble(char isBig1){
			if(verbose==1)printf("Constructor with input char isBig1=%c\n",isBig1);
			bigValue=0;smallValue=0;isBig=isBig1;
			if(isBig1=='n') createDouble();
			else if(BIG_NUM_ENABLED=='N'){
				printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
				exit(-1);
			}
			else createBigNum();

		}
		MyDouble(mpf_t val2){
			if(verbose==1){cout<<"Constructor with input mpf_t val2=";gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, val2);cout<<endl;}
			bigValue=0;smallValue=0;isBig='y';
			createBigNum(val2);	
		}
		MyDouble(double val2){
			if(verbose==1){cout<<"Constructor with input double val2="<<val2<<endl;}
			bigValue=0;smallValue=0;isBig='n';
			createDouble(val2);	
		}
		void createBigNum(){
			if(smallValue!=0){delete smallValue; smallValue = 0;}
			if(BIG_NUM_ENABLED=='N'){
				printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
				exit(-1);
			}
			isBig = 'y';
			if(bigValue==0){
				bigValue = new mpf_t[1];
				mpf_init2(*bigValue,PRECISION);
			}
		}
		void createDouble(){
			double val2 =0.0;
			createDouble(val2);
		}
		void createBigNum(mpf_t val2){
			//smallValue = 0;
			if(smallValue!=0){delete smallValue; smallValue = 0;}
			if(BIG_NUM_ENABLED=='N'){
				printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
				exit(-1);
			}
			isBig = 'y';
			if(bigValue==0){
				if(verbose==1)printf("Allocation mpf_t\n");
				bigValue = new mpf_t[1];
				mpf_init2(*bigValue,PRECISION);//TODO This line was earlier outside this "if"
			}
			//mpf_init2(*bigValue,PRECISION);
			mpf_set(*bigValue, val2);//value=val2;
			//mpf_clear(val2);//Do not un-comment this line, this line will cause errors
		}
		void createDouble(double val2){
			//bigValue = 0;
			if(BIGNUM_ONLY==1){
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, val2);
			       	createBigNum(op2);
				mpf_clear(op2);
				return;
			}
			if(bigValue!=0){ mpf_clear(*bigValue); delete bigValue; bigValue=0;}
			if(smallValue==0){
				if(verbose==1)printf("Allocation double for %f\n",val2);
				smallValue = new double;
			}
			*smallValue = val2;
			isBig='n';
		}
		void deallocate(){
			if(isBig=='y'){ if(bigValue!=0){ if(verbose==1) printf("Deallocation mpf_t\n"); mpf_clear(*bigValue); delete(bigValue); bigValue=0;}isBig='X';}
			else if(isBig=='n'){ if(smallValue!=0){ if(verbose==1) printf("Deallocation double for %f\n",*smallValue); delete(smallValue); smallValue=0;} isBig='X';}
			else if(verbose==1) printf("In MyDouble::deallocate(), Unknown isBig = %c\n", isBig);
			bigValue=0;smallValue=0;isBig='X';
		}
		~MyDouble(){if(verbose==1)printf("Destructor called\n");
			deallocate();	
		}
		bool isInitialized(){
			if(isBig=='y' || isBig=='n') return true;
			return false;
		}
		void reset(){
			if(isBig=='y') mpf_clear(*bigValue);
			else if(isBig=='n') *smallValue = 0;  
			else printf("Unknown isBig = %c\n", isBig);
		}
		void print()const{
			//if(isBig=='y') gmp_printf("fixed point mpf %.*Ff with %d digits\n", 5, *bigValue, 5);
			if(isBig=='y') gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
			//else if(isBig=='n') printf("double %f", *smallValue);//TODO uncomment it
			else if(isBig=='n') printf("%f", *smallValue);
			else printf("Unknown isBig = %c\n", isBig);
		}
		void printInt()const{
			//if(isBig=='y') gmp_printf("fixed point mpf %.*Ff with %d digits\n", 5, *bigValue, 5);
			if(isBig=='y') gmp_printf("mpf %.*Ff", 1, *bigValue);
			//else if(isBig=='n') printf("double %f", *smallValue);//TODO uncomment it
			else if(isBig=='n') printf("%d", (int)(*smallValue));
			else printf("Unknown isBig = %c\n", isBig);
		}
		void print(FILE* outFile)const{
			//if(isBig=='y') gmp_printf("fixed point mpf %.*Ff with %d digits\n", 5, *bigValue, 5);
			if(isBig=='y') gmp_fprintf(outFile, "%.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
			//else if(isBig=='n') printf("double %f", *smallValue);//TODO uncomment it
			else if(isBig=='n') fprintf(outFile, "%f", *smallValue);
			else fprintf(outFile, "Unknown isBig = %c\n", isBig);
		}
		MyDouble operator*(const MyDouble &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_mul(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
				return res;
			}
			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			else if(this->isBig=='y' && obj1.isBig=='n'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_mul(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_mul(*(res.bigValue),op1, *(obj1.bigValue));
				mpf_clear(op1);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n' && obj1.isBig=='n'){
				double a = (*(this->smallValue)) * (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double multiplication\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_mul(*(res.bigValue),op1, op2);
					mpf_clear(op1);
					mpf_clear(op2);
					//res.print();
					if(verbose==1) printf("successful multiplication\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		MyDouble operator*(const double &obj1_double) const {
			const MyDouble obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_mul(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n'){
				double a = (*(this->smallValue)) * (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double multiplication\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_mul(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);
					mpf_clear(op2);
					if(verbose==1)printf("successful multiplication\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}

		MyDouble operator+(const MyDouble &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_add(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
				return res;
			}
			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			else if(this->isBig=='y' && obj1.isBig=='n'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_add(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_add(*(res.bigValue),op1, *(obj1.bigValue));
				mpf_clear(op1);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n' && obj1.isBig=='n'){if(verbose==1)printf("operator+ MyDouble obj1 this->isBig=='n' && obj1.isBig=='n'\n");
				double a = (*(this->smallValue)) + (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double addition\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_add(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					if(verbose==1)printf("successful addition\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		MyDouble operator+(const double &obj1_double) const {
			const MyDouble obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_add(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n'){
				double a = (*(this->smallValue)) + (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double addition\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_add(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);
					mpf_clear(op2);
					if(verbose==1)printf("successful multiplication\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		MyDouble operator-(const MyDouble &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_sub(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
				return res;
			}
			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			else if(this->isBig=='y' && obj1.isBig=='n'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_sub(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_sub(*(res.bigValue),op1, *(obj1.bigValue));
				mpf_clear(op1);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n' && obj1.isBig=='n'){
				double a = (*(this->smallValue)) - (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double subtraction\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_sub(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					if(verbose==1)printf("successful subtraction\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		MyDouble operator-(const double &obj1_double) const {
			const MyDouble obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_sub(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n'){
				double a = (*(this->smallValue)) - (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double subtraction\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_sub(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					if(verbose==1)printf("successful subtraction\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		MyDouble operator/(const MyDouble &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_div(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
				return res;
			}
			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			else if(this->isBig=='y' && obj1.isBig=='n'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_div(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_div(*(res.bigValue),op1, *(obj1.bigValue));
				mpf_clear(op1);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n' && obj1.isBig=='n'){
				double a = (*(this->smallValue)) / (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double division\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_div(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					if(verbose==1)printf("successful division\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		MyDouble operator/(const double &obj1_double) const {
			const MyDouble obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
				mpf_div(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n'){
				double a = (*(this->smallValue)) / (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					if(verbose==1)printf("Overflow Error occurred while native double division\n");
					if(BIG_NUM_ENABLED == 'N'){
						printf("BIG_NUM_ENABLED is not Enabled, hence exiting\n");
						exit(-1);
					}
					MyDouble res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,PRECISION); mpf_set_d(op2, *(obj1.smallValue));
					mpf_div(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					if(verbose==1)printf("successful division\n");
					return res;
				}
			}
			else{
				 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);	
			}
		}
		int compare(const MyDouble &obj1) const{
			//Function: int mpf_cmp (mpf_t op1, mpf_t op2)
			//Function: int mpf_cmp_d (mpf_t op1, double op2)
			//case 1: this object is bigValue and obj1 is also bigValue 
			int result=100;
			if(this->isBig=='y' && obj1.isBig=='y'){
				//return mpf_cmp(*(this->bigValue), *(obj1.bigValue));
				result = mpf_cmp(*(this->bigValue), *(obj1.bigValue));
			}
			//case 2: this object is bigValue and obj1 is smallValue
			else if(this->isBig=='y' && obj1.isBig=='n'){
				//return mpf_cmp_d(*(this->bigValue), *(obj1.smallValue));
				result =  mpf_cmp_d(*(this->bigValue), *(obj1.smallValue));
			}
			//case 3: this object is smallValue and obj2 is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				//mpf_t minusOne; mpf_init2(minusOne,PRECISION); mpf_set_d(minusOne, -1.0);
				//return mpf_mul(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)), minusOne);
				//return -1*(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)));
				result = -1*(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)));
			}
			//case 4: this object is smallValue and obj2 is smallValue
			else if(this->isBig=='n' && obj1.isBig=='n'){
				//return (*(this->smallValue)) - (*(obj1.smallValue));
				double result1 = (*(this->smallValue)) - (*(obj1.smallValue));
				if(result1==0) result=0;
				else if(result1<0) result=-1;
				else result=+1;
			}
			else{
                                 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);
                        }
			//printf("comparing: ");print();printf(" and ");obj1.print();printf(" and result is %d\n",result);
			return result;
		}
		int compare(const double &obj1) const{
			//Function: int mpf_cmp (mpf_t op1, mpf_t op2)
			//Function: int mpf_cmp_d (mpf_t op1, double op2)
			//case 2: this object is bigValue
			int result = 200;
			if(this->isBig=='y'){
				//return mpf_cmp_d(*(this->bigValue), obj1);
				result = mpf_cmp_d(*(this->bigValue), obj1);
			}
			//case 4: this object is smallValue
			else if(this->isBig=='n'){
				//return (*(this->smallValue)) - obj1;
				double result1 = (*(this->smallValue)) - obj1;
				if(result1==0.0) result=0;
                                else if(result1<0.0) result=-1;
                                else result=+1;

			}
			else{
                                 printf("Unknown isBig = %c\n", isBig);
                        }
			//printf("comparing: ");print();printf(" and %f and result is %d\n",obj1,result);
			return result;
		}
		bool operator==(const MyDouble &obj1) const {
			return compare(obj1)==0;		
		}
		bool operator==(const double &obj1) const {
			return compare(obj1)==0;		
		}
		bool operator!=(const MyDouble &obj1) const {
			return compare(obj1)!=0;		
		}
		bool operator!=(const double &obj1) const {
			return compare(obj1)!=0;		
		}
		bool operator<(const MyDouble &obj1) const {
			return compare(obj1)<0;		
		}
		bool operator<(const double &obj1) const {
			return compare(obj1)<0;		
		}
		bool operator>(const MyDouble &obj1) const {
			return compare(obj1)>0;		
		}
		bool operator>(const double &obj1) const {
			return compare(obj1)>0;		
		}
		bool operator<=(const MyDouble &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator<=(const double &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator>=(const MyDouble &obj1) const {
                        return compare(obj1)>=0;
                }
                bool operator>=(const double &obj1) const {
                        return compare(obj1)>=0;
                }

		MyDouble& operator=(const MyDouble &obj1) {
			if(verbose==1){ cout<<"operator= called: obj1=";obj1.print();cout<<", this=";this->print();cout<<endl;}
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			if(this==&obj1) return *this;
			if(isInitialized())this->deallocate();
			else {bigValue=0; smallValue=0;}
			//printf("successful deallocation\n");
			if(obj1.isBig=='n'){ isBig='n'; createDouble(*(obj1.smallValue));}
			else if(obj1.isBig=='y'){isBig='y'; createBigNum(*(obj1.bigValue));}
			return *this;
		}
		MyDouble& operator=(const double &obj1) {
			if(verbose==1){ cout<<"operator= called: obj1="<<obj1;cout<<", this=";this->print();cout<<endl;}
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			//if(this==&obj1) return *this;
			if(isInitialized())this->deallocate();
			else {bigValue=0; smallValue=0;}
			//printf("successful deallocation\n");
			isBig='n';
			createDouble(obj1);
			return *this;
		}
		MyDouble(const MyDouble &obj1) {
			bigValue=0;smallValue=0;isBig='n';
			if(verbose==1){ cout<<"Copy constructor called: obj1=";obj1.print();cout<<endl;}//cout<<", this=";if(isInitialized())this->print();else cout<<"Uninitialized,";cout<<endl;
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			//if(this==&obj1) return ;//*this;
			//if(isInitialized())this->deallocate();
			//printf("successful deallocation\n");
			if(obj1.isBig=='n'){ isBig='n'; createDouble(*(obj1.smallValue));}
			else if(obj1.isBig=='y'){ isBig='y'; createBigNum(*(obj1.bigValue));}
			//return *this;
		}

};
//#endif
#endif
/*
int main(){
	MyDouble c;
	double a_d = pow(10, 308);
	double b_d = pow(10, 307);
	cout<<"a_d = "<<a_d<<endl;
	cout<<"b_d = "<<b_d<<endl;
	//cout<<"a_d * b_d = "<<a_d * b_d<<endl;
	cout<<"a_d + b_d = "<<5*a_d + 5*b_d<<endl;
	MyDouble a(a_d);
	MyDouble b(b_d);
	//c=a*b;
	//c=5*a+5*b;
	//c=a*5+b*5;
	//c=a*5+b+b+b+b+b;
	//c = MyDouble(500.01);
	c=b+b+b+b+b+a*5;
	//printf("a==c->%d, a==a_d->%d, a.compare(c)->%d, c.compare(a_d)->%d\n", a==c, a==a_d, a.compare(c), c.compare(a_d));
	printf("a!=c->%d, a!=b_d->%d, a.compare(c)->%d, c.compare(a_d)->%d\n", a!=c, a!=b_d, a.compare(c), c.compare(a_d));
	c.print();
	printf("\n");
	c.deallocate();
	printf("Default constructor checking\n");
	MyDouble d(1.0); d.print();printf("\n");
	return 0;
}
*/
/*
MyDouble f1(MyDouble a){
	MyDouble b(5.0);
	MyDouble c = a+b;
	return c;
}
MyDouble f2(MyDouble a){
	return f1(a)+MyDouble(3.0);
	//return a+MyDouble(3.0);
}
int main(){
	MyDouble a(4.0);
	a = f2(a);
	//a = a+MyDouble(3.0);
	a.print();
	cout<<endl;
	return 0;
}*/
