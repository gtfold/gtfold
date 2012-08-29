#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "gmp.h"
using namespace std;
const int PRECISION = 64;
static char BIG_NUM_ENABLED = 'Y';
const int PRINT_DIGITS_AFTER_DECIMAL = 10;
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
		MyDouble(){
			//this((double)0.0);
			 createDouble();			
		}
		MyDouble(char isBig1){
			if(isBig1=='n') createDouble();
			else if(BIG_NUM_ENABLED=='N'){
				printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
				exit(-1);
			}
			else createBigNum();

		}
		MyDouble(mpf_t val2){
			createBigNum(val2);	
		}
		MyDouble(double val2){
			createDouble(val2);	
		}
		void createBigNum(){
			smallValue = 0;
			if(BIG_NUM_ENABLED=='N'){
				printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
				exit(-1);
			}
			isBig = 'y';
			bigValue = new mpf_t[1];
			mpf_init2(*bigValue,PRECISION);
		}
		void createDouble(){
			double val2 =0.0;
			createDouble(val2);
		}
		void createBigNum(mpf_t val2){
			smallValue = 0;
			if(BIG_NUM_ENABLED=='N'){
				printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
				exit(-1);
			}
			isBig = 'y';
			bigValue = new mpf_t[1];
			mpf_init2(*bigValue,PRECISION);
			mpf_set(*bigValue, val2);//value=val2;
		}
		void createDouble(double val2){
			bigValue = 0;
			smallValue = new double;
			*smallValue = val2;
			isBig='n';
		}
		void deallocate(){
			if(isBig=='y'){ delete(bigValue); isBig='X';}
                        else if(isBig=='n'){ delete(smallValue); isBig='X';}
			else printf("Unknown isBig = %c\n", isBig);
		}
		~MyDouble(){
			//deallocate();	
		}
		void reset(){
			if(isBig=='y') mpf_clear(*bigValue);
			else if(isBig=='n') *smallValue = 0;  
			else printf("Unknown isBig = %c\n", isBig);
		}
		void print()const{
			//if(isBig=='y') gmp_printf("fixed point mpf %.*Ff with %d digits\n", 5, *bigValue, 5);
			if(isBig=='y') gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
			else if(isBig=='n') printf("double %f", *smallValue);
			else printf("Unknown isBig = %c\n", isBig);
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
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_mul(*(res.bigValue),op1, *(obj1.bigValue));
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
					printf("Overflow Error occurred while native double multiplication\n");
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
					printf("successful multiplication\n");
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
					printf("Overflow Error occurred while native double multiplication\n");
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
					printf("successful multiplication\n");
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
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_add(*(res.bigValue),op1, *(obj1.bigValue));
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else if(this->isBig=='n' && obj1.isBig=='n'){
				double a = (*(this->smallValue)) + (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					MyDouble res;
					res.createDouble(a);
					return res;
				}
				else{
					printf("Overflow Error occurred while native double addition\n");
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
					printf("successful multiplication\n");
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
					printf("Overflow Error occurred while native double multiplication\n");
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
					printf("successful multiplication\n");
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
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_sub(*(res.bigValue),op1, *(obj1.bigValue));
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
					printf("Overflow Error occurred while native double addition\n");
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
					printf("successful multiplication\n");
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
					printf("Overflow Error occurred while native double multiplication\n");
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
					printf("successful multiplication\n");
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
				return res;
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				MyDouble res;
				res.createBigNum();
				mpf_t op1; mpf_init2(op1,PRECISION); mpf_set_d(op1, *(this->smallValue));
				mpf_div(*(res.bigValue),op1, *(obj1.bigValue));
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
					printf("Overflow Error occurred while native double addition\n");
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
					printf("successful multiplication\n");
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
					printf("Overflow Error occurred while native double multiplication\n");
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
					printf("successful multiplication\n");
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
			if(this->isBig=='y' && obj1.isBig=='y'){
				return mpf_cmp(*(this->bigValue), *(obj1.bigValue));
			}
			//case 2: this object is bigValue and obj1 is smallValue
			else if(this->isBig=='y' && obj1.isBig=='n'){
				return mpf_cmp_d(*(this->bigValue), *(obj1.smallValue));
			}
			//case 3: this object is smallValue and obj2 is bigValue
			else if(this->isBig=='n' && obj1.isBig=='y'){
				//mpf_t minusOne; mpf_init2(minusOne,PRECISION); mpf_set_d(minusOne, -1.0);
				//return mpf_mul(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)), minusOne);
				return -1*(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)));
			}
			//case 4: this object is smallValue and obj2 is smallValue
			else if(this->isBig=='n' && obj1.isBig=='n'){
				return (*(this->smallValue)) - (*(obj1.smallValue));
			}
			else{
                                 printf("Unknown isBig = %c, obj1.isBig = %c\n", isBig, obj1.isBig);
                        }
		}
		int compare(const double &obj1) const{
			//Function: int mpf_cmp (mpf_t op1, mpf_t op2)
			//Function: int mpf_cmp_d (mpf_t op1, double op2)
			//case 2: this object is bigValue
			if(this->isBig=='y'){
				return mpf_cmp_d(*(this->bigValue), obj1);
			}
			//case 4: this object is smallValue
			else if(this->isBig=='n'){
				return (*(this->smallValue)) - obj1;
			}
			else{
                                 printf("Unknown isBig = %c\n", isBig);
                        }
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

		MyDouble& operator=(const MyDouble &obj1) {
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			if(this==&obj1) return *this;
			this->deallocate();
			//printf("successful deallocation\n");
			if(obj1.isBig=='n') createDouble(*(obj1.smallValue));
			else if(obj1.isBig=='y') createBigNum(*(obj1.bigValue));
			return *this;
		}
};
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
