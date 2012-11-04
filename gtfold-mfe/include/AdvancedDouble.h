#ifndef _ADVANCED_DOUBLE_H_
#define _ADVANCED_DOUBLE_H_

#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "gmp.h"
using namespace std;

extern int g_bignumprecision;
const int PRINT_DIGITS_AFTER_DECIMAL = 20;

class AdvancedDouble_Native{
	private:
		double value;
	public:
		AdvancedDouble_Native(){
			value=0.0;
		}
		AdvancedDouble_Native(double val){
			value=val;
		}
		AdvancedDouble_Native(const AdvancedDouble_Native &obj1) {
			value=obj1.value;
		}
		void init(){
			value=0.0;
		}
		void deallocate(){

		}
		~AdvancedDouble_Native(){
			//deallocate();
		}
		/*bool isInitialized(){
		      return true;
		}
		void reset(){
			value=0.0;
		}*/
		void print()const{
			//printf("%f", value);
			printf("%.20f", value);
		}
		void printInt()const{
			printf("%d", (int)(value));
		}
		void print(FILE* outFile)const{
			//fprintf(outFile, "%f", value);
			fprintf(outFile, "%.20f", value);
		}
		AdvancedDouble_Native operator*(const AdvancedDouble_Native &obj1) const {
			return (AdvancedDouble_Native)(value*obj1.value);
		}
		AdvancedDouble_Native operator*(const double &obj1_double) const {
			return (AdvancedDouble_Native)(value*obj1_double);
		}
		AdvancedDouble_Native operator+(const AdvancedDouble_Native &obj1) const {
			return (AdvancedDouble_Native)(value+obj1.value);
		}
		AdvancedDouble_Native operator+(const double &obj1_double) const {
			return (AdvancedDouble_Native)(value+obj1_double);
		}
		AdvancedDouble_Native operator-(const AdvancedDouble_Native &obj1) const {
			return (AdvancedDouble_Native)(value-obj1.value);
		}
		AdvancedDouble_Native operator-(const double &obj1_double) const {
			return (AdvancedDouble_Native)(value-obj1_double);
		}
		AdvancedDouble_Native operator/(const AdvancedDouble_Native &obj1) const {
			return (AdvancedDouble_Native)(value/obj1.value);
		}
		AdvancedDouble_Native operator/(const double &obj1_double) const {
			return (AdvancedDouble_Native)(value/obj1_double);
		}
		/*int compare(const AdvancedDouble_Native &obj1) const{
			double cmp = value-obj1.value;
			if(cmp>0) return 1;
			else if(cmp<0) return -1;
			else return 0;
		}
		int compare(const double &obj1) const{
			double cmp = value-obj1;
			if(cmp>0) return 1;
                        else if(cmp<0) return -1;
                        else return 0;

		}
		bool operator==(const AdvancedDouble_Native &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator==(const double &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator!=(const AdvancedDouble_Native &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator!=(const double &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator<(const AdvancedDouble_Native &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator<(const double &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator>(const AdvancedDouble_Native &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator>(const double &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator<=(const AdvancedDouble_Native &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator<=(const double &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator>=(const AdvancedDouble_Native &obj1) const {
                        return compare(obj1)>=0;
                }
                bool operator>=(const double &obj1) const {
                        return compare(obj1)>=0;
                }*/

		bool operator==(const AdvancedDouble_Native &obj1) const {
                        return value==obj1.value;
                }
                bool operator==(const double &obj1) const {
                        return value==obj1;
                }
                bool operator!=(const AdvancedDouble_Native &obj1) const {
                        return value!=obj1.value;
                }
                bool operator!=(const double &obj1) const {
                        return value!=obj1;
                }
                bool operator<(const AdvancedDouble_Native &obj1) const {
                        return value<obj1.value;
                }
                bool operator<(const double &obj1) const {
                        return value<obj1;
                }
                bool operator>(const AdvancedDouble_Native &obj1) const {
                        return value>obj1.value;
                }
                bool operator>(const double &obj1) const {
                        return value>obj1;
                }
                bool operator<=(const AdvancedDouble_Native &obj1) const {
                        return value<=obj1.value;
                }
                bool operator<=(const double &obj1) const {
                        return value<=obj1;
                }
                bool operator>=(const AdvancedDouble_Native &obj1) const {
                        return value>=obj1.value;
                }
                bool operator>=(const double &obj1) const {
                        return value>=obj1;
                }

		AdvancedDouble_Native& operator=(const AdvancedDouble_Native &obj1) {
			if(this==&obj1) return *this;
			//if(isInitialized())this->deallocate();//TODO
			//else {bigValue=0;}
			value=obj1.value;
			return *this;
		}
		AdvancedDouble_Native& operator=(const double &obj1) {
			//if(isInitialized())this->deallocate();//TODO
			//else {bigValue=0; smallValue=0;}
			value=obj1;
			return *this;
		}
};

class AdvancedDouble_BigNum{
	private:
		mpf_t* bigValue;
	public:
		AdvancedDouble_BigNum(){
			bigValue=0;
			createBigNum();
		}
		AdvancedDouble_BigNum(mpf_t val2){
			bigValue=0;
			createBigNum(val2);
		}
		AdvancedDouble_BigNum(double val){
			bigValue=0;
			createBigNum(val);
		}
		AdvancedDouble_BigNum(const AdvancedDouble_BigNum &obj1) {
			bigValue=0;
			createBigNum(*(obj1.bigValue));
		}
		void init(){
			bigValue=0;
			createBigNum();
		}
		void createBigNum(){
			if(bigValue==0){
				bigValue = new mpf_t[1];
				mpf_init2(*bigValue,g_bignumprecision);
			}
		}
		void createBigNum(mpf_t val2){
			createBigNum();
			mpf_set(*bigValue, val2);
		}
		void createBigNum(double val2){
			mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, val2);
			createBigNum(op2);
			mpf_clear(op2);
		}
		void deallocate(){
			if(bigValue!=0){mpf_clear(*bigValue); delete(bigValue); bigValue=0;}
		}
		~AdvancedDouble_BigNum(){
			deallocate();
		}
		bool isInitialized(){
			if(bigValue!=0) return true;
			return false;
		}
		void reset(){
			mpf_clear(*bigValue);
		}
		void print()const{
			if(bigValue!=0) gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
		}
		void printInt()const{
			if(bigValue!=0) gmp_printf("mpf %.*Ff", 1, *bigValue);
		}
		void print(FILE* outFile)const{
			if(bigValue!=0) gmp_fprintf(outFile, "%.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
		}
		AdvancedDouble_BigNum operator*(const AdvancedDouble_BigNum &obj1) const {
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_mul(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator*(const double &obj1_double) const {
			const AdvancedDouble_BigNum obj1(obj1_double);
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_mul(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator+(const AdvancedDouble_BigNum &obj1) const {
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_add(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator+(const double &obj1_double) const {
			const AdvancedDouble_BigNum obj1(obj1_double);
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_add(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator-(const AdvancedDouble_BigNum &obj1) const {
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_sub(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator-(const double &obj1_double) const {
			const AdvancedDouble_BigNum obj1(obj1_double);
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_sub(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator/(const AdvancedDouble_BigNum &obj1) const {
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_div(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		AdvancedDouble_BigNum operator/(const double &obj1_double) const {
			const AdvancedDouble_BigNum obj1(obj1_double);
			AdvancedDouble_BigNum res;
			res.createBigNum();
			mpf_div(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
			return res;
		}
		int compare(const AdvancedDouble_BigNum &obj1) const{
			return mpf_cmp(*(this->bigValue), *(obj1.bigValue));
		}
		int compare(const double &obj1) const{
			return mpf_cmp_d(*(this->bigValue), obj1);
		}
		bool operator==(const AdvancedDouble_BigNum &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator==(const double &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator!=(const AdvancedDouble_BigNum &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator!=(const double &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator<(const AdvancedDouble_BigNum &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator<(const double &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator>(const AdvancedDouble_BigNum &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator>(const double &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator<=(const AdvancedDouble_BigNum &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator<=(const double &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator>=(const AdvancedDouble_BigNum &obj1) const {
                        return compare(obj1)>=0;
                }
                bool operator>=(const double &obj1) const {
                        return compare(obj1)>=0;
                }
		AdvancedDouble_BigNum& operator=(const AdvancedDouble_BigNum &obj1) {
			if(this==&obj1) return *this;
			//if(isInitialized())this->deallocate();//TODO
			//else {bigValue=0;}
			createBigNum(*(obj1.bigValue));
			return *this;
		}
		AdvancedDouble_BigNum& operator=(const double &obj1) {
			//if(isInitialized())this->deallocate();//TODO
			//else {bigValue=0; smallValue=0;}
			createBigNum(obj1);
			return *this;
		}
};

class AdvancedDouble_BigNumOptimized{
	private:
		mpf_t bigValue;
	public:
		AdvancedDouble_BigNumOptimized(){
			mpf_init2(bigValue,g_bignumprecision);
		}
		AdvancedDouble_BigNumOptimized(mpf_t val2){
			mpf_init_set(bigValue, val2);
		}
		AdvancedDouble_BigNumOptimized(double val){
			mpf_init_set_d(bigValue, val);
		}
		AdvancedDouble_BigNumOptimized(const AdvancedDouble_BigNumOptimized &obj1) {
			mpf_init_set(bigValue, obj1.bigValue);
		}
		void init(){
			mpf_init2(bigValue,g_bignumprecision);
		}
		/*
		void createBigNum(){
			//if(bigValue==0){
			//	bigValue = new mpf_t[1];
				mpf_init2(bigValue,g_bignumprecision);
			//}
		}
		void createBigNum(mpf_t val2){
			createBigNum();
			mpf_set(bigValue, val2);
		}
		void createBigNum(double val2){
			mpf_init2(bigValue,g_bignumprecision); mpf_set_d(bigValue, val2);
		}*/
		void deallocate(){
			mpf_clear(bigValue);
		}
		~AdvancedDouble_BigNumOptimized(){
			//deallocate();
			mpf_clear(bigValue);
		}
		/*
		bool isInitialized(){
			if(bigValue!=0) return true;
			return false;
		}
		void reset(){
			mpf_clear(bigValue);
		}*/
		void print()const{
			//if(bigValue!=0) gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
			gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, bigValue);
		}
		void printInt()const{
			//if(bigValue!=0) gmp_printf("mpf %.*Ff", 1, *bigValue);
			gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, bigValue);
		}
		void print(FILE* outFile)const{
			//if(bigValue!=0) gmp_fprintf(outFile, "%.*Ff", PRINT_DIGITS_AFTER_DECIMAL, *bigValue);
			gmp_fprintf(outFile, "%.*Ff", PRINT_DIGITS_AFTER_DECIMAL, bigValue);
		}
		AdvancedDouble_BigNumOptimized operator*(const AdvancedDouble_BigNumOptimized &obj1) const {
			AdvancedDouble_BigNumOptimized res;
			mpf_mul(res.bigValue,this->bigValue, obj1.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator*(const double &obj1_double) const {
			AdvancedDouble_BigNumOptimized res(obj1_double);
			mpf_mul(res.bigValue,this->bigValue, res.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator+(const AdvancedDouble_BigNumOptimized &obj1) const {
			AdvancedDouble_BigNumOptimized res;
			mpf_add(res.bigValue,this->bigValue, obj1.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator+(const double &obj1_double) const {
			AdvancedDouble_BigNumOptimized res(obj1_double);
			mpf_add(res.bigValue,this->bigValue, res.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator-(const AdvancedDouble_BigNumOptimized &obj1) const {
			AdvancedDouble_BigNumOptimized res;
			mpf_sub(res.bigValue,this->bigValue, obj1.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator-(const double &obj1_double) const {
			AdvancedDouble_BigNumOptimized res(obj1_double);
			mpf_sub(res.bigValue,this->bigValue, res.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator/(const AdvancedDouble_BigNumOptimized &obj1) const {
			AdvancedDouble_BigNumOptimized res;
			mpf_div(res.bigValue,this->bigValue, obj1.bigValue);
			return res;
		}
		AdvancedDouble_BigNumOptimized operator/(const double &obj1_double) const {
			AdvancedDouble_BigNumOptimized res(obj1_double);
			mpf_div(res.bigValue,this->bigValue, res.bigValue);
			return res;
		}

		/*int compare(const AdvancedDouble_BigNumOptimized &obj1) const{
			return mpf_cmp(this->bigValue, obj1.bigValue);
		}
		int compare(const double &obj1) const{
			return mpf_cmp_d(this->bigValue, obj1);
		}
		bool operator==(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator==(const double &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator!=(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator!=(const double &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator<(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator<(const double &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator>(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator>(const double &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator<=(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator<=(const double &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator>=(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return compare(obj1)>=0;
                }
                bool operator>=(const double &obj1) const {
                        return compare(obj1)>=0;
                }*/

		bool operator==(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return mpf_cmp(this->bigValue, obj1.bigValue)==0;
                }
                bool operator==(const double &obj1) const {
                        return mpf_cmp_d(this->bigValue, obj1)==0;
                }
                bool operator!=(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return mpf_cmp(this->bigValue, obj1.bigValue)!=0;
                }
                bool operator!=(const double &obj1) const {
                        return mpf_cmp_d(this->bigValue, obj1)!=0;
                }
                bool operator<(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return mpf_cmp(this->bigValue, obj1.bigValue)<0;
                }
                bool operator<(const double &obj1) const {
                        return mpf_cmp_d(this->bigValue, obj1)<0;
                }
                bool operator>(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return mpf_cmp(this->bigValue, obj1.bigValue)>0;
                }
                bool operator>(const double &obj1) const {
                        return mpf_cmp_d(this->bigValue, obj1)>0;
                }
                bool operator<=(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return mpf_cmp(this->bigValue, obj1.bigValue)<=0;
                }
                bool operator<=(const double &obj1) const {
                        return mpf_cmp_d(this->bigValue, obj1)<=0;
                }
                bool operator>=(const AdvancedDouble_BigNumOptimized &obj1) const {
                        return mpf_cmp(this->bigValue, obj1.bigValue)>=0;
                }
                bool operator>=(const double &obj1) const {
                        return mpf_cmp_d(this->bigValue, obj1)>=0;
                }

		AdvancedDouble_BigNumOptimized& operator=(const AdvancedDouble_BigNumOptimized &obj1) {
			if(this==&obj1) return *this;
			//if(isInitialized())this->deallocate();//TODO
			//else {bigValue=0;}
			//deallocate();//TODO
			//createBigNum(*(obj1.bigValue));
			mpf_set(bigValue, obj1.bigValue);
			return *this;
		}
		AdvancedDouble_BigNumOptimized& operator=(const double &obj1) {
			//if(isInitialized())this->deallocate();//TODO
			//else {bigValue=0; smallValue=0;}
			//createBigNum(obj1);
			mpf_set_d(bigValue, obj1);
			return *this;
		}
};


//static int BIGNUM_ONLY=0;
//static int DOUBLE_ONLY=0;
//static int verbose=0;
class AdvancedDouble_Hybrid{
	private: 
		mpf_t* bigValue;
		double* smallValue;
		char isBig;//'y' means it is already BigNum, 'n' means it is native double
	public:
		AdvancedDouble_Hybrid(){//if(verbose==1)printf("Default constructor called\n");
			//this((double)0.0);
			bigValue=0;smallValue=0;isBig='n';
			createDouble();			
		}
		void init(){
			bigValue=0;smallValue=0;isBig='n';
			createDouble();
		}
		AdvancedDouble_Hybrid(char isBig1){
			//if(verbose==1)printf("Constructor with input char isBig1=%c\n",isBig1);
			bigValue=0;smallValue=0;isBig=isBig1;
			if(isBig1=='n') createDouble();
			/*else if(BIG_NUM_ENABLED=='N'){
			  printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
			  exit(-1);
			  }*/
			else createBigNum();

		}
		AdvancedDouble_Hybrid(mpf_t val2){
			//if(verbose==1){cout<<"Constructor with input mpf_t val2=";gmp_printf("mpf %.*Ff", PRINT_DIGITS_AFTER_DECIMAL, val2);cout<<endl;}
			bigValue=0;smallValue=0;isBig='y';
			createBigNum(val2);	
		}
		AdvancedDouble_Hybrid(double val2){
			//if(verbose==1){cout<<"Constructor with input double val2="<<val2<<endl;}
			bigValue=0;smallValue=0;isBig='n';
			createDouble(val2);	
		}
		void createBigNum(){
			if(smallValue!=0){delete smallValue; smallValue = 0;}
			/*if(BIG_NUM_ENABLED=='N'){
			  printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
			  exit(-1);
			  }*/
			isBig = 'y';
			if(bigValue==0){
				bigValue = new mpf_t[1];
				mpf_init2(*bigValue,g_bignumprecision);
			}
		}
		void createDouble(){
			double val2 =0.0;
			createDouble(val2);
		}
		void createBigNum(mpf_t val2){
			//smallValue = 0;
			if(smallValue!=0){delete smallValue; smallValue = 0;}
			/*if(BIG_NUM_ENABLED=='N'){
			  printf("Error in creating mpf_t object as BIG_NUM_ENABLED is %c\n",BIG_NUM_ENABLED);
			  exit(-1);
			  }*/
			isBig = 'y';
			if(bigValue==0){
				//if(verbose==1)printf("Allocation mpf_t\n");
				bigValue = new mpf_t[1];
				mpf_init2(*bigValue,g_bignumprecision);//TODO This line was earlier outside this "if"
			}
			//mpf_init2(*bigValue,g_bignumprecision);
			mpf_set(*bigValue, val2);//value=val2;
			//mpf_clear(val2);//Do not un-comment this line, this line will cause errors
		}
		void createDouble(double val2){
			//bigValue = 0;
			/*if(BIGNUM_ONLY==1){
			  mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, val2);
			  createBigNum(op2);
			  mpf_clear(op2);
			  return;
			  }*/
			if(bigValue!=0){ mpf_clear(*bigValue); delete bigValue; bigValue=0;}
			if(smallValue==0){
				//if(verbose==1)printf("Allocation double for %f\n",val2);
				smallValue = new double;
			}
			*smallValue = val2;
			isBig='n';
		}
		void deallocate(){
			if(isBig=='y'){ if(bigValue!=0){ 
				//if(verbose==1) printf("Deallocation mpf_t\n"); 
				mpf_clear(*bigValue); delete(bigValue); bigValue=0;}isBig='X';
			}
			else if(isBig=='n'){ 
				if(smallValue!=0){ 
					//if(verbose==1) printf("Deallocation double for %f\n",*smallValue); 
					delete(smallValue); smallValue=0;
				}
				isBig='X';}
			else{
				//if(verbose==1) printf("In AdvancedDouble_Hybrid::deallocate(), Unknown isBig = %c\n", isBig);
			}	
			bigValue=0;smallValue=0;isBig='X';
		}
		~AdvancedDouble_Hybrid(){//if(verbose==1)printf("Destructor called\n");
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
		AdvancedDouble_Hybrid operator*(const AdvancedDouble_Hybrid &obj1) const {
			if(this->isBig=='y'){
				//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_mul(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
					return res;
				}
				//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_mul(*(res.bigValue),*(this->bigValue), op2);
					mpf_clear(op2);
					return res;
				}
			}
			else{	
				//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_mul(*(res.bigValue),op1, *(obj1.bigValue));
					mpf_clear(op1);
					return res;
				}
				//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
				else {
					double a = (*(this->smallValue)) * (*(obj1.smallValue));
					if(isfinite(a)){
						AdvancedDouble_Hybrid res;
						res.createDouble(a);
						return res;
					}
					else{
						AdvancedDouble_Hybrid res;
						res.createBigNum();
						mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
						mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
						mpf_mul(*(res.bigValue),op1, op2);
						mpf_clear(op1);
						mpf_clear(op2);
						return res;
					}
				}
			}
		}
		AdvancedDouble_Hybrid operator*(const double &obj1_double) const {
			const AdvancedDouble_Hybrid obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				AdvancedDouble_Hybrid res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
				mpf_mul(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else {
				double a = (*(this->smallValue)) * (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					AdvancedDouble_Hybrid res;
					res.createDouble(a);
					return res;
				}
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_mul(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);
					mpf_clear(op2);
					//if(verbose==1)printf("successful multiplication\n");
					return res;
				}
			}
		}
		AdvancedDouble_Hybrid operator+(const AdvancedDouble_Hybrid &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y'){
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_add(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
					return res;
				}
				//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_add(*(res.bigValue),*(this->bigValue), op2);
					mpf_clear(op2);
					return res;
				}
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else{
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_add(*(res.bigValue),op1, *(obj1.bigValue));
					mpf_clear(op1);
					return res;
				}
				//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
				else{//if(verbose==1)printf("operator+ AdvancedDouble_Hybrid obj1 this->isBig=='n' && obj1.isBig=='n'\n");
					double a = (*(this->smallValue)) + (*(obj1.smallValue));
					//check if a is finite
					if(isfinite(a)){
						AdvancedDouble_Hybrid res;
						res.createDouble(a);
						return res;
					}
					else{
						AdvancedDouble_Hybrid res;
						res.createBigNum();
						mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
						mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
						mpf_add(*(res.bigValue),op1, op2);
						//res.print();
						mpf_clear(op1);mpf_clear(op2);
						//if(verbose==1)printf("successful addition\n");
						return res;
					}
				}
			}
		}

		AdvancedDouble_Hybrid operator+(const double &obj1_double) const {
			const AdvancedDouble_Hybrid obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				AdvancedDouble_Hybrid res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
				mpf_add(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else {
				double a = (*(this->smallValue)) + (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					AdvancedDouble_Hybrid res;
					res.createDouble(a);
					return res;
				}
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_add(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);
					mpf_clear(op2);
					//if(verbose==1)printf("successful multiplication\n");
					return res;
				}
			}
		}
		AdvancedDouble_Hybrid operator-(const AdvancedDouble_Hybrid &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y'){
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_sub(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
					return res;
				}
				//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_sub(*(res.bigValue),*(this->bigValue), op2);
					mpf_clear(op2);
					return res;
				}
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else{
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_sub(*(res.bigValue),op1, *(obj1.bigValue));
					mpf_clear(op1);
					return res;
				}
				//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
				else{
					double a = (*(this->smallValue)) - (*(obj1.smallValue));
					//check if a is finite
					if(isfinite(a)){
						AdvancedDouble_Hybrid res;
						res.createDouble(a);
						return res;
					}
					else{
						AdvancedDouble_Hybrid res;
						res.createBigNum();
						mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
						mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
						mpf_sub(*(res.bigValue),op1, op2);
						//res.print();
						mpf_clear(op1);mpf_clear(op2);
						//if(verbose==1)printf("successful subtraction\n");
						return res;
					}
				}
			}
		}
		AdvancedDouble_Hybrid operator-(const double &obj1_double) const {
			const AdvancedDouble_Hybrid obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				AdvancedDouble_Hybrid res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
				mpf_sub(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else{
				double a = (*(this->smallValue)) - (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					AdvancedDouble_Hybrid res;
					res.createDouble(a);
					return res;
				}
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_sub(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					//if(verbose==1)printf("successful subtraction\n");
					return res;
				}
			}
		}
		AdvancedDouble_Hybrid operator/(const AdvancedDouble_Hybrid &obj1) const {
			//case 1: this object is bigValue and obj1 is also bigValue -- result is bigValue
			if(this->isBig=='y'){
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_div(*(res.bigValue),*(this->bigValue), *(obj1.bigValue));
					return res;
				}
				//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_div(*(res.bigValue),*(this->bigValue), op2);
					mpf_clear(op2);
					return res;
				}
			}
			//case 3: this object is smallValue and obj2 is bigValue -- result is bigValue
			else{
				if(obj1.isBig=='y'){
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_div(*(res.bigValue),op1, *(obj1.bigValue));
					mpf_clear(op1);
					return res;
				}
				//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
				else{
					double a = (*(this->smallValue)) / (*(obj1.smallValue));
					//check if a is finite
					if(isfinite(a)){
						AdvancedDouble_Hybrid res;
						res.createDouble(a);
						return res;
					}
					else{
						AdvancedDouble_Hybrid res;
						res.createBigNum();
						mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
						mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
						mpf_div(*(res.bigValue),op1, op2);
						//res.print();
						mpf_clear(op1);mpf_clear(op2);
						//if(verbose==1)printf("successful division\n");
						return res;
					}
				}
			}
		}
		AdvancedDouble_Hybrid operator/(const double &obj1_double) const {
			const AdvancedDouble_Hybrid obj1(obj1_double);

			//case 2: this object is bigValue and obj1 is smallValue -- result is bigValue
			if(this->isBig=='y'){
				AdvancedDouble_Hybrid res;
				res.createBigNum();
				mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
				mpf_div(*(res.bigValue),*(this->bigValue), op2);
				mpf_clear(op2);
				return res;
			}
			//case 4: this object is smallValue and obj2 is smallValue -- result can be smallValue or bigValue, we need to check
			else{
				double a = (*(this->smallValue)) / (*(obj1.smallValue));
				//check if a is finite
				if(isfinite(a)){
					AdvancedDouble_Hybrid res;
					res.createDouble(a);
					return res;
				}
				else{
					AdvancedDouble_Hybrid res;
					res.createBigNum();
					mpf_t op1; mpf_init2(op1,g_bignumprecision); mpf_set_d(op1, *(this->smallValue));
					mpf_t op2; mpf_init2(op2,g_bignumprecision); mpf_set_d(op2, *(obj1.smallValue));
					mpf_div(*(res.bigValue),op1, op2);
					//res.print();
					mpf_clear(op1);mpf_clear(op2);
					//if(verbose==1)printf("successful division\n");
					return res;
				}
			}
		}
		int compare(const AdvancedDouble_Hybrid &obj1) const{
			//Function: int mpf_cmp (mpf_t op1, mpf_t op2)
			//Function: int mpf_cmp_d (mpf_t op1, double op2)
			//case 1: this object is bigValue and obj1 is also bigValue 
			int result=100;
			if(this->isBig=='y'){
				if(obj1.isBig=='y'){
					//return mpf_cmp(*(this->bigValue), *(obj1.bigValue));
					result = mpf_cmp(*(this->bigValue), *(obj1.bigValue));
				}
				//case 2: this object is bigValue and obj1 is smallValue
				else{
					//return mpf_cmp_d(*(this->bigValue), *(obj1.smallValue));
					result =  mpf_cmp_d(*(this->bigValue), *(obj1.smallValue));
				}
			}
			//case 3: this object is smallValue and obj2 is bigValue
			else{
				if(obj1.isBig=='y'){
					//mpf_t minusOne; mpf_init2(minusOne,g_bignumprecision); mpf_set_d(minusOne, -1.0);
					//return mpf_mul(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)), minusOne);
					//return -1*(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)));
					result = -1*(mpf_cmp_d(*(obj1.bigValue), *(this->smallValue)));
				}
				//case 4: this object is smallValue and obj2 is smallValue
				else{
					//return (*(this->smallValue)) - (*(obj1.smallValue));
					double result1 = (*(this->smallValue)) - (*(obj1.smallValue));
					if(result1==0) result=0;
					else if(result1<0) result=-1;
					else result=+1;
				}
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
			else{
				//return (*(this->smallValue)) - obj1;
				double result1 = (*(this->smallValue)) - obj1;
				if(result1==0.0) result=0;
				else if(result1<0.0) result=-1;
				else result=+1;

			}
			//printf("comparing: ");print();printf(" and %f and result is %d\n",obj1,result);
			return result;
		}
		bool operator==(const AdvancedDouble_Hybrid &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator==(const double &obj1) const {
                        return compare(obj1)==0;
                }
                bool operator!=(const AdvancedDouble_Hybrid &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator!=(const double &obj1) const {
                        return compare(obj1)!=0;
                }
                bool operator<(const AdvancedDouble_Hybrid &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator<(const double &obj1) const {
                        return compare(obj1)<0;
                }
                bool operator>(const AdvancedDouble_Hybrid &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator>(const double &obj1) const {
                        return compare(obj1)>0;
                }
                bool operator<=(const AdvancedDouble_Hybrid &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator<=(const double &obj1) const {
                        return compare(obj1)<=0;
                }
                bool operator>=(const AdvancedDouble_Hybrid &obj1) const {
                        return compare(obj1)>=0;
                }
                bool operator>=(const double &obj1) const {
                        return compare(obj1)>=0;
                }

		AdvancedDouble_Hybrid& operator=(const AdvancedDouble_Hybrid &obj1) {
			//if(verbose==1){ cout<<"operator= called: obj1=";obj1.print();cout<<", this=";this->print();cout<<endl;}
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			if(this==&obj1) return *this;
			if(isInitialized())this->deallocate();
			else {bigValue=0; smallValue=0;}
			//printf("successful deallocation\n");
			if(obj1.isBig=='n'){ isBig='n'; createDouble(*(obj1.smallValue));}
			else if(obj1.isBig=='y'){isBig='y'; createBigNum(*(obj1.bigValue));}
			return *this;
		}
		AdvancedDouble_Hybrid& operator=(const double &obj1) {
			//if(verbose==1){ cout<<"operator= called: obj1="<<obj1;cout<<", this=";this->print();cout<<endl;}
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			//if(this==&obj1) return *this;
			if(isInitialized())this->deallocate();
			else {bigValue=0; smallValue=0;}
			//printf("successful deallocation\n");
			isBig='n';
			createDouble(obj1);
			return *this;
		}
		AdvancedDouble_Hybrid(const AdvancedDouble_Hybrid &obj1) {
			bigValue=0;smallValue=0;isBig='n';
			//if(verbose==1){ cout<<"Copy constructor called: obj1=";obj1.print();cout<<endl;}//cout<<", this=";if(isInitialized())this->print();else cout<<"Uninitialized,";cout<<endl;
			//printf("operator overload= starts\n");obj1.print();//printf("operator overload= ends\n");
			//if(this==&obj1) return ;//*this;
			//if(isInitialized())this->deallocate();
			//printf("successful deallocation\n");
			if(obj1.isBig=='n'){ isBig='n'; createDouble(*(obj1.smallValue));}
			else if(obj1.isBig=='y'){ isBig='y'; createBigNum(*(obj1.bigValue));}
			//return *this;
		}

};
#endif




