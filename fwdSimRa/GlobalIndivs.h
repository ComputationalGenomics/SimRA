/*

SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 08/11/2016

This is a class where each all variables which is used globally and which doesn't change with each generation is kept. 
*/
#include<iostream> 
#include<vector>
#include<algorithm>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<cmath>
#include "Events.h"
#include "ChrInfo.h"
using namespace std; 


//These are declared as extern because they are used by all classes and the main method. 
extern string diploids[];
extern vector<string> diploidVec;
extern string rsid;
extern string bases[];
extern int diploidSize;
extern string chromid; 
extern vector<string> baseVec; 
extern map <string, ChromosomeInfo > AllChromRecords;

//returns a random value of type double between the min and max boundary		
double doubleRand(double min, double max);
//returns the product of all values in a vector
double ComputeProduct(vector<double> values);
//returns the sum of all values in a vector		
double ComputeSum(vector<double> values);
//divide two floating point integers
double divide(double a,double b);
//compute square of a floating point integer		
double computeSquare(double x);
//check if the value is an integer 
bool IsInteger(float k);

//return the vector of pairs from the chromosome container    
vector<string> GetPairs(int n, vector<pair<pair<string,string>, pair<string,string> > > CC);

