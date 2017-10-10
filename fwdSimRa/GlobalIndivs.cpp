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
#include "GlobalIndivs.h"
#include "ChrInfo.h"
#include "Events.h"
#define FITNESS_MAX 0.001
#define FITNESS_MIN -0.001

using namespace std; 
//This array is the definition of the diploids, with which we start
string diploids[] = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};
//converting them to a vector 
vector<string> diploidVec(diploids,diploids+10);
//this is for naming the SNPs
string rsid = "rs";
//Declaring the bases array
string bases[] = {"A","T","G","C"};
//Converting the bases array to a vector
vector<string> baseVec(bases,bases+4); 
int diploidSize = diploidVec.size();
//This is for naming the chromosomes 
string chromid = "chrom";
//This map stores the key-value pairs of the chromosome ID followed by their information of past-present-future
map <string, ChromosomeInfo> AllChromRecords;
double  doubleRand(double min, double max) {
        double val = (double)rand()/RAND_MAX;
        return min + val*(max-min);
}

double  ComputeProduct(vector<double> values){
	double mult = 1.0; 
	for (vector<double>::iterator it = values.begin(); it != values.end(); ++it)
		mult = mult*(*it);
		
	return mult; 
}

double  ComputeSum(vector<double> values){
        double sum = 0.0;
        for (vector<double>::iterator it = values.begin(); it != values.end(); ++it)
                sum = sum+(*it);
        return sum;
}

double  divide(double a,double b) {return a/b;}

double  computeSquare (double x) { return x*x; }

bool IsInteger(float k)
    {
        if( k == (int) k) return true;
        return false;
    }

vector<string>  GetPairs(int n, vector<pair<pair<string,string>, pair<string,string> > > CC){
			
	pair<pair<string,string>, pair<string,string> > chrompair = CC[n];
	vector<string> EachPair;
		
	EachPair.push_back(chrompair.first.first);
	EachPair.push_back(chrompair.first.second);
	EachPair.push_back(chrompair.second.first);
	EachPair.push_back(chrompair.second.second);
	return EachPair;
}