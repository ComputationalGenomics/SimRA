/*

fwdSimRA: A framework for selection forward-in-time with Recombination 
Author: Aritra Bose 
Last Update: 11/04/2018

This is a class where each all variables which is used globally and which doesn't change with each generation is kept. 
*/
#include<iostream> 
#include<vector>
#include<algorithm>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<cmath>
#include<map>
#include<omp.h>
#include<gsl/gsl_rng.h>
#include "Events.h"
#include "ChrInfo.h"
//using namespace std; 


//These are declared as extern because they are used by all classes and the main method. 
extern string diploids[];
extern std::vector<string> diploidVec;
extern string rsid;
extern string bases[];
extern int diploidSize;
extern string chromid; 
extern int RecombCount;
extern int MutCount; 
extern int flag;
extern int *ps;
extern int *numsnp;
extern int *gn;
extern int numberofSNPs; 
extern int epiSNPs; 
extern int nonepiSNPs;
extern double EpiFit1; 
extern double EpiFit2; 
extern double delta;
extern int newnumSNPs; 
extern std::vector<double>PopulationFitnessTable;
extern std::vector<double>EpiFitnessTable; 
extern std::vector<int> epistatus;  
extern std::vector<std::vector<int> > newSNPs;
extern std::vector<std::pair<std::pair<string,string>, std::pair <string,string> > > ExtantHaploids;
extern std::vector<string> baseVec; 
extern std::vector<int> RandNumbers;
extern std::map<int, std::vector<std::pair<std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > >, std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > > > >  AllChromRecords;
extern double FITNESS_MAX;
extern double FITNESS_MIN;
extern double FITNESS;
extern gsl_rng** threadvec;
extern std::map< int, std::vector<std::pair<std::pair<string,std::vector<string> >, std::pair<string,std::vector<string> > > > > AllHapRecords;
extern std::map<int, std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > > MutMap;
extern std::map<int, std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > > RecombMap;
extern int flagmut;
extern int SelectedSNPID; 
double RandU(double min, double max);
extern int NUMRUN;
double GetExpTime(double L, int K);
//returns a random value of type double between the min and max boundary		
double doubleRand(double min, double max);
//returns the product of all values in a std::vector
double ComputeProduct(std::vector<double> values);
//returns the sum of all values in a std::vector		
double ComputeSum(std::vector<double> values);
//divide two floating point integers
double divide(double a,double b);
//compute square of a floating point integer		
double computeSquare(double x);

//check if the value is an integer 
bool IsInteger(float k);
double Dcalc(double x);
std::vector<double> GetLineages(std::vector<std::pair<string,std::vector<string> > > Haps, std::vector<string> ChromID);
std::pair<double,double> GetDiversity(std::vector<std::pair<string,std::vector<string> > > Haps, std::vector<string> ChromID);
void InitializeThread();
void FreeThread();
std::pair<int, std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > > RetrieveListElement(int);

//return the std::vector of std::pairs from the chromosome container    
std::vector<string> Getpairs(int n, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC);

std::vector<double> GetL (std::vector<std::pair<string,std::vector<string> > > Haps, std::vector<string> ChromID);

string firstElement( const std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > &p );

bool nnz(int i); 
