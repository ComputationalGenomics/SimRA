#include<iostream>
#include<string.h>
#include<vector>
#include<set>
#include<algorithm>
#include<fstream>
#include<cstdlib>
#include<time.h>
#include<iterator>
#include<sys/time.h>
#include<stdio.h>
#include<numeric>
#include <stdexcept>
#include<math.h>
#include<limits>
#include<random>
#include<ctime>
#include<map>
#include<omp.h>
#include <gsl/gsl_rng.h>
#include<pthread.h>
#include<chrono>

extern int GenNum; 
extern std::vector<int> ArgPop;
//double PopulationFitnessTable = new double[numberofSNPs*diploidSize];
extern int popSize; 
extern double RecombRate; 
extern double MutationRate; 
extern int m; 
extern int epiflag; 
extern int ChromLen; 