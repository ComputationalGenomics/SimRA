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
#include<omp.h>
#include<map>
#include<time.h>
#include<random>
#include<set>
#include<cmath>
#include<gsl/gsl_rng.h>
#include "GlobalIndivs.h"
#include "ChrInfo.h"
#include "Events.h"


//using namespace std; 
//This array is the definition of the diploids, with which we start
string diploids[] = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};
//converting them to a std::vector 
std::vector<string> diploidVec(diploids,diploids+10);
//this is for naming the SNPs
string rsid = "rs";
//Declaring the bases array
string bases[] = {"A","T","G","C"};
//Converting the bases array to a std::vector
std::vector<string> baseVec(bases,bases+4); 
int diploidSize = diploidVec.size();
//gsl_rng ** threadvec = new gsl_rng*[omp_get_num_threads()];
gsl_rng** threadvec = new gsl_rng*[omp_get_max_threads()];
//This is for naming the chromosomes 
string chromid = "chrom";
std::vector<int> epistatus; 
//This std::map stores the key-value std::pairs of the chromosome ID followed by their information of past-present-future
//std::map <string, ChromosomeInfo> AllChromRecords;
int NUMRUN;
int RecombCount = 0;
int MutCount = 0; 
double FITNESS;
double EpiFit1;
double EpiFit2;
double delta; 
int flag = 0;
int numberofSNPs,newnumSNPs; 
int epiSNPs; 
int nonepiSNPs;
//int GenNum; 
int flagmut = 0;
int SelectedSNPID = -1; 
std::map<int, std::vector<std::pair<std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > >, std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > > > >  AllChromRecords;
std::map< int, std::vector<std::pair<std::pair<string,std::vector<string> >, std::pair<string,std::vector<string> > > > > AllHapRecords;
std::map<int, std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > > MutMap;
std::map<int, std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > > RecombMap;
void InitializeThread(){
	
	gsl_rng_env_setup();		
	for (int b = 0; b < omp_get_max_threads(); b++){
		threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
		gsl_rng_set(threadvec[b],b*clock());
	}
	
	std::cout << "number of threads : " << omp_get_max_threads() << std::endl;
	
}
void FreeThread(){
		for (int b = 0; b < omp_get_max_threads(); b++)
			gsl_rng_free(threadvec[b]);
		
		//free(threadvec);
		//delete[] threadvec; 
}
double RandU(double min, double max) {
    thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}		
double  doubleRand(double min, double max) {
        double val = (double)rand()/RAND_MAX;
        return min + val*(max-min);
}

double  ComputeProduct(std::vector<double> values){
	double mult = 1.0; 
	for (std::vector<double>::iterator it = values.begin(); it != values.end(); ++it)
		mult = mult*(*it);
		
	return mult; 
}

double  ComputeSum(std::vector<double> values){
        double sum = 0.0;
        for (std::vector<double>::iterator it = values.begin(); it != values.end(); ++it)
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

double Dcalc(double x){
	 return x*x;
	}
	
std::pair<double,double> GetDiversity(std::vector<std::pair<string,std::vector<string> > > Haps, std::vector<string> ChromID){
	std::vector<string> SelectedSNPs(ChromID.size()); 
	//std::string allHaps = "";
	for (int bb = 0; bb < ChromID.size(); bb++){
			string SelChrID = ChromID[bb];
			for (std::vector<std::pair<string,std::vector<string> > >::iterator it1 = Haps.begin(); it1 != Haps.end(); ++it1){
				if (SelChrID.compare(it1->first) == 0){
					string strhaps = "";
					std::vector<string> tmphaps = it1->second; 
					for (int z = 0; z < tmphaps.size(); z++)
							strhaps += tmphaps[z];
							//std::cout << std::endl << strhaps << std::endl;
					SelectedSNPs[bb] = strhaps;
					tmphaps.clear();
				}else
					continue;
			}
		//allHaps += SelectedSNPs[bb]; 
		std::cout << std::endl << "Selected SNP was: " << SelectedSNPs[bb] << std::endl;

	}
	int selall = 0;
	std::vector<double> SNPsums(numberofSNPs);
	for( int i = 0; i < numberofSNPs; i++){
		
		int sumAT = 0, sumCG = 0,maxsum; 
		double pj;
		for( int j = 0; j < ChromID.size(); j++){
			string strhaps = SelectedSNPs[j];
			if(strhaps[i] == 'A' || strhaps[i] == 'T')
				sumAT++; 
			if(strhaps[i] == 'C'|| strhaps[i] == 'G')
				sumCG++;
			if(i == SelectedSNPID){
				if(strhaps[i] == 'T')
					selall++;
			}
		}
		if(sumAT > sumCG)
			maxsum = sumAT; 
		else
			maxsum = sumCG; 
		std::cout << "Maxsum is: " << maxsum << std::endl; 
		pj = maxsum/(double)ChromID.size(); 
		std::cout << "pj is: " << pj<< std::endl;
		SNPsums[i] = 2*pj*(1-pj);
	}
	std::cout << std::endl << "Selected SNP id is " << SelectedSNPID << " || Selected allele is in " << selall << " extant indivs" << std::endl;
	double pundersel = (double)selall/(double)ChromID.size();
	double sum_of_elems = 0.0;
	for(std::vector<double>::iterator it = SNPsums.begin(); it != SNPsums.end(); ++it){
		sum_of_elems += *it;
		std::cout << *it << " "; 
	}
	std::cout << std::endl;

	return std::make_pair(sum_of_elems,pundersel); 

	//for (int i = 0; ChromID.size(); i++)
		//allHaps += SelectedSNPs[i];

	/*double accu = 0;
	for (int i = 0; i < allHaps.size(); i++){
		double n = std::count(allHaps.begin(), allHaps.end(), allHaps[i]);
		double ct = n/allHaps.size();
		accu += ((ct/allHaps.size())*(1-(ct/allHaps.size())));
	}
	//double normacc = (accu)/allHaps.size();
	return accu;*/
}
std::vector<double> GetL (std::vector<std::pair<string,std::vector<string> > > Haps, std::vector<string> ChromID){
	std::vector<string> SelectedSNPs(ChromID.size()); 
	for (int bb = 0; bb < ChromID.size(); bb++){
			string SelChrID = ChromID[bb];
			for (std::vector<std::pair<string,std::vector<string> > >::iterator it1 = Haps.begin(); it1 != Haps.end(); ++it1){
				if (SelChrID.compare(it1->first) == 0){
					string strhaps = "";
					std::vector<string> tmphaps = it1->second; 
					for (int z = 0; z < tmphaps.size(); z++)
							strhaps += tmphaps[z];
						//std::cout << std::endl << strhaps << std::endl;
					SelectedSNPs[bb] = strhaps;
					tmphaps.clear();
				}else
					continue;
			}
	}
	//if ( flag == 1 ){
	//for (int ib = 0; ib < SelectedSNPs.size(); ib++){
	//	std::cout << "vector Element: " << ib << " is " << SelectedSNPs[ib] << std::endl;
	//}
	//}
	std::vector<string> ualvec;
	std::set<string> uals(SelectedSNPs.begin(),SelectedSNPs.end());
	ualvec.assign(uals.begin(),uals.end());
	std::vector<double> ExtCount;
	for(std::vector<string>::iterator ii = ualvec.begin(); ii != ualvec.end(); ++ii){
		double myc = std::count(SelectedSNPs.begin(),SelectedSNPs.end(),*ii);
		ExtCount.push_back(myc);
	}
	return ExtCount;
}
std::vector <double> GetLineages(std::vector<std::pair<string,std::vector<string> > > Haps, std::vector<string> ChromID){
	std::vector<string> SelectedSNPs; 
	for (std::vector<string>::iterator it = ChromID.begin(); it != ChromID.end(); ++it){
				string SelChrID = *it;
				for (std::vector<std::pair<string,std::vector<string> > >::iterator it1 = Haps.begin(); it1 != Haps.end(); ++it1){
				if (SelChrID.compare(it1->first) == 0){
						string strhaps = "";
						std::vector<string> tmphaps = it1->second; 
						for (int z = 0; z < tmphaps.size(); z++)
								strhaps += tmphaps[z];
						//std::cout << std::endl << strhaps << std::endl;
						SelectedSNPs.push_back(strhaps);
						tmphaps.clear();
				}else
						continue;
				}
			}
	std::vector<double> ExtCount;
	std::vector<string> AlreadySeen;
	std::vector<string> ualvec;
	std::vector<string>  SelSNPs = SelectedSNPs;
	for (int ii = 0; ii < SelSNPs.size(); ++ii){
	//	std::cout << "std::vector Element: " << ii << " is " << SelectedSNPs[ii] << std::endl; 
		int lincount = 1;
		if(find(ualvec.begin(),ualvec.end(),SelectedSNPs[ii]) == ualvec.end()){
			for (int jj = ii+1; jj < SelectedSNPs.size(); ++jj){
				if(SelectedSNPs[jj].compare(SelSNPs[ii]) == 0){
					lincount++;
					AlreadySeen.push_back(SelSNPs[ii]);
				}
				else
					continue;
			}
		}
		else
			continue;
		std::set <string> uals(AlreadySeen.begin(),AlreadySeen.end());
		ualvec.assign(uals.begin(),uals.end());
		ExtCount.push_back(lincount);
	}
	return ExtCount;
}
double GetExpTime(double NL, int popSize){
	return (4*popSize)*(1 - (1/NL));
}
std::vector<string>  Getpairs(int n, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC){
			
	std::pair<std::pair<string,string>, std::pair<string,string> > chrompair = CC[n];
	std::vector<string> Eachpair;
		
	Eachpair.push_back(chrompair.first.first);
	Eachpair.push_back(chrompair.first.second);
	Eachpair.push_back(chrompair.second.first);
	Eachpair.push_back(chrompair.second.second);
	return Eachpair;
}

string firstElement( const std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > &p ) {
    return p.first;
}
bool nnz(int i) {return (i != 0); }
