/*

SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 10/28/2017

*/
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
//Including required libraries
#include "Individuals.h"
#include "GlobalIndivs.h"
#include "Events.h"
#include "ChrInfo.h"
//Defining the boundaries of the fitness
//using namespace std;

//Defining function prototype for creating ARGs
//std::vector<std::vector<std::pair<string,int> > > createARGs(std::pair<AlleleInfo,AlleleInfo>,std::map <string, ChromosomeInfo >);
void traceARG(std::vector< std::pair<string,int> > ,int );	
int GenNum;
double FITNESS_MAX;
double FITNESS_MIN;
std::vector <std::pair<string,int> > UpdatedChroms;
//std::vector<int> ArgPop {20,50,80,120};
std::vector<int> ArgPop {4,6,8,10};
int main(int argc, char *argv[])
{
	if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " Number of SNPs (any integer greater than 2) --- Population Size (any integer greater than 1) --- Number of Generations (Usually 8*number of populations) ---  Length of Chromosome(any integer greater than 1 and more than the number of SNPs) --- Recombination Rate (any decimal less than 1, preferably lower than 0.1)  --- Mutation Rate (any decimal less than 1, preferably lower than 0.1)  --- Extant Population Size (any number between 1 to less than the population size) --- Fitness Flag (boolean, either 0 or 1) " << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }
	//Takes the number of SNPs from user
	//int numberofSNPs = atoi(argv[1]);
	
	
	//Takes the number of populations from user, considering male and female of equal number, hence 2*N.
    int popSize = atoi(argv[1]);

   	//Takes the number of Generations from user
	int GenNum = popSize*2*2*2;

	//Takes in the rate of recombination, the input has to be a probability, hence between 0 and 1
	double RecombRate = atof(argv[2]); 
	if( RecombRate < 0 || RecombRate > 1 ) {
		std::cerr << "Recombination Rate should be within 0 and 1" << std::endl; 
	}

	//Takes in the rate of mutation, the input has to be a probability, hence between 0 and 1
	double MutationRate = atof(argv[3]);

	if( MutationRate < 0 || MutationRate > 1 ){ 
		std::cerr << "Mutation Rate should be within 0 and 1" << std::endl; 

	}
	
	//Takes the length of Chromosome from user, It checks whether a valid input is used. 
	double ChromLen = atof(argv[4]);
 	
	if(IsInteger(ChromLen) != true) {
		std::cerr << "Chromosome Length should be an integer greater than 1" << std::endl; 
	}
		
	int numberofSNPs = GenNum*MutationRate*ChromLen;
	//int numberofSNPs = 30; 
	std::cout << std::endl << "Number of SNPs being used is: " << numberofSNPs << std::endl; 

 
	//Takes the number of populations to build the ARG, it has to be less than twice the number of populations 
	//int ArgPop = atoi(argv[5]);

	//if( ArgPop > popSize*2 || ArgPop < 0 ){ 
		//std::cout << "The number of populations should be strictly less than 4 times the population size" << std::endl; 
	//}
	//end Input
    int start_s = clock();
	srand((unsigned int)time(NULL));

	FITNESS = atof(argv[5]);
	std::cout << std::endl << "FITNESS is set to: " << FITNESS << std::endl; 
	
		InitializeThread();
	//std::cout << std::endl << "Thread initialized" << std::endl;
	clock_t begin = omp_get_wtime();
	//Define an array of objects of the Individuals class. Each location will contain the information pertaining to that generation. 
	//This is done each for the male and female. 
	Individuals *NewGenMale = new Individuals[GenNum];
	Individuals *NewGenFemale = new Individuals[GenNum];
	
	//Reserving space for faster allocations of the Chromosome Information. This is initialized for the base (0th) generation
	NewGenMale[0].ReserveSpace(numberofSNPs*popSize);
	NewGenFemale[0].ReserveSpace(numberofSNPs*popSize);
	
	
	//Defining types which is to be used later
	typedef std::map <int, std::vector<std::pair<std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > >, std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > > > >::iterator mapIterator;
	typedef std::pair <mapIterator,mapIterator> RangeIterator;
	typedef std::vector<std::pair<std::pair<int,bool>, std::pair<string,int> > >::iterator ChromInfoListIterator;
	typedef std::pair<string, std::pair<std::pair<string,int>,std::pair<string,int> > > ChromInfo;
	typedef std::vector<std::pair<ChromInfo,ChromInfo> > Vecofpairs;
	typedef std::map< int, std::vector<std::pair<std::pair<string,std::vector<string> >, std::pair<string,std::vector<string> > > > >::iterator HapIterator;
	typedef std::map<int, std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > >::iterator MutMapIterator;
	typedef std::map<int, std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > >::iterator RecombMapIterator;
	
	//The Parent ID for the base generation is defined to be -1, to get rid of ambiguations
	int Pid = -1;
	//std::setting up the base chromosome where there's no previous contributions
	std::pair<std::pair<string,int>,std::pair<string,int> > BaseChrom = std::make_pair(std::make_pair("0",0),std::make_pair("0",0));
	//Declaring the object of te class ChrInfo. Through this object we will store information about each chromosome 
	ChrInfo ChrObject;
	//This flag tracks the gender of the individual
	bool sexflag; 
	Vecofpairs AllChromVec(popSize*2);
	//In the following loop we std::set up the populations on the base generation. We populate all the individuals randomly with SNPs and Chromosome IDs. 
	//From 0 to Population Size, we assign Males (std::set the flag to 1) and from Population Size+1 to 2*Population Size, we assign Female . 

	for(int i = 0; i < popSize*2 ; i++){
		if(i < popSize){
			sexflag = 1;
			string chromname1 = NewGenMale[0].AddandGetChromName(0);
			string chromname2 = NewGenMale[0].AddandGetChromName(0);
			//ChromosomeInfo BoxofChrom = ChrObject.PopulateAndGetChrInfo(std::make_pair(Pid,0), 0, i, sexflag, BaseChrom);
			AllChromVec[i] = std::make_pair(std::make_pair(chromname1,BaseChrom), std::make_pair(chromname2,BaseChrom)); 	
					//std::cout << std::endl << "val of i: " << i << std::endl;
			//std::cout << std::endl << "Size of ACV: " << AllChromVec.size() << std::endl;		
			for(int j = 0; j < numberofSNPs; j++) {					
				std::pair <string,string> chrompair1 = std::make_pair(chromname1,bases[rand()%4]);
				std::pair <string,string> chrompair2 = std::make_pair(chromname2,bases[rand()%4]);
				//std::cout << std::endl << "val of j: " << j << std::endl;
				NewGenMale[0].Addpairs(std::make_pair(chrompair1,chrompair2),i*numberofSNPs+j);
				//std::cout << std::endl << "crossed male ?" << std::endl;
			}
		}
		else{
			sexflag = 0; 
			string chromname1 = NewGenFemale[0].AddandGetChromName(0);
			string chromname2 = NewGenFemale[0].AddandGetChromName(0);
			AllChromVec[i] = std::make_pair(std::make_pair(chromname1,BaseChrom), std::make_pair(chromname2,BaseChrom)); 	

			for(int j = 0; j < numberofSNPs; j++) {
				std::pair <string,string> chrompair3 = std::make_pair(chromname1,bases[rand()%4]);
				std::pair <string,string> chrompair4 = std::make_pair(chromname2,bases[rand()%4]);

				NewGenFemale[0].Addpairs(std::make_pair(chrompair3,chrompair4),(i-popSize)*numberofSNPs+j);
				
			}	
		}
	}
	AllChromRecords.insert(std::make_pair(0,AllChromVec));
	AllChromVec.clear();
	//Saving the base generation chromosome information to the containers, to be used by the next generation
	std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > M_ChromContainer = NewGenMale[0].getChromContainer();
	std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > F_ChromContainer = NewGenFemale[0].getChromContainer();
	
	

	
	//We get the fitness tables. This will remain constant throughout the analysis. If we need a new fitness table every generation, that 
	//is an easy fix.
	double *PopulationFitnessTable = NewGenMale[0].getFitnessTable(numberofSNPs, FITNESS_MAX, FITNESS_MIN);

	omp_set_dynamic(0);
	int snpmutflag = 0;

	 clock_t b2 = omp_get_wtime();
	//This marks the start of building the entire network across generations. 
	for (int g = 1; g < GenNum; g++){
		 clock_t b1 = omp_get_wtime();
		//Reserving Space for each generation's chromosome container. 
		NewGenMale[g].ReserveSpace(numberofSNPs*popSize);
		NewGenFemale[g].ReserveSpace(numberofSNPs*popSize);
		//Declaring std::vectors to calculate the probability of each individual selecting a parent
		
			
		std::vector<double> StoreProductFitness_F;
		std::vector<double> StoreProductFitness_M; 
		std::vector<double> indivprobs_M;
		std::vector<double> indivprobs_F; 
		std::vector<double> fitnessVals_M;
		std::vector<double> fitnessVals_F;
		std::vector<string> RetrievedpairFemale;
		std::vector<string> RetrievedpairMale;
		int pos_allele_F, pos_allele_M;
		double fitness_F,fitness_M;
		double revFit_F, revFit_M;
		
		//This piece of code retrieves the fitnesses based on the chromosome containers and calculates the probabilities for each individuals
		for (int i = 0; i < popSize*2; i++){
			if ( i < popSize){ 
				for ( int k = 0; k < numberofSNPs; k++){
					RetrievedpairMale = Getpairs(i*numberofSNPs + k,M_ChromContainer);
					pos_allele_M = NewGenMale[g-1].getPosAllele(RetrievedpairMale);
					fitness_M = PopulationFitnessTable[k*diploidSize + pos_allele_M];
					revFit_M = log(1.0 + fitness_M);
					fitnessVals_M.push_back(revFit_M);
				}
				StoreProductFitness_M.push_back(exp(ComputeSum(fitnessVals_M))); 
			}
			 else{
				 
				for ( int k = 0; k < numberofSNPs; k++){
					RetrievedpairFemale = Getpairs((i-popSize)*numberofSNPs + k,F_ChromContainer); 
					pos_allele_F = NewGenFemale[g-1].getPosAllele(RetrievedpairFemale);
					fitness_F = PopulationFitnessTable[k*diploidSize + pos_allele_F];
					revFit_F = log(1.0 + fitness_F);
					fitnessVals_F.push_back(revFit_F);
				}
				StoreProductFitness_F.push_back(exp(ComputeSum(fitnessVals_F))); 
			}
			
			fitnessVals_F.clear();
			fitnessVals_M.clear();
			RetrievedpairFemale.clear();
			RetrievedpairMale.clear();	
		}
		double sumprobs_F = ComputeSum(StoreProductFitness_F);
		double sumprobs_M = ComputeSum(StoreProductFitness_M);
		//std::cout << std::endl << "Sumprobs Female is " << sumprobs_F << " and Sumprobs Male is " << sumprobs_M << std::endl; 
		for (std::vector<double>::iterator it = StoreProductFitness_F.begin(); it != StoreProductFitness_F.end(); ++it){
			//std::cout << *it << " | " ;
			indivprobs_F.push_back(divide(*it,sumprobs_F));
		}
		//std::cout << std::endl << "Size of indivprobs_F: " << indivprobs_F.size() << std::endl; 
		//std::cout << std::endl;
		for (std::vector<double>::iterator it = StoreProductFitness_M.begin(); it != StoreProductFitness_M.end(); ++it){
				//std::cout << *it << " | " ;
				indivprobs_M.push_back(divide(*it,sumprobs_M));
			}
		//std::cout << std::endl << "Size of indivprobs_M: " << indivprobs_F.size() << std::endl; 
		std::vector<double>oldindivprobs_M(indivprobs_M);
		std::vector<double>oldindivprobs_F(indivprobs_F);
		//std::cout << std::endl; 
		//std::cout << "====================================================";
	
		
		transform(indivprobs_M.begin(), indivprobs_M.end(), indivprobs_M.begin(), computeSquare);	
		transform(indivprobs_F.begin(), indivprobs_F.end(), indivprobs_F.begin(), computeSquare);
					
		partial_sum ( oldindivprobs_M.begin(), oldindivprobs_M.end(), oldindivprobs_M.begin());
		partial_sum ( oldindivprobs_F.begin(), oldindivprobs_F.end(), oldindivprobs_F.begin());
	
		std::vector<int> locationinf_M = NewGenMale[g-1].GetLocSegment(numberofSNPs,ChromLen);
		std::vector<int> locationinf_F = NewGenFemale[g-1].GetLocSegment(numberofSNPs,ChromLen);
		
		std::vector < std::pair<int,int> > Parents(popSize*2);
		std::vector<std::pair<std::pair<string,std::vector<string> >, std::pair<string,std::vector<string> > > > ExtantHaploids(popSize*2);
		std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > ChrMutInfo(popSize*2);
		std::vector<std::pair<std::pair<string,int>, std::pair<string,int> > > ChrRecombInfo(popSize*2);
		std::vector<pair<int,int> > SNPid(popSize*2);
		Vecofpairs AllChromVec(popSize*2);
		//Selecting parent IDs for each individual and getting chromosome and crossover information from each parent 
		
		#pragma omp parallel for 
		for (int i = 0; i < popSize*2; i++){

			AlleleInfo fromFather,fromMother;
			std::pair <int,int>  ParentIdx;
			string chromid_F, chromid_M;
			std::vector<string> haploids_F(numberofSNPs);
			std::vector<string> haploids_M(numberofSNPs);
			//Selecting Father
			
			//double randval_M = gsl_rng_uniform(threadvec[omp_get_thread_num()]);
			double randval_M = RandU(0,1);
			int FatherIndex;
			int MotherIndex;
		
			for( std::vector<double>::iterator itt = oldindivprobs_M.begin(); itt != oldindivprobs_M.end()-1; ++itt ){
				std::vector<double>::iterator nxt1 = itt;
				++nxt1;
				if(*itt <= randval_M && *nxt1 > randval_M) {
					FatherIndex = distance(oldindivprobs_M.begin(),nxt1);
					//break;
				}
				if( randval_M < *(oldindivprobs_M.begin()) ) { 
					FatherIndex = 0;
					//break;
				}
			}	
			//double randval_F = ((double) rand() / (RAND_MAX));

			//Selecting Mother
			//double randval_F = gsl_rng_uniform(threadvec[omp_get_thread_num()]);
			double randval_F = RandU(0,1);
			//std::cout << std::endl << "Random number father is " << randval_M << " and random number mother is " << randval_F << std::endl;
			for( std::vector<double>::iterator it = oldindivprobs_F.begin(); it != oldindivprobs_F.end()-1; ++it ){
				std::vector<double>::iterator nxt = it;
				++nxt;
				if(*it <= randval_F && *nxt > randval_F){
					MotherIndex = distance(oldindivprobs_F.begin(),nxt);
					//break;
				}
				if( randval_F < *(oldindivprobs_F.begin()) ) {
					MotherIndex = 0; 
					//break;
				}
			}
			//gettimeofday(&tp, NULL);
			ParentIdx = std::make_pair(FatherIndex,MotherIndex);
			Parents[i] = ParentIdx;
			//This method returns the information from the parents, whether there was a crossover and mutation. We do this for both mother and father.
			//The contents of type AlleleInfo is discused in detail in the Events class

				fromFather = NewGenMale[g-1].AskChromosome(FatherIndex,M_ChromContainer,locationinf_M,RecombRate, MutationRate, ChromLen,g);
				fromMother = NewGenFemale[g-1].AskChromosome(MotherIndex,F_ChromContainer,locationinf_F,RecombRate, MutationRate, ChromLen,g);
			//clock_t e4 = omp_get_wtime();
			//gettimeofday(&tp, NULL);
			//std::cout << std::endl << "SNPid was " << fromFather.snpid << " and " << fromMother.snpid << std::endl; 
			SNPid[i] = std::make_pair(fromFather.snpid,fromMother.snpid); 
			haploids_M = fromFather.ToTheChild.second;
			haploids_F = fromMother.ToTheChild.second;
			chromid_F =  fromMother.ToTheChild.first;
			chromid_M =  fromFather.ToTheChild.first;
 		
			ExtantHaploids[i] = std::make_pair(std::make_pair(chromid_M,haploids_M), std::make_pair(chromid_F,haploids_F));
			AllChromVec[i] = std::make_pair(std::make_pair(chromid_M,fromFather.ContChrom),std::make_pair(chromid_F,fromMother.ContChrom));
			ChrMutInfo[i] = std::make_pair(std::make_pair(chromid_M,fromFather.MutCount),std::make_pair(chromid_F,fromMother.MutCount));
			ChrRecombInfo[i] = std::make_pair(std::make_pair(chromid_M,fromFather.RecombCount),std::make_pair(chromid_F,fromMother.RecombCount));

			if (i < popSize){
				//#pragma omp parallel for
				for (int j = 0; j < numberofSNPs; j++){
					//std::cout << std::endl << "Val of j: " << j << " Haploids_M is " <<  haploids_M[j] << " and Haploids_F is " << haploids_F[j] << std::endl;
					NewGenMale[g].Addpairs(std::make_pair(std::make_pair(fromFather.ToTheChild.first,haploids_M[j]),std::make_pair(fromMother.ToTheChild.first,haploids_F[j])),i*numberofSNPs+j);
				}
			}
			else{
				//#pragma omp parallel for
				for (int j = 0; j < numberofSNPs; j++){
					//std::cout << std::endl << "Val of j: " << j << " Haploids_M is " << haploids_M[j] << " and Haploids_F is " << haploids_F[j] << std::endl;
					NewGenFemale[g].Addpairs(std::make_pair(std::make_pair(fromFather.ToTheChild.first,haploids_M[j]),std::make_pair(fromMother.ToTheChild.first,haploids_F[j])),(i - popSize)*numberofSNPs+j);	
				}
			}
			haploids_F.clear();
			haploids_M.clear();
		}
		
		if (snpmutflag == 0){
			for (int ll = 0; ll < popSize*2; ll++){
				std::pair<int,int> idp = SNPid[ll]; 
				if (idp.first > 0){
					for (int mm = 0; mm < diploidSize ; mm++){ 
					std::string targetdiploid = diploidVec[mm];
					std::size_t ct = std::count(targetdiploid.begin(),targetdiploid.end(),'T');
					if (ct == 2)
						PopulationFitnessTable[idp.first*diploidSize+mm] = 2*FITNESS;
					else if (ct == 1)
						PopulationFitnessTable[idp.first*diploidSize+mm] = FITNESS; 
					else 
						PopulationFitnessTable[idp.first*diploidSize+mm] = 0;	
					}
					SelectedSNPID = -1; 
					snpmutflag = 1;
					break;
				}
				if (idp.second > 0){
					for (int mm = 0; mm < diploidSize ; mm++){ 
					std::string targetdiploid = diploidVec[mm];
					std::size_t ct = std::count(targetdiploid.begin(),targetdiploid.end(),'T');
					if (ct == 2)
						PopulationFitnessTable[idp.second*diploidSize+mm] = 2*FITNESS;
					else if (ct == 1)
						PopulationFitnessTable[idp.second*diploidSize+mm] = FITNESS; 
					else 
						PopulationFitnessTable[idp.second*diploidSize+mm] = 0;	
					}
					SelectedSNPID = -1; 
					snpmutflag = 1;
					break;
				}
			}
		}
		
		std::cout<<std::endl;
		SNPid.clear();
		M_ChromContainer.clear();
		F_ChromContainer.clear();
		M_ChromContainer.resize(numberofSNPs*popSize);
		F_ChromContainer.resize(numberofSNPs*popSize);
		M_ChromContainer = NewGenMale[g].getChromContainer();
		F_ChromContainer = NewGenFemale[g].getChromContainer();
		StoreProductFitness_M.clear();
		StoreProductFitness_F.clear();
		indivprobs_F.clear();
		indivprobs_M.clear();
		AllChromRecords.insert(std::make_pair(g,AllChromVec));
		AllHapRecords.insert(std::make_pair(g,ExtantHaploids));
		ExtantHaploids.clear();
		ExtantHaploids.resize(popSize*2);
		MutMap.insert(std::make_pair(g,ChrMutInfo));
		ChrMutInfo.clear();
		ChrMutInfo.resize(popSize*2);
		AllChromVec.clear();
		AllChromVec.resize(popSize*2);
		clock_t e1 = omp_get_wtime();
		double et2 = double(e1-b1);
		std::cout << std::endl << "Time taken to compute Generation # " << g << " is " << et2 << " secs" << std::endl; 
	
	}
	clock_t e2 = omp_get_wtime();
	FreeThread();

	double et1 = double(e2 - b2);
	std::cout << std::endl << "Time taken to run the Book of Populations: " << et1 <<  " secs" << std::endl;
	printf("\n -------------Total Number of Recombinations-------------\n %d\n\n",RecombCount);
	printf("\n ----------------Total Number of Mutations--------------\n %d\n\n",MutCount);


	std::string strPop = std::to_string(popSize*4);
	std::string strChrLen = std::to_string(int(ChromLen));
	std::string reco = std::to_string(RecombRate);
	std::string mutu = std::to_string(MutationRate);
	std::string fit = std::to_string(FITNESS); 

	for (int cd = 0; cd < ArgPop.size(); cd++){
		//std::cout << std::endl << "****** Starting Simulation " << ab << "******" << std::endl;
		vector<double> Div(NUMRUN);
		vector<int> Rec(NUMRUN);
		vector<int> Mut(NUMRUN);
		vector<int> HT(NUMRUN);
		vector<double> TtoG(NUMRUN);
		std::vector<double> ExtCount;
		std::string strARGpop = std::to_string(ArgPop[cd]);
		string FName = "N"+strPop+"_g"+strChrLen+"r_"+reco+"mu_"+mutu+"_sel"+"_m"+strARGpop+"fit_"+fit+".txt"; 
		HapIterator iter;
		double NumLineages;
		
		std::vector<std::pair<std::pair<string,std::vector<string> > , std::pair<string,std::vector<string> > > > ExtantHaps(popSize*2);
		std::vector<std::pair<string,std::vector<string> > > haploo;
		for (int ab = 0; ab < NUMRUN; ab++){
			int ArgRecomb = 0;
			int ArgMut = 0;
			int ht=0; 		
			iter = AllHapRecords.find(GenNum-1);
			ExtantHaps = iter->second; 
			//std::vector <double> TimeVec;
			std::vector<string> ExHap(popSize*2);
			for (std::vector<std::pair<std::pair<string,std::vector<string> >,std::pair<string,std::vector<string> > > >::iterator it = ExtantHaps.begin(); it != ExtantHaps.end(); ++it){
				haploo.push_back(it->first);
				haploo.push_back(it->second);
			}
			int *RandNumbers = new int[popSize*2];
			std::vector<int> SelectedPopforArg;
			for (int i = 0; i < popSize*2; i++) 
				RandNumbers[i] = i; 
			//random_shuffle(RandNumbers,RandNumbers+(popSize*2));
			std::random_device rd; 
			std::mt19937 g(rd());
			shuffle(RandNumbers,RandNumbers+(popSize*2), g);
				//std::cout << "RandNumbers: " << std::endl;
			for (int i = 0; i < ArgPop[cd]; i++){
				SelectedPopforArg.push_back(RandNumbers[i]);
				//std::cout << RandNumbers[i] << " "; 
			}
			

			MutMapIterator xter = MutMap.find(GenNum-1);
			std::vector<std::pair<std::pair<string,int>,std::pair<string,int> > > MutChr = xter->second;
			std::vector<std::pair<string,int> > MC; 
			for(std::vector<std::pair<std::pair<string,int>,std::pair<string,int> > >::iterator itm = MutChr.begin(); itm != MutChr.end(); ++itm){
				MC.push_back(itm->first);
				MC.push_back(itm->second);
			}
	
			mapIterator itr = AllChromRecords.find(GenNum-1);
			std::vector<std::pair<ChromInfo,ChromInfo> > LeafVec;
			//int c = 0;
			if(itr != AllChromRecords.end()){
				LeafVec = itr->second;
				
			}
		
			std::vector<std::pair<ChromInfo,ChromInfo> > SelectedChroms;
			for (int i = 0; i < SelectedPopforArg.size(); i++){
		
				SelectedChroms.push_back(LeafVec[SelectedPopforArg[i]]);	
			}
			std::vector<ChromInfo> NewLeafVec;
			
			std::vector<string> ChromIds;
			std::vector<string> ExtantChromIDs;
			double rnd;
					
			for (std::vector<std::pair<ChromInfo,ChromInfo> >::iterator it = SelectedChroms.begin(); it != SelectedChroms.end(); ++it){
				rnd = ((double) rand()/(double) (RAND_MAX));
				if (rnd <= 0.5){
					ExtantChromIDs.push_back(it->first.first);
					if(it->first.second.first.second != 0){
						if (it->first.second.first.second != numberofSNPs)
							ArgRecomb++;
						ChromIds.push_back(it->first.second.first.first);
					}
					if(it->first.second.second.second != 0){
						if(it->first.second.second.second != numberofSNPs)
							ArgRecomb++;
						ChromIds.push_back(it->first.second.second.first);	
					}
				}
				else {
					ExtantChromIDs.push_back(it->second.first);
					if(it->second.second.first.second != 0){
						if (it->second.second.first.second != numberofSNPs)
							ArgRecomb++;
						ChromIds.push_back(it->second.second.first.first);	
					}
					if(it->second.second.second.second != 0){
						if (it->second.second.second.second != numberofSNPs)
							ArgRecomb++;
						ChromIds.push_back(it->second.second.second.first);	
					}
				}
			}
			
			for (std::vector<string>::iterator it = ExtantChromIDs.begin(); it != ExtantChromIDs.end(); ++it){
				string SelChrID = *it;
				//std::cout << std::endl << "Selected Chromosome IDs are: " << SelChrID << std::endl;
				for (std::vector<std::pair<string,int> >::iterator it1 = MC.begin(); it1 != MC.end(); ++it1){
					if (SelChrID.compare(it1->first) == 0)
						ArgMut += it1->second;
				}
		
			}
			
			double D = GetDiversity(haploo, ExtantChromIDs);

			ExtCount.clear();
			haploo.clear();
			ExtantHaps.clear();
			std::set <string> Chrset(ChromIds.begin(),ChromIds.end());
			ChromIds.clear();
			ChromIds.assign(Chrset.begin(),Chrset.end());
			Chrset.clear();
			
			if (ChromIds.size() > 1){
				flag = 0;
			}
			else{
				std::cout << "GMRCA Reached and Depth of ARG is 1" << std::endl;
				flag = 1;
			}
			if (flag == 0){
				std::vector<string> FirstChrIds;
				std::vector<string> ChrIds;
				std::vector<string> NewChromIds;
				for (int g = GenNum-2; g > 0; g--){
					//std::cout << std::endl << " ---- Starting Generation " << g << " ----" << std::endl;
					LeafVec.clear();
					NewLeafVec.clear();
					MC.clear();
					//std::cout << std::endl << "Ht: " << ht << std::endl;
					//Rect.clear();
					MutChr.clear();
					mapIterator itr = AllChromRecords.find(g);
					MutMapIterator mtr = MutMap.find(g);
					MutChr = mtr->second;		
					for(std::vector<std::pair<std::pair<string,int>,std::pair<string,int> > >::iterator itm = MutChr.begin(); itm != MutChr.end(); ++itm){
						MC.push_back(itm->first);
						MC.push_back(itm->second);
					}
				
					iter = AllHapRecords.find(g);
					ExtantHaps = iter->second; 
					for (std::vector<std::pair<std::pair<string,std::vector<string> >,std::pair<string,std::vector<string> > > >::iterator itc = ExtantHaps.begin(); itc != ExtantHaps.end(); ++itc){
						haploo.push_back(itc->first);
						haploo.push_back(itc->second);
					}
					
					LeafVec = itr->second;
						for (std::vector<std::pair<ChromInfo,ChromInfo> >::iterator it = LeafVec.begin(); it != LeafVec.end(); ++it){
								NewLeafVec.push_back(it->first);
								NewLeafVec.push_back(it->second);
							}
					transform(NewLeafVec.begin(), NewLeafVec.end(), back_inserter(FirstChrIds), firstElement);
					
					for (int k = 0; k < ChromIds.size(); k++){
						int fidx = find(FirstChrIds.begin(),FirstChrIds.end(), ChromIds[k]) - FirstChrIds.begin();
						 //std::cout << std::endl << "Main Element: " << fidx << std::endl;
							if (NewLeafVec[fidx].second.first.second != 0){
								if(NewLeafVec[fidx].second.first.second != numberofSNPs)
									ArgRecomb++;
								//std::cout << std::endl << NewLeafVec[fidx].second.first.second << std::endl;
								NewChromIds.push_back(NewLeafVec[fidx].second.first.first);
								for (std::vector<std::pair<string,int> >::iterator it1 = MC.begin(); it1 != MC.end(); ++it1){
									if ((NewLeafVec[fidx].second.first.first).compare(it1->first) == 0)
										ArgMut += it1->second;
								}
															
							}
							if(NewLeafVec[fidx].second.second.second != 0){
								if(NewLeafVec[fidx].second.second.second != numberofSNPs)
									ArgRecomb++;
								//std::cout << std::endl << NewLeafVec[fidx].second.second.second << std::endl;
								NewChromIds.push_back(NewLeafVec[fidx].second.second.first);
								for (std::vector<std::pair<string,int> >::iterator it1 = MC.begin(); it1 != MC.end(); ++it1){
									if ((NewLeafVec[fidx].second.second.first).compare(it1->first) == 0)
										ArgMut += it1->second;
								}
						
						}
					}
					
					ChromIds.clear();
					std::set <string> Chrset(NewChromIds.begin(),NewChromIds.end());
					NewChromIds.clear();
					ChromIds.assign(Chrset.begin(),Chrset.end());
					
					FirstChrIds.clear();
					//std::cout << ChromIds.size() << std::endl;
					if(ChromIds.size() > 1){
					
					}
					else{
						
						ht = (GenNum-g+1);
						flag = 1; 
						break;
					}	
				haploo.clear();
				ExtCount.clear();
				ExtantHaps.clear();
				}
			}
	
			if (flag == 0) 
			{
				std::cout << std::endl << "-----------GMRCA not reached-----------" << std::endl;
					ht = GenNum;
			}
			
			SelectedChroms.clear();
			ExtantChromIDs.clear();
			ChromIds.clear();
			SelectedPopforArg.clear();
			MC.clear();
			LeafVec.clear();
			NewLeafVec.clear();
			Div[ab] = D;
			Rec[ab] = ArgRecomb;
			Mut[ab] = ArgMut;
			HT[ab] = ht;
			TtoG[ab] = (double)ht/(4*popSize);
		}
	
		std::ofstream myfile1(FName);
		if(myfile1.is_open()){
			myfile1 << "ARG Height" << "\t" << "# of Recombinations" << "\t" << "# of Mutations" << "\t" << "Time to GMRCA" << "\t" << "Diversity" << std::endl;
			for (int x = 0; x < NUMRUN; x++)
				myfile1<< HT[x] << "\t" << Rec[x] << "\t" << Mut[x] << "\t" << TtoG[x] << "\t" << Div[x] << std::endl;			
			myfile1.close();
		}
		else cout << "Unable to open file";
		Div.clear();
		Rec.clear();
		Mut.clear();
		HT.clear();
		TtoG.clear();
		//RR.clear();
	}

		
	clock_t end = omp_get_wtime();
		
	double elapsed_time = double(end - begin);
	std::cout << std::endl << "Time taken to run: " << elapsed_time <<  " secs" << std::endl;
}	



	
