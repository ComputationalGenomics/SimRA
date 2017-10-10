/*

SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 08/11/2016

This is the main file which executes the entire simulation.

1. Ask for Input from the user for input parameters Number of SNPs, population Size, Number of Generations (preferable 4*population Size),Length of 
a chromosome, Rate of Recombination and Rate of Mutation and the number of populations on which the Ancestral Recombination Graph is built. 

2. It first sets up a Fitness table which remains constant across all generations and populate the base generation with randomly assigned SNPs and 
corresponding fitnesses from the lookup table. 

3. Then, for each generation, it computes the effective population size, the parents from the previous generation and the chromosome information that
gets passed along each generation. It does this until it reaches the number of generations fed in by the user. 

4. It then randomly selects the number of populations (input by user) out of all the leaf nodes and creates the ARG for each chromosome. Each individual
has two chromosomes and it traces them back across the generations to reach the ancestor for each chromosome. It keeps a track of the recombinations,
as it goes back. 

5. It then matches the links for each chromosome with each other to find the Greatest Most Recent Common Ancestor (GMRCA) and notes the depth. 

***********************************************
This is a code still in progress, hence you will come across many commented sections which are used for debugging purposes and I am still working 
on points 4 and 5.
***********************************************
*/
#include<iostream>
#include<string.h>
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<time.h>
#include<iterator>
#include<stdio.h>
#include<numeric>
#include <stdexcept>
#include<math.h>
#include<limits>
#include<map>
//Including required libraries
#include "Individuals.h"
#include "GlobalIndivs.h"
#include "Events.h"
#include "ChrInfo.h"
//Defining the boundaries of the fitness
#define FITNESS_MAX 0.001
#define FITNESS_MIN -0.001
using namespace std;

//Defining function prototype for creating ARGs
vector<vector<pair<string,int> > > createARGs(pair<AlleleInfo,AlleleInfo>,map <string, ChromosomeInfo >);
	
int main()
{
	//Takes the number of SNPs from user
	int numberofSNPs; 
    cout << "Enter the number of SNPs: ";
    cin >> numberofSNPs;
    cout << endl;
	
	//Takes the number of populations from user, 2*population Size is the total population size, for equal number of male and females. 
    int popSize;
    cout << "Enter the population size for a generation: ";
    cin >> popSize;
    cout << endl;

	//Takes the number of Generations from user
	int GenNum; 
	cout << "Enter the Number of Generations: ";
	cin >> GenNum;
	cout << endl;
	
	//Takes the length of Chromosome from user, It checks whether a valid input is used. 
	float ChromLen;
    cout << "Enter length of the chromosome: ";
    cin >> ChromLen;
	while(IsInteger(ChromLen) != true) {
		cout << "Chromosome Length should be an integer greater than 1" << endl; 
		cout << "Enter length of chromosome again: "; 
		cin >> ChromLen;
	}
    cout << endl;
	
	//Takes in the rate of recombination, the input has to be a probability, hence between 0 and 1
	double RecombRate; 
	cout << "Enter probability of recombination: ";
	cin >> RecombRate;
	while( RecombRate < 0 || RecombRate > 1 ) {
		cout << "Recombination Rate should be within 0 and 1" << endl; 
		cout << "Enter Probability of recombination again: "; 
		cin >> RecombRate; 
	}
	cout << endl;

	//Takes in the rate of mutation, the input has to be a probability, hence between 0 and 1
	double MutationRate;
	cout << "Enter Mutation Rate: "; 
	cin >> MutationRate; 
	while( MutationRate < 0 || MutationRate > 1 ){ 
		cout << "Mutation Rate should be within 0 and 1" << endl; 
		cout << "Enter Mutation Rate again: ";
		cin >> MutationRate;
	}
	cout << endl; 
 
	//Takes the number of populations to build the ARG, it has to be less than twice the number of populations 
	int ArgPop;
    cout << "Enter the number of populations to build the ARG: ";
    cin >> ArgPop;
	while( ArgPop > popSize*2 || ArgPop < 0 ){ 
		cout << "The number of populations should be strictly less than twice the population size" << endl; 
		cin >> ArgPop;
	}
    cout << endl;
	//end Input
    int start_s = clock();
	srand((unsigned int)time(NULL));
	
	//Define an array of objects of the Individuals class. Each location will contain the information pertaining to that generation. 
	//This is done each for the male and female. 
	Individuals *NewGenMale = new Individuals[GenNum];
	Individuals *NewGenFemale = new Individuals[GenNum];
	
	//Reserving space for faster allocations of the Chromosome Information. This is initialized for the base (0th) generation
	NewGenMale[0].ReserveSpace(numberofSNPs*popSize);
	NewGenFemale[0].ReserveSpace(numberofSNPs*popSize);
	
	//Defining types which is to be used later
	typedef map<string, ChromosomeInfo>::iterator MapIterator;
	typedef vector<pair<pair<int,bool>, pair<string,int> > >::iterator ChromInfoListIterator;
	
	//The Parent ID for the base generation is defined to be -1, to get rid of ambiguations
	int Pid = -1;
	//Setting up the base chromosome where there's no previous contributions
	pair<pair<string,int>,pair<string,int> > BaseChrom = make_pair(make_pair("0",0),make_pair("0",0));
	//Declaring the object of te class ChrInfo. Through this object we will store information about each chromosome 
	ChrInfo ChrObject;
	//This flag tracks the gender of the individual
	bool sexflag; 
	//In the following loop we set up the populations on the base generation. We populate all the individuals randomly with SNPs and Chromosome IDs. 
	//From 0 to Population Size, we assign Males (set the flag to 1) and from Population Size+1 to 2*Population Size, we assign Female . 
	for(int i = 0; i < popSize*2 ; i++){
		if(i < popSize){
			sexflag = 1;
			string chromname1 = NewGenMale[0].AddandGetChromName(0);
			string chromname2 = NewGenMale[0].AddandGetChromName(0);
			ChromosomeInfo BoxofChrom = ChrObject.PopulateAndGetChrInfo(make_pair(Pid,0), 0, i, sexflag, BaseChrom);
			AllChromRecords.insert(make_pair(chromname1,BoxofChrom));
			AllChromRecords.insert(make_pair(chromname2,BoxofChrom));		
			for(int j = 0; j < numberofSNPs; j++) {					
				pair <string,string> chrompair1 = make_pair(chromname1,bases[rand()%4]);
				pair <string,string> chrompair2 = make_pair(chromname2,bases[rand()%4]);
				NewGenMale[0].Addpairs(make_pair(chrompair1,chrompair2));
			}
		}
		else{
			sexflag = 0; 
			string chromname1 = NewGenFemale[0].AddandGetChromName(0);
			string chromname2 = NewGenFemale[0].AddandGetChromName(0);
			ChromosomeInfo BoxofChrom = ChrObject.PopulateAndGetChrInfo(make_pair(Pid,1), 0, i-popSize, sexflag, BaseChrom);
			AllChromRecords.insert(make_pair(chromname1,BoxofChrom));
			AllChromRecords.insert(make_pair(chromname2,BoxofChrom));			
			for(int j = 0; j < numberofSNPs; j++) {
				pair <string,string> chrompair3 = make_pair(chromname1,bases[rand()%4]);
				pair <string,string> chrompair4 = make_pair(chromname2,bases[rand()%4]);
				NewGenFemale[0].Addpairs(make_pair(chrompair3,chrompair4));
			}	
		}
	}
	//Saving the base generation chromosome information to the containers, to be used by the next generation
	vector<pair<pair<string,string>, pair<string,string> > > M_ChromContainer = NewGenMale[0].getChromContainer();
	vector<pair<pair<string,string>, pair<string,string> > > F_ChromContainer = NewGenFemale[0].getChromContainer();
	

	//SNP labels, this is not relevant at this point, hence, commented. 
	/*
	for (int i = 0; i < numberofSNPs ; i++){
		
		NewGenMale[0].Addlabel();
		NewGenFemale[0].Addlabel();
	}*/
	
	//We get the fitness tables. This will remain constant throughout the analysis. If we need a new fitness table every generation, that 
	//is an easy fix.
	double *PopulationFitnessTable = NewGenMale[0].getFitnessTable(numberofSNPs);
	
	//This commented out piece of code is just for printing the chromosome containers and for debugging purposes.
	/*
	NewGenMale[0].DisplayTable(PopulationFitnessTable,numberofSNPs);
	cout << "------------------------------------------------------------" <<endl;
	NewGenFemale[0].DisplayTable(PopulationFitnessTable,numberofSNPs);

	cout <<endl << endl << "-----------------------------  MALE -------------------------------------" << endl;
	NewGenMale[0].DisplayContents();
	cout << endl;
	NewGenMale[0].DisplayLabels();
	cout << endl;
	cout << "----------------------------- FEMALE -------------------------------------" << endl;
    NewGenFemale[0].DisplayContents();
	NewGenFemale[0].DisplayLabels();
	cout << endl; 
	*/

	//This marks the start of building the entire network across generations. 
	for (int g = 1; g < GenNum; g++){
		//Reserving Space for each generation's chromosome container. 
		NewGenMale[g].ReserveSpace(numberofSNPs*popSize);
		NewGenFemale[g].ReserveSpace(numberofSNPs*popSize);
		//Declaring vectors to calculate the probability of each individual selecting a parent
		vector<double> StoreProductFitness_F;
		vector<double> StoreProductFitness_M; 
		vector<double> indivprobs_M;
		vector<double> indivprobs_F; 
		//This piece of code retrieves the fitnesses based on the chromosome containers and calculates the probabilities for each individuals
		for (int i = 0; i < popSize; i++){
		
			vector<double> fitnessVals_M;
			vector<double> fitnessVals_F;	
			for ( int k = 0; k < numberofSNPs; k++){
			
				vector<string> RetrievedPairMale = GetPairs(i*numberofSNPs + k,M_ChromContainer);
				vector<string> RetrievedPairFemale = GetPairs(i*numberofSNPs + k,F_ChromContainer); 
				int pos_allele_M = NewGenMale[g-1].getPosAllele(RetrievedPairMale);
				int pos_allele_F = NewGenFemale[g-1].getPosAllele(RetrievedPairFemale);
				
				double fitness_M = PopulationFitnessTable[k*diploidSize + pos_allele_M];
				double fitness_F = PopulationFitnessTable[k*diploidSize + pos_allele_F];
		
				double revFit_M = log(1.0 + fitness_M);
				double revFit_F = log(1.0 + fitness_F);
			
				fitnessVals_M.push_back(revFit_M);
				fitnessVals_F.push_back(revFit_F);
			}

			//This stores the product of S_ij's 
			StoreProductFitness_M.push_back(exp(ComputeSum(fitnessVals_M)));
			StoreProductFitness_F.push_back(exp(ComputeSum(fitnessVals_F))); 
		
		}
	
	
		double sumprobs_F = ComputeSum(StoreProductFitness_F);
		double sumprobs_M = ComputeSum(StoreProductFitness_M);
	
		for (vector<double>::iterator it = StoreProductFitness_F.begin(); it != StoreProductFitness_F.end(); ++it){
			//cout << *it << " | " ;
			indivprobs_F.push_back(divide(*it,sumprobs_F));
		}
		cout << endl;
		for (vector<double>::iterator it = StoreProductFitness_M.begin(); it != StoreProductFitness_M.end(); ++it){
				//cout << *it << " | " ;
				indivprobs_M.push_back(divide(*it,sumprobs_M));
			}

	
		vector<double>oldindivprobs_M(indivprobs_M);
		vector<double>oldindivprobs_F(indivprobs_F);
		//cout << endl; 
		//cout << "====================================================";
	
		//cout << endl;
		
		transform(indivprobs_M.begin(), indivprobs_M.end(), indivprobs_M.begin(), computeSquare);	
		transform(indivprobs_F.begin(), indivprobs_F.end(), indivprobs_F.begin(), computeSquare);

        	
		cout << " Effective population size: " << divide(1,ComputeSum(indivprobs_M)) << endl;
		cout << endl;
		cout <<" Effective population size: " << divide(1,ComputeSum(indivprobs_F)) << endl;	
		
		partial_sum ( oldindivprobs_M.begin(), oldindivprobs_M.end(), oldindivprobs_M.begin());
		partial_sum ( oldindivprobs_F.begin(), oldindivprobs_F.end(), oldindivprobs_F.begin());
	

		vector<int> locationinf_M = NewGenMale[g-1].GetLocSegment(numberofSNPs,ChromLen);
		vector<int> locationinf_F = NewGenFemale[g-1].GetLocSegment(numberofSNPs,ChromLen);
		//for( vector<double>::iterator it = oldindivprobs_M.begin(); it != oldindivprobs_M.end(); ++it )
			//	cout << "Indivprobs: " << *it << endl;
		//cout << endl; 
		vector < pair<int,int> > Parents;
		//Selecting parent IDs for each individual and getting chromosome and crossover information from each parent 
		for ( int i = 0; i < popSize*2; i++){
		
			//Selecting Father
			double randval_M = ((double) rand() / (RAND_MAX));
			int FatherIndex;
			int MotherIndex;
			for( vector<double>::iterator it = oldindivprobs_M.begin(); it != oldindivprobs_M.end()-1; ++it ){
				vector<double>::iterator nxt = it;
				++nxt;
				if(*it <= randval_M && *nxt > randval_M) {
					FatherIndex = distance(oldindivprobs_M.begin(),nxt);
					break;
				}
				if( randval_M < *(oldindivprobs_M.begin()) ) { 
					FatherIndex = 0;
					break;
				}
			}	
			//Selecting Mother
			double randval_F = ((double) rand()/ (RAND_MAX));
			for( vector<double>::iterator it = oldindivprobs_F.begin(); it != oldindivprobs_F.end()-1; ++it ){
				vector<double>::iterator nxt = it;
				++nxt;
				if(*it <= randval_F && *nxt > randval_F){
					MotherIndex = distance(oldindivprobs_F.begin(),nxt);
					break;
				}
				if( randval_F < *(oldindivprobs_F.begin()) ) {
					MotherIndex = 0; 
					break;
				}
			}
		
			pair <int,int> ParentIdx = make_pair(FatherIndex,MotherIndex);
			Parents.push_back(ParentIdx);
			
			//This method returns the information from the parents, whether there was a crossover and mutation. We do this for both mother and father.
			//The contents of type AlleleInfo is discused in detail in the Events class
			AlleleInfo fromFather = NewGenMale[g-1].AskChromosome(FatherIndex,M_ChromContainer,locationinf_M,RecombRate, MutationRate, ChromLen,g);
			AlleleInfo fromMother = NewGenFemale[g-1].AskChromosome(MotherIndex,F_ChromContainer,locationinf_F,RecombRate, MutationRate, ChromLen,g);
			
			vector<string> haploids_M = fromFather.ToTheChild.second;
			vector<string> haploids_F = fromMother.ToTheChild.second;
			string currentChromID;
			ChrInfo ChrObject;
			string currentParentChromID;
			pair<pair<string,int>,pair<string,int> > Parentchroms;
			ChromosomeInfo BoxofChrom;
			pair<pair<int,bool>,pair<string,int> > GivingToChildren;
			//For each chromosome ID,this creates a key-value pair, where the key is the chromosome ID and the value is of type ChromosomeInfo, which 
			//contains information about the chromosomes that created it and it's own contribution to it's children. We also push this info to a list 
			//specific to each generation which contains all the individuals parent chromosome info
			if(i < popSize){
				sexflag = 1;
				currentChromID = fromFather.ToTheChild.first;
				Parentchroms = fromFather.ContChrom;
				BoxofChrom = ChrObject.PopulateAndGetChrInfo(make_pair(fromFather.ParentID,0), g, i, sexflag, Parentchroms);
				AllChromRecords.insert(make_pair(currentChromID,BoxofChrom));		
			
				if(fromFather.ContChrom.first.second != 0){ 
					string currentParentChromID = fromFather.ContChrom.first.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i,0),make_pair(currentChromID,fromFather.ContChrom.first.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				if(fromFather.ContChrom.second.second != 0){ 
					string currentParentChromID = fromFather.ContChrom.second.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i,0),make_pair(currentChromID,fromFather.ContChrom.second.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				
				currentChromID = fromMother.ToTheChild.first;
				Parentchroms = fromMother.ContChrom;	
				BoxofChrom = ChrObject.PopulateAndGetChrInfo(make_pair(fromMother.ParentID,1), g, i, sexflag, Parentchroms);
				AllChromRecords.insert(make_pair(currentChromID,BoxofChrom));		
			
				if(fromMother.ContChrom.first.second != 0){ 
					currentParentChromID = fromMother.ContChrom.first.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i,0),make_pair(currentChromID,fromMother.ContChrom.first.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				if(fromMother.ContChrom.second.second != 0){ 
					currentParentChromID = fromMother.ContChrom.second.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i,0),make_pair(currentChromID,fromMother.ContChrom.second.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}			
				
				NewGenMale[g].PushToList(make_pair(sexflag,make_pair(fromFather,fromMother)));
			}
			else{
				sexflag = 0;
				currentChromID = fromFather.ToTheChild.first;
				Parentchroms = fromFather.ContChrom;
				BoxofChrom = ChrObject.PopulateAndGetChrInfo(make_pair(fromFather.ParentID,0), g, i-popSize, sexflag, Parentchroms);
				AllChromRecords.insert(make_pair(currentChromID,BoxofChrom));		
			
				if(fromFather.ContChrom.first.second != 0){ 
					currentParentChromID = fromFather.ContChrom.first.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i-popSize,1),make_pair(currentChromID,fromFather.ContChrom.first.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				if(fromFather.ContChrom.second.second != 0){ 
					currentParentChromID = fromFather.ContChrom.second.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i-popSize,1),make_pair(currentChromID,fromFather.ContChrom.second.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				
				currentChromID = fromMother.ToTheChild.first;
				Parentchroms = fromMother.ContChrom;	
				BoxofChrom = ChrObject.PopulateAndGetChrInfo(make_pair(fromMother.ParentID,1), g, i-popSize, sexflag, Parentchroms);
				AllChromRecords.insert(make_pair(currentChromID,BoxofChrom));		
			
				if(fromMother.ContChrom.first.second != 0){ 
					currentParentChromID = fromMother.ContChrom.first.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i-popSize,1),make_pair(currentChromID,fromMother.ContChrom.first.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				if(fromMother.ContChrom.second.second != 0){ 
					currentParentChromID = fromMother.ContChrom.second.first; 
					map<string,ChromosomeInfo>::iterator it = AllChromRecords.find(currentParentChromID);
					if(it != AllChromRecords.end()){
						GivingToChildren = make_pair(make_pair(i-popSize,1),make_pair(currentChromID,fromMother.ContChrom.second.second));
						it->second.AllChildrenInfo.push_back(GivingToChildren);
					}
				}
				NewGenFemale[g].PushToList(make_pair(sexflag,make_pair(fromFather,fromMother)));
			}
			//cout << endl; 
			//This pushes back and creates the chromosome container to be used for the next generation.
			for (int j = 0; j < numberofSNPs; j++){
				pair <string,string> chrompair1 = make_pair(fromFather.ToTheChild.first,haploids_M[j]);
				pair <string,string> chrompair2 = make_pair(fromMother.ToTheChild.first,haploids_F[j]);
				if ( i < popSize){
					NewGenMale[g].Addpairs(make_pair(chrompair1,chrompair2));
					//NewGenMale[g].AddChromName(g);
					//NewGenMale[g].AddChromName(g);
				}
				else{
					NewGenFemale[g].Addpairs(make_pair(chrompair1,chrompair2));
					//NewGenFemale[g].AddChromName(g);
					//NewGenFemale[g].AddChromName(g);	
				}
			}
		}
		//Clears all the chromosome container holders and data structures which are specific to each generation, so that they can be used again. 
	
		M_ChromContainer.clear();
		F_ChromContainer.clear();
		
		M_ChromContainer = NewGenMale[g].getChromContainer();
		
		F_ChromContainer = NewGenFemale[g].getChromContainer();

		
		StoreProductFitness_M.clear();
		StoreProductFitness_F.clear();
		indivprobs_F.clear();
		indivprobs_M.clear();
	}
	
//This was done to print the containers of the key-value pairs showing how and what is stored for each chromosome ID. 
	/*ChromosomeInfo access;
	for(MapIterator it = AllChromRecords.begin(); it != AllChromRecords.end();++it){
		cout << endl << "ChromID: " << it->first << endl;
		ChromosomeInfo BoxofChrom = it->second; 
		if(BoxofChrom.MyLocation.first == 0)
			cout << "It is a Boy!" << endl; 
 		else
			cout << "It is a Girl!" << endl; 
		cout << "Location of Chromosome is on Generation " << BoxofChrom.MyLocation.second.first << " and Child ID " << BoxofChrom.MyLocation.second.second << endl;
		cout << "Created by, ParentID: " << BoxofChrom.CreatingParentID << " "; 
		if (BoxofChrom.ParentSex == 0)
			cout << "This is the father" << endl; 
		else 
			cout << "This is the Mother" << endl;
		cout << endl << "Parent Chromosome: " << BoxofChrom.Parentchroms.first.first << " Contributing " <<  BoxofChrom.Parentchroms.first.second << " out of " << numberofSNPs << endl;
		cout << endl << "Parent Chromosome: " << BoxofChrom.Parentchroms.second.first << " Contributing " <<  BoxofChrom.Parentchroms.second.second << " out of " << numberofSNPs << endl;
		for(ChromInfoListIterator iter = BoxofChrom.AllChildrenInfo.begin(); iter != BoxofChrom.AllChildrenInfo.end(); ++iter){
			cout << endl << "Created child ID: " << iter->first.first << endl; 
			if(iter->first.second == 0)
				cout << "Who is a male." << endl; 
			else 
			cout << "Who is a female." << endl; 
			cout << "Contributed to Chromosome ID: " << iter->second.first << " " << " with " << iter->second.second << " number of SNPs among " << numberofSNPs << "total SNPs" << endl; 
			cout << endl << "**********************************************" << endl; 
			}
	   cout << endl << "-	----------------------------------ONE CHROMOSOME ENDS--------------------------------------------" << endl;
	}	 */
	
//Selecting k out of m populations, as specified by the user, for the number of populations on which the ARG will be built
//It randomly shuffles the indices and selects the first k out of m individuals. 
	int *RandNumbers = new int[popSize*2];
	vector<int> SelectedPopforArg;
	for (int i = 0; i < popSize*2; i++) 
		RandNumbers[i] = i; 
	random_shuffle(RandNumbers,RandNumbers+(popSize*2));
	for (int i = 0; i < ArgPop; i++){
		SelectedPopforArg.push_back(RandNumbers[i]); 
	}
//The following piece of code creates the ARG for each chromosome IDs. 
//*****************************************THIS IS A WORK IN PROGRESS*************************************************************
//Hence, I am commenting the following out. 
/*
	pair<bool,pair<AlleleInfo,AlleleInfo> > fromParent;
	vector <vector < pair<string, int> > > ArgList;
	map <string, vector<vector< pair<string,int> > > > ARGRepo;
	string KeyName;
	char numstr[10];
	string male = "Male";
	string female = "Female";
	for (vector<int>::iterator it = SelectedPopforArg.begin(); it != SelectedPopforArg.end(); ++it){
		cout << "pos: " << *it << endl;
		if (*it < popSize) {
			fromParent = NewGenMale[GenNum-1].RetrieveListElement(*it);
			//cout << endl << "in loop" << endl;
			snprintf(numstr,10,"%d",*it);
			KeyName = male + numstr;
			//ArgList = createARGs(fromParent.second);
		

		}
		else{
			fromParent = NewGenFemale[GenNum-1].RetrieveListElement(*it-popSize);
			//cout << endl << "in loop" << endl; 
			snprintf(numstr,10,"%d",*it);
			KeyName = female + numstr;
			//ArgList = createARGs(fromParent.second);
		}
		
	}*/
	delete[] NewGenFemale;
	delete[] NewGenMale;
	int stop_s=clock();
	cout << "time: " << ((stop_s-start_s)/double(CLOCKS_PER_SEC)*1000)/60000 << "mins" << endl;

}
//This function creates the ARGs for each individuals, 2 individuals meaning 4 chromosomes.
/* 	
vector<vector< pair< string,int> > > createARGs(pair<AlleleInfo,AlleleInfo> fromParent){
	AlleleInfo fromFather = fromParent.first; 
	AlleleInfo fromMother = fromParent.second;
	pair<string,int> GrandFatherPat;
	pair<string,int> GrandMotherPat;
	pair<string,int> GrandMotherMat;
	pair<string,int> GrandFatherMat;
	vector <pair<string,int> > ChromLinks;
	if(fromFather.ContChrom.first.second != 0){
		GrandFatherPat = fromFather.ContChrom.first;
		findARG(GrandFatherPat,ChromLinks);//baba-r baba-r ta dekhi 
	}
	if(fromFather.ContChrom.second.second != 0){
		GrandMotherPat = fromFather.ContChrom.second;
		findARG(GrandMotherPat);//baba-r ma er ta dekhi 
	}
	if(fromMother.ContChrom.second.first != 0){ 
		GrandFatherMat = fromFather.ContChrom.first;
		findARG(GrandFatherPat);//ma-r baba-r ta dekhi 
	}
	if(fromMother.ContChrom.second.second != 0){
		GrandMotherMat = fromFather.ContChrom.second;
		findARG(GrandMotherPat);//ma-r ma er ta dekhi 
	}
}
//This function creates the ARG for each run. This is unfinished code and work needs to be done. 
<vector<pair<string,int> >, vector<int> > findARG(pair<string,int> ChromQuant,vector< pair<string, int> > ChromLinks){
		ChromosomeInfo ChrBox = AllChromRecords.find(ChromQuant.first)->second; 
		if(ChromQuant.second <= ChrBox.Parentchroms.first.second && ChrBox.Parentchroms.second.second == 0){
			ChromLinks.push_back(ChrBox.Parentchroms.first);
			findARG(ChrBox.Parentchroms.first,ChromLinks);
		}		
		else if (ChrBox.Parentchroms.first.second == 0 && ChromQuant.second >= ChrBox.Parentchroms.second.second){
			ChromLinks.push_back(ChrBox.Parentchroms.second);
			findARG(ChrBox.Parentchroms.second,ChromLinks);			
		}
		else if (ChromQuant.second < 0 && )
			
}
*/