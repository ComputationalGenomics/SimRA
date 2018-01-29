/*

SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 08/11/2016

This is a class where everything that relates to one run of a generation is generated and kept
It contains the Chromosome Containers for male and female which is referred to while selecting the parents for each individuals
in a generation. 
It contains the Chromsome IDs of each chromosomes and also the chromosome information list for each individual. 
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
#include<gsl/gsl_rng.h>
#include<map>
#include "GlobalIndivs.h"
#include "Events.h"
#ifndef indivs
#define indivs
//using namespace std;

class Individuals
{
	//These are private data structures as they are specific to each generation and no other generation can use them. 
	private:
		std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > ChromContainer;
		std::vector<string> SNPlabels;
		std::vector<string> ChromLabels; 
		std::vector<std::pair<bool,std::pair<AlleleInfo,AlleleInfo> > > ParentInfoList; 
	//All the method used to populate the above data structures are declared public as they are the standard way to get these information
	public:
		Individuals()
		{}
		//This method reserves space for the chromosome containers
		void ReserveSpace(int x);
		//Adding std::pairs to the Chromosome Containers
		void Addpairs(std::pair<std::pair<string,string>,std::pair<string,string> > x, int y); 
		//This just fetches and returns the Chromosome Container
		std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > getChromContainer();
		//This is a destructor for ChromContainer
		void DestroyCC ();
		//Displays the content of the chromcontainer
		void DisplayContents();
		//This is for adding labels to SNPs (not relevant at this point)
		void Addlabel();
		//This is for adding labels to Chromosomes (very relevant)
		string AddandGetChromName(int x);
		//Displaying the Chromosome and SNP labels
		void DisplayLabels();
		//This gets the location segment of each SNP inside the length of the chromosome
		std::vector<int> GetLocSegment(int x,int y); 
		//This fetches the fitness table for a particular generation and the fitness table remains constant throughout the simulation
		double *getFitnessTable(int x, double y, double z);	
		//This is for displaying the fitness table (for debugging purposes)
		void DisplayTable(double* x, int y);
		//This gets the position of the allele from the fitness table with it's respective fitness value
		int getPosAllele(std::vector<string> x);
		//This pushes the parent chromosome info for each individual in the last generation
		void PushToList(std::pair<bool, std::pair <AlleleInfo,AlleleInfo> > x);
		//Displaying the List of Parents for each individual in a given generation (not relevant)
		void DisplayParentList(); 
		void FitnessON(int x);
		//This method mutates the SNPs according to mutation rate
		std::pair<std::vector<string>, int > Mutate(std::vector<string> x,double y, int L);
		//This fetches the List element containing parent chromosome info for each individual in the last generation
		//std::pair<bool,std::pair<AlleleInfo,AlleleInfo> > RetrieveListElement(int);
		//This method gets alleles for the child, when there's no crossover
		AlleleInfo GetAllelesforChild(int x, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC, std::vector<int> SNPlocs, double MutationRate, int L, int GenNum, int chromnum);
		//This method gets alleles for the child, when there's a crossover 
		AlleleInfo GetCrossOverAlleles(int Pidx, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC, std::vector<int> SNPlocs, double MutationRate, int L, int GenNum, int chromnum);
		//This method is invoked when the child asks for chromosome from it's parents and depending on the crossover result either of the above two methods is invoked.
		AlleleInfo AskChromosome (int Pidx, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC, std::vector<int> SNPlocs, double RecombRate, double MutationRate, int L, int GenNum);
};

#endif