/*
SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 08/11/2016

This is a class where each Chromosome information is stored. This stores information of each chromosome's present, past and future states. 
This chromosome information records the parent's chromosome IDs and the number of SNP they contributed to the present chromosome. The present
chromosome's record of it's location, it's gender is stored, so that it can be tracked by generation and then individual. 
It also keeps a track of all the children it contributed to and the amount of contributions. 

*/
#include<iostream> 
#include<vector>
#include<algorithm>
#include<stdio.h>
#include<cstdlib>
#include<string.h> 
#include<math.h>
#include<map>
#include<list> 
#include<vector>
#ifndef CHRINFO_H
#define CHRINFO_H

using namespace std;

struct ChromosomeInfo{
	pair<pair<string,int>,pair<string,int> > Parentchroms; 
	int CreatingParentID; 
	bool ParentSex; 
	pair<bool, pair<int,int> > MyLocation;
	vector<pair<pair<int,bool>, pair<string,int> > > AllChildrenInfo; //This is a vector of of pairs, recording the location, gender and the 
																	 //children chromosomes in the next generation and what it contributed to
};

class ChrInfo{
	
	public: 
		ChromosomeInfo access;
		//this method populates the ChromosomeInfo structure with information from each chromosome. They are eventually sent back to the map 
		//listing all chromosome and their subsequent information. This is required while building the ARG 
		ChromosomeInfo PopulateAndGetChrInfo(pair<int,bool> ParentInfo, int GenNum, int CurrentID, bool sexflag, pair<pair<string,int>,pair<string,int> > Parentchroms);
		//void AddContributions(pair<int, pair<string,int> > ChromContributions);
		vector<pair<pair<int,bool>, pair<string,int> > > GetAllChildrenInfo();
};

#endif