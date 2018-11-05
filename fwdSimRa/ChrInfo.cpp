/*
fwdSimRA: A framework for selection forward-in-time with Recombination 
Author: Aritra Bose 
Last Update: 11/04/2018

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
#include<math.h>
#include<map>
#include<list>
#include<vector> 
#include "Events.h"
#include "ChrInfo.h"

//using namespace std; 
//this method populates the ChromosomeInfo structure with information from each chromosome. They are eventually sent back to the map 
//listing all chromosome and their subsequent information. This is required while building the ARG 
ChromosomeInfo ChrInfo :: PopulateAndGetChrInfo(std::pair<int,bool> ParentInfo, int GenNum, int CurrentID, bool sexflag, std::pair<std::pair<string,int>,std::pair<string,int> > ChromOfParent){
	ChromosomeInfo access;
	//access.CreatingParentID = ParentInfo.first;
	//access.ParentSex = ParentInfo.second;
	access.Parentchroms = ChromOfParent;
	//access.MyLocation = std::make_pair(sexflag,std::make_pair(GenNum,CurrentID));
	//access.AllChildrenInfo = GetAllChildrenInfo();
	
	return access;
}

//This method fetches the vector which containts information about where a chromosome has contributed to in the next generations 
//The first pair relates to the individual ID, where that chromosome belongs to and it's sex. The other pair contain's the chromosome ID and the amount of
//contribution from the present chromosome. When we build the segment trees in forward manner, this will be needed. 
/*vector<pair<pair<int,bool>, pair<string,int> > > ChrInfo :: GetAllChildrenInfo(){
	
	ChromosomeInfo access;
	if (access.AllChildrenInfo.empty())
		return vector<pair<pair<int,bool>, pair<string,int> > >();
	else
		return access.AllChildrenInfo;
}*/
