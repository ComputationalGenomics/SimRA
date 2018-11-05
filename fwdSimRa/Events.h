/*

fwdSimRA: A framework for selection forward-in-time with Recombination 
Author: Aritra Bose 
Last Update: 11/04/2018

This is a class where each event information is stored. 
An event is recorded when the child asks for the chromosomes from the parents and it records the position of the child, the IDs of the parents. 
The contributing chromosomes with the chromosome ID and the how many haploids they are contributing to the child. 
It also keeps a track of the state before Mutation and also after mutation, the haploids along with the chrom ID tag is sent to the child.

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
#ifndef EVENTS_H
#define EVENTS_H

using namespace std;

struct AlleleInfo{
			int position;
			int ParentID;  
			pair<pair<string,int>,pair<string,int> > ContChrom;
			int MutCount;
			int RecombCount; 
			pair <string,vector<string> >  ToTheChild; 
		};
		
class Events{ 

	public: 
		//this is an object to access the structure AlleleInfo
		AlleleInfo access;
		//this method returns a structure of AlleleInfo with populating all the fields of the structure
		AlleleInfo getandPopulateStruct(int pos, int pid, string chromname, pair<pair<string,int>,pair<string,int> > ContChrom, vector<string> haploids, int localmut, int localrecomb);
		
		

};
#endif
