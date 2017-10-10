/*

SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 08/11/2016

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
#include<math.h>
#include<map>
#include<list> 
#include "Events.h"

using namespace std; 
//It gets the input from the crossover method in the Individuals.h and send populates the structure
AlleleInfo Events :: getandPopulateStruct(int pos, int pid, string chromname, pair<pair<string,int>,pair<string,int> > ContChrom, vector<string> haploids, vector<pair<string,int> > &MutPosInfo){
	
	AlleleInfo access;
	access.position = pos;
	access.ParentID = pid; 
	access.ContChrom = ContChrom;
	access.BeforeMutation = MutPosInfo;
	pair <string, vector<string> > chrompair = make_pair(chromname,haploids);
	access.ToTheChild = chrompair; 
	
	return access; 
}

