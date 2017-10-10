/*

SSimRA: A framework for selection in coalescence with Recombination 
Author: Aritra Bose 
Last Update: 08/11/2016

This is the implementation of the Individuals class where everything that relates to one run of a generation is generated and kept
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
#include<map>
#include "Individuals.h"
#include "GlobalIndivs.h"
#include "Events.h"
#define FITNESS_MAX 0.001
#define FITNESS_MIN -0.001
using namespace std;


void Individuals :: ReserveSpace(int totalSize){
	ChromContainer.reserve(totalSize);	
	}
//Chromosome Container is the main building block for every generation. It stores information about what are the diploids for each individual across the number of SNPs. 
void Individuals :: Addpairs(pair<pair<string,string>,pair<string,string> > apair){		
	ChromContainer.push_back(apair);
	}

vector<pair<pair<string,string>, pair<string,string> > > Individuals :: getChromContainer(){
	return ChromContainer;
	}
void Individuals :: DestroyCC (){
	ChromContainer.clear();
	}
				 
void Individuals :: DisplayContents(){
	for( vector<pair<pair<string,string>, pair<string,string> > >::iterator it = ChromContainer.begin(); it != ChromContainer.end(); ++it){
			
		cout << " " << it->first.first << "," << it->first.second << " | ";
		cout << " " << it->second.first << "," << it->second.second << endl;
		}
	cout << endl;
	}

void Individuals :: Addlabel(){
	int num = rand()%9000 + 1000;
	char numstr[10];
	snprintf(numstr,10,"%d",num);
	string label = 	rsid + 	numstr;
	SNPlabels.push_back(label);
	}
//It picks a random number of 4 digits, converts it to string and follows the nomenclature of Chrom_GenNum_RandomNumber to name the chromosome ID
string Individuals :: AddandGetChromName(int GenNum){
	int num = rand()%9000+1000;
	char numstr[10];
	snprintf(numstr,10,"%d",num);
	char numstr1[10];
	snprintf(numstr1,10,"%d",GenNum);
	string label = chromid + "_" + numstr1 + "_" + numstr; 
	ChromLabels.push_back(label);
	return label; 
}
	
void Individuals :: DisplayLabels(){
	for (vector<string>:: iterator it = SNPlabels.begin(); it != SNPlabels.end(); ++it){			
		cout << *it << " " ;
		}
	}

//This gets the c=location segment randomly for each SNP across the length of chromosome, L is always greater than NumberofSNPs
vector<int> Individuals :: GetLocSegment(int NumberofSNPs,int L){
	vector <int> ChromSegLoc; 
	int *RandNumbers = new int[L];
	int uniqueflag;
	int randloc;	
	for (int i = 0; i < L; i++) 
		RandNumbers[i] = i; 
	random_shuffle(RandNumbers,RandNumbers+L);
	for (int i = 0; i < NumberofSNPs; i++)
		ChromSegLoc.push_back(RandNumbers[i]);
							
	return ChromSegLoc; 
	}

//This defines and fetches the fitness table for each individual and this remains constant throughout the simulation	
double* Individuals :: getFitnessTable(int NumberofSNPs){	
	double *PopulationFitnessTable = new double[NumberofSNPs*diploidSize];	
	for ( int i = 0; i <  NumberofSNPs; i++ ){
		vector<int> borderIndices;
		borderIndices.push_back(0);
				
		for (int j = 0; j < diploidSize ; j++){
			double valinsert = doubleRand(FITNESS_MIN,FITNESS_MAX); 
			PopulationFitnessTable[i*diploidSize + j] = valinsert;	
			}	
		int baseIndex = rand()%10;
		PopulationFitnessTable[i*diploidSize+baseIndex] = 0.00; 
	}
				 
	return PopulationFitnessTable;			
	}

		
void Individuals :: DisplayTable(double* PopulationFitnessTable, int NumberofSNPs){
	for (int i = 0; i < NumberofSNPs; i++){	
		for (int j = 0; j < diploidSize; j++){
				cout << PopulationFitnessTable[i*diploidSize + j] << " | " ;
			}
			cout << endl;
		}
	}

//This gets the position of the allele by finding where it lies in the vector of diploids which is in GlobalIndivs class 
int Individuals :: getPosAllele(vector<string> RetrievedPair){
	
	string allelefront = RetrievedPair[1] + RetrievedPair[3];
    string alleleend = RetrievedPair[3] + RetrievedPair[1];
	int pos_allele;

    if (find (diploidVec.begin(), diploidVec.end(), allelefront) != diploidVec.end())

        pos_allele = distance(diploidVec.begin(), find (diploidVec.begin(), diploidVec.end(), allelefront));
    else
        pos_allele = distance(diploidVec.begin(), find (diploidVec.begin(), diploidVec.end(), alleleend));

	return pos_allele; 
	}
pair<bool,pair<AlleleInfo,AlleleInfo> > Individuals :: RetrieveListElement(int x){
	//cout << endl << "eg" << endl;
	return ParentInfoList[x]; 
}
//This method mutates the SNPs depending on the mutation rate specified by the user
//It creates 4 bins, one with no mutation and other into three parts where possible mutation might take place to the other three bases. 
//This is done everytime the parents sends a set of SNPs to the child 
//This also saves the SNPs before and after mutations, helping us track while building the ARGs
pair<vector<string>, vector<pair<string,int> > > Individuals :: Mutate(vector<string> haploids, double MutationRate){
	int numSNP = haploids.size();
	int mutidx; 
	double randmut;
	string targetallele;
	vector<string> tmpbaseVec;
	vector<string> mutatedhaploids; 
	vector<pair<string,int> > mutatedpos;
	pair<vector<string>, vector<pair<string,int> > > MutAlleleContainer;
	int index; 
	for (int i = 0; i < numSNP; i++){
		tmpbaseVec = baseVec;
		randmut = ((double) rand()/(double) (RAND_MAX));
		targetallele = haploids[i];
		//cout << endl << "Target: " << targetallele << endl;
		if(find (baseVec.begin(), baseVec.end(), targetallele) != baseVec.end())
			index = distance(baseVec.begin(), find (baseVec.begin(), baseVec.end(), targetallele));
		tmpbaseVec.erase(tmpbaseVec.begin() + index);
		//cout << endl << "Remaining bases" << endl; 
		//for(vector<string>::iterator it1 = tmpbaseVec.begin();it1 != tmpbaseVec.end(); ++it1)
		//	cout << *it1 << endl;
		if (randmut > 0 && randmut < (1-MutationRate)) mutidx = 1;
		else if (randmut >= (1-MutationRate) && randmut < ((3-(2*MutationRate))/3)) mutidx = 2;
		else if (randmut >= ((3-(2*MutationRate))/3) && randmut < ((3-MutationRate)/3)) mutidx = 3;
		else if (randmut >= ((3-MutationRate)/3) && randmut < 1) mutidx = 4; 
		
		switch(mutidx){
			
			case 1:	mutatedhaploids.push_back(haploids[i]);
					//cout << endl << "at pos: " << i << " " << haploids[i] << " was not changed" << endl; 
					mutatedpos.push_back(make_pair("",i));
					break;					
			case 2:	mutatedhaploids.push_back(tmpbaseVec[0]); 
					//cout << endl << "at pos: " << i << " " << haploids[i] << " was changed to: " << tmpbaseVec[0] << endl;
					mutatedpos.push_back(make_pair(haploids[i],i));
					break; 
			case 3:	mutatedhaploids.push_back(tmpbaseVec[1]);
					//cout << endl << "at pos: " << i << " " << haploids[i] << " was changed to: " << tmpbaseVec[1] << endl;
					mutatedpos.push_back(make_pair(haploids[i],i));
					break; 
			case 4:	mutatedhaploids.push_back(tmpbaseVec[2]); 
					//cout << endl << "at pos: " << i << " " << haploids[i] << " was changed to: " << tmpbaseVec[2] << endl;
					mutatedpos.push_back(make_pair(haploids[i],i));
					break; 
		}
	}
	//for (vector<pair<string,int> >::iterator it = mutatedpos.begin(); it != mutatedpos.end(); ++it){
		//cout << it->first << " | " << it->second << endl; 
	//}
	//cout << endl << "------- " << endl; 
	//for (vector<string>::iterator it = mutatedhaploids.begin(); it != mutatedhaploids.end(); ++it){
		//cout << *it << endl;
	//}
	MutAlleleContainer = make_pair(mutatedhaploids,mutatedpos);
	return MutAlleleContainer; 
} 

//This method gets the alleles for the child when there's no recombination. This can happen two ways, either father sends all of his father's SNPs to child, or,
//the father sends all of his mother's. After it's decided which, the information is kept in the container called Contributing Chromosomes.
//The chromosome which is contributing nothing to the child has 0 in its position, indicating it's contribution. This is helpful while building the ARG. 
//It collapses the edges when there's no recombination.
AlleleInfo Individuals :: GetAllelesforChild(int Pidx, vector<pair<pair<string,string>, pair<string,string> > > CC, vector<int> SNPlocs, double MutationRate, int GenNum, int chromnum){
	//cout << endl << "NO CROSSOVER" << endl;
	int numSNP = SNPlocs.size();
	pair<pair<string,int>,pair<string,int> > ContributingChroms; 
	vector<string> haploids; 
	string chromID,theotherID;	
	pair<string,int> c1,c2; 
	for (int i = 0; i < numSNP; i++){
		pair<pair<string,string>, pair<string,string> > chrompair = ChromContainer[Pidx*numSNP + i];
		if (chromnum == 1){
			haploids.push_back(chrompair.first.second);
			if (i==0){
				chromID = chrompair.first.first; 
				theotherID = chrompair.second.first; 
			}
		}
		else{	
			haploids.push_back(chrompair.second.second);
			if (i==0){
				chromID = chrompair.second.first;
				theotherID = chrompair.first.first; 
			}	
		}
	}
	c1 = make_pair(chromID,numSNP);
	c2 = make_pair(theotherID,0);
	//They pack the chromosome contributions at this stage
	ContributingChroms = make_pair(c1,c2);
	string chromname = AddandGetChromName(GenNum);
	//Now they mutate the SNPs that was sent 
	pair<vector<string>, vector<pair<string,int> > > MutBox = Mutate(haploids,MutationRate);
	vector<string> MutatedHaploids = MutBox.first; 
	vector<pair<string,int> > MutatedposInfo = MutBox.second; 

	Events event; 
	int pos = -1; //might change afterwards, don't worry about this
	//at this point they populate the structure to record the event and return the struct to the child. 
	AlleleInfo toChild = event.getandPopulateStruct(pos, Pidx, chromname, ContributingChroms, MutatedHaploids, MutatedposInfo);
		
		return toChild;
}
//This builds the list for each individual in a given generation with all their parent info in it. Easier to access any given individual's info
//while building the ARG (although, at this stage I am not using this, but it's an useful data structure to have)
void Individuals :: PushToList(pair<bool, pair <AlleleInfo,AlleleInfo> > structpairs){
	this->ParentInfoList.push_back(structpairs);
}

//This displays the aforementioned list for each individual
void Individuals :: DisplayParentList(){
	for(vector<pair<bool, pair<AlleleInfo,AlleleInfo> > >::iterator it = ParentInfoList.begin(); it != ParentInfoList.end(); ++it){
	
		cout << endl << " | " << "FatherID: "  << it->second.first.ParentID << " | " << "Chromosome ID: " << it->second.first.ToTheChild.first << " | " << endl;
		cout << endl <<  "From Father ------ Chromosome: " << it->second.first.ContChrom.first.first << " Contributing: " << it->second.first.ContChrom.first.second << " and Chromosome: " << it->second.first.ContChrom.second.first << " Contributing: " << it->second.first.ContChrom.second.second << endl; 

		for(vector<pair<string,int> >::iterator rot = it->second.first.BeforeMutation.begin(); rot != it->second.first.BeforeMutation.end(); ++rot)
			cout << "-----Baba Original haploid: " << rot->first << " --------- " << "Mutation Position: " << rot->second << endl; 
		
		cout << endl << " | " << "MotherID: "  << it->second.second.ParentID << " | " << "Chromosome ID: " << it->second.second.ToTheChild.first << " | " << endl;
		cout << endl <<  "From Mother ------ Chromosome: " << it->second.second.ContChrom.first.first << " Contributing: " << it->second.second.ContChrom.first.second << " and Chromosome: " << it->second.second.ContChrom.second.first << " Contributing: " << it->second.second.ContChrom.second.second << endl; 

		//for(vector<pair<string,int> >::iterator rot = it->second.BeforeMutation.begin(); rot != it->second.BeforeMutation.end(); ++rot)
			//cout << "-----Maa Original haploid: " << rot->first << " --------- " << "Mutation Position: " << rot->second << endl;
	}
}	

/*This method gets the alleles for the child when there's recombination. This can happen two ways, either father sends some (depending on the random 'C' that is chosen in L) of his father's SNPs to child
and the rest is chosen from his mother's, OR, vice versa.  After it's decided which, the information is kept in the container called Contributing Chromosomes, with positive number in the position
indicating forward contribution and the number of SNPs contributed, the negative number in pos indicates backward contribution and the number of SNPs contributed.
There will be no 0 in the position, this time. This is helpful while building the ARG. 
It collapses the edges when there's no recombination. */

AlleleInfo Individuals :: GetCrossOverAlleles(int Pidx, vector<pair<pair<string,string>, pair<string,string> > > CC, vector<int> SNPlocs, double MutationRate, int L, int GenNum, int chromnum){
	//cout << endl << "YES CROSSOVER" << endl; 
	int pos;
	int numSNP = SNPlocs.size();
	vector <string> haploids; 
	vector<int> SNPidx; 
	pair<pair<string,string>,pair<string,string> > chrompair;
	string chromID;
	pair<pair<string,int>,pair<string,int> > ContributingChroms; 
	pair<string,int> c1,c2; 
	//randomly choses C out of L
	//The first part is crossover1 and the second part is crossover2, delineated with the if-else
	int crossoveridx = rand()%(L + 1);
	if (chromnum == 1){
			for (int i = 0; i < numSNP; i++){
				if (SNPlocs[i] <= crossoveridx)
					SNPidx.push_back(i);
			}
			if(SNPidx.size() == 0){
				chrompair=CC[Pidx*numSNP];
				chromID = chrompair.first.first;
			}			
			for(int i = 0; i < SNPidx.size(); i++){ 
				chrompair = CC[Pidx*numSNP + SNPidx[i]];
				haploids.push_back(chrompair.first.second);
				if (i==0)
					chromID = chrompair.first.first;  
			}
			pos = SNPidx.size();
			c1 = make_pair(chromID,pos);
			SNPidx.clear();
			for (int i = 0; i < numSNP; i++){
				if (SNPlocs[i] > crossoveridx)
					SNPidx.push_back(i);
			}
			if(SNPidx.size() == 0){
				chrompair = CC[Pidx*numSNP];
				chromID = chrompair.second.first;
			}
			for (int i = 0; i < SNPidx.size(); i++){
				chrompair = CC[Pidx*numSNP + SNPidx[i]];
				haploids.push_back(chrompair.second.second);
				if (i==0)
					chromID = chrompair.second.first;
				//cout << "at pos: " << SNPidx[i] << " from Chr2: " << chrompair.second.second << endl;
			}
			pos = SNPidx.size();
			c2 = make_pair(chromID,-pos);		
		ContributingChroms = make_pair(c1,c2);
		}
	else {
			for (int i = 0; i < numSNP; i++){
				if (SNPlocs[i] <= crossoveridx)
					SNPidx.push_back(i);
			}
			if(SNPidx.size() == 0){
				chrompair=CC[Pidx*numSNP];
				chromID = chrompair.second.first;
			}			
			for(int i = 0; i < SNPidx.size(); i++){
				chrompair = CC[Pidx*numSNP + SNPidx[i]];
				haploids.push_back(chrompair.second.second);
				if (i==0)
					chromID = chrompair.second.first;
				//cout << "at pos: " << SNPidx[i] << " from Chr2: " << chrompair.first.second << endl; 
			}
			pos = SNPidx.size();
			c1 = make_pair(chromID,pos);
			SNPidx.clear();
			for (int i = 0; i < numSNP; i++){
				if (SNPlocs[i] > crossoveridx)
					SNPidx.push_back(i);
			}
			//cout << endl << SNPidx.size () << endl;
			if(SNPidx.size() == 0){
				chrompair=CC[Pidx*numSNP];
				chromID = chrompair.first.first;
			}
			for (int i = 0; i < SNPidx.size(); i++){
				chrompair = CC[Pidx*numSNP + SNPidx[i]];
				haploids.push_back(chrompair.first.second);
				if (i==0)
					chromID = chrompair.first.first;
			}
			pos = SNPidx.size();
			c2 = make_pair(chromID,-pos);
			//packs everything into contributing chromosome container indicating positions.
		ContributingChroms = make_pair(c1,c2);
 
		}
		//Mutates the SNPs before passing them on to the child as described in the Mutate method
		pair<vector<string>, vector<pair<string,int> > > MutBox = Mutate(haploids,MutationRate);
		vector<string> MutatedHaploids = MutBox.first; 
		vector<pair<string,int> > MutatedposInfo = MutBox.second; 

		string chromname = AddandGetChromName(GenNum);
		Events event; 
		//passed on to the child in the AlleleInfo container
		AlleleInfo toChild = event.getandPopulateStruct(pos, Pidx, chromname, ContributingChroms, MutatedHaploids, MutatedposInfo);
		
		return toChild;
	}
/*This is the method which executes whether there was a recombination or not, depending on the recombination rate fed by the user. This creates 4 equidistant
bins and tosses a coin to decide whether it's crossover or not and whatever the result whether it comes from the father's strand or mother's strand of SNPs.
Then it obtains the AlleleInfo container after invoking the aforementioned methods to get the haploids to the child and passes it on as a return value */
AlleleInfo Individuals :: AskChromosome (int Pidx, vector<pair<pair<string,string>, pair<string,string> > > CC, vector<int> SNPlocs, double RecombRate, double 	MutationRate, int L, int GenNum){ 
		int recomidx; 
		//Coin Toss
		double randrecom = ((double) rand()/(double) (RAND_MAX));
		//4 bins and seeing where the coin lands
		if (randrecom > 0 && randrecom < ((1-RecombRate)/2)) recomidx = 1;
		else if (randrecom >= ((1-RecombRate)/2) && randrecom < (1-RecombRate)) recomidx = 2;
		else if (randrecom >= (1-RecombRate) && randrecom < (1-(RecombRate/2))) recomidx = 3;
		else if (randrecom >= (1-(RecombRate/2)) && randrecom < 1) recomidx = 4; 
			
		//Depending on the bin where it landed it decides what to do and send the result to the child
		AlleleInfo toChild;
		switch(recomidx){
	
			case 1:  toChild = GetAllelesforChild(Pidx,CC,SNPlocs,MutationRate,GenNum,1); break;
			case 2:  toChild = GetAllelesforChild(Pidx,CC,SNPlocs,MutationRate,GenNum,2); break;
			case 3:  toChild = GetCrossOverAlleles(Pidx,CC,SNPlocs,MutationRate,L,GenNum,1); break;
			case 4:  toChild = GetCrossOverAlleles(Pidx,CC,SNPlocs,MutationRate,L,GenNum,2); break;
		}
		return toChild; 
	}
	 
