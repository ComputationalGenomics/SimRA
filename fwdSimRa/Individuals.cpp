/*

fwdSimRA: A framework for selection forward-in-time with Recombination 
Author: Aritra Bose 
Last Update: 11/04/2018

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
#include<omp.h>
#include<map>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<functional>
#include "Individuals.h"
#include "GlobalIndivs.h"
#include "Events.h"
//using namespace std;

void Individuals :: ReserveSpace(int totalSize){
	ChromContainer.resize(totalSize);	
	}
//Chromosome Container is the main building block for every generation. It stores information about what are the diploids for each individual across the number of SNPs. 
void Individuals :: Addpairs(std::pair<std::pair<string,string>,std::pair<string,string> > apair,int loc){		
	//std::cout << std::endl << "Val of Loc is: " << loc << std::endl <<astd::pair.first.first << "  " << astd::pair.first.second << "  " << astd::pair.second.first << "  " << astd::pair.second.second << std::endl;
	ChromContainer[loc] = apair;
	//ChromContainer.insert(ChromContainer.begin()+loc,astd::pair);
	}

std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > Individuals :: getChromContainer(){
	return ChromContainer;
	}
void Individuals :: DestroyCC (){
	ChromContainer.clear();
	}
				 
void Individuals :: DisplayContents(){
	for( std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > >::iterator it = ChromContainer.begin(); it != ChromContainer.end(); ++it){
			
		std::cout << " " << it->first.first << "," << it->first.second << " | ";
		std::cout << " " << it->second.first << "," << it->second.second << std::endl;
		}
	std::cout << std::endl;
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
	//int num = rand()%9000+1000;
	//std::cout << "thread num is: " << omp_get_thread_num() << "  " << threadvec[omp_get_thread_num()] << std::endl; 
	int num = (int)(gsl_rng_uniform(threadvec[omp_get_thread_num()])*1000000)%9000+1000;
	std::string numstr = std::to_string(num);
	std::string numstr1 = std::to_string(GenNum);
	
	string label = chromid + "_" + numstr1 + "_" + numstr; 
	return label; 
}
	
void Individuals :: DisplayLabels(){
	for (std::vector<string>:: iterator it = SNPlabels.begin(); it != SNPlabels.end(); ++it){			
		std::cout << *it << " " ;
		}
	}

//This gets the c=location segment randomly for each SNP across the length of chromosome, L is always greater than NumberofSNPs
std::vector<int> Individuals :: GetLocSegment(int NumberofSNPs,int L){
	std::vector<int> ChromSegLoc(NumberofSNPs);
	int uniqueflag;
	int randloc;
	srand(time(0));
	int *RandNumbers = new int[L];
	int *SelRandNum = new int [NumberofSNPs];
	for (int z = 0; z < L; z++) 
		RandNumbers[z] = z; 
	//std::vector<int> Rndnums(RandNumbers,RandNumbers+L);
	// std::cout << std::endl << "before shuffle!" << std::endl;
	gsl_ran_shuffle(threadvec[omp_get_thread_num()],RandNumbers, L, sizeof (int));
	//std::random_shuffle(Rndnums.begin(),Rndnums.end());
	//gsl_ran_choose(threadvec[omp_get_thread_num()],SelRandNum, NumberofSNPs, RandNumbers, L, sizeof (int));
	for(int z = 0; z < NumberofSNPs; z++)
		ChromSegLoc[z] = RandNumbers[z];
	delete[] RandNumbers;
	delete[] SelRandNum;
	return ChromSegLoc; 
	}

//This defines and fetches the fitness table for each individual and this remains constant throughout the simulation	
/*double* Individuals :: getFitnessTable(int NumberofSNPs, double FITNESS_MAX, double FITNESS_MIN){	
	double *PopulationFitnessTable = new double[NumberofSNPs*diploidSize];	
	std::cout << std::endl << "Fitness MAX is " << FITNESS_MAX << " and " << "Fitness MIN is " << FITNESS_MIN << std::endl;
	for ( int i = 0; i <  NumberofSNPs; i++ ){
		std::vector<int> borderIndices;
		borderIndices.push_back(0);
				
		for (int j = 0; j < diploidSize ; j++){
			double valinsert = doubleRand(FITNESS_MIN,FITNESS_MAX); 
			PopulationFitnessTable[i*diploidSize + j] = valinsert;	
			}	
		//int baseIndex = rand()%10;
		//PopulationFitnessTable[i*diploidSize+baseIndex] = 0.00; 
	}
				 
	return PopulationFitnessTable;			
	}
*/
void Individuals::PopulateEFT(std::vector<int> epiid, double s1, double s2){
	

	std::vector<double> fitarr1(diploidSize); 
	std::vector<double> fitarr2(diploidSize);

	for (int i = 0; i < diploidSize; i++){
		std::string targetdiploid = diploidVec[i]; 
		if(targetdiploid[0] == 'A' || targetdiploid[0] != 'T'){
			if(targetdiploid[1] != 'A'|| targetdiploid[1] != 'T'){
				fitarr1[i] = 2*s1; 
				fitarr2[i] = 0.0; 
			}
			else{
				fitarr1[i] = s1 + s2; 
				fitarr2[i] = s1 + s2; 
			}
		}
		else{
			if(targetdiploid[1] != 'A' || targetdiploid[1] != 'T'){
				fitarr1[i] = s2 + s1; 
				fitarr2[i] = s1 + s2; 
			}
			else{
				fitarr1[i] = 0.0; 
				fitarr2[i] = 2*s2;
			}
		}
	}
	int ct = 0;
	for(std::vector<std::vector<int>>::iterator it = newSNPs.begin(); it != newSNPs.end(); ++it){
		if (it->size() > 1){
			if(*it == epiid){
				ct++;
				break;
			}
			else 
				ct++;
		}
	}
	
	//std::cout << std::endl << " count was " << ct << std::endl; 
	for (int i = 0; i < diploidSize; i++){
		for (int j = 0; j < diploidSize; j++){
				EpiFitnessTable[(ct*diploidSize*diploidSize)+(i*j)] = ((1+fitarr1[i])*(1+fitarr2[j]))+delta; 
		}
	}
	for (int i = 0; i < diploidSize; i++){	
		for (int j = 0; j < diploidSize; j++){
				std::cout << EpiFitnessTable[ct*diploidSize*diploidSize + (i*j)] << " | " ;
			}
			std::cout << std::endl;
		}
}

void Individuals::FitnessON(int SelectedSNPID){
	//double *PopulationFitnessTable = new double[numberofSNPs*diploidSize];	
	std::cout << std::endl << "Putting Fitness ON" << std::endl;
	/*for ( int i = 0; i <  numberofSNPs; i++ ){
		for (int j = 0; j < diploidSize ; j++){ 
			PopulationFitnessTable[i*diploidSize + j] = 0;	
		}
	}*/
	for (int j = 0; j < diploidSize ; j++){ 
		std::string targetdiploid = diploidVec[j];
		std::size_t ct = std::count(targetdiploid.begin(),targetdiploid.end(),'T');
		//std::cout << std::endl << targetdiploid << " SelectedSNPID: " << SelectedSNPID << std::endl;
		if (ct == 2)
			PopulationFitnessTable[SelectedSNPID*diploidSize+j] = 2*FITNESS;
		else if (ct == 1)
			PopulationFitnessTable[SelectedSNPID*diploidSize+j] = FITNESS; 
		else 
			PopulationFitnessTable[SelectedSNPID*diploidSize+j] = 0.0;	
	}
	
		std::cout << std::endl << "Done with fitness" << std::endl;
}
		
/*void Individuals :: DisplayTable(double* PopulationFitnessTable, int NumberofSNPs){
	for (int i = 0; i < NumberofSNPs; i++){	
		for (int j = 0; j < diploidSize; j++){
				std::cout << PopulationFitnessTable[i*diploidSize + j] << " | " ;
			}
			std::cout << std::endl;
		}
	}
*/
//This gets the position of the allele by finding where it lies in the std::vector of diploids which is in GlobalIndivs class 
int Individuals :: getPosAllele(std::vector<string> Retrievedpair){
	
	string allelefront = Retrievedpair[1] + Retrievedpair[3];
    string alleleend = Retrievedpair[3] + Retrievedpair[1];
	int pos_allele;

    if (find (diploidVec.begin(), diploidVec.end(), allelefront) != diploidVec.end())

        pos_allele = distance(diploidVec.begin(), find (diploidVec.begin(), diploidVec.end(), allelefront));
    else
        pos_allele = distance(diploidVec.begin(), find (diploidVec.begin(), diploidVec.end(), alleleend));

	return pos_allele; 
	}
/*std::pair<bool,std::pair<AlleleInfo,AlleleInfo> > Individuals :: RetrieveListElement(int x){
	//std::cout << std::endl << "eg" << std::endl;
	return ParentInfoList[x]; 
}*/
//This method mutates the SNPs depending on the mutation rate specified by the user
//It creates 4 bins, one with no mutation and other into three parts where possible mutation might take place to the other three bases. 
//This is done everytime the parents sends a std::set of SNPs to the child 
//This also saves the SNPs before and after mutations, helping us track while building the ARGs

std::pair<std::vector<string>, int > Individuals :: Mutate(std::vector<string> haploids, double MutationRate, int L){

	int numSNP = haploids.size();
	int mutidx; 
	string targetallele;
	std::vector<string> tmpbaseVec(baseVec.size());
	std::vector<string> mutatedhaploids(numSNP); 
	std::pair<std::vector<string>, int> MutAlleleContainer;
	int index; 

	int localmut;
	if (MutationRate == 0){
		mutatedhaploids = haploids;
		localmut = 0;
	}
	else{
		
		int *RandNumbers = new int[numSNP];
		for (int z = 0; z < numSNP; z++) 
			RandNumbers[z] = z; 
	    double scaledMutRate = MutationRate*(double)L;
		unsigned int SNPsToMutate = gsl_ran_poisson(threadvec[omp_get_thread_num()],scaledMutRate);
		//std::cout << std::endl << "Number of SNPs to Mutate is " << SNPsToMutate << std::endl;
		
		int *MutateLocs= new int[SNPsToMutate];

		tmpbaseVec = baseVec;
		double randmut = gsl_rng_uniform(threadvec[omp_get_thread_num()]);
		gsl_ran_shuffle(threadvec[omp_get_thread_num()],RandNumbers, numSNP, sizeof (int));
		gsl_ran_choose(threadvec[omp_get_thread_num()], MutateLocs, SNPsToMutate, RandNumbers, numSNP, sizeof (int));
		
		mutatedhaploids = haploids;
		localmut = SNPsToMutate;
		/*for (int la = 0; la < SNPsToMutate; la++)
			std::cout  << MutateLocs[la] << " ";
		std::cout << std::endl;*/
		for (int z = 0; z < SNPsToMutate; z++){
			//#pragma omp critical
			if(flagmut == 0){
				if((EpiFit2 + EpiFit1) !=0){
					std::vector<int> MutLocs(MutateLocs,MutateLocs+SNPsToMutate);
					std::vector<int> tmpmutloc = MutLocs; 
					int ct = 0; 
					std::vector<int> tmpid; 
					for(std::vector<std::vector<int>>::iterator it = newSNPs.begin(); it != newSNPs.end(); ++it){					
						if(std::find(it->begin(),it->end(),MutateLocs[z]) != it->end()){
							tmpid = *it;
							tmpmutloc.erase(tmpmutloc.begin()+z);
							break;
						}
					}
					for (int ml = 0; ml < tmpid.size(); ml++)
						std::cout << " " << tmpid[ml]; 				
					std::cout << std::endl; 
					for (int ml = 0; ml < tmpid.size(); ml++){
						if(std::find(MutLocs.begin(), MutLocs.end(), tmpid[ml]) != MutLocs.end()){
							ct++;
						}
					}
					if(tmpid.size() > 1){
						if(ct == 2)
						PopulateEFT(tmpid,EpiFit1,EpiFit2);
						else if(ct == 1){
							double r = gsl_rng_uniform(threadvec[omp_get_thread_num()]);
							if (r < 0.5)
								PopulateEFT(tmpid,EpiFit1,0.0);
							else
								PopulateEFT(tmpid,0.0,EpiFit2);
						}
					} 
					else{
						SelectedSNPID = MutateLocs[0];
						FitnessON(SelectedSNPID);
					}
					//std::cout << std::endl << "Updated at " << SelectedSNPID << std::endl;
					flagmut = 1;
					MutLocs.clear(); 
					tmpid.clear();
					tmpmutloc.clear();
					std::cout << std::endl << "signing off " << std::endl;
				}
				else{
					SelectedSNPID = MutateLocs[z];
					FitnessON(SelectedSNPID);
					//std::cout << std::endl << "Updated at " << SelectedSNPID << std::endl;
					flagmut = 1;
				}
			}
			
			targetallele = mutatedhaploids[RandNumbers[z]];
			for( int kk = 0; kk < baseVec.size(); kk++){
				if(baseVec[kk].compare(targetallele) == 0)
					tmpbaseVec.erase(tmpbaseVec.begin()+kk);
			}
			if (randmut > 0 && randmut < 1/3) {
				mutatedhaploids[RandNumbers[z]] = tmpbaseVec[0]; 
				MutCount++;
			}
			else if (randmut >= 1/3 && randmut < 2/3) {
					mutatedhaploids[RandNumbers[z]] = tmpbaseVec[1];
					MutCount++;
			}
			else if (randmut >= 2/3 && randmut < 1) {
					mutatedhaploids[RandNumbers[z]] = tmpbaseVec[2]; 
					MutCount++;
			}
		}
		delete[] RandNumbers;
		delete[] MutateLocs;
	}
	
	MutAlleleContainer = std::make_pair(mutatedhaploids,localmut);
	//if (flagmut == 1)
	//	SelectedSNPID = -1;
	mutatedhaploids.clear();	
	return MutAlleleContainer; 
}
//This method gets the alleles for the child when there's no recombination. This can happen two ways, either father sends all of his father's SNPs to child, or,
//the father sends all of his mother's. After it's decided which, the information is kept in the container called Contributing Chromosomes.
//The chromosome which is contributing nothing to the child has 0 in its position, indicating it's contribution. This is helpful while building the ARG. 
//It collapses the edges when there's no recombination.
AlleleInfo Individuals :: GetAllelesforChild(int Pidx, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC, std::vector<int> SNPlocs, double MutationRate,int L, int GenNum, int chromnum){
  
 //	std::cout << std::endl << "NO CROSSOVER" << std::endl;
	AlleleInfo toChild; 
	int numSNP = SNPlocs.size();
	std::pair<std::pair<string,int>,std::pair<string,int> > ContributingChroms; 	
	std::pair<string,int> c1,c2; 
	int i;
	std::vector<string> MutatedHaploids;
	int localmut;
			string chromID,theotherID;	
		std::vector<string> haploids(numSNP); 

	if (chromnum == 1){
	for (int z = 0; z < numSNP; z++){
		std::pair<std::pair<string,string>, std::pair<string,string> > chrompair = CC[Pidx*numSNP + z];
			haploids[z] = chrompair.first.second;
			if (z==0){
				chromID = chrompair.first.first; 
				theotherID = chrompair.second.first; 
			}
		}
	}
	else{	
	for (int z = 0; z < numSNP; z++){
		std::pair<std::pair<string,string>, std::pair<string,string> > chrompair = CC[Pidx*numSNP + z];

			haploids[z] = chrompair.second.second;
			if (z==0){
				chromID = chrompair.second.first;
				theotherID = chrompair.first.first; 
			}	
		}
	}
	string chromname; 
	c1 = std::make_pair(chromID,numSNP);
	c2 = std::make_pair(theotherID,0);
	//They pack the chromosome contributions at this stage
	ContributingChroms = std::make_pair(c1,c2);
	chromname = AddandGetChromName(GenNum);
	//}
		std::string s = "";
		for (const auto &piece1 : haploids) s += piece1;

	std::pair<std::vector<string>, int >MutBox;
	
	 MutBox = Mutate(haploids,MutationRate,L);
	
	MutatedHaploids = MutBox.first; 
	localmut = MutBox.second;
	
	//MutatedposInfo = MutBox.second; 
	Events event; 
	int localrecomb = 0;
	int pos = -1; //might change afterwards, don't worry about this
	
	toChild = event.getandPopulateStruct(pos, Pidx, chromname, ContributingChroms, MutatedHaploids, localmut, localrecomb);
	
	haploids.clear();
	MutatedHaploids.clear();

	return toChild;
}

//This builds the list for each individual in a given generation with all their parent info in it. Easier to access any given individual's info
//while building the ARG (although, at this stage I am not using this, but it's an useful data structure to have)
/*void Individuals :: PushToList(std::pair<bool, std::pair <AlleleInfo,AlleleInfo> > structstd::pairs){
	this->ParentInfoList.push_back(structstd::pairs);
}*/

//This displays the aforementioned list for each individual
void Individuals :: DisplayParentList(){
	for(std::vector<std::pair<bool, std::pair<AlleleInfo,AlleleInfo> > >::iterator it = ParentInfoList.begin(); it != ParentInfoList.end(); ++it){
	
		std::cout << std::endl << " | " << "FatherID: "  << it->second.first.ParentID << " | " << "Chromosome ID: " << it->second.first.ToTheChild.first << " | " << std::endl;
		std::cout << std::endl <<  "From Father ------ Chromosome: " << it->second.first.ContChrom.first.first << " Contributing: " << it->second.first.ContChrom.first.second << " and Chromosome: " << it->second.first.ContChrom.second.first << " Contributing: " << it->second.first.ContChrom.second.second << std::endl; 

	
		std::cout << std::endl << " | " << "MotherID: "  << it->second.second.ParentID << " | " << "Chromosome ID: " << it->second.second.ToTheChild.first << " | " << std::endl;
		std::cout << std::endl <<  "From Mother ------ Chromosome: " << it->second.second.ContChrom.first.first << " Contributing: " << it->second.second.ContChrom.first.second << " and Chromosome: " << it->second.second.ContChrom.second.first << " Contributing: " << it->second.second.ContChrom.second.second << std::endl; 

	}
}	

/*This method gets the alleles for the child when there's recombination. This can happen two ways, either father sends some (depending on the random 'C' that is chosen in L) of his father's SNPs to child
and the rest is chosen from his mother's, OR, vice versa.  After it's decided which, the information is kept in the container called Contributing Chromosomes, with positive number in the position
indicating forward contribution and the number of SNPs contributed, the negative number in pos indicates backward contribution and the number of SNPs contributed.
There will be no 0 in the position, this time. This is helpful while building the ARG. 
It collapses the edges when there's no recombination. */

AlleleInfo Individuals :: GetCrossOverAlleles(int Pidx, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC, std::vector<int> SNPlocs, double MutationRate, int L, int GenNum, int chromnum){
//	std::cout << std::endl << "YES CROSSOVER" << std::endl; 
	int pos;
	int numSNP = SNPlocs.size();
	std::vector <string> haploids(numSNP); 
	std::vector<int> SNPidx(numSNP); 
	std::pair<std::pair<string,string>,std::pair<string,string> > chrompair;
	string chromID, chromname;
	std::pair<std::pair<string,int>,std::pair<string,int> > ContributingChroms; 
	std::pair<string,int> c1,c2;
	std::pair<std::vector<string>, int > MutBox;
	std::vector<string> MutatedHaploids;
	Events event;
	std::string s;
	int localrecomb = 0;
	int localmut;
	//std::vector<std::pair<string,int> > MutatedposInfo;
	AlleleInfo toChild; 
	//randomly choses C out of L
	//The first part is crossover1 and the second part is crossover2, delineated with the if-else
	//int crossoveridx = RandU(0,1);
	//int crossoveridx = gsl_rng_uniform(threadvec[omp_get_thread_num()])%(L+1);
	int crossoveridx = gsl_rng_uniform_int(threadvec[omp_get_thread_num()], L);
	
	//int crossoveridx = rand()%(L + 1);
	//std::cout << std::endl << "Crossoveridx was: " << crossoveridx << endl;
	if (chromnum == 1){
		
			int ct = 0;
			for (int i = 0; i < numSNP; i++){
				if (SNPlocs[i] <= crossoveridx)
				ct++;	
			}
			SNPidx.resize(ct);
			if(SNPidx.size() == 0){
				toChild = GetAllelesforChild(Pidx,CC,SNPlocs,MutationRate,L,GenNum,2);

			}
			else{
				localrecomb++;
				for(int z = 0; z < SNPidx.size(); z++){ 
			
					chrompair = CC[Pidx*numSNP + z];
					haploids[z] = chrompair.first.second;
		
					if (z==0)
						chromID = chrompair.first.first;  
				}
				pos = SNPidx.size();
				c1 = std::make_pair(chromID,pos);
				SNPidx.clear();
				SNPidx.resize(numSNP-pos);

				if(SNPidx.size() == 0){
					//RecombCount++;
					chrompair = CC[Pidx*numSNP];
					chromID = chrompair.second.first;
					c2 = std::make_pair(chromID, 0);
				}
			    else{
			
					for (int z = 0; z < SNPidx.size() ; z++){
						chrompair = CC[Pidx*numSNP + pos+z];	
						haploids[pos+z] = chrompair.second.second;
				
						if (z==0)
							chromID = chrompair.second.first;
					}
					pos = SNPidx.size();
					RecombCount++;
					c2 = std::make_pair(chromID,-pos);	
				}
		
				ContributingChroms = std::make_pair(c1,c2);
				MutBox = Mutate(haploids,MutationRate,L);
				MutatedHaploids = MutBox.first; 
				localmut = MutBox.second; 
				
				chromname = AddandGetChromName(GenNum);
		//passed on to the child in the AlleleInfo container
				toChild = event.getandPopulateStruct(pos, Pidx, chromname, ContributingChroms, MutatedHaploids, localmut, localrecomb);
		
				haploids.clear();
				haploids.resize(numSNP);
				SNPidx.clear();
				
			}
		}
	else {
			
			int ct = 0;
			localrecomb++;
			for (int i = 0; i < numSNP; i++){
				if (SNPlocs[i] <= crossoveridx)
				ct++;	
			}
			SNPidx.resize(ct);
			if(SNPidx.size() == 0){
				toChild = GetAllelesforChild(Pidx,CC,SNPlocs,MutationRate,L,GenNum,1);
			}
			else{
		
				for(int z = 0; z < SNPidx.size(); z++){
				
					chrompair = CC[Pidx*numSNP + z];
					haploids[z] = chrompair.second.second;
					if (z==0)
						chromID = chrompair.second.first;
				}
				pos = SNPidx.size();
				c1 = std::make_pair(chromID,pos);
				SNPidx.clear();
				SNPidx.resize(numSNP-pos);
			
				if(SNPidx.size() == 0){
					//RecombCount++;
					chrompair=CC[Pidx*numSNP];
					chromID = chrompair.first.first;
					c2 = std::make_pair(chromID, 0);
				}
				else{
			
					for (int z = 0; z < SNPidx.size(); z++){
						chrompair = CC[Pidx*numSNP + pos+z];
						haploids[pos+z] = chrompair.first.second;
				
						if (z==0)
						chromID = chrompair.first.first;
					}
					pos = SNPidx.size();
					c2 = make_pair(chromID,-pos);
			//packs everything into contributing chromosome container indicating positions.
					RecombCount++;
				}	
				//s = "";
				//for (const auto &piece1 : haploids) s += piece1;
			//	std::cout << std::endl << "haploids: " << s << std::endl;
				ContributingChroms = std::make_pair(c1,c2);
 				MutBox = Mutate(haploids,MutationRate,L);
				MutatedHaploids = MutBox.first; 
				localmut = MutBox.second; 
				chromname = AddandGetChromName(GenNum);
				
				toChild = event.getandPopulateStruct(pos, Pidx, chromname, ContributingChroms, MutatedHaploids, localmut, localrecomb);
				
				haploids.clear();
				haploids.resize(numSNP);
				
				SNPidx.clear();
	
		//Mutates the SNPs before passing them on to the child as described in the Mutate method
			}
	}
	return toChild;
}
/*This is the method which executes whether there was a recombination or not, depending on the recombination rate fed by the user. This creates 4 equidistant
bins and tosses a coin to decide whether it's crossover or not and whatever the result whether it comes from the father's strand or mother's strand of SNPs.
Then it obtains the AlleleInfo container after invoking the aforementioned methods to get the haploids to the child and passes it on as a return value */
AlleleInfo Individuals :: AskChromosome (int Pidx, std::vector<std::pair<std::pair<string,string>, std::pair<string,string> > > CC, std::vector<int> SNPlocs, double rrate, double MutationRate, int L, int GenNum){ 
		int recomidx; 	
		//Coin Toss
		string s;
		//std::cout << std::endl << "lo" << std::endl;
		double randrecom = gsl_rng_uniform(threadvec[omp_get_thread_num()]);
		//std::cout << std::endl << "randrecom is " << randrecom << std::endl; 
		//4 bins and seeing where the coin lands

		double RecombRate = rrate*L;
		//std::cout << std::endl << "RecombRate: " << RecombRate << std::endl;
		//Depending on the bin where it landed it decides what to do and send the result to the child
		AlleleInfo toChild;

	
		if (randrecom >= 0 && randrecom < ((1-RecombRate)/2)) toChild = GetAllelesforChild(Pidx,CC,SNPlocs,MutationRate,L,GenNum,1);
		else if (randrecom >= ((1-RecombRate)/2) && randrecom < (1-RecombRate)) toChild = GetAllelesforChild(Pidx,CC,SNPlocs,MutationRate,L,GenNum,2);
		else if (randrecom >= (1-RecombRate) && randrecom < (1-(RecombRate/2))) toChild = GetCrossOverAlleles(Pidx,CC,SNPlocs,MutationRate,L,GenNum,1);
		else if (randrecom >= (1-(RecombRate/2)) && randrecom < 1) toChild = GetCrossOverAlleles(Pidx,CC,SNPlocs,MutationRate,L,GenNum,2); 
		return toChild; 
	}
	 
