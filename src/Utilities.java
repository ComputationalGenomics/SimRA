/*
  Copyright 2015 IBM Corporation


Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under 
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF 
ANY KIND, either express or implied. 
See the License for the specific language governing permissions and limitations under the License.


*/


import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

/**
 * The class Utilities has many useful functions and procedures that are called and used in the rest of the program
 * 
 * @author Anna Paola Carrieri
 */

public class Utilities {
	

	/**
	 * Given the active lineages l1, l2, .., L the function returns the sum of all the recombination rate
	 * r_l (for l = 1, .., L) associated to each active lineage
	 * @param activeEdges ArrayList of Integer objects representing the ID of the active lineages in the ARG
	 * @param graphEdges Map of containing all edges of the ARG created during the backward simulation. 
	 * The key field is the Integer representing the ID of the edge
	 * @return a real value that is the sum of all the recombination rate
	 * r_l (for l = 1, .., L) fol active lineage l
	 */
	public static double computeSumRates(ArrayList<Integer> activeEdges, Map<Integer, Edge> graphEdges){
		double sum = 0;
		
		for(int i = 0; i < activeEdges.size(); i++) {
			sum = graphEdges.get(activeEdges.get(i)).getRate()+sum;
		}
		
		return sum; 
	}
	/**
	 * This function implements the exponential distribution with lambda as parameter
	 * @param lambda parameter of the exponential distribution
	 * @return the time t of the next event
	 */
	public  static double computeNextTimeEvent(double lambda){
		double random = Math.random();
		
		//x = ln(1-p)/-lamda
		double t = (Math.log(1-random))/(0-lambda);
		
		return t;
	}
	
	
	/** 
	 * This computes the time T = T+t to the next event using the exponential distribution
	 * @param time t computed by exponential distribution time is in [0,1)
	 * @param currentgeneration is the current age of the ARG (not scaled in generations)
	 * @return the time T updated based on to the current level or generation and the time computed by the exponential distribution
	 */
	public static  double computeNextGeneration(double time, double currentgeneration){
		return time+currentgeneration;
	}
	
	/**
	 * This function computes the kind of the next event (coalescent if
	 * index is equal to 0, recombination if index is greater than  0)
	 * @param lambda is the coefficient of the exponential distribution
	 * @param binCoef (L 2) binomial coefficient where L is the number of the active lineages
	 * @param activeEdges list of the active lineages
	 * @param graphEdges map containing the edges in the whole ARG
	 * @return the index indicating the type of the next event
	 */
	public static int computeIndexNextEvent(double lambda, double binCoef, ArrayList<Integer> activeEdges, HashMap<Integer, Edge> graphEdges){
		
		double x = Math.random()*lambda;
		double left = 0;
		double right = binCoef;
		boolean found = false;
		
		//0 is associated with the coalescent event, from 1 to L with recombination events
		int index = 0;
		
		while(!found) {
			
			if(x >= left && x < right) {
				found = true;
			}
			else{
				index++;
				left = right;
				right = right+graphEdges.get(activeEdges.get(index-1)).getRate();
			}
		}
		return index;
	}
	
	/**
	 * Function that returns the binomial coefficient (n r)
	 * @param n integer value = number of elements in the set X
	 * @param r integer value = number of elements in each subset of X
	 * @return binomial coefficient = number of distinct k-elements subsets of X
	 */
	public static double binomialCoefficient(int n, int r) {
	        double t = 1;
	        
	        int m = n - r; // r = Math.max(r, n - r);
	        if (r < m) {
	            r = m;
	        }
	        
	        for (int i = n, j = 1; i > r; i--, j++) {
	            t = t * i / j;
	        }
	        
	        return t;
	    }
	
	/**
	 * This procedures prints on video the number of SNPs on the edges that are involved in the split of one population 
	 * @param pop object of PopulationARG
	 * @see PopulationARG class
	 * @param num_edges_Splitted number of edges involved in the splitting of the population pop
	 */
	public static void printMutationOnEdgesFromSplitting(PopulationARG pop, int num_edges_Splitted){
		
		for(int i = 0; i < num_edges_Splitted; i++){
			System.out.println("Edge "+i);
			System.out.println("Number of mutations: "+pop.getGraphEdges().get(i).getnMut());
		}
		
	}
	
	/**
	 * This procedure prints on video the number and the list of mutation for each leaf of the ARG
	 * @param pop object of PopulationARG
	 * @see PopulationARG class
	 */
	public static void printMutationsForLeaf(PopulationARG pop){
		
		for(int i = 0; i < pop.getExtantUnits(); i++){
			TreeSet<Double> muts = pop.getNodeSet().get(i).getMutation_set();
			System.out.println("EXTANT UNIT "+i);
			System.out.println("Number of mutations : "+muts.size());
			System.out.println("List of mutations :");
			Iterator<Double> muts_it = muts.iterator();
			while(muts_it.hasNext()){
				System.out.println(muts_it.next());
			}
		}
	}
	
	/**
	 * This procedure prints the list of all mutations in one ARG
	 * @param pop object of PopulationARG
	 * @see PopulationARG class
	 */
	public static void printAllMutations(PopulationARG pop){
		Map<Double,Mutation> muts = pop.getMutationSet();
		Iterator<Double> it_muts = muts.keySet().iterator();
		while(it_muts.hasNext()){
			Double id_mut = it_muts.next();
			Mutation mut = muts.get(id_mut);
			mut.printMutation();
		}
		
	}
	
} //end class


//*************************** SOME OLD PROCEDURES USEFUL FOR DEBUGGING *************************** 
	/*import java.io.BufferedWriter;
	import java.io.FileNotFoundException;
	import java.io.FileOutputStream;
	import java.io.IOException;
	import java.io.OutputStreamWriter;
	import java.io.UnsupportedEncodingException;
	import java.io.Writer;*/
	
	/*public static void printOnFileAllSNPsPerPopulation(String fileName, PopulationARG pop){
	
	Writer writer = null;

		
	    try {
			writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), "utf-8"));
			
			writer.write("-------- SNPs POPULATION "+pop.getId_pop()+" -------- \n\n");
			TreeSet<Double> list_SNPs = pop.getSNPpositionsList();
			writer.write("Total number of SNPs : "+list_SNPs.size()+"\n\n");
			writer.write("List of SNPs\n\n");
			
			Iterator<Double> it_snps = list_SNPs.iterator();
			int i = 1;
			while(it_snps.hasNext()){
				writer.write("#"+i+" : "+it_snps.next()+"\n");
				i++;
			}
			
			writer.close();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
}*/

/*	public static void debugMutationsSinglePop(){
		PopulationARG pop1 = new PopulationARG();
		pop1.generateLeafPopulation(1, 1, 10000, 5, 0.1, 0.7, 75);
		pop1.printKind();
		CreatingFilesForStructure.createTxtFile_PopulationStructure("pop1.txt", pop1);
		DecorateMutationsSNP.decoratingWithMutations(pop1.getMu(), pop1);
		System.out.println("ALL MUTATIONS: ");
		Utilities.printAllMutations(pop1);
		System.out.println("ALL LEAVES: ");
		Utilities.printMutationsForLeaf(pop1);
	}
	
public static void debugMutationsCoalescence(){
		
		PopulationARG pop1 = new PopulationARG();
		PopulationARG pop2 = new PopulationARG();
		PopulationARG pop3 = new PopulationARG();
		
		pop1.generateLeafPopulation(1, 0.5, 10000, 4, 0.1, 0.7, 75);
		pop1.printKind();
		CreatingFilesForStructure.createTxtFile_PopulationStructure("pop1.txt", pop1);
		DecorateMutationsSNP.decoratingWithMutations(pop1.getMu(), pop1);
		System.out.println("------------ Mutations of POP1 ------------");
		Utilities.printMutationsForLeaf(pop1);
		System.out.println("-------------------------------------------");
		
		pop2.generateLeafPopulation(2, 0.5, 10000, 4, 0.1, 0.7, 75);
		pop2.printKind();
		CreatingFilesForStructure.createTxtFile_PopulationStructure("pop2.txt", pop2);
		DecorateMutationsSNP.decoratingWithMutations(pop2.getMu(), pop2);
		System.out.println("------------ Mutations of POP2 ------------");
		Utilities.printMutationsForLeaf(pop2);
		System.out.println("-------------------------------------------");
		
		pop3.generateCoalescentPopulation(3, pop1, pop2, 0.5, 2, 10000, 0.1, 0.7, 75);
		pop3.printKind();
		CreatingFilesForStructure.createTxtFile_PopulationStructure("pop3.txt", pop3);
		DecorateMutationsSNP.decoratingWithMutations(pop3.getMu(), pop3);
		
		
		DecorateMutationsSNP.updateMutationsForTheBottomPop(pop3);	
		
		
		System.out.println("------------ Mutations of POP3 ------------");
		System.out.println("extant units pop3 = "+pop3.getExtantUnits());
		Utilities.printMutationsForLeaf(pop3);
		System.out.println("-------------------------------------------");
		
		System.out.println("----------- UPDATED MUTATIONS-------------");
		
		System.out.println("------------ Mutations of POP1 ------------");
		Utilities.printMutationsForLeaf(pop1);
		System.out.println("-------------------------------------------");
		
		System.out.println("------------ Mutations of POP2 ------------");
		Utilities.printMutationsForLeaf(pop2);
		System.out.println("-------------------------------------------");
	}

	public static void debugMutationsSplit(){
		
		PopulationARG pop1 = new PopulationARG();
		
		pop1.generateLeafPopulation(1, 0.5, 10000, 5, 0.1, 0.7, 75);
		pop1.printKind();
		CreatingFilesForStructure.createTxtFile_PopulationStructure("pop1.txt", pop1);
	
		//split population 2 at time t
		ArrayList<ArrayList<Integer>> left_right_activeEdges = pop1.splitPopulation();
		if(left_right_activeEdges.size()==2){
			System.out.println("Splitting population 1");
			
			PopulationARG pop2 = new PopulationARG();
			pop2.generatePopulationFromSplitting(2, pop1, left_right_activeEdges.get(0), 0.5, 2, 10000, 0.1, 0.7, 75);
			CreatingFilesForStructure.createTxtFile_PopulationStructure("pop2.txt", pop2);
			
			
			PopulationARG pop3 = new PopulationARG();
			pop3.generatePopulationFromSplitting(3, pop1, left_right_activeEdges.get(1), 0.5, 2, 10000, 0.1, 0.7, 75);
			CreatingFilesForStructure.createTxtFile_PopulationStructure("pop3.txt", pop3);
			
			
			System.out.println("End Splitting");
			
			
			DecorateMutationsSNP.decoratingWithMutations(pop1.getMu(), pop1);
			System.out.println("------------ Mutations of POP1 ------------");
			Utilities.printMutationsForLeaf(pop1);
			System.out.println("-------------------------------------------");
			
			DecorateMutationsSNP.decoratingWithMutations(pop2.getMu(), pop2);
			System.out.println("------------ Mutations of POP2 ------------");
			Utilities.printMutationsForLeaf(pop2);
			System.out.println("-------------------------------------------");
			
			DecorateMutationsSNP.decoratingWithMutations(pop3.getMu(), pop3);
			System.out.println("------------ Mutations of POP3 ------------");
			Utilities.printMutationsForLeaf(pop3);
			System.out.println("-------------------------------------------");
			
			DecorateMutationsSNP.updateMutationsForTheBottomPop(pop2);	
			DecorateMutationsSNP.updateMutationsForTheBottomPop(pop3);	
			
			
			System.out.println("----------- UPDATED MUTATIONS-------------");
			
			System.out.println("------------ Mutations of POP1 ------------");
			Utilities.printMutationsForLeaf(pop1);
			System.out.println("-------------------------------------------");
				
		}
		else{
			System.out.println("The population has not been splitted. Exit!");
		}	
			
	} */
