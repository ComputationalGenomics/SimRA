/*
  Copyright 2015 IBM Corporation


Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under 
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF 
ANY KIND, either express or implied. 
See the License for the specific language governing permissions and limitations under the License.


*/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

/**
 * This is the main class that implements SimRA algorithm. This class allows to  simulating complex scenarios of related multiple populations
 * with subdivision and admixture. 
 * It provides all the output files that describe the resulting populations in terms of SNPs and STRs and the relevant evolutionary history of them.
 * generate SNPs for admixture populations as specified in the text file in input
 * 
 * @author Anna Paola  Carrieri
 */
public class GenerateAdmixturePopulations {
	
				// ***** PRIVATE STATIC FIELDS *****
	
	//map containing all the populations that form the scaffold					
	private static HashMap<Integer, PopulationARG> populations;
	//list of ordered IDs of the populations based on the ascending order in which they have been created
	private static ArrayList<Integer> ordered_id_populations;
	//list of events that form the scaffold (merge and split of populations)
	private static ArrayList<Event> events;
	//total number of populations in the scaffold
	private static int num_tot_pops = 0;
	//split points (due to recombination events) in the whole scaffold
	private static TreeSet<Double> splitPoints;
	//all SNP mutation positions in the segment that decorate the whole scaffold
	private static TreeSet<Double> allMutations;
	//set of all instances of the class Mutation
	private static Map<Double,Mutation> allmutationSet;
	//number of contemporary populations in the scaffold
	private static int num_actual_pops = 0;
	//length of the segment
	private static int g = 0;
	//recombination rate
	private static double r = 0;
	//mutation rate
	private static double mu = 0;
	//total number of recombination events
	private static int recomb_total_number;
	//total number of coalescent events
	private static int coalesc_total_number;
	//total number of edges that have at least one mutations
	private static int num_edges_more_than_one_mut;
	//total number of edges that have no mutations
	private static int num_edges_without_mut;
	
	
					// ***** GET AND SET FUNCTIONS *****
	/**
	 * The function returns the list of the events that form the scaffold
	 * @return an ArrayList of objects of the class Event 
	 */
	public static ArrayList<Event> getEvents() {
		return events;
	}
	/**
	 * The procedure sets the list of the events (objects of the class Event) that form the scaffold
	 * @param events ArrayList of the events (objects of the class Event) that form the scaffold
	 */
	public static void setEvents(ArrayList<Event> events) {
		GenerateAdmixturePopulations.events = events;
	}
	/**
	 * The functions returns the total number of populations/edges in the scaffold
	 * @return an integer that is the total number of populations/edges in the scaffold
	 */
	public static int getNum_tot_pops() {
		return num_tot_pops;
	}
	/**
	 * The procedure sets the total number of populations in the scaffold 
	 * @param num_tot_pops total number of populations/edges in the scaffold
	 */
	public static void setNum_tot_pops(int num_tot_pops) {
		GenerateAdmixturePopulations.num_tot_pops = num_tot_pops;
	}
	/**
	 * This function returns the number of leaves (contemporary populations) in the scaffold
	 * @return an integer that is the number of leaves in the scaffold
	 */
	public static int getNum_actual_pops() {
		return num_actual_pops;
	}
	/**
	 * This procedure sets the number of leaves (contemporary populations) in the scaffold
	 * @param num_actual_pops integer representing the number of leaves (contemporary populations) in the scaffold
	 */
	public static void setNum_actual_pops(int num_actual_pops) {
		GenerateAdmixturePopulations.num_actual_pops = num_actual_pops;
	}
	/**
	 * The function returns the length of the segment length by the edges and nodes
	 * @return segment length integer representing the segment length carried by the edges and nodes
	 */
	public static int getG() {
		return g;
	}
	/**
	 * The procedure sets the segment length carried by the edges and nodes
	 * @param g segment length integer representing the length of the segment carried by the edges and nodes
	 */
	public static void setG(int g) {
		GenerateAdmixturePopulations.g = g;
	}
	/**
	 * This function returns the recombination rate used during the backward simulation of the whole scaffold
	 * @return r real value representing the recombination rate 
	 */
	public static double getR() {
		return r;
	}
	/**
	 * This procedure sets the recombination rate used during the backward simulation of the whole scaffold
	 * @param r = real value representing the recombination rate 
	 */
	public static void setR(double r) {
		GenerateAdmixturePopulations.r = r;
	}
	
	/**
	 * This function returns the mutation rate used during the backward simulation of the whole scaffold
	 * @return mutation rate real number representing the mutation rate
	 */
	public static double getMu() {
		return mu;
	}
	/**
	 * This procedure sets the mutation rate used during the backward simulation of the whole scaffold
	 * @param mu real number representing the mutation rate
	 */
	public static void setMu(double mu) {
		GenerateAdmixturePopulations.mu = mu;
	}
	
	/**
	 * This function returns a map containing all the populations (objects of the class PopulationARG) of the scaffold, the key field is the population ID
	 * @return a map of pairs where the key field is the ID of the single ARG or population while the objects are the ARGs (objects of the class PopulationARG)
	 * @see PopulationARG
	 */
	public static HashMap<Integer, PopulationARG> getPopulations() {
		return populations;
	}
	/**
	 * This procedure sets the map containing all the populations (objects of the class PopulationARG) of the scaffold, the key field is the population ID
	 * @param pops map of pairs where Integer is the ID of PopulationARG and an object of the class PopulationARG represents a single population/edge in the scaffold
	 * @see PopulationARG
	 */
	public static void setPopulations(HashMap<Integer, PopulationARG> pops) {
		populations = pops;
	}
	/**
	 * This function returns the set of all SNP mutations.
	 * @return the set of all SNP mutations. Each value in the set is an ID of a mutation and corresponds to its
	 * position in the interval (0,1)
	 */
	public static TreeSet<Double> getAllMutations() {
		return allMutations;
	}
	/**
	 * This function sets the set of all SNP mutations in the scaffold
	 * @param allMutations object of the class TreeSet containing all SNP mutations in the scaffold ordered by their IDs
	 */
	public static void setAllMutations(TreeSet<Double> allMutations) {
		GenerateAdmixturePopulations.allMutations = allMutations;
	}
	/**
	 * This function returns the ordered list of population IDs. More precisely the IDs are ordered based on the order in which the populations have been created
	 * @return the ordered list of the population IDs
	 */
	public static ArrayList<Integer> getOrdered_id_populations() {
		return ordered_id_populations;
	}
	
	/**
	 * The procedure sets the ordered list of population IDs. The order depends on the order in which the populations have been generated
	 * @param ordered_id_populations list of Integer, each values correspond to an ID of a population
	 */
	public static void setOrdered_id_populations(
			ArrayList<Integer> ordered_id_populations) {
		GenerateAdmixturePopulations.ordered_id_populations = ordered_id_populations;
	}
	/**
	 * This function returns all the split points generated by recombination events during the backward simulation of the whole scaffold
	 * @return an ordered set of split points generated by recombination events during the backward simulation of the whole scaffold
	 */
	public static TreeSet<Double> getSplitPoints() {
		return splitPoints;
	}
	/**
	 * This procedure sets the ordered set of all the split points generated by recombination events during the backward simulation of the whole scaffold
	 * @param splitPoints a TreeSet of objects of the class Double. Each instance of Double represents a position in the interval (0,1) where a segment has been split by a recombination event
	 */
	public static void setSplitPoints(TreeSet<Double> splitPoints) {
		GenerateAdmixturePopulations.splitPoints = splitPoints;
	}
	/**
	 * This function returns a map, where the key is an object Double representing the ID of a SNP mutation (its position in the interval (0,1)), while the element is an object of the class Mutation
	 * @see Mutation
	 * @return the set of all instances of the class Mutation in the whole scaffold, that contain informations about the SNPs in the scaffold
	 */
	public static Map<Double,Mutation> getAllmutationSet() {
		return allmutationSet;
	}
	/**
	 * This procedure sets the set of all instances of the class Mutation in the whole scaffold, that contain informations about the SNPs in the scaffold
	 * @param allmutationSet map, where the key is an object Double representing the ID of a SNP mutation (its position in the interval (0,1)), while the element is an object of the class Mutation
	 */
	public static void setAllmutationSet(Map<Double,Mutation> allmutationSet) {
		GenerateAdmixturePopulations.allmutationSet = allmutationSet;
	}
	/**
	 * This function returns the total number of recombination events generated during the backward simulation of the scaffold
	 * @return an integer representing the total number of recombination events (that is recombination nodes) in the whole scaffold
	 */
	public static int getRecomb_total_number() {
		return recomb_total_number;
	}
	/**
	 * This procedure sets the total number of recombination events generated during the backward simulation of the scaffold
	 * @param recomb_total_number total number of recombination events (that is recombination nodes) in the whole scaffold
	 */
	public static void setRecomb_total_number(int recomb_total_number) {
		GenerateAdmixturePopulations.recomb_total_number = recomb_total_number;
	}
	/**
	 * This function returns the total number of coqlescent events generated during the backward simulation of the scaffold
	 * @return an integer representing the total number of coalescent events (that is coalescent nodes) in the whole scaffold
	 */
	public static int getCoalesc_total_number() {
		return coalesc_total_number;
	}
	/**
	 * This procedure sets the total number of coalescent events generated during the backward simulation of the scaffold
	 * @param coalesc_total_number total number of coalescent events (that is coalescent nodes) in the whole scaffold
	 */
	public static void setCoalesc_total_number(int coalesc_total_number) {
		GenerateAdmixturePopulations.coalesc_total_number = coalesc_total_number;
	}
	/**
	 * This function returns the total number of edges in the whole scaffold that have at least one mutation
	 * @return total number of edges in the whole scaffold that have at least one mutation
	 */
	public static int getNum_edges_more_than_one_mut() {
		return num_edges_more_than_one_mut;
	}
	/**
	 * This procedure sets the total number of edges in the whole scaffold that have at least one mutation
	 * @param num_edges_more_than_one_mut: total number of edges in the whole scaffold that have at least one mutation
	 */
	public static void setNum_edges_more_than_one_mut(
			int num_edges_more_than_one_mut) {
		GenerateAdmixturePopulations.num_edges_more_than_one_mut = num_edges_more_than_one_mut;
	}
	/**
	 * This function returns the total number of edges in the whole scaffold that have no mutation
	 * @return total number of edges in the whole scaffold that have no mutation
	 */
	public static int getNum_edges_without_mut() {
		return num_edges_without_mut;
	}
	/**
	 * This procedure sets the total number of edges in the whole scaffold that have no mutation
	 * @param num_edges_without_mut total number of edges in the whole scaffold that have no mutation
	 */
	public static void setNum_edges_without_mut(int num_edges_without_mut) {
		GenerateAdmixturePopulations.num_edges_without_mut = num_edges_without_mut;
	}
	
	/**
	 * Main procedure: it takes in input a txt file where there are specified all the parameters together with the description of how and when 
	 * the populations merge and split. The main procedure reads the input file, generates the structure of the scaffold and finally decorates it with SNPs mutations. 
	 * For each actual population this procedure creates a .txt file 
	 * containing the information about the SNPs of that population
	 * @param args list of parameters stored in a vector of Strings 
	 * 			-	args[0] whole path of the directory containing the input .txt file
	 * 			- 	args[1] name of the input .txt file
	 * 			-	args[2] whole path of the directory where the output files will be stored
	 * Additional optional parameters: - STR [num] [s] [muSTRs]
	 *			-	args[3] = STR : to require output files with information about STRs mutations for each contemporary population
	 * 			-	args[4] = num : integer representing the number of STRs loci
	 *			-	args[5] = s : inital state for each STR locus
	 *		    -	args[6] = muSTRs : mutation rate for STRs loci
	 */
	public static void main(String[] args) {
		
		populations =  new HashMap<Integer, PopulationARG>();;
		allMutations = new TreeSet<Double>();
		ordered_id_populations = new ArrayList<Integer>();
		events = new ArrayList<Event>();
		allmutationSet = new HashMap<Double, Mutation>();
		splitPoints = new TreeSet<Double>();
		num_edges_more_than_one_mut = 0;
		//total number of edges that have no mutations
		num_edges_without_mut = 0;
		
		
		Runtime runtime = Runtime.getRuntime();
	    long start = System.currentTimeMillis();
		
		//1A. READ FILE IN INPUT JAR VERSION
		String argPathInput;
		String argFileName_Input;
		String argPathOutput;
		String argFileName_Output;
		
		if(args.length < 3) {
			errorInsertingParameters();
		}
		
		argPathInput = args[0].toString();
		argFileName_Input = args[1].toString();
		argPathOutput = args[2].toString();
		argFileName_Output = args[3].toString();
		
		read_input_from_file(""+argPathInput+argFileName_Input);
		
		//1B. READ FILE IN INPUT (version not for jar file)
		//read_input_from_file("scaffold1.txt");
		//read_input_from_file("scaffold2.txt");
		//read_input_from_file("single_population.txt");
	   // read_input_from_file("scaffold_3admix.txt");
		
		
		
		//2. GENERATION OF POPULATION EVENTS AND STRUCTURE
		generatePopulations_Events();
		
		//3. DECORATION WITH SNPs 
		decorateSNPs(""+argPathOutput,argFileName_Output);
		//decorateSNPs("","scaffold_3admix");
		
		//4. DECORATION WITH STRs [optional]
		if(args.length > 4) {
			
			//Check if the 6th parameter is the string "-STR"
			if(args[4].toString().equalsIgnoreCase("-STR")) {
				//Check if the the number of parameter after "STR" option is correct
				if(args.length == 8) {
					int numSTRs = Integer.parseInt(args[5]);
					int initialeState = Integer.parseInt(args[6]);
					double mutRateSTRs = Double.parseDouble(args[7]);
					decorateSTRs(""+argPathOutput,argFileName_Output,numSTRs,initialeState,mutRateSTRs);
					
				}
				else
				    errorInsertingParameters();
			}
			else
				errorInsertingParameters();	
		}	
		
		//without jar file
		//decorateSTRs("", "scaffold_3admix", 15,40,6.0);
		
		//5. CREATION OF L AND S FILES
		CreatingFilesForStructure.setGlobalID_nodes(populations);
		
		//CreatingFilesForStructure.printMapNodesTimes();
		CreatingFilesForStructure.createLtxt(""+argPathOutput,argFileName_Output,events);
		CreatingFilesForStructure.createStxt(""+argPathOutput,argFileName_Output,GenerateAdmixturePopulations.getOrdered_id_populations(), GenerateAdmixturePopulations.getPopulations());
		CreatingFilesForStructure.CreatefileForCytoscape(""+argPathOutput,argFileName_Output);
		
		//without jar file
		//CreatingFilesForStructure.createLtxt("","scaffold_3admix",events);
		//CreatingFilesForStructure.createStxt("","scaffold_3admix",GenerateAdmixturePopulations.getOrdered_id_populations(), GenerateAdmixturePopulations.getPopulations());
		//CreatingFilesForStructure.CreatefileForCytoscape("","scaffold_3admix");
		
		
		//for stats file
		//computeRecombANDCoalescentNumber();
		//CreatingFilesForStructure.createSTATSfile("");
		 
		System.out.println("\nMemory - Time\n");
		long end = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		System.out.print(formatter.format((end - start) / 1000d));
		int mb = 1024;
		runtime = Runtime.getRuntime();
		System.out.println("\t"+(runtime.totalMemory() - runtime.freeMemory()) / mb);
		
	}
	
/**
 * This procedure decorates the whole scaffold with SNPs mutations and creates the output files describing the SNPs for each actual population
 * @param argPathOutput string that specifies the directory where the output files will be stored
 * @param argFileName_Output name of the output 
 */
public static void decorateSNPs(String argPathOutput, String argFileName_Output){
	
	//Decoration of all populations 
	Iterator<Integer> it_pops = populations.keySet().iterator();
	while(it_pops.hasNext()){
		Integer id_pop = it_pops.next();
		PopulationARG pop = populations.get(id_pop);
		try{
		DecorateMutationsSNP.decoratingWithMutations(pop.getMu(), pop);
		}
		catch(ErrorTimesInputFileException e){
			System.out.println("\n\n !!!!!!\n\n Attention: error in input file! Check time parameters! \n Inconsistency in event times specification!\n\n !!!!! ");
			System.exit(1);
		}
	}
			
	createTxtFile_DescriptionScaffold(""+argPathOutput+""+argFileName_Output+"_Scaffold_Description.txt");
	
	
	//Update the SNPs starting from the last population generated
	//System.out.println("Total Number of ordered populations: "+ordered_id_populations.size());
	//System.out.println(" Number of leaf populations: "+num_actual_pops);
	for(int j = ordered_id_populations.size()-1; j > num_actual_pops-1; j--){
		int id_pop = ordered_id_populations.get(j);
		//System.out.println("Index j = "+j);
		//System.out.println("Id population = "+id_pop);
		DecorateMutationsSNP.updateMutationsForTheBottomPop(populations.get(id_pop));
	}
	
	//Create SNPs output files 
	for(int j = 0; j < num_actual_pops; j++){
		int id_pop = ordered_id_populations.get(j);
		System.out.println(id_pop);
		DecorateMutationsSNP.createSNP_TxtFile(""+argPathOutput+""+argFileName_Output+"_SNP_pop"+id_pop+".txt", populations.get(id_pop));		
	}
}

/**
 * The procedure decorates the all scaffold with STR mutations. For each leaf (contemporary) population it creates a txt file containing the information  about STRs of the leaves 
 * @param argPathOutput whole path of the directory where the output file will be stored
 * @param argFileName_Output file name of the output file
 * @param numSTRs total number of STRs
 * @param initialState initial state value for each STR
 * @param mutRateSTRs mutation rate relative to STRs
 */
public static void decorateSTRs(String argPathOutput, String argFileName_Output, int numSTRs, int initialState, double mutRateSTRs){
		//Decoration of all populations 
		Iterator<Integer> it_pops = populations.keySet().iterator();
		while(it_pops.hasNext()){
			Integer id_pop = it_pops.next();
			PopulationARG pop = populations.get(id_pop);
			DecorateSTRs.inizializeStrs(numSTRs,initialState,mutRateSTRs*Math.pow(10, -3),pop);
			DecorateSTRs.decoratingDeltaStrsEdges(pop);
		}
		
		//Update STRs in order..the first population and then the other
		//System.out.println("Updating STRs");
		//System.out.println("Total Number of ordered populations: "+ordered_id_populations.size());
		//System.out.println(" Number of leaf populations: "+num_actual_pops);
		for(int j = ordered_id_populations.size()-1; j > num_actual_pops-1; j--){
			int id_pop = ordered_id_populations.get(j);
			//System.out.println("Index j = "+j);
			//System.out.println("Id population = "+id_pop);
			///***updateSTRs
			DecorateSTRs.updateSTRsRoots(populations.get(id_pop));
		}
		
		//Create STRs for leaves population and print them in output files 
		for(int j = 0; j < num_actual_pops; j++){
			int id_pop = ordered_id_populations.get(j);
			//System.out.println(id_pop);
			DecorateSTRs.updateStrs(populations.get(id_pop));
			//System.out.println("Creating STRs txt file");
			DecorateSTRs.createSTR_TxtFile(""+argPathOutput+""+argFileName_Output+"_STRs_pop"+id_pop+".txt", populations.get(id_pop));
		}
}

/**
 * This procedure creates the structure of the scaffold that is how the ancestral populations are connected among them and with the leaves
 * For each population (ancestral or leaf) it creates a file containing information about the structure of the ARG that describes the evolution of the population
 */
public static void generatePopulations_Events(){
	
	for(int i = 0; i < events.size(); i++){
		Event e = events.get(i);
		System.out.println("----- Event "+e.getId_event()+" -----");
		System.out.println("Kind : "+e.getKind());
		
		// *** if it's a leaf of the scaffold
		if(e.getKind() == 0) {
			System.out.println("Start time: "+e.getTime_event());
			System.out.println("Creating leaf population "+e.getPops_input()[0]+" with end time "+e.getEnd_times()[0]);
			System.out.println("With N = "+e.getN()[0]+" and m = "+e.getM());
			///Create the leaf population
			//start_time = 0
			int id_leaf = e.getPops_input()[0];
			PopulationARG pop_leaf = new PopulationARG();
			pop_leaf.generateLeafPopulation(id_leaf, e.getEnd_times()[0], e.getN()[0], e.getM(), r, mu, g);
			//CreatingFilesForStructure.createTxtFile_PopulationStructure(""+argPathOutput+"pop"+id_leaf+".txt", pop_leaf);
			populations.put(id_leaf, pop_leaf);
			ordered_id_populations.add(id_leaf);
		}
		
		// ** if it's a merge event
		if(e.getKind() == 1){
			System.out.println("Start time: "+e.getTime_event());
			System.out.println("Merging populations "+e.getPops_input()[0]+" and "+e.getPops_input()[1]);
			System.out.println("Merged population "+e.getPops_output()[0]+" with end time "+e.getEnd_times()[0]+" and N = "+e.getN()[0]);
			int id_popA = e.getPops_input()[0];
			int id_popB = e.getPops_input()[1];
			int id_popM = e.getPops_output()[0];
			PopulationARG pop_merge = new PopulationARG();
			pop_merge.generateCoalescentPopulation(id_popM, populations.get(id_popA), populations.get(id_popB), e.getTime_event(), e.getEnd_times()[0], e.getN()[0], r, mu, g);
			//CreatingFilesForStructure.createTxtFile_PopulationStructure(""+argPathOutput+"pop"+id_popM+".txt", pop_merge);
			populations.put(id_popM, pop_merge);
			ordered_id_populations.add(id_popM);
		}
		
		// ** if it's a splitting event
		if(e.getKind() == 2){
			System.out.println("Start time: "+e.getTime_event());
			System.out.println("Splitting population "+e.getPops_input()[0]);
			System.out.println("Resulting population "+e.getPops_output()[0]+" with end time "+e.getEnd_times()[0]+" and N = "+e.getN()[0]);
			System.out.println("Resulting population "+e.getPops_output()[1]+" with end time "+e.getEnd_times()[1]+" and N = "+e.getN()[1]);
			int id_pop = e.getPops_input()[0];
			int id_popA = e.getPops_output()[0];
			int id_popB = e.getPops_output()[1];
			int N_popA = e.getN()[0];
			int N_popB = e.getN()[1];
			ArrayList<ArrayList<Integer>> left_right_activeEdges = populations.get(id_pop).splitPopulation();
			if(left_right_activeEdges.size()==2){
				System.out.println("Splitting population "+id_pop);
				
				PopulationARG popA = new PopulationARG();
				popA.generatePopulationFromSplitting(id_popA, populations.get(id_pop), left_right_activeEdges.get(0), e.getTime_event(), e.getEnd_times()[0], N_popA, r, mu, g);
				//CreatingFilesForStructure.createTxtFile_PopulationStructure(""+argPathOutput+"pop"+id_popA+".txt", popA);
				
				PopulationARG popB = new PopulationARG();
				popB.generatePopulationFromSplitting(id_popB, populations.get(id_pop), left_right_activeEdges.get(1), e.getTime_event(), e.getEnd_times()[1], N_popB, r, mu, g);
				//CreatingFilesForStructure.createTxtFile_PopulationStructure(""+argPathOutput+"pop"+id_popB+".txt", popB);
				
				System.out.println("End Splitting");
				
				populations.put(id_popA, popA);
				populations.put(id_popB, popB);
				ordered_id_populations.add(id_popA);
				ordered_id_populations.add(id_popB);
			}
			else 
				System.out.println("The population "+id_pop+"cannot be splitted");
				//add exception
		}
	}
}

/**
 * This procedure reads the input .txt file and stores the  read information that will be used for the generation of the scaffold
 * @param input_file string that represents the name of the input .txt file
 */
public static void read_input_from_file(String input_file){
		
	int event_count = 0;
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(input_file));
			String line;
			
			//reading the parameters
			line = br.readLine();
			//System.out.println(line);
			line = line+"\n";
		
			//Read segment length g
			int c = 2;
			int start_g = c; 
			while(line.charAt(c) != ' '){
				c++;
			}	
			String g_s = line.substring(start_g, c);
			g = (int)Integer.parseInt(g_s);
			//System.out.println(g);
			
			//Read recombination rate r
			c = c+3; //' 'r=
			int start_r = c;
			while(line.charAt(c) != ' '){
				c++;
			}	
			String r_s = line.substring(start_r, c);
			r = (double)Double.parseDouble(r_s);
			//System.out.println(r);
			
			//Read mutation rate mu
			c = c+4; //' 'mu=
			int start_mu = c;
			while(line.charAt(c) != '\n'){
				c++;
			}	
			String mu_s = line.substring(start_mu, c);
			mu = (double)Double.parseDouble(mu_s);
			//System.out.println(mu);
			
			// *** Read the information about actual populations ***
			line = br.readLine();
			//System.out.println(line);
			
			while(!line.equals("#begin_actual_populations")){
				//System.out.println(line);
				line = br.readLine();	
			}
			
			line = br.readLine();
			while(!line.equals("#end_actual_populations")) {
				
				//System.out.println(line);
				c = 0;
				while(line.charAt(c)!='='){
					c++;
				}
				
				c++;
				int start_id_pop = c;
				while(line.charAt(c)!=' '){
					c++;		
				}
				String id_pop_s;
				
				id_pop_s = line.substring(start_id_pop, c);
				int id_pop = (int)Integer.parseInt(id_pop_s);		
				//System.out.println(id_pop);
				
				//Reading m
				while(line.charAt(c)!='='){
					c++;
				}
				c++;
				int start_m = c;
				while(line.charAt(c)!=' '){
					c++;		
				}
				String m_pop_s;
				m_pop_s = line.substring(start_m,c);
				//System.out.println(m_pop_s);
				int m_pop = (int)Integer.parseInt(m_pop_s);
				
				//Reading N
				while(line.charAt(c)!='='){
					c++;
				}
				c++;
				int start_N = c;
				while(line.charAt(c)!=' '){
					c++;		
				}
				String N_pop_s = line.substring(start_N,c);
				//System.out.println(N_pop_s);
				int N_pop = (int)Integer.parseInt(N_pop_s);
				
				//Reading end_time
				while(line.charAt(c)!='='){
					c++;
				}
				c++;
				int start_time = c;
				while(line.charAt(c)!=']'){
					c++;
				}
				String end_time_s = line.substring(start_time, c);
				double end_time = (double)Double.parseDouble(end_time_s);
				//System.out.println(end_time);
				
				//Creat leaf event
				Event e = new Event();
				e.setId_event(event_count);
				event_count++;
				e.setTime_event(0.0);
				e.setKind(0);
				e.setM(m_pop);
				e.getN()[0]=N_pop;
				e.getPops_input()[0] = id_pop;
				e.getEnd_times()[0] = end_time;
				
				events.add(e);
				
				num_actual_pops++;
				num_tot_pops++;
				
				line = br.readLine();	
			}
			
			//System.out.println("Number of samples: "+num_actual_pops);
			
			line = br.readLine();
			
			while(!line.equals("#begin_events")){
				//System.out.println(line);
				line = br.readLine();
			}
			line = br.readLine();
			while(!line.equals("#end_events")) {
				
				//System.out.println(line);
				
				//time of the event
				c = 0;
				int start_time_fc = c;
				while(line.charAt(c)!=' '){
					c++;
				}
				String time_event_s = line.substring(start_time_fc, c);
				double time_event = (double)Double.parseDouble(time_event_s);
				//System.out.println(time_event);
				
				c++;
				int kind_start = c;
				while(line.charAt(c)!='('){
					c++;
				}
				
				//event kind: split or merge
				String kind = line.substring(kind_start, c);
				
				// *** MERGE ***
				if(kind.equals("merge")) {
					num_tot_pops++;
					//System.out.println("c="+c);
					c=c+1;
					int id_popA_start  = c;
					while(line.charAt(c)!=','){
						c++;
					}
					String id_popA_s = line.substring(id_popA_start, c);
					//System.out.println(id_pop_s);
					int id_popA = (int)Integer.parseInt(id_popA_s);
					System.out.println(id_popA);
					
					c=c+1;
					int id_popB_start  = c;
					while(line.charAt(c)!=')'){
						c++;
					}
					String id_popB_s = line.substring(id_popB_start, c);
					//System.out.println(id_pop_s);
					int id_popB = (int)Integer.parseInt(id_popB_s);
					//System.out.println(id_popB);
					
					//Reading pop merge
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int id_popM_start  = c;
					while(line.charAt(c)!=' '){
						c++;
					}
					String id_popM_s = line.substring(id_popM_start, c);
					int id_popM = (int)Integer.parseInt(id_popM_s);
					//System.out.println(id_popM);
					
					//Reading N of pop_merge
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int N_popM_start = c;
					while(line.charAt(c)!=' '){
						c++;
					}
					String N_popM_s = line.substring(N_popM_start, c);
					//System.out.println(N_popM_s);
					int N_popM = (int)Integer.parseInt(N_popM_s);
					
					//Reading end_time
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int endtime_popM_start = c;
					while(line.charAt(c)!=']'){
						c++;
					}
					String endtime_popM_s = line.substring(endtime_popM_start, c);
					double endtime_popM = (double)Double.parseDouble(endtime_popM_s);
					//System.out.println(endtime_popM);
					
					Event e = new Event();
					e.setId_event(event_count);
					event_count++;
					e.setTime_event(time_event);
					e.setKind(1);
					e.getN()[0] = N_popM;
					e.getPops_input()[0] = id_popA;
					e.getPops_input()[1] = id_popB;
					e.getPops_output()[0] = id_popM;
					e.getEnd_times()[0] = endtime_popM;
					events.add(e);
				}
				
				// *** SPLIT ***
				if(kind.equals("split")) {
					num_tot_pops=num_tot_pops+2;
					c=c+1;
					int id_pop_start  = c;
					while(line.charAt(c)!=')'){
						c++;
					}
					String id_pop_s = line.substring(id_pop_start, c);
					int id_pop = (int)Integer.parseInt(id_pop_s);
					//System.out.println(id_pop);
					
					//Reading id_popA
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int id_popA_start  = c;
					while(line.charAt(c)!=' '){
						c++;
					}
					
					String id_popA_s = line.substring(id_popA_start, c);
					int id_popA = (int)Integer.parseInt(id_popA_s);
					//System.out.println(id_popA);
					
					//Reading N popA
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int N_popA_start  = c;
					while(line.charAt(c)!=' '){
						c++;
					}
					String N_popA_s = line.substring(N_popA_start, c);
					int N_popA = (int)Integer.parseInt(N_popA_s);
					//System.out.println(N_popA);
					
					//Reading end_time popA
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int endtime_popA_start  = c;
					while(line.charAt(c)!=']'){
						c++;
					}
					String endtime_popA_s = line.substring(endtime_popA_start, c);
					double endtime_popA = (double)Double.parseDouble(endtime_popA_s);
					//System.out.println(endtime_popA);
					
					
					//Reading pop B
					//Reading id_popA
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int id_popB_start  = c;
					while(line.charAt(c)!=' '){
						c++;
					}
					
					String id_popB_s = line.substring(id_popB_start, c);
					int id_popB = (int)Integer.parseInt(id_popB_s);
					//System.out.println(id_popB);
					
					//Reading N popA
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int N_popB_start  = c;
					while(line.charAt(c)!=' '){
						c++;
					}
					String N_popB_s = line.substring(N_popB_start, c);
					int N_popB = (int)Integer.parseInt(N_popB_s);
					//System.out.println(N_popB);
					
					//Reading end_time popA
					while(line.charAt(c)!='='){
						c++;
					}
					c++;
					int endtime_popB_start  = c;
					while(line.charAt(c)!=']'){
						c++;
					}
					String endtime_popB_s = line.substring(endtime_popB_start, c);
					double endtime_popB = (double)Double.parseDouble(endtime_popB_s);
					//System.out.println(endtime_popB);
					
					
					Event e = new Event();
					e.setId_event(event_count);
					event_count++;
					e.setKind(2);
					e.setTime_event(time_event);
					e.getN()[0]=N_popA;
					e.getN()[1]=N_popB;
					e.getPops_input()[0]=id_pop;
					e.getPops_output()[0]=id_popA;
					e.getPops_output()[1]=id_popB;
					e.getEnd_times()[0]=endtime_popA;
					e.getEnd_times()[1]=endtime_popB;
					events.add(e);
				}
				
				line = br.readLine();
			}
			br.close();
		}//end try
		catch (IOException ex) {
			ex.printStackTrace();
		} 	
		catch(Exception e){
			e.printStackTrace();
		} //or write your own exceptions
	}

//***** UTILITIES *****

/**
 * This procedure prints to video an ordered list of ARGs as they have been created
 */
public static void printIdsPopulations(){
	System.out.println("List IDs of Populations");
	for(int i = 0; i < ordered_id_populations.size(); i++){
		System.out.println(ordered_id_populations.get(i));
	}
	
	System.out.println("Populations that are not leaves");
	for(int j = ordered_id_populations.size()-1; j > num_actual_pops-1; j--){
		System.out.println(ordered_id_populations.get(j));
	}
	
	System.out.println("Leaves of the scaffold");
	//Create SNPs output files 
	for(int j = 0; j < num_actual_pops; j++){
		System.out.println(ordered_id_populations.get(j));
	}
	
}

/**
 * This procedure prints to video an ordered list of events that are merge and/or split of ARGs
 */
public static void printEvents(){
	for(int i = 0; i < events.size(); i++){
		Event e = events.get(i);
		System.out.println("----- Event "+e.getId_event()+" -----");
		
		System.out.println("Kind : "+e.getKind());
		if(e.getKind() == 0) {
			System.out.println("Start time: "+e.getTime_event());
			System.out.println("Creating leaf population "+e.getPops_input()[0]+" with end time "+e.getEnd_times()[0]);
		}
		if(e.getKind() == 1){
			System.out.println("Start time: "+e.getTime_event());
			System.out.println("Merging populations "+e.getPops_input()[0]+" and "+e.getPops_input()[1]);
			System.out.println("Merged population "+e.getPops_output()[0]+" with end time "+e.getEnd_times()[0]);
		}
		if(e.getKind() == 2){
			System.out.println("Start time: "+e.getTime_event());
			System.out.println("Splitting population "+e.getPops_input()[0]);
			System.out.println("Resulting population "+e.getPops_output()[0]+" with end time "+e.getEnd_times()[0]);
			System.out.println("Resulting population "+e.getPops_output()[1]+" with end time "+e.getEnd_times()[1]);
		}
	}
	
}
/**
 * This procedure creates a .txt file describing the structure of the scaffold that has been created. In other words this file describes has the different ARGs are connected to each other
 * @param wholePath String representing the wholePath where the output file will be stored (the name of the file is in the wholePath)
 */
public static void createTxtFile_DescriptionScaffold(String wholePath) {
	Writer writer = null;
	File file = new File(""+wholePath);
	
	try {
	    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
	    
	    Iterator<Integer> it_pops = populations.keySet().iterator();
	    
	    writer.write("------------		SCAFFOLD DESCRIPTION		------------ \n\n");
	    writer.write("FIXED PARAMETERS \n");
	    writer.write("	- SNP mutation rate  =  "+mu+"\n");
	    writer.write("	- recombination rate =  "+r+"\n");
	    writer.write("	- segment length   =  "+g+"\n\n");
	    
	    writer.write("Number of edges with more than one SNP mutation =  "+num_edges_more_than_one_mut+"\n");
	    writer.write("Number of edges with no SNP mutation =  "+num_edges_without_mut+"\n\n");
	    
	    while(it_pops.hasNext()){
	    	
	    	PopulationARG pop = populations.get(it_pops.next());
	    	
		    writer.write("-------- POPULATION " +pop.getId_pop()+" --------\n");
    
		    if(pop.getKind() == 0){
		    	writer.write("Population "+pop.getId_pop()+" is a LEAF of the scaffold \n");
		    }
		    
		    if(pop.getKind() == 2){
		    	writer.write("Population "+pop.getId_pop()+"  derives by the MERGING of 2 populations: "+pop.getPopanc()[0]+" and "+pop.getPopanc()[1] +"at time "+populations.get(pop.getPopanc()[0]).getThresholdT()+"\n");
		    	
		    }
		    if(pop.getKind() == 1){
		    	writer.write("Population "+pop.getId_pop()+" derives by the SPLITTING of population: "+pop.getPopanc()[0]+"at time "+populations.get(pop.getPopanc()[0]).getThresholdT()+"\n");
		    }
		    
		    writer.write("Number of extant units : " +pop.getExtantUnits()+"\n");
		    writer.write("Effective population size : " +pop.getN()+"\n");
		    writer.write("Total number of nodes : " +pop.getNodeSet().size()+"\n");
		    writer.write("Total number of edges : " +pop.getGraphEdges().size()+"\n");
		    writer.write("Time threshold : " +pop.getThresholdT()+"\n");
		    writer.write("Number of active lineages at time "+pop.getThresholdT()+" : "+pop.getActiveEdges().size()+"\n");
		    writer.write("Total number of SNP mutations : " +pop.getMutationSet().size()+"\n");
		    writer.write("----------------------------------\n\n");
		
	    }
	   writer.close();   
	}
	catch (IOException ex) {
		ex.printStackTrace();
	} 	
}

/**
 * In case of mistake by the user, this procedure display instruction about the correct order and type of parameters to execute SimRA from command line. 
 */
public static void errorInsertingParameters(){
	System.out.println("Usage:\n\t java -jar SimRa.jar [pathInputDirectory] [FileNameInput] [pathOutputDirectory]");
	System.out.println("\t\t [pathInputDirectory] : whole path of the directory where the input txt file is located");
	System.out.println("\t\t [fileName] : name of the input text file");
	System.out.println("\t\t [pathOutputDirectory] : whole path of the directory where the output files will be stored");
	System.out.print("Additional optional parameters:\n\t");
	System.out.println("-STR [num] [s] [muSTRs]");
	System.out.println("\t\t -STR : to require output files with information about STRs mutations for each contemporary population;");
	System.out.println("\t\t num : integer representing the number of STRs loci;");
	System.out.println("\t\t s : inital state for each STR locus;");
	System.out.println("\t\t muSTRs : mutation rate for STRs loci;");
	System.out.println("\n EXAMPLE 1: \n\t java -jar SimRa.jar ~/input_SimRA/ scaffold1.txt ~/outputSimRA/");
	System.out.println("\n EXAMPLE 2: \n\t java -jar SimRa.jar ~/input_SimRA/ scaffold2.txt ~/outputSimRA/ -STR 15 40 6.0");
	System.exit(1);
}

/**
 * This procedure compute the total number of recombination and coalescence events in the whole scaffold
 */
public static void computeRecombANDCoalescentNumber(){
	
	setRecomb_total_number(0);
	setCoalesc_total_number(0);
	
	for(int j = 0; j < GenerateAdmixturePopulations.getOrdered_id_populations().size(); j++){
		int id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
		PopulationARG pop = GenerateAdmixturePopulations.getPopulations().get(id_pop);
		GenerateAdmixturePopulations.setRecomb_total_number(GenerateAdmixturePopulations.getRecomb_total_number()+pop.getRecombinationsNumber());	
		GenerateAdmixturePopulations.setCoalesc_total_number(GenerateAdmixturePopulations.getCoalesc_total_number()+pop.getCoalescentNumber());	
		}	
}


}
