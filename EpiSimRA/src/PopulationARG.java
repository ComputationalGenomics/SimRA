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
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;


/**
 * The class PopulationARG represents the evolution of a single population or edge in the scaffold. The evolution of a population is modeled by ARG.
 * This class contains all the information and functions to manage the ARG.
 * 
 * @author Anna Paola Carrieri
 */

public class PopulationARG {
	
	//***************************** PRIVATE FIELDS  ***************************************
	
	//Id of population
	private int id_pop;
	//List of the current active lineages 
	private  ArrayList<Integer> activeEdges;
	//Set E of actual edges of the ARG G=(V,E)
	private  Map<Integer,Edge> graphEdges;
	//Set V of actual nodes of the ARG G=(V,E)
	private  Map<Integer,Node> nodeSet;
	
	//kind = 0 if the population is a leaf of the scaffold,
	//kind = 1 if population derives from 1 population (it's from recombination), 
	//kind = 2 if population derives from 2 populations (it's coalescence event)
	//kind = 3 if a root of the scaffold
	private int kind = 0;
	
	//array that contains the id of population/s 
	private int popanc[] = new int[2];
	
	//Map between the id of the root nodes and the id of the corresponding leave node of the ancestral
	//private Map<Integer, Integer> leaves_roots;
	//ArrayList with at most two {id_pop ; Map<Integer, Integer>}
	private HashMap<Integer, HashMap<Integer, Integer>> maps_leaves_roots;

	//Map of SNP Mutations in the ARG
	private  Map<Double,Mutation> mutationSet;
	//Set of SplitPoints due to recombination events
	private  TreeSet<Double> splitPoints;
	//Ordered list of SNP positions
	private  TreeSet<Double> SNPpositionsList;
	//m = number of extant units
	private  int extantUnits=0;
	//L = number of live lineages L inizialized to m
	private  int L=0;
	//g = length of the chromosomal segment
	private  double g=0;
	//N = population size
	private  int N=0;
	//current generation ot time
	private  double curr_generation=0;
	//private static double scaled generation j = TN
	private  double scaledGeneration=0;
	//recombination rate recomb=(g/c)
	private  double recomb = 0;
	//stimated number of mutations
	private  int mutationsNumber = 0;
	//mutation rate
	private  double mu = 0;
	//Name of output file
	private  String filename ="";
	//Total number of recombination events
	private  int recombinationsNumber = 0;
	//Total number of coalescent events
	private  int coalescentNumber = 0;
	//time t in which the population stops growing
	private double thresholdT = 0;
	
	/**
	 * default constructor
	 */
	PopulationARG(){ 
	}
	
	//***************************** GETTERS AND SETTERS ***************************************
	
	/**
	 * This function returns the map that contains the IDs of mutations and the Mutation objects
	 * @return the map containing IDs of mutations and the Mutation objects 
	 */
	public  Map<Double, Mutation> getMutationSet() {
		return mutationSet;
	}
	/**
	 * The function returns the list of IDs of active lineages in the ARG
	 * @return the list of IDs of active lineages in the ARG
	 */
	public ArrayList<Integer> getActiveEdges() {
		return activeEdges;
	}
	/**
	 * The procedure sets the list of IDs of active lineages in the ARG
	 * @param activeEdges list of IDs of active lineages in the ARG
	 */
	public void setActiveEdges(ArrayList<Integer> activeEdges) {
		this.activeEdges = activeEdges;
	}
	
	/**
	 * The function returns the map containing the IDs of the edges and the Edge objects of the ARG
	 * @return the map containing the IDs of the edges (keys) and the Edge objects of the ARG
	 */
	public Map<Integer, Edge> getGraphEdges() {
		return graphEdges;
	}
	/**
	 * The procedure sets the map containing the IDs of the edges (keys) and the Edge objects of the ARG
	 * @param graphEdges Map object containing the IDs of the edges (keys) and the Edge objects of the ARG
	 */
	public void setGraphEdges(Map<Integer, Edge> graphEdges) {
		this.graphEdges = graphEdges;
	}
	/**
	 * The function returns the map containing the IDs of the nodes and the Node objects of the ARG
	 * @return the map containing the IDs of the nodes (keys) and the Node objects of the ARG
	 */
	public Map<Integer, Node> getNodeSet() {
		return nodeSet;
	}
	/**
	 * The procedure sets the map containing the IDs of the nodes (keys) and the Node objects of the ARG
	 * @param nodeSet Map object containing the IDs of the nodes (keys) and the Node objects of the ARG
	 */
	public void setNodeSet(Map<Integer, Node> nodeSet) {
		this.nodeSet = nodeSet;
	}
	/**
	 * The function returns the the order list of split points in the whole ARG that have been generated by recombination events
	 * @return the order list of split points in the whole ARG
	 */
	public TreeSet<Double> getSplitPoints() {
		return splitPoints;
	}
	/**
	 * The procedure sets the order list of split points in the whole ARG
	 * @param splitPoints order list of split points in the whole ARG
	 */
	public void setSplitPoints(TreeSet<Double> splitPoints) {
		this.splitPoints = splitPoints;
	}
	/**
	 * The function returns the sample size
	 * @return sample size (number of leaves in the ARG)
	 */
	public int getExtantUnits() {
		return extantUnits;
	}
	/**
	 * The function sets the size of the sample (number of leaves in the ARG)
	 * @param extantUnits integer representing the sample size 
	 */
	public void setExtantUnits(int extantUnits) {
		this.extantUnits = extantUnits;
	}
	/**
	 * The function returns the number of active lineages in the current ARG
	 * @return the number of active lineages
	 */
	public int getL() {
		return L;
	}
	/**
	 * The procedure sets the number of active lineages in the current ARG
	 * @param l the  number of active lineages
	 */
	public void setL(int l) {
		L = l;
	}
	/**
	 * This function returns the length of the chromosome 
	 * @return the length of the chromosome (of segment of it)
	 */
	public double getG() {
		return g;
	}
	/**
	 * Sets the length of the chromosome (of segment of it)
	 * @param g segment length 
	 */
	public void setG(double g) {
		this.g = g;
	}
	/**
	 * This function return the effective population size
	 * @return the effective population size
	 */
	public int getN() {
		return N;
	}
	/**
	 * Sets the effective population size
	 * @param n effective population size
	 */
	public void setN(int n) {
		N = n;
	}
	/**
	 * This function returns the current age (time) of the ARG (not scaled in generations)
	 * @return the current age (time) of the ARG (not scaled in generation)
	 */
	public double getCurr_generation() {
		return curr_generation;
	}
	/**
	 * This procedure sets the current age (time) of the ARG (not scaled in generations)
	 * @param curr_generation the current age (time) of the ARG (not scaled in generations)
	 */
	public void setCurr_generation(double curr_generation) {
		this.curr_generation = curr_generation;
	}
	/**
	 * This function returns the scaled depth (age) of the ARG in generations
	 * @return the scaled depth (age) of the ARG in generations
	 */
	public double getScaledGeneration() {
		return scaledGeneration;
	}
	/**
	 * sets the scaled depth (age) of the ARG in generations
	 * @param scaledGeneration scaled depth (age) of the ARG in generations
	 */
	public void setScaledGeneration(double scaledGeneration) {
		this.scaledGeneration = scaledGeneration;
	}
	/**
	 * returns the recombination rate
	 * @return the recombination rate
	 */
	public double getRecomb() {
		return recomb;
	}
	/**
	 * sets the recombination rate 
	 * @param recomb recombination rate 
	 */
	public void setRecomb(double recomb) {
		this.recomb = recomb;
	}
	/**
	 * returns the total number of mutations in the whole ARG
	 * @return the total number of mutations in the whole ARG
	 */
	public int getMutationsNumber() {
		return mutationsNumber;
	}
	/**
	 * Sets the total number of SNP mutations in the whole ARG
	 * @param mutationsNumber total number of SNP mutations in the whole ARG
	 */
	public void setMutationsNumber(int mutationsNumber) {
		this.mutationsNumber = mutationsNumber;
	}
	/**
	 * this function retursn the SNPs mutation rate
	 * @return the SNP mutation rate 
	 */
	public double getMu() {
		return mu;
	}
	/**
	 * Sets the SNP mutation rate
	 * @param mu SNP mutation rate
	 */
	public void setMu(double mu) {
		this.mu = mu;
	}
	/**
	 * returns the name of the file
	 * @return the name of the file
	 */
	public String getFilename() {
		return filename;
	}
	/**
	 * Sets the name of the file
	 * @param filename string representing the name of the file
	 */
	public void setFilename(String filename) {
		this.filename = filename;
	}
	/**
	 * returns  the total number of recombinations in the ARG
	 * @return the total number of recombinations in the ARG
	 */
	public  int getRecombinationsNumber() {
		return recombinationsNumber;
	}
	/**
	 * Sets the total recombinations number in the ARG
	 * @param recombinationsNumber  total number of recombination nodes in the ARG
	 */
	public  void setRecombinationsNumber(int recombinationsNumber) {
		this.recombinationsNumber = recombinationsNumber;
	}
	/**
	 * 
	 * @return the total number of coalescence nodes in the whole ARG
	 */
	public  int getCoalescentNumber() {
		return coalescentNumber;
	}
	/**
	 * Set the total number of coalescent nodes in the whole ARG
	 * @param coalescentNumber the total number of coalescent nodes
	 */
	public  void setCoalescentNumber(int coalescentNumber) {
		this.coalescentNumber = coalescentNumber;
	}
	/**
	 *  The function returns the threshold of time for the backward simulation of the ARG
	 * @return threshold of time for the backward simulation of the ARG
	 */
	public double getThresholdT() {
		return thresholdT;
	}
	/**
	 * The procedure sets the threshold of time for the backward simulation of the
	 * @param thresholdT threshold of time for the backward simulation of the ARG
	 */
	public void setThresholdT(double thresholdT) {
		this.thresholdT = thresholdT;
	}
	/**
	 * The function returns the ordered set of SNP mutation positions in the chromosome
	 * @return the ordered set of SNP mutation positions in the chromosome
	 */
	public  TreeSet<Double> getSNPpositionsList() {
		return SNPpositionsList;
	}
	/**
	 * The procedure sets the ordered set of SNP mutation positions in the chromosome
	 * @param sNPpositionsList ordered set of SNP mutation positions in the chromosome
	 */
	public  void setSNPpositionsList(TreeSet<Double> sNPpositionsList) {
		SNPpositionsList = sNPpositionsList;
	}
	/**
	 * The procedure sets the map of SNP mutations in the ARG
	 * @param mutationSet Map containing all the instance of Mutation (SNPs) in the ARG. The key field is the SNP position in the segment (that is also its ID)
	 */
	public void setMutationSet(Map<Double, Mutation> mutationSet) {
		this.mutationSet = mutationSet;
	}
	/**
	 * The function returns the ID of the population (ARG)
	 * @return the ID of the population (ARG)
	 */
	public int getId_pop() {
		return id_pop;
	}
	/**
	 * The procedure sets the ID of the population (ARG)
	 * @param id_pop integer representing the ID of the population
	 */
	public void setId_pop(int id_pop) {
		this.id_pop = id_pop;
	}
	/**
	 * This function returns the kind of the population
	 * @return the kind of population (leaf or deriving from a merge or deriving from a splitting event)
	 * 	kind = 0 if the population is a leaf of the scaffold,
	 * kind = 1 if population derives from 1 population (it's from recombination), 
	 * kind = 2 if population derives from 2 populations (it's coalescence event)
	 * kind = 3 if a root of the scaffold
	 */
	public int getKind() {
		return kind;
	}
	/**
	 * The procedure sets the kind of population (leaf or deriving from a merge or deriving from a splitting event)
	 * kind = 0 if the population is a leaf of the scaffold,
	 * kind = 1 if population derives from 1 population (it's from recombination), 
	 * kind = 2 if population derives from 2 populations (it's coalescence event)
	 * kind = 3 if a root of the scaffold
	 * @param kind type of the population
	 */
	public void setKind(int kind) {
		this.kind = kind;
	}
	/**
	 * The function returns the ID/s of the ancestral population/s. The ancestral populations can be almost two.
	 * @return the ID/s of the ancestral population/s
	 */
	public int[] getPopanc() {
		return popanc;
	}
	/**
	 * The procedure sets the ID/s of the ancestral population/s
	 * @param popanc ID/s of the ancestral population/s. The integer vector can have at most size equal to two
	 */
	public void setPopanc(int[] popanc) {
		this.popanc = popanc;
	}
	/**
	 * The function returns the map between the leaves of this ARG and the roots of the underlying ARG
	 * @return the map between the leaves of this ARG and the roots of the underlying ARG.
	 * The key is the ID of the underlying population and the object is a map between roots of the underlying population (keys) and leaves
	 * of this PopulationARG object
	 */
	public HashMap<Integer, HashMap<Integer, Integer>> getMaps_leaves_roots() {
		return maps_leaves_roots;
	}
	/**
	 * The procedure sets the map between the leaves of the ARG and the roots of the underlying ARG
	 * The key is the ID of the underlying population and the object is a map between the roots of the underlying population (keys) and the leaves
	 * of this PopulationARG object
	 * @param maps_leaves_roots map between this ARG and the underlying ARG. More precisely the leaves of this ARG are mapped to the roots of the underlying ARG
	 */
	public void setMaps_leaves_roots(
			HashMap<Integer, HashMap<Integer, Integer>> maps_leaves_roots) {
		this.maps_leaves_roots = maps_leaves_roots;
	}


	
	//***************************** USEFUL PROCEDURES ***************************************
	/**
	 * Given the active lineages l1, l2, .., L the function return the sum of all the recombination rate
	 * rl (for l = 1, .., L) associated to each active lineage
	 * @return double number sum of recombination rates of the active lineages
	 */
	public  double computeSumRates(){
		double sum = 0;
		
		for(int i = 0; i < activeEdges.size(); i++) {
			sum = graphEdges.get(activeEdges.get(i)).getRate()+sum;
		}
		
		return sum; 
	}

	/**
	 * The function implements the exponential distribution with lambda as parameter
	 * @param lambda parameter of the exponential distribution
	 * @return the time t of the next event
	 */
	public  double computeNextTimeEvent(double lambda){
		double random = Math.random();
		
		//x = ln(1-p)/-lamda
		double t = (Math.log(1-random))/(0-lambda);
		
		return t;
	}
	
	
	/** 
	 * The function computes the time T = T+t to the next event using the exponential distribution
	 * @param time t computed by exponential distribution time is in [0,1)
	 * @return the time T updated based on to the current generation and the next time generated by the exponential distribution
	 */
	public  double computeNextGeneration(double time){
		return time+getCurr_generation();
	}
	
	/**
	 * Function that computes the type of the next event (coalescent if
	 * index is equal to 0, recombination if index is greater 0)
	 * @param lambda parameter of the exponential distribution that generates the next time of the event
	 * @param binCoef (L 2) where L is the # of the active lineages
	 * @return the index indicating the type of the next event
	 */
	public int computeIndexNextEvent(double lambda, double binCoef){
		
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
	public double binomialCoefficient(int n, int r) {
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
	 * The function generates the ARG structure representing the evolution of a single population
	 * @param threshold time in which the population either stops growing or splits in two population or merge with a second population
	 * @param populationSize actual population size
	 * @param m number of extant units
	 * @param r recombination rate 
	 * @param mu mutation rate
	 * @param g segment length of each leaf node
	 * @param id ID of the new leaf population 
	 */
	public void generateLeafPopulation(int id, double threshold, int populationSize, int m, double r, double mu, int g ) {
		

		//***************************** PARAMETER INITIALIZATION  ***************************************
		
		//***** Set id
		this.setId_pop(id);
		
		//***** Threshold time t *****
		setThresholdT(threshold);
		
		//***** Population Size *****
		setN(populationSize);
				
		//***** Extant Units ****
		setExtantUnits(m);
		
		//***** Recombination Rate ****
		setRecomb(r*Math.pow(10, -8));
		
		//**** SNP Mutation Rate ****
		setMu(mu*Math.pow(10, -8));
		
		//**** Segment length g ***
		setG(g*Math.pow(10, 3));
		
		//**** Set current generation ***
		setCurr_generation(0);
		
		setRecombinationsNumber(0);
		
		setCoalescentNumber(0);
		
				
		//***************************** ARG RANDOM GENERATION ALGORITHM   ***********************************
		
		//********* INITIALIZATION PRIVATE FIELD**************/
		
		nodeSet = new HashMap<Integer,Node>();
		graphEdges = new HashMap<Integer,Edge>();
		activeEdges = new ArrayList<Integer>();
		mutationSet = new HashMap<Double,Mutation>();
		SNPpositionsList = new TreeSet<Double>();
		splitPoints = new TreeSet<Double>();
		splitPoints.add(0.0);
		splitPoints.add(1.0);
		//it is a leaf population it could have ancestral population on it
		this.setKind(0);
	
		//********* INITIALIZATION OF LEAVE NODES ****************/
		
		for (int i=0; i<extantUnits; i++)
		{
			//Create a leave node
			Node n = new Node(i,false,getCurr_generation());
			ArrayList<Interval> segments = new ArrayList<Interval>();
			segments.add(new Interval(0,1));
			n.setSegments(segments);
			nodeSet.put(i, n);
			
			//Create a new active lineages that carries a segment of length 1
			Edge e = new Edge(i, n.getId(), -1, segments);
			e.computeLength();
			e.computeDensity();
			e.computeRate(getN(), getRecomb(), getG());
			
			//Add the new lineage to the list of active lineages
			activeEdges.add(e.getIdEDGE());
			graphEdges.put(e.getIdEDGE(), e);
							
		}
			
		setL(activeEdges.size());
		backwardConstructionARG(0,threshold,false);
	}
/**
 * This procedure simulate the merge of two population, say pop1 and pop2, merging in a new ancestral population (instance of the class PopulationARG)
 * @param id ID of the new ancestral population resulting from the merge of two populations
 * @param pop1 ID of the first population 
 * @param pop2 ID of the second population
 * @param startTime time in which the merge between pop1 and pop2 occurs creating the new ancestral population
 * @param endTime threshold of time in which the new population deriving from the merge will stop its growing
 * @param populationSize effective population size of the new population resulting from the merge
 * @param r recombination rate of the new population resulting from the merge
 * @param mu mutation rate of the new population resulting from the merge
 * @param g segment length of the sample in the new population
 */
public void generateCoalescentPopulation(int id, PopulationARG pop1, PopulationARG pop2, double startTime, double endTime, int populationSize, double r, double mu, int g) {
		

		//***************************** PARAMETER INITIALIZATION  ***************************************
		
		//***** Id population *****
		this.setId_pop(id);
		
		//***** Kind of Population *****
		this.setKind(2);
		
		//***** Threshold time t *****
		setThresholdT(endTime);
		
		//***** Population Size *****
		setN(populationSize);
				
		//***** Recombination Rate ****
		setRecomb(r*Math.pow(10, -8));
		
		//**** SNP Mutation Rate ****
		setMu(mu*Math.pow(10, -8));
		
		//**** Segment length g ***
		setG(g*Math.pow(10, 3));
		
		//*** Set current generation 
		this.setCurr_generation(startTime);
		
		nodeSet = new HashMap<Integer,Node>();
		graphEdges = new HashMap<Integer,Edge>();
		activeEdges = new ArrayList<Integer>();
		maps_leaves_roots = new HashMap<Integer, HashMap<Integer, Integer>>();
		mutationSet = new HashMap<Double,Mutation>();
		SNPpositionsList = new TreeSet<Double>();
		splitPoints = new TreeSet<Double>();
		splitPoints.add(0.0);
		splitPoints.add(1.0);
		popanc[0] = pop1.getId_pop();
		popanc[1] = pop2.getId_pop();
		
		setRecombinationsNumber(0);
		
		setCoalescentNumber(0);
		
		//********* INITIALIZATION OF LEAVE NODES and ACTIVE EDGES ************/
	
		//For each active edge of pop 1 take the node son of active edge and 
		//copy in this population the node and active edge and a mapping
		ArrayList<Integer> act_edges_pop1 = pop1.getActiveEdges();
		
		//**** COPY ACTIVE NODES IN POP1 AS LEAVES IN POP ****
		TreeSet<Integer> active_nodes = new TreeSet<Integer>();
		HashMap<Integer, Integer> map_nodes_pop1 = new HashMap<Integer, Integer>();
		int id_leaf = 0;
		Edge act_e;
		for(int i = 0; i < act_edges_pop1.size(); i++) {
			//take the node of pop1 that is son of the active edge and copy it has leaf
			act_e = pop1.getGraphEdges().get(act_edges_pop1.get(i));
			
			if(!active_nodes.contains(act_e.getId_son())){
				active_nodes.add(act_e.getId_son());
				Node n_pop1 = pop1.getNodeSet().get(act_e.getId_son());
				Node n = new Node();
				n.setId(id_leaf);
				n.setLevel(n_pop1.getLevel());
				n.setRecomb(n_pop1.isRecomb());
				n.setSegments(n_pop1.getSegments());
				nodeSet.put(n.getId(), n);
				map_nodes_pop1.put(n_pop1.getId(), n.getId());
				id_leaf++;
			}
		}//fine for insertion nodes
		
		//**** COPY ACTIVE EDGES IN POP1 AS ACTIVE EDGES IN COALESCENT POP ****
		int id_edge = 0;
		for(int i = 0; i < act_edges_pop1.size(); i++) {
			act_e = pop1.getGraphEdges().get(act_edges_pop1.get(i));
			//add the edge at the list of active edges and the graph edge
			Edge e = new Edge();
			e.setIdEDGE(id_edge);
			e.setDensity(act_e.getDensity());
			//mappa nodi
			e.setId_son(map_nodes_pop1.get(act_e.getId_son()));
			e.setId_fath(-1);
			e.setLength(act_e.getLength());
			e.setRate(act_e.getRate());
			e.setSegments(act_e.getSegments());
			e.setDensity(act_e.getDensity());
			//e.setTime(act_e.getTime());
			e.setTime();
			activeEdges.add(id_edge);
			graphEdges.put(id_edge, e);
			id_edge++;
		}
		maps_leaves_roots.put(pop1.getId_pop(), map_nodes_pop1);
		
		
		//For each active edge of pop 2 take the node son of active edge and 
		//copy in this population the node and active edge and a mapping
		ArrayList<Integer> act_edges_pop2 = pop2.getActiveEdges();
		
		//**** COPY ACTIVE NODES IN POP2 AS LEAVES IN POP ****
		active_nodes = new TreeSet<Integer>();
		HashMap<Integer, Integer> map_nodes_pop2 = new HashMap<Integer, Integer>();

		for(int i = 0; i < act_edges_pop2.size(); i++) {
			//take the node of pop1 that is son of the active edge and copy it has leaf
			act_e = pop2.getGraphEdges().get(act_edges_pop2.get(i));
			
			if(!active_nodes.contains(act_e.getId_son())){
				active_nodes.add(act_e.getId_son());
				Node n_pop2 = pop2.getNodeSet().get(act_e.getId_son());
				Node n = new Node();
				n.setId(id_leaf);
				n.setLevel(n_pop2.getLevel());
				n.setRecomb(n_pop2.isRecomb());
				n.setSegments(n_pop2.getSegments());
				nodeSet.put(n.getId(), n);
				map_nodes_pop2.put(n_pop2.getId(), n.getId());
				id_leaf++;
			}
		}//fine for insertion nodes
		
		//**** COPY ACTIVE EDGES IN POP2 AS ACTIVE EDGES IN COALESCENT POP ****
		for(int i = 0; i < act_edges_pop2.size(); i++) {
			act_e = pop2.getGraphEdges().get(act_edges_pop2.get(i));
			//add the edge at the list of active edges and the graph edge
			Edge e = new Edge();
			e.setIdEDGE(id_edge);
			e.setDensity(act_e.getDensity());
			//mappa nodi
			e.setId_son(map_nodes_pop2.get(act_e.getId_son()));
			e.setId_fath(-1);
			e.setLength(act_e.getLength());
			e.setRate(act_e.getRate());
			e.setSegments(act_e.getSegments());
			e.setDensity(act_e.getDensity());
			//e.setTime(act_e.getTime());
			e.setTime();
			activeEdges.add(id_edge);
			graphEdges.put(id_edge, e);
			id_edge++;
		}
		maps_leaves_roots.put(pop2.getId_pop(), map_nodes_pop2);
		
		//***** Extant Units ****
		setExtantUnits(nodeSet.size());
		setL(activeEdges.size());
		backwardConstructionARG(startTime,endTime, false);
	}

/**
 * This procedure simulates the evolution of an ancestral population resulting from the splitting of a more recent underlying population in the scaffold
 * @param id ID of this new ancestral population
 * @param pop underlying population in the scaffold that it is split
 * @param act_lineages_leaves list of the active lineages of the population that is splitting. The nodes connected to these active lineages (roots of the population that is splitting) will become the corresponding leaves of the new ancestral population
 * @param startTime time in which the merge between pop1 and pop2 occurs creating the new ancestral population
 * @param endTime threshold of time in which the new population deriving from the merge will stop its growing
 * @param populationSize effective population size of the new population resulting from the merge
 * @param r recombination rate of the new population resulting from the merge
 * @param mu mutation rate of the new population resulting from the merge
 * @param g segment length of the sample in the new population
 */
public void generatePopulationFromSplitting(int id, PopulationARG pop, ArrayList<Integer> act_lineages_leaves, double startTime, double endTime, int populationSize, double r, double mu, int g) {
	

	//***************************** PARAMETER INITIALIZATION  ***************************************
	
	//***** Id population *****
	this.setId_pop(id);
	
	//***** Kind of Population *****
	this.setKind(1);
	
	//***** Threshold time t *****
	setThresholdT(endTime);
	
	//***** Population Size *****
	setN(populationSize);
	
	//***** Recombination Rate ****
	setRecomb(r*Math.pow(10, -8));
	
	//**** SNP Mutation Rate ****
	setMu(mu*Math.pow(10, -8));
	
	//**** Segment length g ***
	setG(g*Math.pow(10, 3));
	
	//*** Set current generation 
	this.setCurr_generation(startTime);
	
	nodeSet = new HashMap<Integer,Node>();
	graphEdges = new HashMap<Integer,Edge>();
	activeEdges = new ArrayList<Integer>();
	maps_leaves_roots = new HashMap<Integer, HashMap<Integer, Integer>>();
	mutationSet = new HashMap<Double,Mutation>();
	SNPpositionsList = new TreeSet<Double>();
	splitPoints = new TreeSet<Double>();
	splitPoints.add(0.0);
	splitPoints.add(1.0);
	popanc[0] = pop.getId_pop();
	setRecombinationsNumber(0);
	setCoalescentNumber(0);
	
	
	//********* INITIALIZATION OF LEAVE NODES and ACTIVE EDGES ************/

	//For each active edge of pop take the node son of active edge and 
	//copy in this population the node and active edge and a mapping
	
	
	//**** COPY ACTIVE NODES IN POP AS LEAVES IN POP ****
	TreeSet<Integer> active_nodes = new TreeSet<Integer>();
	HashMap<Integer, Integer> map_nodes_pop = new HashMap<Integer, Integer>();
	int id_leaf = 0;
	Edge act_e;
	for(int i = 0; i < act_lineages_leaves.size(); i++) {
		//take the node of pop1 that is son of the active edge and copy it has leaf
		act_e = pop.getGraphEdges().get(act_lineages_leaves.get(i));
		
		if(!active_nodes.contains(act_e.getId_son())){
			active_nodes.add(act_e.getId_son());
			Node n_pop1 = pop.getNodeSet().get(act_e.getId_son());
			Node n = new Node();
			n.setId(id_leaf);
			n.setLevel(n_pop1.getLevel());
			n.setRecomb(n_pop1.isRecomb());
			n.setSegments(n_pop1.getSegments());
			nodeSet.put(n.getId(), n);
			map_nodes_pop.put(n_pop1.getId(), n.getId());
			id_leaf++;
		}
	}//fine for insertion nodes
	
	//**** COPY ACTIVE EDGES IN POP AS ACTIVE EDGES IN COALESCENT POP ****
	int id_edge = 0;
	for(int i = 0; i < act_lineages_leaves.size(); i++) {
		act_e = pop.getGraphEdges().get(act_lineages_leaves.get(i));
		//add the edge at the list of active edges and the graph edge
		Edge e = new Edge();
		e.setIdEDGE(id_edge);
		e.setDensity(act_e.getDensity());
		//mappa nodi
		e.setId_son(map_nodes_pop.get(act_e.getId_son()));
		e.setId_fath(-1);
		e.setLength(act_e.getLength());
		e.setRate(act_e.getRate());
		e.setSegments(act_e.getSegments());
		e.setDensity(act_e.getDensity());
		//e.setTime(act_e.getTime());
		e.setTime();
		activeEdges.add(id_edge);
		graphEdges.put(id_edge, e);
		id_edge++;
	}
	
	maps_leaves_roots.put(pop.getId_pop(), map_nodes_pop);
	
	setExtantUnits(nodeSet.size());
	setL(activeEdges.size());
	backwardConstructionARG(startTime,endTime,true);
}
	
	/**
	 * This procedure represents the core of the backward simulation of the ARG 
	 * @param startTime the time in which the process of backwards reconstruction starts from the leaves of the ARG 
	 * @param endTime threshold of time in which the process ends
	 * @param split set to true if the population to reconstruct derives from the splitting of another population
	 * @return true if at the end of the construction process there is only one active lineage meaning that the GMRCA has been reached. 
	 */
	public boolean backwardConstructionARG(double startTime, double endTime, boolean split){
		
		//*************** Construction of the ARG backwards starting from the leaves ****************/
		
		this.setCurr_generation(startTime);
		double time = startTime;
		double binCoef = 0;
		double SumRates = 0;
		double lamda = 0;
		boolean done = false;
		
		if(startTime == 0 || split)
			done = true;
		
		//Until we have created the GRMCA or we reach a threshold of time 
		while(activeEdges.size() != 1 && getCurr_generation() < endTime)
		{
			
			//COMPUTE THE TIME T TO THE NEXT EVENT USING THE EXPONENTIAL DISTRIBUTION
			
			//1) Compute  binomial coefficient (L 2)
			binCoef = BinomialCoefficient.binomialCoeff(getL(), 2);
		
			//2) Compute THE SUM OF ALL r_l for each active lineage
			SumRates = computeSumRates();
				
			//3) Compute lamda
			lamda = binCoef+SumRates;
			
			if(done) {
			//4) Compute the time to the next event using exponential distribution
			time = computeNextTimeEvent(lamda);
			
			//5) Update the current generation
			setCurr_generation(computeNextGeneration(time));
			}
			else 
				//this turn we use t has time of the next event
				done = true; 
			
			if(getCurr_generation() < endTime) {
				
				//COMPUTE THE TYPE OF THE NEXT EVENT (COALESCENT or RECOMBINATION)
				int index;
				
				if(getRecomb() > 0){
				//ARG with recombination nodes
				 index = computeIndexNextEvent(lamda, binCoef);
				}
				else{
				//Tree only coalescent nodes
				 index = 0;
				}
				
				Random random = new Random();
				
				//***** CASE 1: THE NEXT EVENT IS COALESCENT *****
				if(index == 0) {
					
					this.setCoalescentNumber(this.getCoalescentNumber()+1);
					//Create a new Coalescent node
					Node coal_node = new Node(nodeSet.size(), false, getCurr_generation());
					nodeSet.put(coal_node.getId(), coal_node);
				
					//Pick randomly the first active edge from the set of active edges
					int j = random.nextInt(activeEdges.size());
					int id_first_edge = activeEdges.get(j);
					
					//to understand how many consecutive coalescent and recombination we have (with the same two lineages)
					int son1 = graphEdges.get(id_first_edge).getId_son();
					//Set the father of the active edge selected
					graphEdges.get(id_first_edge).setId_fath(coal_node.getId());
					//Set the time passed for this event (computed by the exponential distribution)
					//graphEdges.get(id_first_edge).setTime(nodeSet.get(coal_node.getId()).getLevel()-nodeSet.get(son1).getLevel());
					
					
					//Get the arraylist of Intervals/Segments of the left edge
					ArrayList<Interval> left_segments = graphEdges.get(activeEdges.get(j)).getSegments();
					
					//Remove the lineage from the set of the active ones
					activeEdges.remove(j);
					
					//Pick randomly the second active edge from the set of active edges
					j = random.nextInt(activeEdges.size());
					int id_second_edge = graphEdges.get(activeEdges.get(j)).getIdEDGE();
					
					int son2 = graphEdges.get(id_second_edge).getId_son();
					
					//Set the father of the selected active lineage as coalescent node
					graphEdges.get(id_second_edge).setId_fath(coal_node.getId());
					//Set the time passed for this event (computed by the exponential distribution)
					//graphEdges.get(id_second_edge).setTime(nodeSet.get(coal_node.getId()).getLevel()-nodeSet.get(son2).getLevel());
					
					//Get the arraylist of Intervals/Segments of the right edge
					ArrayList<Interval> right_segments = graphEdges.get(activeEdges.get(j)).getSegments();
					
					//Remove the lineage from the set of active ones
					activeEdges.remove(j);
					
					//Compute the union of the left and the right intervals
					MergeIntervals merge = new MergeIntervals();
					ArrayList<Interval> union = merge.ConcatTwoIntervalList(left_segments, right_segments);
					union = merge.merge(union);
					
					//Memorize the new set of segments in the coalescent node and the two sons
					coal_node.setSegments(union);
					coal_node.setIDsonsx(son1);
					coal_node.setIDsondx(son2);
				
					//Create a new active edge outgoing from the new coalescent node
					Edge new_edge = new Edge(graphEdges.size(), coal_node.getId(), -1, union);
					new_edge.computeLength();
					new_edge.computeRate(getN(), getRecomb(), getG());
					graphEdges.put(new_edge.getIdEDGE(), new_edge);
					activeEdges.add(new_edge.getIdEDGE());
					
				}
				
				//***** CASE 2: THE NEXT EVENT IS RECOMBINATION *****
				else{
					
					//Increment the number of recombinations events
					setRecombinationsNumber(this.getRecombinationsNumber()+1);
					
					//Create a new recombination node
					Node recomb_node = new Node(nodeSet.size(), true, getCurr_generation());
					nodeSet.put(recomb_node.getId(), recomb_node);
					
					//get the ID of the active edge selected
					int id_edge = activeEdges.get(index-1);
					
					//Set the father of the active edge as the ID of the new recomb node 
					graphEdges.get(id_edge).setId_fath(recomb_node.getId());
					
					//Set the time computed by exponential distribution
					int son = graphEdges.get(id_edge).getId_son();
					//graphEdges.get(id_edge).setTime(nodeSet.get(recomb_node.getId()).getLevel()-nodeSet.get(son).getLevel());
					
					//Set the list of solid intervals
					recomb_node.setSegments(graphEdges.get(activeEdges.get(index-1)).getSegments());
					recomb_node.setIDsonsx(son);
					
					//Remove the picked edge from the active ones
					activeEdges.remove(index-1);
					
					//Create two new active edges that have as son the new recombination node
					/*ArrayList<ArrayList<Interval>> splitted = SplittingIntervals.split(recomb_node.getSegments(), recomb_node.getId(), this);
					Edge left_edge = new Edge(graphEdges.size(), recomb_node.getId(), -1, splitted.get(0));
					left_edge.computeLength();
					left_edge.computeRate(getN(), getRecomb(), getG());
					graphEdges.put(left_edge.getIdEDGE(), left_edge);
					activeEdges.add(left_edge.getIdEDGE());
					
					Edge right_edge = new Edge(graphEdges.size(), recomb_node.getId(), -1, splitted.get(1));
					right_edge.computeLength();
					right_edge.computeRate(getN(), getRecomb(), getG());
					graphEdges.put(right_edge.getIdEDGE(), right_edge);
					activeEdges.add(right_edge.getIdEDGE());*/		
				}
				
				//update L
				setL(activeEdges.size());
				
			} 
		} //end while
		if(activeEdges.size()==1)
			return true;
		else //reached the threshold
			return false;		
	} //end procedure generate structure of ARG
	
	
	/**
	 * The function prints to video the id and the kind of the population 
	 */
	public void printKind(){
		System.out.println("Population "+ getId_pop());
		System.out.println("Kind "+ this.getKind());
	}
	
	/**
	 * This procedure splits the active lineages of the population in two sets that will represent the two new ancestral population to simulate
	 * @return two ArrayLists of Integer objects. Each lists contains a list of active lineage IDs in the current ARG
	 */
	public ArrayList<ArrayList<Integer>> splitPopulation(){
		
		
		//Set of active nodes that need to be split as leaves nodes of two population
		ArrayList<Integer> active_nodes = new ArrayList<Integer>();
		
		//For each active edge memorize in the tree set the son of the active_edge
		for(int i = 0; i < activeEdges.size(); i++){
			
			Edge e = this.getGraphEdges().get(activeEdges.get(i));
			if(!active_nodes.contains(e.getId_son()))
				active_nodes.add(e.getId_son());
		}
		System.out.println("Number of active edges: "+activeEdges.size());
		System.out.println("Number of active nodes: "+active_nodes.size());
		
		ArrayList<ArrayList<Integer>> activeEdges_splitted = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> act_edges_left = new ArrayList<Integer>();
		ArrayList<Integer> act_edges_right = new ArrayList<Integer>();
		
		if(active_nodes.size() > 1){
			System.out.println("Splitting population");
			int half = active_nodes.size()/2;
			Integer id_node;
			//take the first half of nodes and the corresponding active lineages
			for(int i = 0; i < half; i++){
				id_node = active_nodes.get(i);
				for(int j = 0; j < activeEdges.size(); j++){
					if(this.graphEdges.get(activeEdges.get(j)).getId_son()==id_node)
						act_edges_left.add(activeEdges.get(j));	
				}
			}

			for(int i = half; i < active_nodes.size(); i++){
				id_node = active_nodes.get(i);
				for(int j = 0; j < activeEdges.size(); j++){
					if(this.graphEdges.get(activeEdges.get(j)).getId_son()==id_node)
						act_edges_right.add(activeEdges.get(j));	
				}
			}
		
		}
		else{
			System.out.println("Cannot split the population because there is only one active lineage");
		}
		
		activeEdges_splitted.add(act_edges_left);
		activeEdges_splitted.add(act_edges_right);
		return activeEdges_splitted;
	}
}
