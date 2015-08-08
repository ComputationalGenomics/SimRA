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


}
