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

import org.apache.commons.math3.distribution.PoissonDistribution;

/**
 * The class represents an edge in an ARG.
 * An edge has is own ID and the id of the its son node and its father node (instances of the class Node). If the father has value -1 this means that the edge is actually an
 * active lineage. This class has a private field representing the segment carried by the edge. The segment is here represented as a list of Interval instances. 
 * The class Interval represents a solid unit carried by the edge.
 * The density of the segment carried by the edge is the sum of the solid unit lengths. While the length of the segment is computed considering also the eventual holes among the solids.
 * An edge has its rate for recombination event that depends from the length of the segment that is carrying. 
 * When the edge is not active the id of the father node is a positive integer value and the time passed between the father event node and the son event node has to be computed. 
 *
 * @see Interval
 * @see Node
 * 
 */
 
 	
	//id of the edge
	private int idEDGE;
	//id of the object Node that is the son of this directed edge
	private int id_son;
	//id of the object Node that is the father of this directed edge
	private int id_fath;
	//List of solid segments carried by the edges
	private ArrayList<Interval> segments;
	//Recombination rate r_l associated to the lineage l
	private double rate;
	//Length of the whole segment carried by the edge (holes are included)
	private double length;
	//time passed between the event node father and the event node son of the edge = time of the node father - time of the node son
	private double time;
	//density = sum of lengths of solid intervals (holes are not considered)
	private double density;
	//nMut = #of mutations in this edge derived from a Poisson distribution
	private int nMut = 0;
	//for each STR locus we have the value delta 
	private HashMap<Double,Integer> deltaStr;
	//number of mutations/steps for str
	private double num_mut_str;
	
	/**
	 * Constructor with input parameters
	 * @param idEDGE  id of the new edge
	 * @param id_son  id of the son node of the edge
	 * @param id_fath  id of the father node of the edge
	 * @param s  list of intervals representing the segment carried by the edge
	 */
	public Edge(int idEDGE, int id_son, int id_fath, ArrayList<Interval> s) {
		super();
		this.id_son = id_son;
		this.id_fath = id_fath;
		this.idEDGE = idEDGE;
		this.segments = s;
	}
	/**
	 * Constructor without parameters
	 */
	public Edge() {	
	}
	/**
	 * This function returns the ID of the edge
	 * @return the id of the edge
	 */
	public int getIdEDGE() {
		return idEDGE;
	}
	/**
	 * This procedure sets the id of the edge
	 * @param idEDGE integer representing the id of the edge
	 */
	public void setIdEDGE(int idEDGE) {
		this.idEDGE = idEDGE;
	}
	/**
	 * The function returns the id of the node son of the edge
	 * @return the id of the son node
	 */
	public int getId_son() {
		return id_son;
	}
	/**
	 * The procedure sets the id of the node son of the edge
	 * @param id_son  integer representing the id of the son node of the edge
	 */
	public void setId_son(int id_son) {
		this.id_son = id_son;
	}
	/**
	 * The function returns the id of the node that is father of the edge
	 * @return the id of the father node
	 */
	public int getId_fath() {
		return id_fath;
	}
	/**
	 * The procedure sets the id of the father node
	 * @param id_fath  integer representing the id of the father node of the edge
	 */
	public void setId_fath(int id_fath) {
		this.id_fath = id_fath;
	}
	/**
	 * The function returns the list of solid genetic material carried by an edge
	 * @see Interval
	 * @return the list of intervals (solids) carried by the edge
	 */
	public ArrayList<Interval> getSegments() {
		return segments;
	}
	/**
	 * The procedure sets the list of solid genetic material carried by an edge
	 * @param segments list of intervals (solids) carried by the edge
	 */
	public void setSegments(ArrayList<Interval> segments) {
		this.segments = segments;
	}
	/**
	 * The function returns the recombination rate associated to the edge
	 * @return the recombination rate associated to the edge
	 */
	public double getRate() {
		return rate;
	}
	/**
	 * The procedure sets the recombination rate for the edge
	 * @param rate real value representing the recombination rate for the edge
	 */
	public void setRate(double rate) {
		this.rate = rate;
	}
	/**
	 * The function returns the number of the SNP mutations occurred in the interval of time represented by the edge
	 * @return nMut number of the SNP mutations occurred in the interval of time represented by the edge
	 */
	public int getnMut() {
		return nMut;
	}
	/**
	 * The procedure sets the number of SNPs in the interval of time represented by the edge
	 * @param nMut number of SNPs in the interval of time represented by the edge
	 */
	public void setnMut(int nMut) {
		this.nMut = nMut;
	}
	/**
	 * The function returns the density of the genetic material carried by the edge (hole lengths are not considered)
	 * @return the density of the segment 
	 */
	public double getDensity() {
		return density;
	}
	/**
	 * The procedure sets the density of the genetic material carried by the edge
	 * @param density double value representing the sum of the solids carried by the edge (hole lengths are not considered)
	 */
	public void setDensity(double density) {
		this.density = density;
	}
	/**
	 * The function returns the time of the edge that is the length of the edge representing the time of the edge passed between the father node  and the son node
	 * @return the time of the edge
	 */
	public double getTime() {
		return time;
	}	
	/**
	 * The procedure sets set the time of the edge
	 * @param length_edge the length of the edge that represents the time of the edge passed between the father node  and the son node
	 */
	public void setTime(double length_edge) {
		time = length_edge;
	}
	/**
	 * The function returns he length of the whole segment carried by the edge 
	 * @return the length of the whole segment carried by the edge (eventual holes are considered in the length of the segment)
	 */
	public double getLength() {
		return length;
	}
	/**
	 * The procedure sets the length of the segment carried by the edge
	 * @param length = the length of the whole segment carried by the edge (eventual holes are considered in the length of the segment)
	 */
	public void setLength(double length) {
		this.length = length;
	}
	/**
	 * The function returns the total number of STR mutations for the edge
	 * @return total number of STR mutations for the edge
	 */
	public double getNum_mut_str() {
		return num_mut_str;
	}
	/**
	 * The procedure sets the total number of STR mutations for the edge
	 * @param num_mut_str total number of STR mutations for the edge
	 */
	public void setNum_mut_str(double num_mut_str) {
		this.num_mut_str = num_mut_str;
	}
	/**
	 * The function returns the delta value for STR mutations 
	 * @return the delta value for STR mutations 
	 */
	public  HashMap<Double, Integer> getDeltaStr() {
		return deltaStr;
	}
	/**
	 * The procedure sets the delta value for STR mutations
	 * @param deltaStr = delta value for STR mutations
	 */
	public  void setDeltaStr(HashMap<Double, Integer> deltaStr) {
		this.deltaStr = deltaStr;
	}
	/**
	 * The procedure returns the binomial value from a binomial distribution
	 * @param n number of trials
	 * @param p probability of success
	 * @return the binomial value
	 */
	public static int getBinomial(int n, double p) {
		 //TO BE IMPLEMENTED
		}
	/**
	 * The function computes and returns the number of STR mutations for the edge by Poisson Distribution
	 * @param mu_str  STRs mutation rate 
	 * @param N effective population size
	 * @param arg instance of the class PopulationARG
	 * @see PopulationARG
	 * @return an integer representing the total number of STR mutations
	 */
	public int computeNumberMutationsStr(double mu_str, int N, PopulationARG arg) {
		   
		  //TO BE IMPLEMENTED
	}
	/**
	 * The procedure computes the delta value for each STR locus that is present in one of the segments carried by the edge
	 * Create an HashMap of Delta values where the key is the location/id of the STR
	 * @param N effective population size
	 * @param arg instance of the class PopulationARG
	 * @see PopulationARG
	 */
	public void computeDeltaStrs(int N, PopulationARG arg) {
		
	//TO BE IMPLEMENTED
	}
	/**
	 * This procedure computes and sets the recombination rate for this particular edge or lineage
	 * @param N effective population size (integer value)
	 * @param recomb recombination rate (double value parameter of the program given by the user)r
	 * @param g total length of the segment in the extant units (parameter given by the user - integer value)
	 */
	public void computeRate(int N, double recomb, double g){
		rate = N*g*recomb*length;
	}
	/**
	 * The function computes the total density of the segment carried by the edge. That is the density is the sum of the lengths of the solids
	 * @see Interval
	 * @return the density of the segment carried by the edge.
	 */
	public double computeDensity(){
		//TO BE IMPLEMENTED
	}
	/**
	 * The function computes the total length of the segment carried by the edge. The total length considers also holes that can be present among
	 * the solid intervals
	 * @see Interval
	 */
	public void computeLength(){
	//TO BE IMPLEMENTED
	}
	
	/**
	 * The function returns a String representing the solid portions of genetic material carried by the edge 
	 * @return a String representing the portions of genetic material carried by the edge 
	 */
	public String printSegments(){
		//TO BE IMPLEMENTED
	}
 
 
