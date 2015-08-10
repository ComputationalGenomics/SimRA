/*
  Copyright 2015 IBM Corporation


Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under 
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF 
ANY KIND, either express or implied. 
See the License for the specific language governing permissions and limitations under the License.


*/


import java.util.Iterator;
import java.util.TreeSet;
/**
 * The class Mutation represents a Single Neuclotide Polymorphism. An instance of this class stores and manages
 * several information associated with a mutation: 
 * - the ID of the mutation is its position (real value) in the interval [0,1];
 * - the Interval (solid) in which the mutation is places
 * - the ID of the edge in which the mutation appears 
 * - the set of leaves that have that mutation
 * 
 * @author Anna Paola  Carrieri
 */

public class Mutation {
	
	private double positionID;
	private int IDedge;
	private TreeSet<Integer> leaves;
	private Interval segment;
	private int IDson;
	private int IDfather;
	private int ID_original_pop;
	
	/**
	 * This function returns the position of the SNP mutation in [0,1]
	 * @return a real value that is the position of the SNP in [0,1] representing its ID
	 */
	public double getPositionID() {
		return positionID;
	}
	/**
	 * This procedure sets the position of the SNP mutation in [0,1]
	 * @param positionID real value that is the position of the SNP in [0,1] representing its ID
	 */
	public void setPositionID(double positionID) {
		this.positionID = positionID;
	}
	
	/**
	 * This procedure returns the set of IDs of leaves that have this instance of the class Mutation 
	 * @return the set of IDs of leaves that have this instance of the class Mutation 
	 */
	public TreeSet<Integer> getLeaves() {
		return leaves;
	}
	/**
	 * This procedure returns the set of IDs of leaves that have this instance of the class Mutation 
	 * @param leaves the set of IDs of leaves that have this instance of the class Mutation 
	 */
	public void setLeaves(TreeSet<Integer> leaves) {
		this.leaves = leaves;
	}
	/**
	 * This function returns the specific instance of Interval where this object of Mutation is located
	 * @see Interval
	 * @return the object Interval that has the mutation 
	 */
	public Interval getSegment() {
		return segment;
	}
	/**
	 * This function returns the specific instance of Interval where this object of Mutation is located
	 * @see Interval
	 * @param segment the object of the class Interval that has the mutation 
	 */
	public void setSegment(Interval segment) {
		this.segment = segment;
	}
	/**
	 * This function returns the ID of the node son of the edge where this instance of Mutation occurred
	 * @return an integer that is the ID of the node son of the edge where this instance of Mutation occurred
	 */
	public int getIDson() {
		return IDson;
	}
	/**
	 * This procedure sets the ID of the node son of the edge where this instance of Mutation occurred
	 * @param iDson an integer that is the ID of the node son of the edge where this instance of Mutation occurred
	 */
	public void setIDson(int iDson) {
		IDson = iDson;
	}
	/**
	 * This function returns the ID of the node parent of the edge where this instance of Mutation occurred
	 * @return an integer that is the ID of the node parent of the edge where this instance of Mutation occurred
	 */
	public int getIDfather() {
		return IDfather;
	}
	/**
	 * This procedure sets the ID of the node parent of the edge where this instance of Mutation occurred
	 * @param iDfather an integer that is the ID of the node parent of the edge where this instance of Mutation occurred
	 */
	public void setIDfather(int iDfather) {
		IDfather = iDfather;
	}
	/**
	 * This function returns the ID of the edge where this instance of Mutation occurred
	 * @return an integer that is the ID of the edge where this instance of Mutation occurred
	 */
	public int getIDedge() {
		return IDedge;
	}
	/**
	 * This procedure sets the ID of the edge where this instance of Mutation occurred
	 * @param iDedge an integer that is the ID of the edge where this instance of Mutation occurred
	 */
	public void setIDedge(int iDedge) {
		IDedge = iDedge;
	}
	/**
	 * This function returns the ID of the population where the SNP originally occurred
	 * @return an integer representing the ID of the population where the SNP originally occurred
	 */
	public int getID_original_pop() {
		return ID_original_pop;
	}
	/**
	 * This procedure sets the ID of the population where the SNP originally occurred
	 * @param iD_original_pop an integer representing the ID of the population where the SNP originally occurred
	 */
	public void setID_original_pop(int iD_original_pop) {
		ID_original_pop = iD_original_pop;
	}
	/**
	 * This procedure prints the information about the object Mutation
	 */
	public void printMutation(){
		System.out.println("MUTATION "+this.getPositionID());
		System.out.println("Id edge "+this.getIDedge());
		System.out.println("Id node father "+this.getIDfather());
		System.out.println("Id node son "+this.getIDson());
		System.out.println("Set of leaves that has the mutation");
		this.printLeaves();
		
	}
	
	/**
	 * This procedure prints the leaves of the ARG that have that mutation
	 */
	public void printLeaves(){
		Iterator<Integer> it_l = leaves.iterator();
		while(it_l.hasNext()){
			System.out.print(it_l.next()+" ");
		}
		System.out.println();
	}
	
}
