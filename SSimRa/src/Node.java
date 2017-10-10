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
import java.util.TreeSet;
/**
 * The class Node represents a node in an ARG. A node can be either a leaf or a recombination node or a coalescent node.
 * An instance of Node stores and manages several information: its ID, its level (that is the time in which the node has been created), the segment
 * stored in the node that is represented by a list of Interval objects, the set of SNP mutations.  
 * 
 * @author Anna Paola Carrieri
 *
 */
public class Node {
	
	//id of the node
	private int id;
	//true if node is recombination node. In this case the node has a split point
	private boolean recomb;
	//level = time of the node
	private double level;
	//true if the node is a leaf
	private boolean leaf;
	//the split point in case the node has been split for a recombination event
	private double splitPoint;
	//the segment carried by the node is represented by a list of intervals
    private ArrayList<Interval> segments;
    //id of the node son on the left
    private int IDsonsx;
	//id of the node son on the right
    private int IDsondx;
    //global id in the whole scaffold
    private int global_ID;
    //list of Str objects contained in th endoe
    private ArrayList<Str> strs;
	//set of SNP mutation positions
	private TreeSet<Double> mutation_set;
	//this field is true if the node is a root of the ARG
	private boolean root = false;
	
	/**
	 * Constructor of Node without parameters
	 */
	Node(){
		global_ID = -1;
	}
	/**
	 * Constructor for the class Node with parameters
	 * @param id integer representing the ID of the node
	 * @param recomb boolean value that is true if the node is a recombination node
	 * @param level real value representing the time when the node is created
	 */
	public Node(int id, boolean recomb, double level) {
		super();
		this.id = id;
		this.recomb = recomb;
		this.level = level;
	}
	/**
	 * This function returns the ID of the node
	 * @return an integer ID of the node
	 */
	public int getId() {
		return id;
	}
	/**
	 * The procedure sets the ID of the node
	 * @param id integer representing the ID of the node
	 */
	public void setId(int id) {
		this.id = id;
	}
	/**
	 * If the node is a recombination node, the function returns the split point in [0,1]
	 * @return a real value representing the split point in [0,1]
	 */
	public double getSplitPoint() {
		return splitPoint;
	}
	/**
	 * The procedure sets the split point (only if the node is a recombination node)
	 * @param splitPoint real value is split point in [0,1]
	 * 
	 */
	public void setSplitPoint(double splitPoint) {
		this.splitPoint = splitPoint;
	}
	/**
	 * The function returns an ordered list of intervals representing the whole genetic material stored in the node
	 * @see Interval
	 * @return an array list of objects of the class Interval representing the genetic material stored in the node
	 */
	public ArrayList<Interval> getSegments() {
		return segments;
	}
	/**
	 * The procedure sets the list of Interval objects representing the genetic material stored in the node
	 * @see Interval
	 * @param segments ArrayList of objects of the class Interval
	 */
	public void setSegments(ArrayList<Interval> segments) {
		this.segments = segments;
	}
	/**
	 * The function returns the time of the node or its age (the time when the node has been created)
	 * @return a real value representing the age of the node 
	 */
	public double getLevel() {
		return level;
	}
	/**
	 * 
	 * The procedure sets the age of the node (the time when the node has been created)
	 * @param level real value representing the time of the node
	 */
	public void setLevel(double level) {
		this.level = level;
	}
	/**
	 * The function returns true if the node is a recombination node
	 * @return a boolean value that is true if the node is a recombination node
	 */
	public boolean isRecomb() {
		return recomb;
	}
	/**
	 * The function sets the value to true if the node is recombination, false otherwise
	 * @param recomb boolean value that is true if the node is recombination, false otherwise
	 */
	public void setRecomb(boolean recomb) {
		this.recomb = recomb;
	}
	/**
	 * The function returns true if the node is a leaf of the ARG
	 * @return a boolean value that is true if the node is a leaf of the ARG
	 */
	public boolean isLeaf() {
		return leaf;
	}
	/**
	 * The procedure sets the value to true if the node is a leaf of the ARG
	 * @param leaf boolean value that is true if the node is a leaf of the ARG, false otherwise
	 */
	public void setLeaf(boolean leaf) {
		this.leaf = leaf;
	}
	/**
	 * The function returns true if the node is root of the ARG (it has no node father)
	 * @return a boolean value that is true if the node is a root of the ARG, otherwise false
	 */
	public boolean isRoot() {
		return root;
	}
	/**
	 * The procedure sets the value to true if the node is a root of the ARG (it has no node father)
	 * @param root a boolean value that is true if the node is a root of the ARG
	 */
	public void setRoot(boolean root) {
		this.root = root;
	}
	/**
	 * The function returns an ordered set of SNP mutations. The SNPs are in ascending order based on their locations in the segment
	 * @return an ordered set of real values representing the position of the SNPs in the segment stored in the node
	 */
	public TreeSet<Double> getMutation_set() {
		return mutation_set;
	}
	/**
	 * The procedure sets an ordered set of SNP mutations. The SNPs are in ascending order based on their locations in the segment
	 * @param mutation_set set of real values representing the positions of the SNP mutations in the segment stored in the node
	 */
	public void setMutation_set(TreeSet<Double> mutation_set) {
		this.mutation_set = mutation_set;
	}
	/**
	 * The function returns the list of STRs in the node
	 * @return the ArrayList of STR objects in the node
	 */
	public ArrayList<Str> getStrs() {
		return strs;
	}
	/**
	 * The procedure sets the list of STRs in the node
	 * @param strs the ArrayList of the STR objects in the node
	 */
	public void setStrs(ArrayList<Str> strs) {
		this.strs = strs;
	}
	/**
	 * The function returns the global ID of the node with respect to all nodes in the whole scaffold
	 * @return an integer representing the global ID of the node
	 */
	public int getGlobal_ID() {
		return global_ID;
	}
	/**
	 * The procedures sets the global ID of the node with respect to all nodes in the whole scaffold
	 * @param global_ID integer representing the global ID of the node 
	 */
	public void setGlobal_ID(int global_ID) {
		this.global_ID = global_ID;
	}
	/**
	 * The functions returns the ID of the left son of the node 
	 * @return integer representing the id of the left son node
	 */
	public int getIDsonsx() {
		return IDsonsx;
	}
	/**
	 * The procedure sets the ID of the left son of the node
	 * @param iDsonsx integer representing the iD of the left son of the node
	 */
	public void setIDsonsx(int iDsonsx) {
		IDsonsx = iDsonsx;
	}
	/**
	 * The functions returns the ID of the right son of the node 
	 * @return integer representing the id of the right son node
	 */
	public int getIDsondx() {
		return IDsondx;
	}
	/**
	 * The procedure sets the ID of the right son of the node
	 * @param iDsondx integer representing the iD of the right son of the node
	 */
	public void setIDsondx(int iDsondx) {
		IDsondx = iDsondx;
	}
	

}

