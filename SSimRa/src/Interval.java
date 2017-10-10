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

/**
 * 
 * 
 * This class represents a solid carried by an edge or a node in the ARG. In other words this is a solid portion
 * of the whole segment carried by an edge and stored in a node of the graph.
 * Since the whole segment is represented in an interval [0,1], then an instance of Interval has to be a subinterval of [0,1], that is start_position >=0 and end_position <=1. 
 * 
 *
 */
class Interval {
	
	 //starting position of the interval/solid
	 double start;
	 //ending position of the interval/solid
	 double end;
	 //list of all mutations in that specific interval
	 ArrayList<Double> posMutations;
	 
	 /**
	  * Constructor without parameters that creates a [0,0] solid
	  */
	 Interval() { start = 0; end = 0; posMutations = new ArrayList<Double>(); }
	 /**
	  * Constructor with parameters representing the start and end positions of the interval/solid [start, end] 
	  * @param s start position of the interval
	  * @param e end position of the interval
	  */
	 Interval(double s, double e) { start = s; end = e; posMutations = new ArrayList<Double>(); }
	
	/**
	 * The function returns an the list of SNP positions of in the interval
	 * @return an ArrayList of Double objects that represent SNP positions in the interval
	 */
	public ArrayList<Double> getPosMutations() {
		return posMutations;
	}
	/**
	 * The procedure sets the list of SNP positions in the interval
	 * @param an ArrayList of Double objects that represent SNP positions in the interval
	 */
	public void setPosMutations(ArrayList<Double> mutations) {
		this.posMutations = mutations;
	}
	/**
	 * This function returns the start position of the interval
	 * @return a double representing the start position of the interval
	 */
	public double getStart() {
		return start;
	}
	/**
	 * This procedure sets the start position of the interval
	 * @param start double representing the position where the interval begins 0 <= start <=1
	 */
	public void setStart(double start) {
		this.start = start;
	}
	/**
	 * This function returns the end position of the interval
	 * @return a double representing the end position of the interval
	 */
	public double getEnd() {
		return end;
	}
	/**
	 * This procedure sets the end position of the interval
	 * @param end double representing the position where the interval end 0 <= end <=1
	 */
	public void setEnd(double end) {
		this.end = end;
	}
	/**
	 * This procedures print the interval in this form [start-end]
	 */
	public void printInterval(){
		 System.out.println("["+start+" - "+end+ "]");
	 }
 
}
