/*
  Copyright 2015 IBM Corporation


Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under 
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF 
ANY KIND, either express or implied. 
See the License for the specific language governing permissions and limitations under the License.


*/

/**
 * An Event occurs among populations in the scaffold that represents the evolution of multiple populations with admixture.
 * An event can be the merge of two populations in only one ancestral at time t.
 * An event can be the split of one population in two ancestral populations at time t.
 * An event can be simply an actual population (leaf of the scaffold) that grow until something happens.
 *
 *
 */
public class Event {
	
	//id of the event
	private int id_event;
	//kind of event 0 if the event is a leaf population, 1 if is the merge of two populations, 2 if is the split
	//of two populations in one
	private int kind;
	//time in which the event occurs
	private double time_event;
	//input populations: can be either only one if the event is the split, or two if the event is the merge
	private int[] pops_input = new int[2];
	//output populations: can be either two if the even is the split, or one if the event is the merge
	private int[] pops_output = new int[2];
	//time/s in which the resulting population/s stops/stop growing or is/are involved in another event
	private double[] end_times = new double[2];
	private int N[] = new int[2];
	private int m;
	
	/**
	 * The function returns the id of the event
	 * @return an integer representing the id of the event
	 */
	public int getId_event() {
		return id_event;
	}
	/**
	 * The procedure sets the id of the event
	 * @param id_event = integer representing the id of the event
	 */
	public void setId_event(int id_event) {
		this.id_event = id_event;
	}
	/**
	 * The function returns the kind of the event
	 * @return an integer representing the kind of the event:
	 * - 0 if the event corresponds to the simulation of a leaf population
	 * - 1 if the event is the merge of two populations
	 * - 2 if the event is the split of a population in two ancestral populations
	 */
	public int getKind() {
		return kind;
	}
	/**
	 * The procedure sets the kind of the event
	 * @param kind = integer representing the kind of the event:
	 * - 0 if the event corresponds to the simulation of a leaf population
	 * - 1 if the event is the merge of two populations
	 * - 2 if the event is the split of a population in two ancestral populations
	 */
	public void setKind(int kind) {
		this.kind = kind;
	}
	/**
	 * The function returns the time/when the event occurred
	 * @return the time when the event have occurred
	 */
	public double getTime_event() {
		return time_event;
	}
	/**
	 * The procedure sets the time when the event have occurred
	 * @param time_event the time when the event have occurred
	 */
	public void setTime_event(double time_event) {
		this.time_event = time_event;
	}
	/**
	 * The function returns one/two population/s involved in an event that is/are the recent one/s in the event
	 * @return an integer vector that has at most size 2 and that stores the ids of at most two populations (the population/s 
	 * represents/represent the input of the event
	 * - it has size 1 if the event is a splitting of a population. The id of the population is stored in the
	 * first entry of the vector
	 * - it has size 2 if the event is a merging of two populations. The ids of the two populations are stored respectively
	 * in the first and in the second (and last) entries of the vector
	 * - it has size 0 if the event is the simple simulation of a leaf (extant) population
	 */
	public int[] getPops_input() {
		return pops_input;
	}
	/**
	 * The procedure sets the integer vector containing the id/ids of the population/populations that is/are involved in the event
	 * @param pops_input an integer vector containing the id/ids of the population/populations that is/are involved in the event
	 */
	public void setPops_input(int[] pops_input) {
		this.pops_input = pops_input;
	}
	/**
	 * The function returns one/two population/s involved in an event that is/are the recent one/s in the event
	 @return an integer vector that has at most size 2 and that stores the ids of at most two populations the population/s 
	 * represents/represent the output of the event
	 * - it has size 1 if the event is a merging of two populations. The id of the population is stored in the
	 * first entry of the vector and represent the population obtained after the merging
	 * - it has size 2 if the event is a splitting of two populations. The ids of the two populations are stored respectively
	 * in the first and in the second (and last) entries of the vector and represent the two populations resulting from the split 
	 * of an input population
	 * - it has size 0 if the event is the simple simulation of a leaf (extant) population
	 */
	public int[] getPops_output() {
		return pops_output;
	}
	/**
	 * The procedure sets  the vector containing the id/ids of the population/populations that is/are involved in the event and are the ancestral one/s
	 * @see getPops_output
	 * @param pops_output an integer vector that has at most size 2 and that stores the ids of at most two populations (the population/s 
	 * represents/represent the output of the event
	 */
	public void setPops_output(int[] pops_output) {
		this.pops_output = pops_output;
	}
	/**
	 * The function returns  the time/s when the population/populations in output stop growing because another event involve it/them
	 * @return a double vector storing  the times when the population/populations in output stop growing because another event involve it/them
	 */
	public double[] getEnd_times() {
		return end_times;
	}
	
	/**
	 * The procedure sets  the values of a double vector representing times when the population/populations generated by the event stop growing
	 * @param end_times a double vector that contains the values of a double vector representing times when the population/populations generated by the event stop growing
	 */
	public void setEnd_times(double[] end_times) {
		this.end_times = end_times;
	}
	/**
	 * The function returns containing the information about the population size of the population/s
	 * @return a vector containing the information about the population size of the population/s
	 */
	public int[] getN() {
		return N;
	}
	/**
	 * The procedure sets the information about the population size of the population/s
	 * @param n vector containing the information about the population size of the population/s
	 */
	public void setN(int[] n) {
		N = n;
	}
	/**
	 * The function returns the number of samples
	 * @return the number of samples
	 */
	public int getM() {
		return m;
	}
	/**
	 * The procedures sets the number of samples
	 * @param m number of samples
	 */
	public void setM(int m) {
		this.m = m;
	}
}
