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
 * This class represents a Short Tandem Repeat
 *  
 * @author Anna Paola Carrieri
 *
 */
public class Str {
	//number of repeats for the specific STR locus
	private int num_repeats;
	//STR locus
	private double location;
	//STR presence flag in a given node where it is stored
	private boolean presence;
	//mutation rate for the STR
	private double mu_rate;
	
	/**
	 * default constructor of the class STR
	 */
	public Str() {
	
	}
	
	/**
	 * Constructor with parameters
	 * @param num_repeats  number of repeats for the specific STR locus
	 * @param location  STR locus
	 * @param presence  STR presence flag in a given node where it is stored
	 * @param mu_rate  mutation rate for the STR
	 */
	public Str(int num_repeats, double location, boolean presence,
			double mu_rate) {
		super();
		this.num_repeats = num_repeats;
		this.location = location;
		this.presence = presence;
		this.mu_rate = mu_rate;
	}
	/**
	 * Function that returns the number of repeats for the given locus
	 * @return integer number of repeats for the given locus
	 */
	public int getNum_repeats() {
		return num_repeats;
	}
	/**
	 * Procedure that sets the number of repeats for the given locus
	 * @param num_repeats integer number of repeats for the given locus
	 */
	public void setNum_repeats(int num_repeats) {
		this.num_repeats = num_repeats;
	}
	/**
	 * Function that returns the STR locus
	 * @return real value representing the STR locus in the segment [0,1]
	 */
	public double getLocation() {
		return location;
	}
	/**
	 * Procedure that sets the STR locus
	 * @param location real value representing the STR locus in the segment [0,1]
	 */
	public void setLocation(double location) {
		this.location = location;
	}
	/**
	 * Function that returns true if the STR is present in the given node based on the fact that the node carries genetic material containing the given STR locus
	 * @return a boolean value that is true if the STR is present in the given node based on the fact that the node carries genetic material containing the given STR locus
	 */
	public boolean isPresence() {
		return presence;
	}
	/**
	 * Function that sets the presence of the STR in the given node
	 * @param presence boolean value that is true if the STR locus is present in the given node
	 */
	public void setPresence(boolean presence) {
		this.presence = presence;
	}
	/**
	 * Function that returns the mutation rate for the STR locus
	 * @return a real value representing mutation rate for the STR locus
	 */
	public double getMu_rate() {
		return mu_rate;
	}
	/**
	 * Procedure that sets the mutation rate for the STR locus
	 * @param mu_rate real value representing the mutation rate for the STR locus
	 */
	public void setMu_rate(double mu_rate) {
		this.mu_rate = mu_rate;
	}
	
}
