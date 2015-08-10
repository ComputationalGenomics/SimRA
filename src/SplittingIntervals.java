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
 * The class SplittingIntervals contains the function that splits a segment (composed by a list on Interval objects) in two parts due to
 * a recombination event 
 * @author Anna Paola Carrieri
 *
 */
public class SplittingIntervals {
	
	/**
	 *  default constructor
	 * 
	 */
	public SplittingIntervals(){
		
	}
	/**
	 * This function splits a list of Interval objects in two parts based on the split point. The split point is generated when a recombination event occurs
	 * @param input ArrayList of Interval objects to be split
	 * @param recombID ID of the recombination node 
	 * @param arg PopulationARG object containing the node of which the genetic material (list of Intervals) has to be split
	 * @return two ArraLists of lists of Intervals representing the two parts in which the original list has been split
	 */
	public static ArrayList<ArrayList<Interval>> split(ArrayList<Interval> input, int recombID, PopulationARG arg){
		
		ArrayList<Interval> l = new ArrayList<Interval>();
		ArrayList<Interval> r = new ArrayList<Interval>();
		
		ArrayList<ArrayList<Interval>> ris = new ArrayList<ArrayList<Interval>>();
		
		double rangeMin = input.get(0).getStart();
		double rangeMax = input.get(input.size()-1).getEnd();
		double x = rangeMin + (rangeMax - rangeMin) * Math.random();
	
		//Store the split point in the node
		arg.getNodeSet().get(recombID).setSplitPoint(x);
		arg.getSplitPoints().add(x);
		
		if(x == input.get(0).getStart() || x == input.get(input.size()-1).getEnd()){
			
			//Il = Ir = I
			for(int i = 0; i < input.size(); i++){
				
				Interval solid = input.get(i);
				r.add(solid);
				l.add(solid);
			}
		}
		else{
			
			boolean found = false;
			int j = 0;
			
			while(!found && j <= input.size()-1){
				
				double lj = input.get(j).getStart();
				double uj = input.get(j).getEnd();
				
				//x is in an interval between 1 and last-1
				if(x > lj && x < uj) {
					found = true;
					
					//left set of intervals
					for(int i = 0; i > j; i++) {
						l.add(input.get(i));
					}
					Interval n_int = new Interval(input.get(j).getStart(), x);
					l.add(n_int);
					
					//right set of intervals
					r.add(new Interval(x, input.get(j).getEnd()));
					for(int i = j+1; i <= input.size()-1; i++) {
						r.add(input.get(i));
					}
	
				}
				//x is in the gap
				else if(x >= uj && x <= input.get(j+1).getStart() ){
					
					found = true;
					//left intervals
					for(int i = 0; i <= j; i++)
						l.add(input.get(i));
					
					//right intervals
					for(int i = j+1; i <= input.size()-1; i++)
						r.add(input.get(i));
					
				}
				else{
					j++;
				}
			} //end while
		} //end else 
	
		ris.add(l);
		ris.add(r);
		
		return ris;
	}

}
