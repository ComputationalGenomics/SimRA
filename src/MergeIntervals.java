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
import java.util.Collections;
import java.util.Comparator;

/**
 * The class MergeIntervals is useful to do some operations on objects of the class Interval, such as the merge of two objects of Intervals
 * @see Interval
 */
public class MergeIntervals{
/**
 * 	
 */
public MergeIntervals(){
	
}
	
/**
 * This function takes in input a list of objects of the class Interval and merge them ordering them in case there are holes between intervals
 * @param intervals a list of object of the class Interval that has to be merged and ordered
 * @return an ArrayList of objects of the class Interval that have been merged and ordered
 */
public ArrayList<Interval> merge(ArrayList<Interval> intervals) {
        
        if(intervals.size() == 0)
            return intervals;
        if(intervals.size() == 1)
            return intervals;
        
        Comparator<Interval> IntervalsComp = new IntervalComparator();
        // sort intervals by using self-defined Comparator
        Collections.sort(intervals, IntervalsComp);
        
        Interval first = intervals.get(0);
        double start = first.start;
        double end = first.end;
        
        ArrayList<Interval> result = new ArrayList<Interval>();
        
        for(int i = 1; i < intervals.size(); i++){
            Interval current = intervals.get(i);
            if(current.start <= end){
                end = Math.max(current.end, end);
            }else{
                result.add(new Interval(start, end));
                start = current.start;
                end = current.end;
            }
            
        }
        
        result.add(new Interval(start, end));
        
        return result;
        
    }
/**
 * This procedure prints the list of objects of the class Interval 
 * @param l ArrayList of objects of the class Interval to be print
 */
public static void printListIntervals(ArrayList<Interval> l){
	
	for(int i = 0; i < l.size(); i++ ) {
		l.get(i).printInterval();
	}
	
}

/**
 * This procedure creates and returns a String representing a list of Interval objects
 * @param l is an ArrayList of objects of the class Interval
 * @return a String object representing the list of Interval objects
 */
public static String intervalList(ArrayList<Interval> l){
	
	String segms = new String("");
	Interval intv;
	
	for(int i = 0; i < l.size()-1; i++ ) {
		intv = l.get(i);
		segms.concat("["+intv.getStart()+" - "+intv.getEnd()+"] , ");
		
		
	}
	
	intv = l.get(l.size()-1);
	segms.concat("["+intv.getStart()+" - "+intv.getEnd()+ "]");
	
	return segms;
	
}

/**
 * 
 *
 */
class IntervalComparator implements Comparator<Interval>{
    public int compare(Interval o1, Interval o2){
        int r = 0;        
        if(o1.start < o2.start)
        	r = -1;
        else if(o1.start > o2.start)
        	r = 1;
        return r;
    }
}

/**
 * This function connects two list of Interval objects of the class Interval
 * @param l ArrayList of Interval objects that will be precede the other list 
 * @param r ArrayList of Interval objects that will be follow the other list 
 * @return the ArrayList of Interval object resulting from the concatenation of two lists
 */
public ArrayList<Interval> ConcatTwoIntervalList(ArrayList<Interval> l, ArrayList<Interval> r) {
	
	ArrayList<Interval> ris = new ArrayList<Interval>();
	
	for(int i = 0; i < l.size(); i++) {
		ris.add(l.get(i));
	}
	
	for(int j = 0; j < r.size(); j++) {
		ris.add(r.get(j));
	}
	
	return ris;
}

}

