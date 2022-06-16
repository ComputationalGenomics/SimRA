import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class MergeIntervals{
	
public MergeIntervals(){
	
}
	
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

public static void printListIntervals(ArrayList<Interval> l){
	
	for(int i = 0; i < l.size(); i++ ) {
		l.get(i).printInterval();
	}
	
}

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


