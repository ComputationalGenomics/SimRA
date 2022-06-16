import java.util.ArrayList;


public class SplittingIntervals {

	private static double rec_pos;
	
	public SplittingIntervals(){
		
	}
	
	public static ArrayList<ArrayList<Interval>> split(ArrayList<Interval> input, int recombID){
		
		ArrayList<Interval> l = new ArrayList<Interval>();
		ArrayList<Interval> r = new ArrayList<Interval>();
		
		ArrayList<ArrayList<Interval>> ris = new ArrayList<ArrayList<Interval>>();
		
		double rangeMin = input.get(0).getStart();
		double rangeMax = input.get(input.size()-1).getEnd();
		double x = rangeMin + (rangeMax - rangeMin) * Math.random();
		
		rec_pos=x;
		
		//Store the split point in the node
		GenerateARG.getNodeSet().get(recombID).setSplitPoint(x);
		GenerateARG.getSplitPoints().add(x);
		
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

	public static double getRec_pos() {
		return rec_pos;
	}
	
	
}
