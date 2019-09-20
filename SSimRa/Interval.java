import java.util.ArrayList;


class Interval {
	 
	 double start;
	 double end;
	 ArrayList<Double> posMutations;
	 Interval() { start = 0; end = 0; posMutations = new ArrayList<Double>(); }
	 Interval(double s, double e) { start = s; end = e; posMutations = new ArrayList<Double>(); }
	 
	public ArrayList<Double> getPosMutations() {
		return posMutations;
	}
	public void setPosMutations(ArrayList<Double> mutations) {
		this.posMutations = mutations;
	}
	public double getStart() {
		return start;
	}
	public void setStart(double start) {
		this.start = start;
	}
	public double getEnd() {
		return end;
	}
	public void setEnd(double end) {
		this.end = end;
	}
	public void printInterval(){
		 System.out.println("["+start+" - "+end+ "]");
	 }
 
}