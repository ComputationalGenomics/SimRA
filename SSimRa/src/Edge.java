import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.distribution.PoissonDistribution;


public class Edge {

	private int id_son;
	private int id_fath;
	private int idEDGE;
	private ArrayList<Interval> segments;
	
	//Recombination rate r_l associated to the lineage l
	private double rate;
	//Total length of segments
	private double g;
	//Lenght of the whole segment carried by the edge (holes are included)
	private double length;
	//time t to the next event (length of the edge)
	private double time;
	//density = sum of lenght of solid intervals
	private double density;
	//nMut = #of mutations derived from a Poisson distribution
	private int nMut = 0;
	//for each STR locus we have the value delta 
	private HashMap<Double,Integer> deltaStr;
	//number of mutations/steps for str
	private double num_mut_str;
	private double Nedge;
	
	//
	private boolean underSelection;
	

	public double getNum_mut_str() {
		return num_mut_str;
	}

	public void setNum_mut_str(double num_mut_str) {
		this.num_mut_str = num_mut_str;
	}

	public  HashMap<Double, Integer> getDeltaStr() {
		return deltaStr;
	}

	public  void setDeltaStr(HashMap<Double, Integer> deltaStr) {
		this.deltaStr = deltaStr;
	}

	public Edge(int idEDGE, int id_son, int id_fath, ArrayList<Interval> s, double l) {
		super();
		this.id_son = id_son;
		this.id_fath = id_fath;
		this.idEDGE = idEDGE;
		this.segments = s;
		this.length = l;
		this.underSelection = false;
	}
	
	public Edge(int idEDGE, int id_son, int id_fath, ArrayList<Interval> s, double l, boolean selection, double Nedge) {
		super();
		this.id_son = id_son;
		this.id_fath = id_fath;
		this.idEDGE = idEDGE;
		this.segments = s;
		this.length = l;
		this.underSelection = selection;
		this.Nedge = Nedge;
	}


	public Edge(int idEDGE, int id_son, int id_fath, ArrayList<Interval> s) {
		super();
		this.id_son = id_son;
		this.id_fath = id_fath;
		this.idEDGE = idEDGE;
		this.segments = s;
		this.underSelection=false;
	}
	
	public Edge(int idEDGE, int id_son, int id_fath, ArrayList<Interval> s, boolean selection, double Nedge) {
		super();
		this.id_son = id_son;
		this.id_fath = id_fath;
		this.idEDGE = idEDGE;
		this.segments = s;
		this.underSelection = selection;
		this.Nedge = Nedge;
	}
	
	public Edge() {	
	}
	
	public int getId_son() {
		return id_son;
	}
	public void setId_son(int id_son) {
		this.id_son = id_son;
	}
	public int getId_fath() {
		return id_fath;
	}
	public void setId_fath(int id_fath) {
		this.id_fath = id_fath;
	}
	public ArrayList<Interval> getSegments() {
		return segments;
	}
	public void setSegments(ArrayList<Interval> segments) {
		this.segments = segments;
	}
	public double getRate() {
		return rate;
	}
	public void setRate(double rate) {
		this.rate = rate;
	}
	
	public double getG() {
		return g;
	}
	public void setG(double g) {
		this.g = g;
	}
	public int getIdEDGE() {
		return idEDGE;
	}
	public void setIdEDGE(int idEDGE) {
		this.idEDGE = idEDGE;
	}
	public void updateG(){
		
	}
	public int getnMut() {
		return nMut;
	}

	public void setnMut(int nMut) {
		this.nMut = nMut;
	}

	public double getDensity() {
		return density;
	}

	public void setDensity(double density) {
		this.density = density;
	}

	public double getTime() {
		return time;
	}

	public void setTime() {
		this.time = GenerateARG.getNodeSet().get(this.getId_fath()).getLevel()-GenerateARG.getNodeSet().get(this.getId_son()).getLevel();
	}

	public double getLength() {
		return length;
	}

	public void setLength(double length) {
		this.length = length;
	}
	
	public double getNedge()
	{
		return Nedge;
	}
	
	public static int getBinomial(int n, double p) {
		  int x = 0;
		  for(int i = 0; i < n; i++) {
		    if(Math.random() < p)
		      x++;
		  }
		  return x;
		}
	
	public int computeNumberMutationsStr(double mu_str, int N) {
		   
		   double lamda = mu_str*this.getTime()*GenerateARG.getN();
		   PoissonDistribution Prand = new PoissonDistribution(lamda);
		   int y = Prand.sample();
		   
		   return y;
	}
	/**
	 * Compute the delta value for each STR that is present in one of the segments carried by the edge
	 * Create an HashMap of Delta values where the key is the location/id of the STR
	 * @param N
	 */
	public void computeDeltaStrs(int N) {
		
		//Get the list of STR of the parent node
		ArrayList<Str> strs_parent = GenerateARG.getNodeSet().get(this.getId_fath()).getStrs();
		this.setDeltaStr(new HashMap<Double, Integer>());
		
		
		for(int i = 0; i < strs_parent.size(); i++) {
			
			
			boolean found = false;
			
			for(int j = 0; j < segments.size() && !found; j++) {
				if(strs_parent.get(i).getLocation() >= segments.get(j).getStart() && strs_parent.get(i).getLocation() <= segments.get(j).getEnd()){
					found = true;
				}
			}
			
			if(found == false){
				GenerateARG.getNodeSet().get(this.getId_fath()).getStrs().get(i).setPresence(false);
				
			}
			else {
			
				//1 compute the number of mutations 
				int y = computeNumberMutationsStr(strs_parent.get(i).getMu_rate(), N);
				
				//2 extract number of successes from binomial distribution
				int x = getBinomial(y, 0.5);
				
				//3 compute the delta for the strs
				int delta = 2*x-y;
				
				//4 store the delta in the map carried by the edge with the ID of the edge
				deltaStr.put(strs_parent.get(i).getLocation(),delta);
			}
			
		}
	}
	
	public void computeRate(int N, double recomb, double g){
		
		rate = N*g*recomb*length;
		
		//Hudson version
		//rate=N*g*recomb;
	}
	
	public double computeDensity(){
		double sumD = 0;
		for(int i = 0; i < segments.size(); i++) {
			sumD = sumD + (segments.get(i).getEnd()-segments.get(i).getStart());
		}
		density = sumD;
		return sumD;
	}
	
	public void computeLength(){
		ArrayList<Interval> seg = this.getSegments();
		length = seg.get(seg.size()-1).getEnd() - seg.get(0).getStart();
	}

	
	public boolean isUnderSelection() {
		return underSelection;
	}

	@Override
	public boolean equals(Object obj) {
		// TODO Auto-generated method stub
		//return super.equals(obj);
		
		if (obj == null) return false;
	    if (obj == this) return true;
	    if (!(obj instanceof Edge))return false;
	    Edge otherMyClass = (Edge)obj;
	    
	    if (this.idEDGE==otherMyClass.getIdEDGE())
	    	return true;
	    else
	    	return false;
	    
	}


	
	
	
	
	
}
