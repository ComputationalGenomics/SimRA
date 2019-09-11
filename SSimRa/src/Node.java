
import java.util.ArrayList;


public class Node {
	
	Node(){
	}
	
	private int id;
	
	//true if node is recombination node
	private boolean recomb;
	private double level;
	private double splitPoint;
    private ArrayList<Interval> segments;
	private ArrayList<Str> strs;
	private int numVisitedIncomingE = 0;
	private int IDfathsx;
	private int IDfathdx;
	private int IDsonsx;
	private int IDsondx;
	
	
	public int getIDfathsx() {
		return IDfathsx;
	}
	public void setIDfathsx(int iDfathsx) {
		IDfathsx = iDfathsx;
	}
	public int getIDfathdx() {
		return IDfathdx;
	}
	public void setIDfathdx(int iDfathdx) {
		IDfathdx = iDfathdx;
	}
	public int getIDsonsx() {
		return IDsonsx;
	}
	public void setIDsonsx(int iDsonsx) {
		IDsonsx = iDsonsx;
	}
	public int getIDsondx() {
		return IDsondx;
	}
	public void setIDsondx(int iDsondx) {
		IDsondx = iDsondx;
	}
	public double getSplitPoint() {
		return splitPoint;
	}
	public void setSplitPoint(double splitPoint) {
		this.splitPoint = splitPoint;
	}
	public int getNumVisitedIncomingE() {
		return numVisitedIncomingE;
	}
	public void setNumVisitedIncomingE(int numVisitedIncomingE) {
		this.numVisitedIncomingE = numVisitedIncomingE;
	}
	public ArrayList<Str> getStrs() {
		return strs;
	}
	public void setStrs(ArrayList<Str> strs) {
		this.strs = strs;
	}
	public ArrayList<Interval> getSegments() {
		return segments;
	}
	public void setSegments(ArrayList<Interval> segments) {
		this.segments = segments;
	}
	public double getLevel() {
		return level;
	}
	public void setLevel(double level) {
		this.level = level;
	}
	public boolean isRecomb() {
		return recomb;
	}
	public void setRecomb(boolean recomb) {
		this.recomb = recomb;
	}
	public Node(int id, boolean recomb, double level) {
		super();
		this.id = id;
		this.recomb = recomb;
		this.level = level;
		this.IDfathdx = -1;
		this.IDfathsx = -1;
		this.IDsonsx = -1;
		this.IDsondx = -1;
	}
	public int getId() {
		return id;
	}
	public void setId(int id) {
		this.id = id;
	}
}

