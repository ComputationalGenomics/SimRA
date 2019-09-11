import java.util.TreeSet;


public class Mutation {
	
	private double positionID;
	private int IDson;
	private int IDfather;
	private int IDedge;
	private TreeSet<Integer> leaves;
	private Interval segment;
	public double getPositionID() {
		return positionID;
	}
	public void setPositionID(double positionID) {
		this.positionID = positionID;
	}
	public int getIDson() {
		return IDson;
	}
	public void setIDson(int iDson) {
		IDson = iDson;
	}
	public int getIDfather() {
		return IDfather;
	}
	public void setIDfather(int iDfather) {
		IDfather = iDfather;
	}
	public int getIDedge() {
		return IDedge;
	}
	public void setIDedge(int iDedge) {
		IDedge = iDedge;
	}
	public TreeSet<Integer> getLeaves() {
		return leaves;
	}
	public void setLeaves(TreeSet<Integer> leaves) {
		this.leaves = leaves;
	}
	public Interval getSegment() {
		return segment;
	}
	public void setSegment(Interval segment) {
		this.segment = segment;
	}
	
	
}
