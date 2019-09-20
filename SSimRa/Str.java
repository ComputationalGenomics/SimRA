
public class Str {

	private int num_repeats;
	private double location;
	private boolean presence;
	private double mu_rate;
	
	public Str() {
	
	}
	
	public Str(int num_repeats, double location, boolean presence,
			double mu_rate) {
		super();
		this.num_repeats = num_repeats;
		this.location = location;
		this.presence = presence;
		this.mu_rate = mu_rate;
	}
	
	public int getNum_repeats() {
		return num_repeats;
	}
	public void setNum_repeats(int num_repeats) {
		this.num_repeats = num_repeats;
	}
	public double getLocation() {
		return location;
	}
	public void setLocation(double location) {
		this.location = location;
	}
	public boolean isPresence() {
		return presence;
	}
	public void setPresence(boolean presence) {
		this.presence = presence;
	}
	public double getMu_rate() {
		return mu_rate;
	}
	public void setMu_rate(double mu_rate) {
		this.mu_rate = mu_rate;
	}
	
}
