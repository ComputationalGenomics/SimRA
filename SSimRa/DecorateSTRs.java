import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;


public class DecorateSTRs {
	
	
	
	DecorateSTRs(){}
	
	
	public static void inizializeStrs(int number, int initialState, double rate) {
		
		//Inizialize rates
		double[] rates = new double[number];
		for(int i = 0; i < rates.length; i++) {
			rates[i] = rate;
		}
		
		//Initialize locations
		double[] locations = new double[number];
		for(int i = 0; i < locations.length; i++) {
			locations[i] = Math.random();
		}
		
		//Initialize states to 0
		int[] states = new int[number];
		for(int i = 0; i < states.length; i++) {
			states[i] = initialState;
		}
		
		setLocationsRatesStatesStrs(locations, rates, states);
		
	}
	
	public static void setLocationsRatesStatesStrs(double[] locations, double[] rates, int[] states){
		
		Iterator<Integer> it_nodes = GenerateARG.getNodeSet().keySet().iterator();  
		
		while (it_nodes.hasNext()) {  
			
			Integer keyN = it_nodes.next();
			ArrayList<Str> strs = new ArrayList<Str>();
					
			for(int i = 0; i < locations.length; i++) {
				
				Str str_i = new Str();
				str_i.setMu_rate(rates[i]);
				str_i.setLocation(locations[i]);
				str_i.setNum_repeats(states[i]);
				
				strs.add(str_i);
			}

			GenerateARG.getNodeSet().get(keyN).setStrs(strs);
	    }
		
		//Set the presence of the root to true
		ArrayList<Str> s_root = GenerateARG.getNodeSet().get(GenerateARG.getNodeSet().size()-1).getStrs();
		for(int i = 0; i < s_root.size(); i++) {
			s_root.get(i).setPresence(true);
		}
	}
	
	/**
	 * For each edge of the ARG and for each STR locus, it computes the delta variation of number of repeat between father and son
	 */
	public static void decoratingDeltaStrsEdges() {
		
		Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator();  
		while (it_edges.hasNext()) {  
			Integer keyE = it_edges.next();
			if(GenerateARG.getGraphEdges().get(keyE).getId_fath() != -1) 
				GenerateARG.getGraphEdges().get(keyE).computeDeltaStrs(GenerateARG.getN());	
		}		
	}
	
	public static void printStrs() {
	 
		Iterator<Integer> it_nodes = GenerateARG.getNodeSet().keySet().iterator();  
		
		
		while (it_nodes.hasNext()) {  
			
			Integer keyN = it_nodes.next();
			System.out.println("################# STRs NODE "+ keyN +"########################");
			ArrayList<Str> strs = GenerateARG.getNodeSet().get(keyN).getStrs();
			System.out.println("List solids: ");
			MergeIntervals.printListIntervals(GenerateARG.getNodeSet().get(keyN).getSegments());
			
			for(int i = 0; i < strs.size(); i++) {
				System.out.println("-------------- STR locus "+ strs.get(i).getLocation() +"----------------");
				System.out.println("Present: "+strs.get(i).isPresence());
				System.out.println("Mutation rate:"+strs.get(i).getMu_rate());
				System.out.println("Number repeats:"+strs.get(i).getNum_repeats());
				System.out.println("--------------------------------------------------------");
			}
			
			System.out.println("######################################################");
		}
					
	}
	
	public static void updateStrs(int fath, int son, int num){
		
		Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator();  
		int found = 0;
		
		 //finding the edge/s to get delta 
		 while (it_edges.hasNext() && found != num) {  
				   
				   Integer keyE = it_edges.next();
				   
				   if(GenerateARG.getGraphEdges().get(keyE).getId_fath() != -1){
					   if(GenerateARG.getGraphEdges().get(keyE).getId_fath() == fath && GenerateARG.getGraphEdges().get(keyE).getId_son() == son){
						   found++;
						  
						   //Get the delta of this edge
						   HashMap<Double,Integer> deltaSTRs = GenerateARG.getGraphEdges().get(keyE).getDeltaStr();		   
						   updatingUsingDelta(fath, son, deltaSTRs);
					   }
				   }
		 	}
	}
	
	public static void updateStrs(){
		
		int id_node = GenerateARG.getNodeSet().size()-1;
		
		
		id_node--;
		Iterator<Integer> it_edges;
		
		while(id_node >= 0) {
			
			it_edges = GenerateARG.getGraphEdges().keySet().iterator(); 
			int num_incomingEdges;
			int found = 0;
			
			if(GenerateARG.getNodeSet().get(id_node).isRecomb())
				num_incomingEdges = 2;
			else 
				num_incomingEdges = 1;
			
			while (it_edges.hasNext() && found != num_incomingEdges) {
				Integer keyE = it_edges.next();
				if(GenerateARG.getGraphEdges().get(keyE).getId_son()==id_node) {
					//This is the edge with all the delta values
					found++;
					HashMap<Double,Integer> deltaSTRs = GenerateARG.getGraphEdges().get(keyE).getDeltaStr();		   
					updatingUsingDelta(GenerateARG.getGraphEdges().get(keyE).getId_fath(), id_node, deltaSTRs);
				}	
			}
			id_node--;	
		}	 
	}
	
	
	public static void updatingUsingDelta(int fath, int son, HashMap<Double,Integer> deltaSTRs){
		
	ArrayList<Str> strs_f = GenerateARG.getNodeSet().get(fath).getStrs();
	
	for(int i = 0; i < strs_f.size(); i++) {
		
		//if the father and the son has the same str
		if(deltaSTRs.containsKey(strs_f.get(i).getLocation())) {
			
			//Get delta value
			int d = deltaSTRs.get(strs_f.get(i).getLocation());
			GenerateARG.getNodeSet().get(son).getStrs().get(i).setNum_repeats(strs_f.get(i).getNum_repeats()+d);
			GenerateARG.getNodeSet().get(son).getStrs().get(i).setPresence(true);
			
		}
		else {}
	}
}
	
 public static void printSTRsLeaves(){
	 
	 System.out.println("@@@@@@@@@@@@@@@@ Printing STRs of Leaves @@@@@@@@@@@@@@@@@@@");
	 
	 for(int i = 0; i < GenerateARG.getExtantUnits(); i++) {
		 System.out.println("----- Strs of leave "+i+ " -----");
		 ArrayList<Str> strs = GenerateARG.getNodeSet().get(i).getStrs();
		 for(int j = 0; j < strs.size(); j++) {
				System.out.println("-------------- STR locus"+ strs.get(j).getLocation() +"----------------");
				System.out.println("Present: "+strs.get(j).isPresence());
				System.out.println("Mutation rate:"+strs.get(j).getMu_rate());
				System.out.println("Number repeats:"+strs.get(j).getNum_repeats());
				System.out.println("--------------------------------------------------------");
			}
		 System.out.println("------------------------------");
	 }
	 System.out.println("@@@@@@@@@@@@@@@@ Printing STRs of Leaves @@@@@@@@@@@@@@@@@@@@@");
 }

 
 public static void createSTR_TxtFile(String wholePath, String fileName, Map<Integer,Node> nodeSet) {
		Writer writer = null;
		File file = new File(""+wholePath+fileName);
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    
		    for(int i = 0; i < GenerateARG.getExtantUnits(); i++) {
				
		    	//get STRs of leave i
		    	ArrayList<Str> strs = nodeSet.get(i).getStrs();
		    	//print STRs of leave i in one row
		    	for(int j = 0; j < strs.size(); j++){
		    		writer.write(strs.get(j).getNum_repeats()+"\t");
		    	}
		    	writer.write("\n");
		    }
		    
		   writer.close();   
		}//fine try 
		
		catch (IOException ex) {
		// report
		} 	
	}//fine procedure
 
}
