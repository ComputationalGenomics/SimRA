/*
  Copyright 2015 IBM Corporation


Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under 
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF 
ANY KIND, either express or implied. 
See the License for the specific language governing permissions and limitations under the License.


*/

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

/**
 * This class contains the main procedures to decorate an Ancestral Recombination Graph with STRs
 * 
 * @author Anna Paola Carrieri
 *
 */
public class DecorateSTRs {
	
	
	/**
	 * default constructor
	 */
	DecorateSTRs(){}
	
	/**
	 * The procedure creates the data structures containing the information about STRs
	 * @param number total number of STR loci 
	 * @param initialState initial state of each STR locus
	 * @param rate STR mutation rate 
	 * @param arg ancestral recombination graph to decorate with STRs
	 */
	public static void inizializeStrs(int number, int initialState, double rate, PopulationARG arg) {
		
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
		
		ArrayList<Str> strs = new ArrayList<Str>();
		
		for(int i = 0; i < locations.length; i++) {
			
			Str str_i = new Str();
			str_i.setMu_rate(rates[i]);
			str_i.setLocation(locations[i]);
			str_i.setNum_repeats(states[i]);
			
			strs.add(str_i);
		}
		
		setLocationsRatesStatesStrs(strs, arg);
		
	}
	
	/**
	 * This procedure stores in each node of the ARG the set of initialized STR loci and sets the presence of each
	 * STR locus in the root
	 * @param strs list of STRs
	 * @param arg ancestral recombination graph to be decorated with STRs
	 * @see Str
	 * @see PopulationARG
	 */
	public static void setLocationsRatesStatesStrs(ArrayList<Str> strs, PopulationARG arg){
		
		Iterator<Integer> it_nodes = arg.getNodeSet().keySet().iterator();  
		
		while (it_nodes.hasNext()) {  
			
			Integer keyN = it_nodes.next();
			
			arg.getNodeSet().get(keyN).setStrs(strs);
	    }
		
	
		//Set the presence of the str root to true
		ArrayList<Str> s_root = arg.getNodeSet().get(arg.getNodeSet().size()-1).getStrs();
		for(int i = 0; i < s_root.size(); i++) {
			s_root.get(i).setPresence(true);
		}
	}
	
	/**
	 * For each edge of the ARG and for each STR locus, this procedure computes the delta variation of number of repeat between father and son nodes of the edge
	 * @param arg instance of the class PopulationARG containing all the information about a single population modeled by an ARG
	 */
	public static void decoratingDeltaStrsEdges(PopulationARG arg) {
		
		Iterator<Integer> it_edges = arg.getGraphEdges().keySet().iterator();  
		while (it_edges.hasNext()) {  
			Integer keyE = it_edges.next();
			if(arg.getGraphEdges().get(keyE).getId_fath() != -1) 
				arg.getGraphEdges().get(keyE).computeDeltaStrs(arg.getN(),arg);	
		}		
	}
	
	/**
	 * The procedure prints the information about STR loci for each node of the ARG
	 @param arg instance of the class PopulationARG containing all the information about a single population modeled by an ARG
	 */
	public static void printStrs(PopulationARG arg) {
	 
		Iterator<Integer> it_nodes = arg.getNodeSet().keySet().iterator();  
		
		
		while (it_nodes.hasNext()) {  
			
			Integer keyN = it_nodes.next();
			System.out.println("################# STRs NODE "+ keyN +"########################");
			ArrayList<Str> strs = arg.getNodeSet().get(keyN).getStrs();
			System.out.println("List solids: ");
			MergeIntervals.printListIntervals(arg.getNodeSet().get(keyN).getSegments());
			
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
	
	/**
	 * For each node of the ARG the procedure updates the state of each STR locus based on delta values stored in each edge of the ARG
	 * @param arg instance of the class PopulationARG containing all the information about a single population modeled by an ARG
	 * @see PopulationARG
	 */
	public static void updateStrs(PopulationARG arg){
		
		int id_node = arg.getNodeSet().size()-1;
		
		
		id_node--;
		Iterator<Integer> it_edges;
		
		while(id_node >= 0) {
			
			it_edges = arg.getGraphEdges().keySet().iterator(); 
			int num_incomingEdges;
			int found = 0;
			
			if(arg.getNodeSet().get(id_node).isRecomb())
				num_incomingEdges = 2;
			else 
				num_incomingEdges = 1;
			
			while (it_edges.hasNext() && found != num_incomingEdges) {
				Integer keyE = it_edges.next();
				if(arg.getGraphEdges().get(keyE).getId_son()==id_node) {
					//if the edge has a node father
					if(arg.getGraphEdges().get(keyE).getId_fath() != -1){
						//This is the edge with all the delta values
						found++;
						HashMap<Double,Integer> deltaSTRs = arg.getGraphEdges().get(keyE).getDeltaStr();		   
						updatingUsingDelta(arg.getGraphEdges().get(keyE).getId_fath(), id_node, deltaSTRs, arg);
					}
				}	
			}
			id_node--;	
		}	 
	}
	
	/**
	 * The procedure updates the state of each STR locus in a node, given the status of the father locus and the delta value on the edge
	 * @param fath node father where there are stored the states of the STRs useful for updated the states in the son node
	 * @param son node in which the procedure updates the states of the STRs
	 * @param deltaSTRs list of delta values for each STR locus
	 * @param arg instance of the class PopulationARG containing all the information about a single population modeled by an ARG
	 */
	public static void updatingUsingDelta(int fath, int son, HashMap<Double,Integer> deltaSTRs, PopulationARG arg){
		
	ArrayList<Str> strs_f = arg.getNodeSet().get(fath).getStrs();
	
	for(int i = 0; i < strs_f.size(); i++) {
		
		//if the father and the son has the same str
		if(deltaSTRs.containsKey(strs_f.get(i).getLocation())) {
			
			//Get delta value
			int d = deltaSTRs.get(strs_f.get(i).getLocation());
			arg.getNodeSet().get(son).getStrs().get(i).setNum_repeats(strs_f.get(i).getNum_repeats()+d);
			arg.getNodeSet().get(son).getStrs().get(i).setPresence(true);
			
		}
		else {}
	}
}
/**
 * The procedure prints to video the information about the STRs in the leaves of the ARG
 * @param arg instance of the class PopulationARG containing all the information about a single population modeled by an ARG
 */
 public static void printSTRsLeaves(PopulationARG arg){
	 
	 System.out.println("@@@@@@@@@@@@@@@@ Printing STRs of Leaves @@@@@@@@@@@@@@@@@@@");
	 
	 for(int i = 0; i < arg.getExtantUnits(); i++) {
		 System.out.println("----- Strs of leave "+i+ " -----");
		 ArrayList<Str> strs = arg.getNodeSet().get(i).getStrs();
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

 /**
  * The procedure creates a text file containing the information of all STRs in the leaves of the ARG
  * @param wholePath whole path of the directory where the output .txt file will be stored
  * @param arg instance of the class PopulationARG containing all the information about a single population modeled by an ARG. 
  * The leaves of the ARG contain the information about STRs
  */
 public static void createSTR_TxtFile(String wholePath, PopulationARG arg) {
	 
	 Writer writer = null;
	 File file = new File(""+wholePath);
		
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    Map<Integer,Node> nodeSet = arg.getNodeSet();
		    
		    for(int i = 0; i < arg.getExtantUnits(); i++) {
				
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
 
 /**
  * This procedure updates the STRs of the population on the top and the eventual one or two populations attached to that due to previous
  *  merging or splitting events
  * @param pop_up instance of the class PopulationARG containing all the information about a single population modeled by an ARG
  * @see PopulationARG 
  */
 public static void updateSTRsRoots(PopulationARG pop_up){
	 
	 if(pop_up.getKind()==1 || pop_up.getKind()==2){
		 
		 	//Decorate the 
		 	updateStrs(pop_up);
			
			//Add all mutations to the mutation set
			Iterator<Integer> it_map = pop_up.getMaps_leaves_roots().keySet().iterator();
			//System.out.println("**** UPDATING STRs FROM POPULATION  "+pop_up.getId_pop()+" ****");
			
			while(it_map.hasNext()){
				
				Integer key_pop = it_map.next();
				PopulationARG pop_buttom = GenerateAdmixturePopulations.getPopulations().get(key_pop);
				//System.out.println("------- Mapping with Population "+key_pop+" --------- \n");
				HashMap<Integer, Integer> map_nodes = pop_up.getMaps_leaves_roots().get(key_pop);
				Iterator<Integer> it_mapping_nodes = map_nodes.keySet().iterator();
				
				while(it_mapping_nodes.hasNext()) {
					Integer node_buttom = it_mapping_nodes.next();
					Integer node_up = map_nodes.get(node_buttom);
					//System.out.println("node_up - node_buttom\n");
					//System.out.println(node_up+" - "+node_buttom+"\n");
					//Copy STRs node up in STRs node 
					pop_buttom.getNodeSet().get(node_buttom).setStrs(pop_up.getNodeSet().get(node_up).getStrs());
				}
				//System.out.println("------------------------ End mapping ------------------------------- \n");
			}
		}
 }
 
}
