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
import java.util.TreeSet;


/**
 * This class creates a file .txt containing the information about the ARG that describes the evolution 
 * of a population
 * 
 * @author Anna Paola Carrieri
 * 
 * @see PopulationARG
 */
public class CreatingFilesForStructure {
	
	
	private static HashMap<Double, ArrayList<Integer>> map_times_nodePop;
	private static TreeSet<Double> ordered_times;
	private static HashMap<Double, GlobalNode> map_global_nodes;
	
	/**
	 * This function returns a map where the key value is the age of a node in a population and the object is list of two elements: the ID of the population of the node and ID of the node
	 * @return the map where the key value is the age of a node in a population and the object is list of two elements: the ID of the population of the node and ID of the node
	 */
	public static HashMap<Double, ArrayList<Integer>> getMap_times_nodePop() {
		return map_times_nodePop;
	}

	/**
	 * This procedure sets the map where the key value is the age of a node in a population and the object is list of two elements: the ID of the population of the node and ID of the node
	 * @param map between the age of a node and the list of two elements: ID of population where the node belongs and ID of the node
	 */
	public static void setMap_times_nodePop(HashMap<Double, ArrayList<Integer>> map) {
		map_times_nodePop = map;
	}
	/**
	 * This function returns an ordered list of times/ages of the nodes in the ARGs of the scaffold
	 * @return an ordered list of times associated with each single event node
	 */
	public static TreeSet<Double> getOrdered_times() {
		return ordered_times;
	}
	
	/**
	 * Set the ordered list of times/ages of the nodes in the ARGs of the scaffold
	 * @param ordered_times an ordered list of times/ages of the nodes in the ARGs of the scaffold
	 */
	public static void setOrdered_times(TreeSet<Double> ordered_times) {
		CreatingFilesForStructure.ordered_times = ordered_times;
	}
	/**
	 * This function returns a map where the key is the age of a single node in the whole scaffold and the object is the corresponding instance of the class GlobalNode containing the information about
	 * the node that is global in the scaffold
	 * @return the map where the key is the age of a single node in the whole scaffold and the object is the corresponding instance of the class GlobalNode containing the information about
	 * the node that is global in scaffold 
	 * @see GlobalNode 
	 */
	public static HashMap<Double, GlobalNode> getMap_global_nodes() {
		return map_global_nodes;
	}
	/**
	 * Set the map between the time/age of the node and its global information, that is the key is the age of a single node in the whole scaffold and the object is the corresponding instance of the class GlobalNode containing the information about
	 * the node that is global in the scaffold
	 * @param map_global_nodes the map between the time/age of the node and its global information 
	 * @see GlobalNode
	 */
	public static void setMap_global_nodes(HashMap<Double, GlobalNode> map_global_nodes) {
		CreatingFilesForStructure.map_global_nodes = map_global_nodes;
	}
	
/**
 * This procedure creates a file .txt containing the information about the ARG that describes the evolution 
 * of a population
 * @see PopulationARG
 * @param pop an instance of the class PopulationARG that stores all the information of the ARG associated with a population
 * @param fileName string representing the output file name 
 */
 public static void createTxtFile_PopulationStructure(String fileName, PopulationARG pop) {
	 	
	 	Writer writer = null;
		File file = new File(""+fileName);
	
			try {
				writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
			    
			    writer.write("####### Information about Population #######\n");
			    writer.write("Population ID :" +pop.getId_pop()+"\n");
			    writer.write("Number of nodes :" +pop.getNodeSet().size()+"\n");
			    writer.write("Number of edges :" +pop.getGraphEdges().size()+"\n");
			    writer.write("Number of activeEdges :" +pop.getActiveEdges().size()+"\n");
			    writer.write("Number of extantUnits :" +pop.getExtantUnits()+"\n");
			    writer.write("Population Kind :" +pop.getKind()+"\n");
			    
			    
			    writer.write("\n\n####### Information about edges and connectivity #######\n");
			    Integer keyE;
				Edge e;
			    Iterator<Integer> it_edges = pop.getGraphEdges().keySet().iterator(); 
				
				while (it_edges.hasNext()) {
			    
					keyE = it_edges.next();
					e = pop.getGraphEdges().get(keyE);
				    	
				    writer.write(e.getId_fath()+"-->"+e.getId_son()+"\t");
				    
				  }
				
			    if(pop.getKind() == 2){
			    	writer.write("Population derives by coalescence of 2 populations:"+pop.getPopanc()[0]+" and "+pop.getPopanc()[1] +"\n");
			    	
			    }
			    if(pop.getKind() == 1){
			    	writer.write("Population derives by recombination of 1 population:"+pop.getPopanc()[0]+"\n");
			    }
			    
			    if(pop.getKind()==1 || pop.getKind()==2){
			    	writer.write("\n\n####### Information about mapping between root nodes and leaves #######\n");
					Iterator<Integer> it_map = pop.getMaps_leaves_roots().keySet().iterator();
					while(it_map.hasNext()){
						Integer key_pop = it_map.next();
						writer.write("------- Mapping with Population "+key_pop+" --------- \n");
						HashMap<Integer, Integer> map_nodes = pop.getMaps_leaves_roots().get(key_pop);
						Iterator<Integer> it_mapping_nodes = map_nodes.keySet().iterator();
						while(it_mapping_nodes.hasNext()) {
							Integer key_node = it_mapping_nodes.next();
							writer.write(key_node+" - "+map_nodes.get(key_node)+"\n");
						}
						writer.write("-------------------------------------------------------- \n");
					}
			    }
			   
				
				writer.write("\n\n####### Information about times in nodes and lengths of edges #######\n");
				it_edges = pop.getGraphEdges().keySet().iterator(); 
				
				while (it_edges.hasNext()) {
					 keyE = it_edges.next();
					 e = pop.getGraphEdges().get(keyE);
					 if(e.getId_fath() != -1){
						 writer.write("------------ EDGE "+ keyE +" ------------\n");
						 writer.write("father "+e.getId_fath()+"\n");
						 writer.write("time of node father "+pop.getNodeSet().get(e.getId_fath()).getLevel()+"\n");
						 writer.write("son "+e.getId_son()+"\n");
						 writer.write("time of node son "+pop.getNodeSet().get(e.getId_son()).getLevel()+"\n");
						 writer.write("lenght of the edge "+e.getId_fath()+" -> "+e.getId_son()+" "+e.getTime()+"\n");
						 writer.write("-------------------------------------------\n");
					}
				 } 
				
			   writer.close();   
			}//fine try 
			
			catch (IOException ex) {
				ex.printStackTrace();
			} 	
		}//fine procedure
/**
 * The procedure creates a .txt file containing detailed information about the SNP mutations and their location in the non-mixing 
 * segments in the extant units of the contemporary populations
 * @param argPathOutput full path of the directory where the output .txt file will be stored
 * @param argFileName_Output name to assign to the input file
 * @param ordered_pops list of IDs of populations ordered by the time they have been created
 * @param pops map containing all the population created
 * @see PopulationARG
 */
 public static void createStxt(String argPathOutput, String argFileName_Output, ArrayList<Integer> ordered_pops, HashMap<Integer, PopulationARG> pops) {
	 
	    Writer writer = null;
	    String filename = ""+argPathOutput+""+argFileName_Output+"_S.txt";
		File file = new File(filename);
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    //First rows
		    
		    writer.write("paramfile: params\n");
		    //for each leaf population
		    for(int p = 0; p < GenerateAdmixturePopulations.getNum_actual_pops(); p++){
		    	
		    	PopulationARG pop = pops.get(ordered_pops.get(p));
		    	writer.write("Pop "+pop.getId_pop()+" "+pop.getExtantUnits()+"\n");
		    	}
		    
		    writer.write("params: length "+GenerateAdmixturePopulations.getG()+" mu "+GenerateAdmixturePopulations.getMu()+"\n");
		    writer.write("L "+GenerateAdmixturePopulations.getG()+"\n");
		    writer.write("  pos \t"+"anc1 anc2 \t"+"freq \t"+"nodes...\n");
		    
		    //create a TreeSet<Integer> containing all the split points of all the populations
		    TreeSet<Double> total_split_points = new TreeSet<Double>();
		    for(int i = 0; i < ordered_pops.size(); i++){
		    	total_split_points.addAll(pops.get(ordered_pops.get(i)).getSplitPoints());
		    }
		    
		   
		    /*//for each actual population
		    for(int p = 0; p < GenerateAdmixturePopulations.getNum_actual_pops(); p++){
		    	
		    	PopulationARG pop = pops.get(ordered_pops.get(p));
		    	
		    	System.out.println("Id actual population: "+pop.getId_pop());
		    	
		    	//Map<Integer,Edge> graphEdges = pop.getGraphEdges();
			    //Map<Integer,Node> nodeSet = pop.getNodeSet();
			    Iterator<Double> iterator = pop.getSplitPoints().iterator();*/
		    	
	    	Iterator<Double> iterator = total_split_points.iterator();
		    double[] splits = new double[total_split_points.size()];
		    int i = 0;
		    while(iterator.hasNext()) {
		    	splits[i] = (double)iterator.next();
		    	i++;
		    }
		   
		    double extSx;
		    double extDx;
		    int numMut;
	    	TreeSet<Double> allMutations = GenerateAdmixturePopulations.getAllMutations();	
	    	Iterator<Double> it;
			Map<Double, Mutation> allMutationsSet = GenerateAdmixturePopulations.getAllmutationSet();
			Mutation mut;
			int id_pop;
			PopulationARG pop;
		    
		    for(int j = 0; j < splits.length-1; j++) {
		    	
		    	extSx = splits[j];
		    	extDx = splits[j+1]-splits[j];
		    	numMut = computeNumberMutationsSegment(splits[j], splits[j+1]);
		    	
		    	writer.write("> ["+extSx+","+extDx+"]  time 1  E[muts] = #mut ("+numMut+")\n");	
		    	it = allMutations.iterator();
		    	Double nextM;
		    	
				while(it.hasNext()) {
					nextM = it.next();
					if(nextM > splits[j] && nextM <= splits[j+1]){
						
						mut = allMutationsSet.get(nextM);
						id_pop = mut.getID_original_pop();
						pop = GenerateAdmixturePopulations.getPopulations().get(id_pop);
						TreeSet<Integer> leaves_set = compute_leaves_have_mutations(mut);
						writer.write("M "+nextM+" \t"+pop.getNodeSet().get(mut.getIDfather()).getGlobal_ID()+" "+pop.getNodeSet().get(mut.getIDson()).getGlobal_ID()+" \t"+leaves_set.size()+"  \t");
						
						
						Iterator<Integer> it_leaves = leaves_set.iterator();
						while(it_leaves.hasNext()) {
							writer.write(""+it_leaves.next()+" ");
						}
						writer.write("\n");
					}
				}	
		    }
		    //}    
		    writer.close();  
		}
		catch (IOException ex) {
			// report
		} 	
	}
 
 /**
  * Given an non-mixing segment this function returns the number of mutations in that non-mixing segments (based on their locations in the chromosome)
  * @param start start position of the non-mixing segment
  * @param end start position of the non-mixing segment
  * @return the number total number of mutations in the non-mixing segment
  */
 public static int computeNumberMutationsSegment(double start, double end){
		
	 	TreeSet<Double> allMutations = GenerateAdmixturePopulations.getAllMutations();
		
		int count = 0;
		
		Iterator<Double> it = allMutations.iterator();
		
		while(it.hasNext()) {
			
			Double nextM = it.next();
			if(nextM > start && nextM <= end){
				count++;
			}
		}
		return count;	
}
 
/**
 * This procedure creates a .txt file that describes the structure of the scaffold and can be given to Cytoscape in order to visualize the detailed structure of scaffold
 * @param argPathOutput full path of the output directory where the output file for Cytoscape will be stored
 * @param argFileName_Output name for the output file
 */
public static void CreatefileForCytoscape(String argPathOutput, String argFileName_Output){
	Writer writer = null;
	File file = new File(""+argPathOutput+""+argFileName_Output+"_Cytoscape.txt");
	try {
	    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
	    //for each population
	    for(int i = 0; i < GenerateAdmixturePopulations.getOrdered_id_populations().size(); i++){
	    	
	    	PopulationARG pop = GenerateAdmixturePopulations.getPopulations().get(GenerateAdmixturePopulations.getOrdered_id_populations().get(i));
	    	//For each edge print the information about the edge
	    	 Map<Integer,Edge> graphEdges = pop.getGraphEdges();
	    	 Iterator<Integer> it_edges = graphEdges.keySet().iterator();
	    	 Integer it;
	    	 Edge e;
	    	 
	    	 while(it_edges.hasNext()){
	    		 it = it_edges.next();
	    		 e = graphEdges.get(it);
	    		 if(e.getId_fath() != -1){
	    			 writer.write(e.getId_fath()+" "+e.printSegments()+" "+e.getId_son()+"\n");
	    		 }		
	    	 }
	    }
	    writer.close();  
	}
	catch (IOException ex) {
		// report
	}
}
 
/**
 * This function returns a ordered list of extant unit IDs that have the SNP mutation given in input
 * @param mut object of the class Mutation containing the information about a single SNP mutation
 * @see Mutation
 * @return a ordered list of extant units IDs that have the mutations in input
 */
 public static TreeSet<Integer> compute_leaves_have_mutations(Mutation mut){
	 TreeSet<Integer> leaves = new TreeSet<Integer>();
	 for(int p = 0; p < GenerateAdmixturePopulations.getNum_actual_pops(); p++){
	    	
	    	PopulationARG pop = GenerateAdmixturePopulations.getPopulations().get(GenerateAdmixturePopulations.getOrdered_id_populations().get(p));
	    	//System.out.println("Id actual population: "+pop.getId_pop());
	    	//for each leaf of the actual population pop
	    	for(int i = 0; i < pop.getExtantUnits(); i++){
	    		//if the leaf has the mutation then 
	    		if(pop.getNodeSet().get(i).getMutation_set().contains(mut.getPositionID()))
	    			leaves.add(pop.getNodeSet().get(i).getGlobal_ID());
	    	}	   
	 }
	 return leaves;
 }
 
 
 /**
  * Based on the age/time of each node, this procedure assigns a global ID to each node in the whole scaffold. More precisely,
  * the nodes (so the IDs) are globally ordered in increasing order of time in which they have been created (based on their age, from the most recent (time 0), 
  * to the oldest (time of the GMRCA) of the all scaffold for instance). 
  * @param populations HashMap containing all the populations (extant and ancestral) that form the scaffold, the key is the ID of the population and the object is the corresponding instance of the class PopulationARG
  * @see PopulationARG
  */
 public static void setGlobalID_nodes(HashMap<Integer, PopulationARG> populations){ 
	 
	 map_times_nodePop = new HashMap<Double, ArrayList<Integer>>();
	 ordered_times = new TreeSet<Double>();
	 Map<Integer, Node> nodes_set;
	 int id_node_global = 0;
	 ArrayList<Integer> ordered_ids_pops = GenerateAdmixturePopulations.getOrdered_id_populations();
	
	 for(int i = ordered_ids_pops.size()-1; i >= 0; i--){
		 PopulationARG pop = populations.get(ordered_ids_pops.get(i));
			
		 //for each node in the population pop
		 nodes_set = pop.getNodeSet();
		 Iterator<Integer> it_nodes = nodes_set.keySet().iterator();
			
		 while(it_nodes.hasNext()){
			
			Node n = pop.getNodeSet().get(it_nodes.next());
		
			//if the current Node n is not a leaf node in the current PopulationARG pop
			if(n.getId() >= pop.getExtantUnits()){
				
				ArrayList<Integer> ids_pop_node = new ArrayList<Integer>();
				ids_pop_node.add(pop.getId_pop());
				ids_pop_node.add(n.getId());
				if(!ordered_times.contains(n.getLevel())) {
					ordered_times.add(n.getLevel());
					map_times_nodePop.put(n.getLevel(), ids_pop_node);
				}
			}
			else {
				//if instead it is a leaf node we have to check that it is not a leaf of a leaf population
				if(n.getLevel() != 0){
					
					//n is a leaf of an ancestral population and has a corresponding root node in the 
					//bottom population
					n.setRoot(true);
					ArrayList<Integer> ids_pop_node = new ArrayList<Integer>();
					ids_pop_node.add(pop.getId_pop());
					ids_pop_node.add(n.getId());
					if(!ordered_times.contains(n.getLevel())) {
						ordered_times.add(n.getLevel());
						map_times_nodePop.put(n.getLevel(), ids_pop_node);
					}
					
					/*HashMap<Integer, HashMap<Integer, Integer>> map = pop.getMaps_leaves_roots();
					Iterator<Integer> it_map = map.keySet().iterator();
					while(it_map.hasNext()){
						
						Integer pop_buttom = it_map.next();
						HashMap<Integer, Integer> nodeDown_nodeUp = map.get(pop_buttom);
						Iterator<Integer> it_mapping_nodes = nodeDown_nodeUp.keySet().iterator();
						while(it_mapping_nodes.hasNext()) {
							Integer node_buttom = it_mapping_nodes.next();
							Integer node_up = nodeDown_nodeUp.get(node_buttom);
						}
					}*/
				} //end if it is not a leaf of a current population
			} //end else the Node n is a leaf of an ancestral population
		}
	} //end for all populations
	 
	 //For each leaf populations set the global IDs only leaves of the population
	 for(int j = 0; j < GenerateAdmixturePopulations.getNum_actual_pops(); j++){
			int id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
			PopulationARG pop = populations.get(id_pop);
			//For each leaf node in pop..
			for(int i = 0; i < pop.getExtantUnits(); i++){
				pop.getNodeSet().get(i).setGlobal_ID(id_node_global);
				id_node_global++;
			}	
	 }
	 
	 
	 //A partire dal primo elemento del TreeSet fino a finire all'ultimo prendo il tempo,
	 //recupero il nodo della corretta popolazione e gli assegno l'indice incrementato di 1 
	 //il primo indice associato al primo tempo nel TreeSet diverso da 0, sara' X+1
	 
	 	Iterator<Double> iterator = ordered_times.iterator();
		//System.out.print("Tree set data - global IDs: ");
	 
		// Displaying the Tree set data
		while (iterator.hasNext()) {
			
			double time = iterator.next();
			//System.out.println(time + " ");
			ArrayList<Integer> id_pop_node = map_times_nodePop.get(time);
			populations.get(id_pop_node.get(0)).getNodeSet().get(id_pop_node.get(1)).setGlobal_ID(id_node_global);
			id_node_global++;	
		}
 }
 
/**
 * This procedure prints local and global IDs  of all nodes in the whole scaffold. It has been used mainly for debugging
 */
public static void printMapNodesTimes(){
	
	System.out.println("********* ORDERED NODES OF WHOLE SCAFFOLD *********\n");
	Iterator<Double> iterator = ordered_times.iterator();
	
	int count = 0;
	
	//PRINT LEAF NODES
	for(int j = 0; j < GenerateAdmixturePopulations.getNum_actual_pops(); j++){
		int id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
		PopulationARG pop = GenerateAdmixturePopulations.getPopulations().get(id_pop);
		
		//For each leaf node in pop..
		for(int i = 0; i < pop.getExtantUnits(); i++){
			System.out.print("#"+count+" Age/Time : "+ pop.getNodeSet().get(i).getLevel() + "\t Pop: "+ pop.getId_pop() + "\t Id Node: "+pop.getNodeSet().get(i).getId()+ "\t Id global: "+ pop.getNodeSet().get(i).getGlobal_ID() + "\n");
			count++;			
		}	
	}
	
	//PRINT ALL THE OTHER NODES IN ORDER
	// Displaying the Tree set data
	while (iterator.hasNext()) {
		
		double time = iterator.next();
		ArrayList<Integer> id_pop_node = map_times_nodePop.get(time);
		System.out.print("#"+count+" Age/Time: "+ time + "\t Pop: "+ GenerateAdmixturePopulations.getPopulations().get(id_pop_node.get(0)).getId_pop() + "\t Id Node: "+ GenerateAdmixturePopulations.getPopulations().get(id_pop_node.get(0)).getNodeSet().get(id_pop_node.get(1)).getId() + "\t Id global: "+GenerateAdmixturePopulations.getPopulations().get(id_pop_node.get(0)).getNodeSet().get(id_pop_node.get(1)).getGlobal_ID() +"\n");
		count++;		
	}
	
	System.out.println("********* END LIST *********\n");
	
}



/**
 * The procedure creates a map called map_global_nodes (static field of this class) that associated to each node of each population in the scaffold a global ID.
 * All nodes of the scaffold have IDs ordered based on increasing orderer of time in which the nodes have been created during the simulation process. 
 */
public static void createMapGlobalNodes(){
	
	map_global_nodes = new HashMap<Double, GlobalNode>();
	
	//For all the population starting from the most ancestral one
	for(int j = GenerateAdmixturePopulations.getOrdered_id_populations().size()-1; j >= 0; j--){
		
		int id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
		PopulationARG pop = GenerateAdmixturePopulations.getPopulations().get(id_pop);
		
		//For each node except for the leaf nodes of the contemporary populations
		
		Iterator<Integer> it_nodes = pop.getNodeSet().keySet().iterator();
		while(it_nodes.hasNext()){
			
			int id_n = it_nodes.next();
			Node n = pop.getNodeSet().get(id_n);
			
			// if n is not a leaf of a contemporary population and it is not a root of a population in the middle 
			if(n.getLevel() != 0 && n.getGlobal_ID() != -1){
			
				GlobalNode gNode = new GlobalNode();
				gNode.setGlobalID(n.getGlobal_ID());
				gNode.setRecomb(n.isRecomb());
				gNode.setTime(n.getLevel());
				
				//find sons and parents ID_global
				int sons[] = new int[2];
				int parents[] = new int[2];
				Iterator<Integer>  it_edges;
				
				//if is not a leaf of the current population
				if(n.getId() >= pop.getExtantUnits()){
					
					gNode.setId_pop(pop.getId_pop());
					
					//if n is coalescent
					if(!n.isRecomb()) {
						
						
						//set global_id_sonSx
						sons[0] = pop.getNodeSet().get(n.getIDsonsx()).getGlobal_ID();
						sons[1] = pop.getNodeSet().get(n.getIDsondx()).getGlobal_ID();
						//find the single parent
						boolean found = false;
						it_edges = pop.getGraphEdges().keySet().iterator();  
						 while (it_edges.hasNext() && !found) {  
							Integer keyE = it_edges.next();
							if(pop.getGraphEdges().get(keyE).getId_son()==n.getId()){
									found = true;
									int father = pop.getGraphEdges().get(keyE).getId_fath();
									if(father != -1)
										parents[0] = pop.getNodeSet().get(father).getGlobal_ID();
									else
										parents[0] = -1;
							}	   
						  }
					} //end if the node is coalescent
					
					else{//if the node is recomb
						
						gNode.setSplitPoint(n.getSplitPoint());
						
						sons[0] = pop.getNodeSet().get(n.getIDsonsx()).getGlobal_ID();
						int temp_parents[] = new int[2];
						//search the two parents and checks who is the left and who is the right
						int found = 0;
						it_edges = pop.getGraphEdges().keySet().iterator(); 
						while (it_edges.hasNext() && found < 2) {  
							Integer keyE = it_edges.next();
						    if(pop.getGraphEdges().get(keyE).getId_son()==n.getId()){
						    		temp_parents[found] = pop.getGraphEdges().get(keyE).getId_fath();
									found++;
						    }   
						 }
						//if both parents are different from -1, understand who is the right and who is the left
						if(temp_parents[0] != -1 && temp_parents[1] != -1){
							//Select the parents of sx and dx of the split point
							if(pop.getNodeSet().get(temp_parents[0]).getSegments().get(0).getStart() < pop.getNodeSet().get(temp_parents[1]).getSegments().get(0).getStart()) {
							    	parents[0] = pop.getNodeSet().get(temp_parents[0]).getGlobal_ID();
							    	parents[1] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
							 }else{
								 	parents[0] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
								 	parents[1] = pop.getNodeSet().get(temp_parents[0]).getGlobal_ID();
							 }
						}
						else if(temp_parents[0] == -1 && temp_parents[1] != -1){
							parents[0] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
							parents[1] = -1;
						}
						else if(temp_parents[1] == -1 && temp_parents[0] != -1){
							parents[0] = pop.getNodeSet().get(temp_parents[0]).getGlobal_ID();
							parents[1] = -1;
						}
						else{
							parents[0] = -1;
							parents[1] = -1;
						}
					}//end else the node is recomb
				
				}//end if node n is not a leaf of the current population
				
				else{ //if the node n is a leaf of the current population
					
					//search for the corresponding root in the population down 
					//find the two sons in the population down
					int pop_node_down[] = find_PopulationNode_buttom(pop,n);
					PopulationARG pop_down = GenerateAdmixturePopulations.getPopulations().get(pop_node_down[0]);
					Node node_down = pop_down.getNodeSet().get(pop_node_down[1]);
					
					gNode.setId_pop(pop_down.getId_pop());
					
					//if n is a coalescent node
					if(!n.isRecomb()){
						
						//find the only parent in the current pop
						it_edges = pop.getGraphEdges().keySet().iterator();
						int father = -1;
						boolean found = false;
						
						while (it_edges.hasNext() && !found) {  
							Integer keyE = it_edges.next();
						    if(pop.getGraphEdges().get(keyE).getId_son()==n.getId()){
						    		father = pop.getGraphEdges().get(keyE).getId_fath();
									found = true;
						    }   
						 }
						if(found)
							parents[0] = pop.getNodeSet().get(father).getGlobal_ID();
						else 
							parents[0] = -1;
						
						//find the two sons in the population down
						sons[0] = pop_down.getNodeSet().get(node_down.getIDsonsx()).getGlobal_ID();
						sons[1] = pop_down.getNodeSet().get(node_down.getIDsondx()).getGlobal_ID();
					}//end if the node is coalescent
					
					else{//the node is recomb
						
						gNode.setSplitPoint(node_down.getSplitPoint());
						//search for the son
						sons[0] = pop_down.getNodeSet().get(node_down.getIDsonsx()).getGlobal_ID();
						
						
						int temp_parents[] = new int[2];
						
						//search for one of the parents in the population down
				
						boolean found = false;
						it_edges = pop_down.getGraphEdges().keySet().iterator(); 
						while (it_edges.hasNext() && !found) {  
							Integer keyE = it_edges.next();
						    if(pop_down.getGraphEdges().get(keyE).getId_son()==n.getId()){
						    		int temp_fath = pop_down.getGraphEdges().get(keyE).getId_fath();
									if(temp_fath != -1){
										found = true;
										temp_parents[0] = temp_fath;
									}
						    }   
						 }
						//if found one parent in the population down then find the other parent in the current population  
						if(found){
							boolean found2 = false;
							it_edges = pop.getGraphEdges().keySet().iterator(); 
							while (it_edges.hasNext() && !found2) {  
								Integer keyE = it_edges.next();
							    if(pop.getGraphEdges().get(keyE).getId_son()==n.getId()){
							    		temp_parents[1] = pop.getGraphEdges().get(keyE).getId_fath();
										found2 = true;
							    }   
							 }
							//if both parents are different from -1, understand who is the right and who is the left
							if(temp_parents[0] != -1 && temp_parents[1] != -1){
								//Select the parents of sx and dx of the split point
								if(pop_down.getNodeSet().get(temp_parents[0]).getSegments().get(0).getStart() < pop.getNodeSet().get(temp_parents[1]).getSegments().get(0).getStart()) {
								    	parents[0] = pop_down.getNodeSet().get(temp_parents[0]).getGlobal_ID();
								    	parents[1] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
								 }else{
									 	parents[0] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
									 	parents[1] = pop_down.getNodeSet().get(temp_parents[0]).getGlobal_ID();
								 }
							}
							else if(temp_parents[0] == -1 && temp_parents[1] != -1){
								parents[0] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
								parents[1] = -1;
							}
							else if(temp_parents[1] == -1 && temp_parents[0] != -1){
								parents[0] = pop_down.getNodeSet().get(temp_parents[0]).getGlobal_ID();
								parents[1] = -1;
							}
							else{
								parents[0] = -1;
								parents[1] = -1;
							}	
						}//end if one parent is in pop down and the other in population current
						else{
							//search for the two parents in the same population (that is the current)
							int gotIt = 0;
							it_edges = pop.getGraphEdges().keySet().iterator(); 
							while (it_edges.hasNext() && gotIt < 2) {  
								Integer keyE = it_edges.next();
							    if(pop.getGraphEdges().get(keyE).getId_son()==n.getId()){
							    		temp_parents[gotIt] = pop.getGraphEdges().get(keyE).getId_fath();
							    		gotIt++;
							    }   
							 }
							//if both parents are different from -1, understand who is the right and who is the left
							if(temp_parents[0] != -1 && temp_parents[1] != -1){
								//Select the parents of sx and dx of the split point
								if(pop.getNodeSet().get(temp_parents[0]).getSegments().get(0).getStart() < pop.getNodeSet().get(temp_parents[1]).getSegments().get(0).getStart()) {
								    	parents[0] = pop.getNodeSet().get(temp_parents[0]).getGlobal_ID();
								    	parents[1] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
								 }else{
									 	parents[0] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
									 	parents[1] = pop.getNodeSet().get(temp_parents[0]).getGlobal_ID();
								 }
							}
							else if(temp_parents[0] == -1 && temp_parents[1] != -1){
								parents[0] = pop.getNodeSet().get(temp_parents[1]).getGlobal_ID();
								parents[1] = -1;
							}
							else if(temp_parents[1] == -1 && temp_parents[0] != -1){
								parents[0] = pop.getNodeSet().get(temp_parents[0]).getGlobal_ID();
								parents[1] = -1;
							}
							else{
								parents[0] = -1;
								parents[1] = -1;
							}
						}//end else both parents are in the current populations
					}//end the node is recomb and is a leaf node 
				}//end the node is a leaf in the current population (that is not contemporary)
				gNode.setParents(parents);
				gNode.setSons(sons);
				map_global_nodes.put(gNode.getTime(), gNode);
			}//end if it is not a leaf of a contemporary populations
		}//end while for each node of the current population
	}//end for all populations starting from the most ancestral one
}//end procedure

/**
 * @param pop an ARG, instance of the class PopulationARG
 * @param n node of the ARG pop, instance of the class Node
 * @see Node
 * @see PopulationARG
 * @return an integer vector of size 2, where the first element is an ID of the population in the scaffold that has the root node corresponding to the leaf node n in the ARG pop. The second
 * element is the ID of the root node corresponding to the node n in pop
 */
public static int[] find_PopulationNode_buttom(PopulationARG pop, Node n){
	
	int ids_pop_node[] = new int[2];
	boolean found = false;
	
	Iterator<Integer> it_map = pop.getMaps_leaves_roots().keySet().iterator();

	while(it_map.hasNext() && !found){
		int pop_buttom = it_map.next();
		HashMap<Integer, Integer> map_nodes = pop.getMaps_leaves_roots().get(pop_buttom);
		Iterator<Integer> it_mapping_nodes = map_nodes.keySet().iterator();
		
		while(it_mapping_nodes.hasNext() && !found) {
			Integer node_buttom = it_mapping_nodes.next();
			Integer node_up = map_nodes.get(node_buttom);
			if(node_up == n.getId()){
				found = true;
				ids_pop_node[0] = pop_buttom;
				ids_pop_node[1] = node_buttom;
			} //end if found
		}//end while I don't find the node
	} //end while I don't find the pop
	return ids_pop_node;
}

/**
 *  The procedures creates a .txt file containing detailed information about the structure of the scaffold generated during the backwards simulation process.
 * @param argPathOutput full path of the directory where the output file will be stored
 * @param argFileName_Output name of the output file 
 * @param events list of events that formed the scaffold
 */
public static void createLtxt(String argPathOutput, String argFileName_Output, ArrayList<Event> events) {
	
	createMapGlobalNodes();
	
	Writer writer = null;
	File file = new File(""+argPathOutput+""+argFileName_Output+"_L.txt");
	try {
	    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
	    
	   
	    TreeSet<Double> times_events = new TreeSet<Double>();
	    
	    for(int i = 0; i < events.size(); i++){
	    	//if the event is not a leaf population
	    	if(events.get(i).getKind() != 0) {
	    		times_events.add(events.get(i).getTime_event());
	    	}
	    }
	    
	    PopulationARG pop;
	    int id_pop;
	    
	    // PRINT THE FIRST ROWS ABOUT LEAF POPULATIONS
	    //for all contemporary populations
		for(int j = 0; j < GenerateAdmixturePopulations.getNum_actual_pops(); j++){
			id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
			writer.write(""+0.000000+"\t"+"create_pop"+"\t"+"pop: "+id_pop+"\t"+"size: 0\n");
		}
		for(int j = 0; j < GenerateAdmixturePopulations.getNum_actual_pops(); j++){
			id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
			pop = GenerateAdmixturePopulations.getPopulations().get(id_pop);
		    writer.write(""+0.000000+"\t"+"change_size"+"\t"+"pop: "+id_pop+"\t"+"size: "+pop.getN()+"\n");
		}
	
		//PRINT ALL LEAF NODES OF ALL CONTEMPORARY POPULATIONS
	    
		//for all contemporary populations print all the leaves that have time 0
		for(int j = 0; j < GenerateAdmixturePopulations.getNum_actual_pops(); j++){
			id_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(j);
			pop = GenerateAdmixturePopulations.getPopulations().get(id_pop);
			
			//For each leaf node in pop..
			for(int i = 0; i < pop.getExtantUnits(); i++){
				Node n = pop.getNodeSet().get(i);
		    	writer.write(""+n.getLevel()+"\t"+"ADD"+"\t"+"node: "+n.getGlobal_ID()+" "+"pop: "+pop.getId_pop()+"\n");
							
			}	
		}
	    Iterator<Double> iterator = ordered_times.iterator();
		
		//PRINT ALL THE OTHER NODES IN ORDER OF TIME CONSIDERING ALSO THE OTHER EVENTS
		// Displaying the Tree set data
		while (iterator.hasNext()) {
			
			double time = iterator.next();
		
			if(!times_events.isEmpty() && time > times_events.first()){
				//write in the file the type of event
				
				//Search the specific event in the arraylist of the events
				for(int i = 0; i < events.size(); i++){
					if(times_events.first() == events.get(i).getTime_event()) {
						
						if(events.get(i).getKind() == 1){
							//merge
							writer.write(""+events.get(i).getTime_event()+"\t H \t merge pop "+events.get(i).getPops_input()[0]+" and pop "+events.get(i).getPops_input()[1]+" in pop "+events.get(i).getPops_output()[0]+"\n");	
						}
						else if(events.get(i).getKind() == 2){
							//split
							writer.write(""+events.get(i).getTime_event()+"\t H \t split pop "+events.get(i).getPops_input()[0]+" in "+events.get(i).getPops_output()[0]+" and pop "+events.get(i).getPops_output()[1]+"\n");
						}
						else{} //leaf population 
						
					}
				}
				//delete the time from times_event
				times_events.pollFirst();
			
			} //end if I have to write the event first
			
			GlobalNode gNode = map_global_nodes.get(time);
			if(gNode.isRecomb()){
				
				writer.write(""+gNode.getTime()+"\t"+"Y"+"\t"+gNode.getSons()[0]+" -> "+gNode.getGlobalID()+" pop: "+gNode.getId_pop()+"\n");
				writer.write(""+gNode.getTime()+"\t"+"R"+"\t"+gNode.getGlobalID()+" -> "+gNode.getParents()[0]+" "+gNode.getParents()[1]+" "+gNode.getId_pop()+" "+gNode.getSplitPoint()+"\n");	
				
			}
			else{ //gNode is coalescent
				writer.write(""+gNode.getTime()+"\t"+"C"+"\t"+gNode.getSons()[0]+" "+gNode.getSons()[1]+" -> "+gNode.getGlobalID()+" pop: "+gNode.getId_pop()+"\n");
			}
		}//end while iterator has next
	    
	  writer.close();  
	}
	catch (IOException ex) {
		// report
		} 	
} 

/*
public static void createSTATSfile(String argPathOutput){
	
		Writer writer = null;
		File file = new File(""+argPathOutput+"output_STATS.txt");
		try {
		   
	
		writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
			    
		//System.out.println("*** PROPERTY (26):  mutation rate ~ (Y/ gNT ) ***");
		int id_last_pop = GenerateAdmixturePopulations.getOrdered_id_populations().get(GenerateAdmixturePopulations.getOrdered_id_populations().size()-1);
		PopulationARG last_pop = GenerateAdmixturePopulations.getPopulations().get(id_last_pop);
		double computed_mu = GenerateAdmixturePopulations.getAllMutations().size()/(last_pop.getG()*last_pop.getN()*last_pop.getCurr_generation());
		double computed_r = GenerateAdmixturePopulations.getRecomb_total_number()/(last_pop.getN()*last_pop.getG()*last_pop.getCurr_generation());
		double scaled_generations = last_pop.getCurr_generation()*last_pop.getN();
		
		// PRINTING FILE STATS
		writer.write("Time of GMRCA \t "+last_pop.getCurr_generation()+"\n");
		writer.write("Time in generations \t "+scaled_generations+"\n");
		writer.write("Number of recombination nodes \t "+GenerateAdmixturePopulations.getRecomb_total_number()+"\n");
		writer.write("Number of coalescent nodes \t "+GenerateAdmixturePopulations.getCoalesc_total_number()+"\n");
		writer.write("Total number of SNPs \t "+GenerateAdmixturePopulations.getAllMutations().size()+"\n");
		writer.write("Recombination Rate \t "+computed_r+"\n");
		writer.write("Mutation Rate \t "+computed_mu+"\n");
		
		writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
}*/

}
