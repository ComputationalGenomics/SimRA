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
import java.util.LinkedList;
import java.util.Queue;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;


/**
 * DecorateMutationsSNP decorates with SNPs mutations each edge of the ARG with mutations and that are transmitted to the extant units (leave).
 * The class can decorate both a simple single population and admixture populations. Also it has a procedure that creates a .txt file containing the information about the SNP mutations in the extant units
 * of an actual populations. In case of admixture populations it there will be txt file for each actual population considered. 
 * 
 * @author Anna Paola Carrieri
 */

public class DecorateMutationsSNP {
	
	/**
	 * The procedure computes for each edge of the ARG the number of mutations occurred in that period of time (length of the edge) by means of Poisson and Normal distributions
	 * @param mu real value that is SNP mutation rate
	 * @param arg instance of the class PopulationARG representing the ARG of a single population
	 * @throws ErrorTimesInputFileException this exception is thrown if there is an inconsistency in time parameters setting
	 * @see PopulationARG
	 */
	public static void decoratingWithMutations(double mu, PopulationARG arg) throws ErrorTimesInputFileException{
		
		//For each leaf of arg create a set of mutations in the node
		for(int j = 0; j < arg.getExtantUnits(); j++){
			arg.getNodeSet().get(j).setLeaf(true);
			arg.getNodeSet().get(j).setMutation_set(new TreeSet<Double>());
		}
		
		
		Iterator<Integer> it_edges = arg.getGraphEdges().keySet().iterator();  
		double lamda = 0;
		
		//System.out.println("------- SNPs Decoration of Population "+arg.getId_pop()+" ------------ ");
		
		//For each edge in the ARG, compute the number of mutations and decorate the solids
	    while (it_edges.hasNext()) {  
			   
			   Integer keyE = it_edges.next();
			   
			   if(arg.getGraphEdges().get(keyE).getId_fath() != -1){
				   
				   
				   double density = arg.getGraphEdges().get(keyE).computeDensity();
				   double time = arg.getGraphEdges().get(keyE).getTime();
				   
				   //*** if time is negative then throw the exception that 
				   //there is probably an error in how the time parameters are specified in the input file ***//
				   if(time <= 0){
					   throw new ErrorTimesInputFileException();
				   }
				   
				   //System.out.println("Time: "+arg.getGraphEdges().get(keyE).getTime());
				   double n = arg.getG()*density;
				   double p = mu*arg.getN()*time;
				   
				   //******* Compute the # of mutations by Poisson Distribution ******
				   lamda = n*p;
				  
				   PoissonDistribution Prand = new PoissonDistribution(lamda);
				   int randomPois = Prand.sample();
				 
				   arg.getGraphEdges().get(keyE).setnMut(randomPois);
				   
				   if(randomPois == 0)
					   GenerateAdmixturePopulations.setNum_edges_without_mut(GenerateAdmixturePopulations.getNum_edges_without_mut()+1);
				   else
					   GenerateAdmixturePopulations.setNum_edges_more_than_one_mut(GenerateAdmixturePopulations.getNum_edges_more_than_one_mut()+1);
				   
				   //For each mutation we assign a position
				   for(int i = 0; i < arg.getGraphEdges().get(keyE).getnMut(); i++) {
					   
					  //maps the real solids in (0,1) in an interval (0,densityTotal)
					  //by creating a new list of segments.
					   
					  double start = 0;
					  
					  ArrayList<Interval> tempSegs = new ArrayList<Interval>();
					  
					  for(int j = 0; j < arg.getGraphEdges().get(keyE).getSegments().size(); j++) {  
						  
						  double lengthInt =  arg.getGraphEdges().get(keyE).getSegments().get(j).getEnd()-arg.getGraphEdges().get(keyE).getSegments().get(j).getStart();
						  tempSegs.add(new Interval(start, start+lengthInt));
						  start = start+lengthInt;
					  }
					  
					  double x = Math.random()*arg.getGraphEdges().get(keyE).getDensity();
					  boolean found = false;
					  int index = 0;
					  
					  while(!found && index < tempSegs.size()){
						  if(x >= tempSegs.get(index).getStart() && x < tempSegs.get(index).getEnd()){
							  found = true;
						  }
						  else{
							 index++;
						  }
					  }
					  
					  double offset = x-tempSegs.get(index).getStart();
					  double posMutation =  arg.getGraphEdges().get(keyE).getSegments().get(index).getStart()+offset;
					  arg.getGraphEdges().get(keyE).getSegments().get(index).getPosMutations().add(posMutation);
					  
					  
					  //add the mutation to the set of the all mutations
					  GenerateAdmixturePopulations.getAllMutations().add(new Double(posMutation));
					  
					  //Create the mutation object
					  Mutation mut = new Mutation();
					  mut.setPositionID(posMutation);
					  mut.setIDfather(arg.getGraphEdges().get(keyE).getId_fath());
					  mut.setIDson(arg.getGraphEdges().get(keyE).getId_son());
					  mut.setSegment(arg.getGraphEdges().get(keyE).getSegments().get(index));
					  mut.setLeaves(computeLeaves(arg.getGraphEdges().get(keyE).getId_son(), arg));
					  mut.setID_original_pop(arg.getId_pop());
					  
					  arg.getSNPpositionsList().add(posMutation);
					  arg.getMutationSet().put(posMutation, mut);
					  GenerateAdmixturePopulations.getAllmutationSet().put(mut.getPositionID(), mut);
					  
					  //For each leaf (in the set), that has the mutation mut, add mut to the leaf 
					  TreeSet<Integer> leaves = mut.getLeaves();
					  Iterator<Integer> it_leaves = leaves.iterator();
					  
					  while(it_leaves.hasNext()){
						int l = it_leaves.next();  
						arg.getNodeSet().get(l).getMutation_set().add(mut.getPositionID());
					  }  
				   } //end of computing the exact position of each mutation 
			   } //end if the edge has a father	
			   
	       } //end while for each edge of the ARG	
	    //System.out.println("------- End of SNPs Decoration of Population "+arg.getId_pop()+" ------------ ");
	}//end of the procedure
		
	
	/**
	 * This function computes and returns the total number of mutations occurred in the generated Ancestral Recombination Graph
	 * @param arg instance of the class PopulationARG that represents the ARG of a single population
	 * @see PopulationARG
	 * @return sum_nMut integer representing the total number of mutations in the ARG
	 */
	public static int computeTotalMutationsNumber(PopulationARG arg){
		
		int sum_nMut = 0;
		Iterator<Integer> it_edges = arg.getGraphEdges().keySet().iterator();  
	    Edge e;
	    
	    while (it_edges.hasNext()) {  
			
			   Integer keyE = it_edges.next();
			   e = arg.getGraphEdges().get((int)keyE);
			   
			   if(e.getId_fath() != -1){
				   sum_nMut = sum_nMut+e.getnMut();
			   }
	    }
		return sum_nMut;
	}
	
	/**
	 * The procedure creates a .txt file containing the information about the SNP mutations in the extant units
	 * The first row represents the position of each mutation that is a double value in the interval (0,1)
	 * The other rows represent the leave nodes as a matrix of 0 and 1 indicating the absence/presence of the corresponding mutation
	 * @param wholePath path of the directory where the output file is stored
	 * @param arg instance of the class PopulationARG that represents the ARG of a single population (in this case and actual population)
	 * @see PopulationARG
	 */
	public static void createSNP_TxtFile(String wholePath, PopulationARG arg) {
		Writer writer = null;
		File file = new File(""+wholePath);
		
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    TreeSet<Double> posMutations = GenerateAdmixturePopulations.getAllMutations();
		    //System.out.println("Creating file SNPs for population "+arg.getId_pop());
		    //System.out.println("# of all mutations in all populations "+posMutations.size());
		    
		    writer.write("Population "+arg.getId_pop()+"\n");
		    writer.write("Number of extant units (rows): "+arg.getExtantUnits()+"\n");
		    writer.write("Number of SNPs for each extant unit (columns): "+posMutations.size()+"\n\n");
		    
		    Iterator<Double> it_posMuts = posMutations.iterator();
		    Double position;
		    while(it_posMuts.hasNext()){
		    	position = it_posMuts.next();
		    	String troncato = String.format ("%.4f", position);
		    	writer.write(troncato+" ");
		    }  
		    writer.write("\n\n");
		    
		    //For each leave print a row representing the absence or presence of SNP mutations
		    for(int i = 0; i < arg.getExtantUnits(); i++) {
		    	
		    	it_posMuts = posMutations.iterator();
		    	
			    while(it_posMuts.hasNext()) {
			    	
			    	position = it_posMuts.next();
			    	
			    	//check if the arg has the mutation
			    	if(arg.getMutationSet().containsKey(position)){
			    		
			    		//if yes then check if the leaf has the mutation
			    		//Mutation mut = arg.getMutationSet().get(position);
			    		TreeSet<Double> muts_set = arg.getNodeSet().get(i).getMutation_set();
				    	
			    		if(muts_set.contains(position))
				    		writer.write("1 ");
				    	else
				    		writer.write("0 ");
			    	}
			    	else
			    		writer.write("0 ");
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
	 * The procedure updates the SNPs of all populations that are directly connected to the ancestral population in input.
	 * There can be at most two populations directly connected to an ancestral one. The procedure updates their SNPs based on the SNPs of the ancestral population that are transmitted 
	 * @param pop_up instance of the class PopulationARG that represents the ARG of a single population (in this case an ancestral population)
	 * @see PopulationARG
	 */
	public static void updateMutationsForTheBottomPop(PopulationARG pop_up){
		
		if(pop_up.getKind()==1 || pop_up.getKind()==2){
			
			//Add all mutations to the mutation set
			Iterator<Integer> it_map = pop_up.getMaps_leaves_roots().keySet().iterator();
			//System.out.println("**** UPDATING MUTATIONS FROM POPULATION  "+pop_up.getId_pop()+" ****");
			
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
					TreeSet<Integer> leaves_from_node_buttom = computeLeaves(node_buttom, pop_buttom);
					
					//for each mutation in node_up
					TreeSet<Double> muts = pop_up.getNodeSet().get(node_up).getMutation_set();
					Iterator<Double> it_muts = muts.iterator();
					
					while(it_muts.hasNext()){
						//Prendo la mutazione
						Double mut_id = it_muts.next();
						Mutation mut_up = pop_up.getMutationSet().get(mut_id);
						Mutation mut = new Mutation();
						mut.setPositionID(mut_up.getPositionID());
						mut.setSegment(mut_up.getSegment());
						mut.setLeaves(leaves_from_node_buttom);
						mut.setIDfather(mut_up.getIDfather());
						mut.setIDson(mut_up.getIDson());
						mut.setID_original_pop(mut_up.getID_original_pop());
						  
						pop_buttom.getSNPpositionsList().add(mut.getPositionID());
						pop_buttom.getMutationSet().put(mut.getPositionID(), mut);
						  
						//For each leaf (in the treeset) that has the mutation mut, add mut to the leaf 
						TreeSet<Integer> leaves = mut.getLeaves();
						Iterator<Integer> it_leaves = leaves.iterator();
						  
						  while(it_leaves.hasNext()){
							int l = it_leaves.next();  
							pop_buttom.getNodeSet().get(l).getMutation_set().add(mut.getPositionID());
						  }  
					}
				}
				//System.out.println("------------------------ End mapping ------------------------------- \n");
			}
		}	
	}
	
	/**
	 * This function returns the set of leaves in the ARG that are reachable from a node give in input 
	 * @param root integer value representing the ID of the node in the arg
	 * @param arg  object of the class PopulationARG
	 * @see PopulationARG
	 * @see Node
	 * @return the set of integer values representing the leaves in the ARG that are reachable from a node give in input (root)
	 */
	public static TreeSet<Integer> computeLeaves(int root, PopulationARG arg) {
		
		TreeSet<Integer> leaves = new TreeSet<Integer>();
		Queue<Integer> q = new LinkedList<Integer>();
		boolean visited[] = new boolean[arg.getNodeSet().size()];
		/*System.out.println("ID Population = "+arg.getId_pop());
		System.out.println("Population Kind = "+arg.getKind());
		System.out.println("Extant units = "+arg.getExtantUnits());
		System.out.println("Node set size = "+arg.getNodeSet().size());*/
		
		for(int i = 0; i < visited.length; i++) {
			visited[i] = false;
		}
		
		//System.out.println("Root = "+root);
		//If root is a leaf I have to add just that leaf in leaves, no reason to do the visit
		if(root >= 0 && root < arg.getExtantUnits()) {
			leaves.add(new Integer(root));
		}
		//Else bfv to search for all the leaves of the subtree with root root
		else { 
			visited[root] = true;
			q.add(new Integer(root));
			
			while (!q.isEmpty())
			{ 
				// remove a labeled vertex from the queue
				int w = ((Integer) q.remove()).intValue(); 
				
				// mark unreached vertices adjacent from w
				Iterator<Integer> it_edges = arg.getGraphEdges().keySet().iterator();  
				while (it_edges.hasNext()) {  
					Integer keyE = it_edges.next();
					Edge e = arg.getGraphEdges().get(keyE);
					//if I have found an edge with w has father
					if (e.getId_fath()==w) { 
						int u = e.getId_son();
						//System.out.println("Visiting edge "+w+"-"+u);
						if(!visited[u]){
				          q.add(new Integer(u));
				          visited[u] = true;
				          //System.out.println("Visiting node :"+u);
				          if(u >= 0 && u < arg.getExtantUnits()) {
				        	  //System.out.println("u is a leaf :"+u);
				        	  
				        	  //u is a leaf node 
				        	  leaves.add(u);
				          }
					 	}
				  	  }	//end if w is a father of the edge
					}//for all edge..select the one that has w has father and take u as son
				}//end while  
		} // end else root is not a leaf
		return leaves;
	}
	
} //end class

