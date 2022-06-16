
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;

import java.io.IOException;
import java.io.OutputStreamWriter;

import java.io.Writer;

import java.util.Iterator;

import java.util.Map;

import java.util.TreeSet;


public class CreatingFilesForStructure {
	
	public static void createLtxt(String wholePath, String fileName, Map<Integer,Edge> graphEdges, Map<Integer,Node> nodeSet) {
		Writer writer = null;
		File file = new File(""+wholePath+fileName);
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    writer.write(""+0.000000+"\t"+"create_pop"+"\t"+"pop: 1"+"\t"+"size: 0\n");
		    writer.write(""+0.000000+"\t"+"change_size"+"\t"+"pop: 1"+"\t"+"size: "+GenerateARG.getN()+"\n");
		    
		    //Store the leaves 
		    for(int i = 0; i < GenerateARG.getExtantUnits(); i++){
		    	
		    	Node n = GenerateARG.getNodeSet().get(i);
		    	writer.write(""+n.getLevel()+"\t"+"ADD"+"\t"+"node: "+n.getId()+" "+"pop: 1"+"\n");
		    	
		    }
		    
		    //Store all the others nodes
		    for(int j = GenerateARG.getExtantUnits(); j < GenerateARG.getNodeSet().size(); j++){
		    	Node n = GenerateARG.getNodeSet().get(j);
		    	
		    	//If n is a recombination node
		    	if(n.isRecomb()) {
		    		
		    		//I have to scroll all the edges and find the two parents
		    		int found = 0;
		    		int parents[] = new int[2];
		    		int parSx;
		    		int parDx;
		    		int son = -1;
		    		
		    		//Searchng for the two parents
				    Iterator<Integer> it_edges = graphEdges.keySet().iterator();  
				    while (it_edges.hasNext() && found < 2) {  
						
						   Integer keyE = it_edges.next();
						   if(graphEdges.get(keyE).getId_son()==n.getId()){
							   parents[found] = graphEdges.get(keyE).getId_fath();
							   found++;
						   }
						   
				    }
				    
				    //Searching for the son
				    found = 0;
				    it_edges = graphEdges.keySet().iterator();  
				    while (it_edges.hasNext() && found == 0) {  
						
						   Integer keyE = it_edges.next();
						   if(graphEdges.get(keyE).getId_fath()==n.getId()){
							   son = graphEdges.get(keyE).getId_son();
							   found++;
						   }	   
				    }
		    		
				    //Select the parents of sx and dx of the split point
				   
				    if(nodeSet.get(parents[0]).getSegments().get(0).getStart() < nodeSet.get(parents[1]).getSegments().get(0).getStart()) {
				    	parSx = parents[0];
				    	parDx = parents[1];
				    }else{
				    	parSx = parents[1];
				    	parDx = parents[0];
				    }
				    
				    //Print in the file
				    writer.write(""+n.getLevel()+"\t"+"Y"+"\t"+son+" -> "+n.getId()+" pop: 1"+"\n");
				    writer.write(""+n.getLevel()+"\t"+"R"+"\t"+n.getId()+" -> "+parSx+" "+parDx+" 1 "+n.getSplitPoint()+"\n");
		    	}
		    	//if it's coalescent node
		    	else{
		    		writer.write(""+n.getLevel()+"\t"+"C"+"\t"+n.getIDsonsx()+" "+n.getIDsondx()+" -> "+n.getId()+" pop: 1"+"\n");
		    	}
		    }
		    
		  writer.close();  
		}
		catch (IOException ex) {
			// report
			} 	
	}
	
	public static void createStxt(String wholePath, String fileName, Map<Integer,Edge> graphEdges, Map<Integer,Node> nodeSet) {
		Writer writer = null;
		File file = new File(""+wholePath+fileName);
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    //Intestazione
		    writer.write("paramfile: params\n");
		    writer.write("A 1 "+GenerateARG.getExtantUnits()+"\n");
		    writer.write("params: length "+GenerateARG.getG()+" mu "+GenerateARG.getMu()+"\n");
		    writer.write("L "+GenerateARG.getG()+"\n");
		    writer.write("  pos \t"+"anc1 anc2 \t"+"freq \t"+"nodes...\n");
		    
		    Iterator<Double> iterator = GenerateARG.getSplitPoints().iterator();
		    double[] splits = new double[GenerateARG.getSplitPoints().size()];
		    int i = 0;
		    while(iterator.hasNext()) {
		    	splits[i] = (double)iterator.next();
		    	i++;
		    }
		   
		   
		    double extSx;
		    double extDx;
		    for(int j = 0; j < splits.length-1; j++) {
		    	
		    	extSx = splits[j];
		    	extDx = splits[j+1]-splits[j];
		    	
		    	int numMut = computeNumerMutationsSegment(splits[j], splits[j+1]);
		    	
		    	writer.write("> ["+extSx+","+extDx+"]  time 1  E[muts] = 12345 ("+numMut+")\n");	
		    	
		    	//Each row for each mutations
		    	Map<Double,Mutation> mutSet = GenerateARG.getMutationSet();
				
				Iterator<Double> it = mutSet.keySet().iterator();
				
				
				while(it.hasNext()) {
					Double nextM = it.next();
					if(nextM > splits[j] && nextM <= splits[j+1]){
						
						writer.write("M "+nextM+" \t"+mutSet.get(nextM).getIDfather()+" "+mutSet.get(nextM).getIDson()+" \t"+mutSet.get(nextM).getLeaves().size()+"  \t");
						TreeSet<Integer> leaves = mutSet.get(nextM).getLeaves();
						Iterator<Integer> it_leave = leaves.iterator();
						while(it_leave.hasNext()) {
							writer.write(""+it_leave.next()+" ");
						}
						writer.write("\n");
					}
				}	
		    }
		    
		    
		    writer.close();  
		}
		catch (IOException ex) {
			// report
			} 	
	}
	
public static int computeNumerMutationsSegment(double start, double end){
		
		Map<Double,Mutation> mutSet = GenerateARG.getMutationSet();
		int count = 0;
		
		Iterator<Double> it = mutSet.keySet().iterator();
		
		while(it.hasNext()) {
			
			Double nextM = it.next();
			if(nextM > start && nextM <= end){
				count++;
			}
		}
		
		return count;
		
	}
	
	 public static void createTxtFile_timesNodeConnectivity(String fileName, DiGraph graph) {
			Writer writer = null;

			try {
			    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), "utf-8"));
			    
			    writer.write("####### Connectivity of nodes in the tree #######\n");
			    int i, j;
			    writer.write("Number of vertices :" +graph.getAdjMatrix().length+"\n");
				
				for (i = 0; i < graph.getAdjMatrix().length; i++) {
					  
				    for (j = 0; j < graph.getAdjMatrix()[0].length; j++) {
				    	if(graph.getAdjMatrix()[i][j]==1)
				    		 writer.write(i+"-->"+j+"\t");
				    		 //writer.write("\n");
				    	if(graph.getAdjMatrix()[i][j]==2) {
				    		writer.write(i+"-->"+j+"\t");
				    		//writer.write("\n");
				    		writer.write(i+"-->"+j+"\t");
				    		
				    	}
				    }
				    
				  }
				  
				 writer.write("\n####### Information about times in nodes and lengths of edges #######\n");
				 Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator(); 
				 while (it_edges.hasNext()) {
					 Integer keyE = it_edges.next();
					 Edge e = GenerateARG.getGraphEdges().get(keyE);
					 if(e.getId_fath() != -1){
						 writer.write("------------ EDGE "+ keyE +" ------------\n");
						 writer.write("father "+e.getId_fath()+"\n");
						 writer.write("time of node father "+GenerateARG.getNodeSet().get(e.getId_fath()).getLevel()+"\n");
						 writer.write("son "+e.getId_son()+"\n");
						 writer.write("time of node son "+GenerateARG.getNodeSet().get(e.getId_son()).getLevel()+"\n");
						 writer.write("lenght of the edge "+e.getId_fath()+" -> "+e.getId_son()+" "+e.getTime()+"\n");
						 writer.write("-------------------------------------------\n");
					}
				 }
				 writer.write("\n##########################################################\n");
				
			   writer.close();   
			}//fine try 
			
			catch (IOException ex) {
			// report
			} 	
		}//fine procedure
	
}