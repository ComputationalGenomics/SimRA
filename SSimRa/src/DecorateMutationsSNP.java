import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;


/**
 * DecorateMutationsSNP decorates with SNPs mutations each edge of the ARG with mutations and then each extant units (leave) 
 * @author AnnaPaola Carrieri
 * @version July 16 2014
 */

public class DecorateMutationsSNP {
	
	public static void decoratingWithMapMutations(double mu, DiGraph graph, List<Edge> MutMap){
		Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator();  
		double lamda = 0;
		
		//For each edge in the ARG, compute the number of mutations and decorate the solids
	    while (it_edges.hasNext()) {  
			   
			   Integer keyE = it_edges.next();
			   
			   if(GenerateARG.getGraphEdges().get(keyE).getId_fath() != -1){
				   int flag = 0;
				   double density = GenerateARG.getGraphEdges().get(keyE).computeDensity();
				   double time = GenerateARG.getGraphEdges().get(keyE).getTime();
				   //System.out.println("Density: "+density);
				    //System.out.println("Time: "+time);
				   double n = GenerateARG.getG()*density;
				   double p = mu*time; //* GenerateARG.getGraphEdges().get(keyE).getNedge();//GenerateARG.getN()
				  // System.out.println("g: "+ GenerateARG.getG()+" density: "+density+" mu: "+mu+" N: "+GenerateARG.getGraphEdges().get(keyE).getNedge()+" time: "+time);
				   //******* Compute the # of mutations by Poisson Distribution ******
				   lamda = n*p;
				  // System.out.println("Lamda: "+lamda);
				   PoissonDistribution Prand = new PoissonDistribution(lamda);
				   int randomPois = Prand.sample();
				   for (int ij = 0; ij < MutMap.size(); ij++) {
					   if (MutMap.get(ij).getIdEDGE() == keyE) {
						   flag = 1; 
						   break;
					   }
				   }
				   if(flag == 1 && randomPois == 0)
					   randomPois = 1;
				   
				  // System.out.println("Random integer from Poisson: "+randomPois);
				   GenerateARG.getGraphEdges().get(keyE).setnMut(randomPois);
				   //  System.out.println("Number of mutations: "+ GenerateARG.getGraphEdges().get(keyE).getnMut());
				  
				  /* System.out.println("List of segments carried in the edge "+GenerateARG.getGraphEdges().get(keyE).getIdEDGE());
				   MergeIntervals.printListIntervals(GenerateARG.getGraphEdges().get(keyE).getSegments());
				   System.out.println();*/
				  
				   //For each mutation we assign a position
				   for(int i = 0; i < GenerateARG.getGraphEdges().get(keyE).getnMut(); i++) {
					   
					  //maps the real solids in (0,1) in an interval (0,densityTotal)
					  //by creating a new list of segments.
					   
					  double start = 0;
					  
					  ArrayList<Interval> tempSegs = new ArrayList<Interval>();
					  
					  for(int j = 0; j < GenerateARG.getGraphEdges().get(keyE).getSegments().size(); j++) {  
						  
						  double lengthInt =  GenerateARG.getGraphEdges().get(keyE).getSegments().get(j).getEnd()-GenerateARG.getGraphEdges().get(keyE).getSegments().get(j).getStart();
						  tempSegs.add(new Interval(start, start+lengthInt));
						  start = start+lengthInt;
					  }
					  
					  
					  
					 /*System.out.println("List of temporary segments:");
					 MergeIntervals.printListIntervals(tempSegs);
					 System.out.println();*/
					  
					  //generate random number
					  double x = Math.random()*GenerateARG.getGraphEdges().get(keyE).getDensity();
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
					  double posMutation =  GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getStart()+offset;
					  GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getPosMutations().add(posMutation);
					  //System.out.println("Solid "+index+" ["+ GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getStart()+","+GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getEnd()+"]");
					  System.out.println("Position in solid "+index+" : "+posMutation);
					  
					  //Create the mutation object
					  Mutation mut = new Mutation();
					  mut.setPositionID(posMutation);
					  mut.setIDfather(GenerateARG.getGraphEdges().get(keyE).getId_fath());
					  mut.setIDson(GenerateARG.getGraphEdges().get(keyE).getId_son());
					  mut.setSegment(GenerateARG.getGraphEdges().get(keyE).getSegments().get(index));
					  mut.setLeaves(graph.computeLeaves(GenerateARG.getGraphEdges().get(keyE).getId_son()));
					  
					  GenerateARG.getSNPpositionsList().add(posMutation);
					  GenerateARG.getMutationSet().put(posMutation, mut);
					  
				   } //end of computing the exact position of each mutation 
			   } //end if the edge has a father	   
	       } //end while for each edge of the ARG	
	    //System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%% END DECORATING WITH MUTAIONS %%%%%%%%%%%%%%%%%%%%%%%%%%");
	}//end of the procedure
		  
	   	/**
	 * For each edge of the ARG computes the number of mutations occurred in that arc of time by means of Poisson and Normal distributions
	 * @param mu = SNP mutation rate
	 */
	public static void decoratingWithMutations(double mu, DiGraph graph){
		
		Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator();  
		double lamda = 0;
		
		//For each edge in the ARG, compute the number of mutations and decorate the solids
	    while (it_edges.hasNext()) {  
			   
			   Integer keyE = it_edges.next();
			   
			  // System.out.println("*** EDGE "+GenerateARG.getGraphEdges().get(keyE).getIdEDGE()+ " ***");
			   
			   if(GenerateARG.getGraphEdges().get(keyE).getId_fath() != -1){
				   
				   double density = GenerateARG.getGraphEdges().get(keyE).computeDensity();
				   double time = GenerateARG.getGraphEdges().get(keyE).getTime();
				  // System.out.println("Density: "+density);
				  //  System.out.println("Time: "+time);
				   double n = GenerateARG.getG()*density;
				   double p = mu*time; //* GenerateARG.getGraphEdges().get(keyE).getNedge();//GenerateARG.getN()
				   //System.out.println("g: "+ GenerateARG.getG()+" density: "+density+" mu: "+mu+" N: "+GenerateARG.getGraphEdges().get(keyE).getNedge()+" time: "+time);
				   //******* Compute the # of mutations by Poisson Distribution ******
				   lamda = n*p;
				  // System.out.println("Lamda: "+lamda);
				   PoissonDistribution Prand = new PoissonDistribution(lamda);
				   int randomPois = Prand.sample();
				  // System.out.println("Random integer from Poisson: "+randomPois);
				   GenerateARG.getGraphEdges().get(keyE).setnMut(randomPois);
				   //  System.out.println("Number of mutations: "+ GenerateARG.getGraphEdges().get(keyE).getnMut());
				  
				  /* System.out.println("List of segments carried in the edge "+GenerateARG.getGraphEdges().get(keyE).getIdEDGE());
				   MergeIntervals.printListIntervals(GenerateARG.getGraphEdges().get(keyE).getSegments());
				   System.out.println();*/
				  
				   //For each mutation we assign a position
				   for(int i = 0; i < GenerateARG.getGraphEdges().get(keyE).getnMut(); i++) {
					   
					  //maps the real solids in (0,1) in an interval (0,densityTotal)
					  //by creating a new list of segments.
					   
					  double start = 0;
					  
					  ArrayList<Interval> tempSegs = new ArrayList<Interval>();
					  
					  for(int j = 0; j < GenerateARG.getGraphEdges().get(keyE).getSegments().size(); j++) {  
						  
						  double lengthInt =  GenerateARG.getGraphEdges().get(keyE).getSegments().get(j).getEnd()-GenerateARG.getGraphEdges().get(keyE).getSegments().get(j).getStart();
						  tempSegs.add(new Interval(start, start+lengthInt));
						  start = start+lengthInt;
					  }
					  
					  
					  
					 /*System.out.println("List of temporary segments:");
					 MergeIntervals.printListIntervals(tempSegs);
					 System.out.println();*/
					  
					  //generate random number
					  double x = Math.random()*GenerateARG.getGraphEdges().get(keyE).getDensity();
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
					  double posMutation =  GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getStart()+offset;
					  GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getPosMutations().add(posMutation);
					  //System.out.println("Solid "+index+" ["+ GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getStart()+","+GenerateARG.getGraphEdges().get(keyE).getSegments().get(index).getEnd()+"]");
					  //System.out.println("Position in solid "+index+" : "+posMutation);
					  
					  //Create the mutation object
					  Mutation mut = new Mutation();
					  mut.setPositionID(posMutation);
					  mut.setIDfather(GenerateARG.getGraphEdges().get(keyE).getId_fath());
					  mut.setIDson(GenerateARG.getGraphEdges().get(keyE).getId_son());
					  mut.setSegment(GenerateARG.getGraphEdges().get(keyE).getSegments().get(index));
					  mut.setLeaves(graph.computeLeaves(GenerateARG.getGraphEdges().get(keyE).getId_son()));
					  
					  GenerateARG.getSNPpositionsList().add(posMutation);
					  GenerateARG.getMutationSet().put(posMutation, mut);
					  
				   } //end of computing the exact position of each mutation 
			   } //end if the edge has a father	   
	       } //end while for each edge of the ARG	
	    //System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%% END DECORATING WITH MUTAIONS %%%%%%%%%%%%%%%%%%%%%%%%%%");
	}//end of the procedure
		
	
	/**
	 * Compute the total number of mutations occurred in the generated ARG
	 * @return sum_nMut : integer representing the total number of mutations in the ARG
	 */
	public static int computeTotalMutationsNumber(){
		
		int sum_nMut = 0;
		Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator();  
	    Edge e;
	    
	    while (it_edges.hasNext()) {  
			
			   Integer keyE = it_edges.next();
			   e = GenerateARG.getGraphEdges().get((int)keyE);
			   
			   if(e.getId_fath() != -1){
				   sum_nMut = sum_nMut+e.getnMut();
			   }
	    }
		return sum_nMut;
	}
		
	/**
	 * The procedure creates a .txt file contaning the information about the SNP mutations in the extant units
	 * The first row represents the position of each mutation that is a double in the interval (0,1)
	 * The other rows represent the leave nodes as a matrix of 0 and 1 indicating the absence/presence of the corresponding mutation
	 * @param wholePath path of the directory where the output file is stored
	 * @param fileName of the output file
	 */
	public static void createSNP_TxtFile(String wholePath, String fileName) {
		Writer writer = null;
		File file = new File(""+wholePath+fileName);
		try {
		    writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
		    
		    TreeSet<Double> posMutations = GenerateARG.getSNPpositionsList();
		    
		    /*Iterator<Double> it_posMuts1 = posMutations.iterator();
		    while(it_posMuts1.hasNext()) {
		    	writer.write(it_posMuts1.next()+"\t");
		    	
		    }
		    writer.write("\n");*/
		    
		    Map<Double,Mutation> muts = GenerateARG.getMutationSet();
		    //For each leave print a row representing the absence or presence of SNP mutations
		    for(int i = 0; i < GenerateARG.getExtantUnits()+GenerateARG.getExtantUnitsunderSelection(); i++) {
		    	Iterator<Double> it_posMuts = posMutations.iterator();
			    while(it_posMuts.hasNext()) {
			    	Double position = it_posMuts.next();
			    	Mutation mut = muts.get(position);
			    	if(mut.getLeaves().contains(i))
			    		writer.write("1\t");
			    	else
			    		writer.write("0\t");
			    }
		    	
		    	writer.write("\n");
		    }
		    
		   writer.close();   
		}//fine try 
		
		catch (IOException ex) {
		// report
		} 
	}//fine procedure
}//end class

