import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Scanner; 
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Iterator; 




/**
 * Main class of SimRa.
 * Generates the ARG randomly and then decorates it with SNP mutations and (optionally) STRs mutations.
 * 
 * @author Anna Paola Carrieri
 * @author Filippo Utro
 * @author Aritra Bose
 * @version August 1, 2018.
 */

public class GenerateARG {

	//***************************** STATIC VARIABLES  ***************************************

	//total number of run for each set up
	private static int totalRuns;
	//List of the current active lineages 
	private static List<Edge> activeEdges;

	private static Map<String, List<Edge>> MapActiveEdge = new HashMap<String, List<Edge>>(); 

	//Set E of actual edges of the ARG G=(V,E)
	private static Map<Integer,Edge> graphEdges;
	//Set V of actual nodes of the ARG G=(V,E)
	private static Map<Integer,Node> nodeSet;
	//Map of SNP Mutations in the ARG
	private static Map<Double,Mutation> mutationSet;
	//Ordered list of SNP positions
	private static TreeSet<Double> SNPpositionsList;
	//Set of SplitPoints due to recombination events
	private static TreeSet<Double> splitPoints;
	//m = number of extant units
	private static int extantUnits=0;

	//m'
	private static int extantUnitsunderSelection=0;


	//L = number of live lineages L inizialized to m
	private static int L=0;
	private static int Ls = 0;
	private static int Lsbar = 0;
	//g = length of the chromosomal segment
	private static double g=0;
	//c = lenght of the whole chromosome
	private static double c=0;
	//N = population size
	private static int N=0;
	private static int Ns=0; //population for selection

	private static Map<String,Integer> SelN = new HashMap<String,Integer>(); 
	private static Map<String, Integer> SelL = new HashMap<String, Integer>(); 

	//current generation at time
	private static double curr_generation=0;
	private static Map<String, Double>  currGeneration = new HashMap<String,Double>(); 
	//private static double scaled generation j = TN
	private static double scaledGeneration=0;
	//recombination rate recomb=(g/c)
	private static double recomb = 0;
	//estimated number of mutations
	private static int mutationsNumber = 0;
	//mutation rate
	private static double mu = 0;
	//Name of output file
	private static String filename ="";
	//Total number of recombination events
	private static int recombinationsNumber = 0;
	//Total number of coalescent events
	private static int coalescentNumber = 0;
	//Total number of consecutive recombination and coalescent events involving the same 2 edges
	private static int numberRECandCoal = 0;

	private static double sumDenom = 0;


	//***************************** GETTERS AND SETTERS ***************************************
	public static Map<Double, Mutation> getMutationSet() {
		return mutationSet;
	}
	public static TreeSet<Double> getSplitPoints() {
		return splitPoints;
	}
	public static void setSplitPoints(TreeSet<Double> splitPoints) {
		GenerateARG.splitPoints = splitPoints;
	}
	public static void setMutationSet(Map<Double, Mutation> mutationSet) {
		GenerateARG.mutationSet = mutationSet;
	}
	public static int getCoalescentNumber() {
		return coalescentNumber;
	}
	public static void setCoalescentNumber(int coalescentNumber) {
		GenerateARG.coalescentNumber = coalescentNumber;
	}
	public static int getRecombinationsNumber() {
		return recombinationsNumber;
	}
	public static void setRecombinationsNumber(int recombinationNumber) {
		GenerateARG.recombinationsNumber = recombinationNumber;
	}
	public static double getMu() {
		return mu;
	}
	public static void setMu(double mu) {
		GenerateARG.mu = mu;
	}
	public static int getMutationsNumber() {
		return mutationsNumber;
	}
	public static void setMutationsNumber(int mutationsNumber) {
		GenerateARG.mutationsNumber = mutationsNumber;
	}
	public static double getRecomb() {
		return recomb;
	}
	public static void setRecomb(double recomb) {
		GenerateARG.recomb = recomb;
	}
	public static double getG() {
		return g;
	}
	public static void setG(double g) {
		GenerateARG.g = g;
	}
	public static int getExtantUnits() {
		return extantUnits;
	}
	public static void setExtantUnits(int extantUnits) {
		GenerateARG.extantUnits = extantUnits;
	}
	public static double getC() {
		return c;
	}
	public static void setC(double c) {
		GenerateARG.c = c;
	}
	public static int getL() {
		return L;
	}
	public static void setL(int l) {
		L = l;
	}

	public static void SetL(String str, int val) {
		SelL.put(str, val);
	}
	public static int GetL(String str) {
		return SelL.get(str);
	}
	public static int getLs() {
		return Ls;
	}
	public static void setLs(int ls) {
		Ls = ls;
	}
	public static int getLsbar() {
		return Lsbar;
	}
	public static void setLsbar(int lsbar) {
		Lsbar = lsbar;
	}
	public static int getN() {
		return N;
	}
	public static void setN(int n) {
		N = n;
	}


	public static int getNs() {
		return Ns;
	}
	public static void setNs(int ns) {
		Ns = ns;
	}
	public static void SetEffPop(String str, int val) {
		SelN.put(str, 2*val);
	}


	private static String getRndIt(List<String> items)
	{
		String s=null;


		double prod = 1.0, ValtoPut = 0.0;
		double sum = 0.0; 
		for(String  i: items) {
			if (!i.equals("99"))
				sum += (double)SelN.get("["+i+"]")/(double)(2*N);

			//prod *= FunS(val);
		}
		ArrayList<Double> tmp = new ArrayList<Double>();
		double newsum = 0.0;
		for(String  i: items){
			if (!i.equals("99")){
				newsum += SelN.get("["+i+"]")/(2*N*sum);
				tmp.add(newsum);
			}
			else
				tmp.add(newsum);
		}
		int flag = 0;
		double r1 = Math.random();
		//System.out.println("Random number is: "+r1);
		for(int i=0; i < tmp.size(); i++) {

			if(r1 < tmp.get(i)) {
				flag = 1; 
				return items.get(i);
			}
		}
		if (flag == 0) {
			for(String  i: items) {
				if ((SelN.get("["+i+"]")/(2*N)) > ValtoPut)
					return i;
			}
		}
		return s;
	}


	public static int GetEffPop(String str) {
		if (str.contains(","))
		{
			List<String> items = Arrays.asList(str.replace("[", "").replace("]","").replace(" ","").split(","));
			String s = getRndIt(items);
			//Collections.shuffle(items);
			//System.out.println(SelN.keySet());
			//System.out.println("["+items.get(0)+"]");
			//String s = items.get(0);
			return SelN.get("["+s+"]");
		}
		else
			return SelN.get(str);
	}
	public static double getcurrGeneration(String key) {
		return currGeneration.get(key);
	}
	public static void setGeneration(String key, double g) {
		GenerateARG.currGeneration.put(key,g);
	}


	public static double getGeneration() {
		return curr_generation;
	}
	public static void setGeneration(double g) {
		curr_generation = g;
	}

	public static double getScaledGeneration() {
		return scaledGeneration;
	}
	public static void setScaledGeneration(double scaledGeneration) {
		GenerateARG.scaledGeneration = scaledGeneration;
	}
	public static void SetMapActiveEdge(String str, List<Edge> e) {
		MapActiveEdge.put(str,e);
	}
	public static List<Edge> GetMapActiveEdge(String str){
		if (MapActiveEdge.get(str) == null) {
			return null;
		}
		else 
			return MapActiveEdge.get(str);
	}
	public static List<Edge> getActiveEdges() {
		return activeEdges;
	}
	public static void setActiveEdges(List<Edge> activeEdges) {
		GenerateARG.activeEdges = activeEdges;
	}
	public static Map<Integer, Edge> getGraphEdges() {
		return graphEdges;
	}
	public static void setGraphEdges(Map<Integer, Edge> graphEdges) {
		GenerateARG.graphEdges = graphEdges;
	}
	public static Map<Integer, Node> getNodeSet() {
		return nodeSet;
	}
	public static void setNodeSet(Map<Integer, Node> nodeSet) {
		GenerateARG.nodeSet = nodeSet;
	}
	public static double getCurr_generation() {
		return curr_generation;
	}
	public static void setCurr_generation(double curr_generation) {
		GenerateARG.curr_generation = curr_generation;
	}
	public static String getFilename() {
		return filename;
	}
	public static void setFilename(String filename) {
		GenerateARG.filename = filename;
	}
	public static int getNumberRECandCoal() {
		return numberRECandCoal;
	}
	public static void setNumberRECandCoal(int numberRECandCoal) {
		GenerateARG.numberRECandCoal = numberRECandCoal;
	}
	public static double getSumDenom() {
		return sumDenom;
	}
	public static void setSumDenom(double sumDenom) {
		GenerateARG.sumDenom = sumDenom;
	}
	public static int getTotalRuns() {
		return totalRuns;
	}
	public static void setTotalRuns(int totalRuns) {
		GenerateARG.totalRuns = totalRuns;
	}
	public static TreeSet<Double> getSNPpositionsList() {
		return SNPpositionsList;
	}
	public static void setSNPpositionsList(TreeSet<Double> sNPpositionsList) {
		SNPpositionsList = sNPpositionsList;
	}
	public static int StringToInt(String tmpstr) {
		tmpstr = tmpstr.replace("[", "");
		tmpstr = tmpstr.replace(",", "");
		tmpstr = tmpstr.replace("]","");
		int comb = Integer.parseInt(tmpstr);
		return comb;
	}
	public static String IntToString(int comb) {
		String tmpstr = Integer.toString(comb);
		String newstr =""; 
		for (int i=0; i<tmpstr.length(); i++) {
			if (i < tmpstr.length()-1)
				newstr = newstr+tmpstr.charAt(i)+",";
			else
				newstr = newstr+tmpstr.charAt(i)+"]";
		}
		return newstr;

	}
	public static String StrTrun(String str) {
		String tmpstr = str.replace("[", "");
		tmpstr = tmpstr.replace("]","");
		tmpstr = tmpstr.replace(",","");
		tmpstr = tmpstr.replace(" ","");
		return tmpstr; 
	}
	public static int getExtantUnitsunderSelection() {
		return extantUnitsunderSelection;
	}
	public static void setExtantUnitsunderSelection(int extantUnitsunderSelection) {
		GenerateARG.extantUnitsunderSelection = extantUnitsunderSelection;
	}
	/**
	 * Given the active lineages l1, l2, .., L the function return the sum of all the recombination rate
	 * rl (for l = 1, .., L) associated to each active lineage
	 * @return double number: sum of recombination rates of the active lineages
	 */
	public static double computeSumRates(){
		double sum = 0;

		for(int i = 0; i < activeEdges.size(); i++) {
			sum = activeEdges.get(i).getRate()+sum;
		}

		return sum; 
	}

	/**
	 * @param forSelection it is true if the value is computed for selection lineages, false otherwise.forSelection
	 * Given the active lineages l1, l2, .., L the function return the sum of all the recombination rate
	 * rl (for l = 1, .., L) associated to each active lineage
	 * @return double number: sum of recombination rates of the active lineages either in neutral or under selection
	 */
	public static double computeSumRates(boolean forSelection){
		double sum = 0;

		if (forSelection)
		{
			for(int i = 0; i < activeEdges.size(); i++) {
				if (activeEdges.get(i).isUnderSelection() )
					sum = activeEdges.get(i).getRate()+sum;
			}
		}
		else
		{
			for(int i = 0; i < activeEdges.size(); i++) {
				if (!activeEdges.get(i).isUnderSelection())
					sum = activeEdges.get(i).getRate()+sum;
			}
		}	
		return sum; 
	}

	public static double ComputeSumRates(List<Edge> ActiveEdges) {
		double sum = 0; 
		for (int i=0; i < ActiveEdges.size(); i++)
			sum = ActiveEdges.get(i).getRate()+sum;

		return sum;
	}
	/**
	 * Implement the exponential distribution with lambda as parameter
	 * @param lambda parameter of the exponential distribution
	 * @return the time t of the next event
	 */
	public static double computeNextTimeEvent(double lambda){
		double random = Math.random();

		//x = ln(1-p)/-lamda
		double t = (Math.log(1-random))/(0-lambda);

		return t;
	}

	/**
	 * Implement the exponential distribution with lambda and N as parameters
	 * @param lambda parameter of the exponential distribution
	 * @return the time t of the next event
	 */
	public static double computeNextTimeEvent(double lambda, int N){
		double random = Math.random();

		//x = ln(1-p)/-lamda
		double t = (Math.log(N*(1-random)))/(0-lambda);

		return t;
	}



	/** 
	 * Compute the time T = T+t to the next event using the exponential distribution
	 * @param time t computed by exponential distribution time is in [0,1)
	 * @return the time T updated basing on to the current level or generation
	 */
	public static double computeNextGeneration(double time){
		return time+getGeneration();
	}

	public static double computeNextGeneration(String key, double time){
		return time+getcurrGeneration(key);
	}

	/**
	 * ComputeIndexNextEvent function: compute the type of the next event (coalescent if
	 * index == 0, recombination if index > 0)
	 * @param lambda
	 * @param binCoef (L 2) where L is the # of the active lineages
	 * @return the index indicating the type of the next event
	 */
	public static int computeIndexNextEvent(double lambda, double binCoef){

		double x = Math.random()*lambda;
		double left = 0;
		double right = binCoef;
		boolean found = false;
		System.out.println("X is "+x+" -- bincoef is "+binCoef);

		//0 is associated with the coalescent event, from 1 to L with recombination events
		int index = 0;

		while(!found) {

			if(x >= left && x < right) {
				found = true;
			}
			else{
				index++;
				left = right;
				right = right+activeEdges.get(index-1).getRate();
			}
		}
		return index;
	}

	/**
	 * ComputeIndexNextEvent function: compute the type of the next event (coalescent if
	 * index == 0, recombination if index > 0)
	 * @param lambda
	 * @param binCoef (L 2) where L is the # of the active lineages
	 * @param inselection, true if the event has to be in the lineage under selection, false otherwise
	 * @return the index indicating the type of the next event
	 */
	public static int computeIndexNextEvent(double lambda, double binCoef, ArrayList<Integer> edges){

		double x = Math.random()*lambda;
		double left = 0;
		double right = binCoef;
		boolean found = false;
		//System.out.println("X is "+x+" -- bincoef is "+binCoef);
		//0 is associated with the coalescent event, from 1 to L with recombination events
		int index = 0;


		while(!found) {
			if(x >= left && x < right) {
				found = true;
			}
			else{
				index++;
				left = right;
				if (index>=edges.size()); //MARK: comment if 
				index=edges.size();	  //MARK: comment
				int idedge=edges.get(index-1);
				Edge e = graphEdges.get(idedge);
				int pos = activeEdges.indexOf(e);

				right = right+activeEdges.get(pos).getRate();

			}
		}
		return index;
	}

	/**
	 * Function that returns the binomial coefficient (n r)
	 * @param n integer value = number of elements in the set X
	 * @param r integer value = number of elements in each subset of X
	 * @return binomial coefficient = number of distinct k-elements subsets of X
	 */
	public static double binomialCoefficient(int n, int r) {
		double t = 1;

		int m = n - r; // r = Math.max(r, n - r);
		if (r < m) {
			r = m;
		}

		for (int i = n, j = 1; i > r; i--, j++) {
			t = t * i / j;
		}

		return t;
	}

	public static double ArrSum(double[] arr) {
		double sum = 0;
		for (int i = 0; i < arr.length; i++) 
			sum += arr[i];
		return sum; 
	}

	public static void findCombination(ArrayList<Double> fitvect, int arr[], int tmp[], int start, int end, int idx, int k, ArrayList<Double> order) {
		double tmpprod = 1.0; 

		if(idx == k) {
			// Current combination is ready to be printed, print it
			for (int j=0; j<k; j++) 
				tmpprod *= fitvect.get(tmp[j]); 
			order.add(tmpprod);
			return;	
		}

		for(int i = start; i <= end && end-i+1 >= k-idx; i++) {
			//System.out.print(i + " || " + idx + " --- ");
			tmp[idx] = arr[i]; 

			findCombination(fitvect, arr, tmp, i+1, end, idx+1, k, order);

		
		}

	}
	public static void findcombidx(double fitvect[], int arr[], int tmp[], int start, int end, int idx, int k, ArrayList<ArrayList<Integer>> order) {
		ArrayList<Integer> tmporder = new ArrayList<Integer>();
		if(idx == k) {
			// Current combination is ready to be printed, print it
			for (int j=0; j<k; j++)  
				tmporder.add(tmp[j]);
		
			order.add(tmporder);
			return;	
		}
		tmporder.clear();
		for(int i = start; i <= end && end-i+1 >= k-idx; i++) {
			//System.out.print(i + " || " + idx + " --- ");
			tmp[idx] = arr[i]; 

			findcombidx(fitvect, arr, tmp, i+1, end, idx+1, k, order);
		}

	}
	public static double findCombinationSum(ArrayList<Double> fitvect, int arr[], int n, int k) {
		int[] tmp = new int[k]; 
		double tmpsum = 0.0;
		ArrayList<Double> order = new ArrayList<Double>();
		Arrays.sort(arr);

		findCombination(fitvect,arr,tmp,0,n-1,0,k,order);

		for (int i = 0; i < order.size(); i++) {
			tmpsum += order.get(i); 
		}
		
		order.clear();
		return tmpsum;
	}

	public static ArrayList<ArrayList<Integer>> findCombinationIdx(double fitvect[],int arr[], int n, int k) {
		int[] tmp = new int[k];
		ArrayList<ArrayList<Integer>> order = new ArrayList<ArrayList<Integer>>(); 

		findcombidx(fitvect, arr, tmp, 0, n-1, 0, k, order);

		/*System.out.println("---- Print Vals -----");
		for(int i = 0; i < order.size(); i++)
			System.out.print(order.get(i)+" ");
		System.out.println("");*/
		return order; 

	}
	public static double ComputeSstar(double[] eBase) {
		double numer = 0.0; 
		double denom = 0.0; 
		for (int i = 0; i < eBase.length; i++) {
			numer += ((Math.pow(2,((eBase.length-1)-i-1)) - 1)*eBase[i]); 
			denom += ((1 - (Math.pow(2,((eBase.length-1) -i))))*eBase[i]); 
		}
		double sstar = (2*numer)/denom; 

		return sstar; 
	}

	public static double FunS(double val) {
		return ((1+val)/(2+val));
	}
	/*public static void errorInsertingParameters(){
		System.out.println("Usage:\n\t java -jar SimRa.jar [N] [m] [r] [mu] [g] [outputDirectory] [fileName]");
		System.out.println("\t\t N : integer representing the population size;");
		System.out.println("\t\t m : integer representing the sample size;");
		System.out.println("\t\t r : double representing the recombination rate in cM/Mb/gen;");
		System.out.println("\t\t mu : double representing the mutation rate in  mut/bp/gen x 10^(-8);");
		System.out.println("\t\t g : integer representing the segment length in Kb");
		System.out.println("\t\t [outputDirectory] : whole path of the directory where the output files will be stored");
		System.out.println("\t\t [fileName] : name for the output file");
		System.out.print("Additional optional parameters:\n\t");
		System.out.println("-STR [num] [s] [muSTRs]");
		System.out.println("\t\t -STR : to require an output file with information about STRs mutations;");
		System.out.println("\t\t num : integer representing the number of STRs for each node;");
		System.out.println("\t\t s : inital state for each STR locus;");
		System.out.println("\t\t muSTRs : mutation rate for STRs loci;");
		System.out.println("\n EXAMPLE 1: \n\t java -jar SimRa.jar 10000 200 0.1 0.7 75 output_SimRa/ output1");
		System.out.println("\n EXAMPLE 2: \n\t java -jar SimRa.jar 10000 200 0.1 0.7 75 output_SimRa/ output2 -STR 15 40 2.0");
		System.exit(1);
	}*/

	public static void errorInsertingParameters(){
		System.out.println("Usage:\n\t java -jar SimRa.jar [N] [r] [mu] [g] [iter] [eflag] [m] [s]");
		System.out.println("\t\t N : integer representing the population size;");
		System.out.println("\t\t r : double representing the recombination rate in cM/Mb/gen;");
		System.out.println("\t\t mu : double representing the mutation rate in  mut/bp/gen x 10^(-8);");
		System.out.println("\t\t g : integer representing the segment length in Kb");
		System.out.println("\t\t iter : Number of iterations;");
		System.out.println("\t\t eflag : Epistatic interaction flag;"); 
		System.out.println("\t\t m : array of integers less than N/2 representing the extant sample size;");
		System.out.println("\t\t s : array of doubles representing selection coefficient at each loci;");
		System.out.println("\n EXAMPLE: \n\t java -jar SimRa.jar -N 10000 -r 0 -mu 1 -g 25 -m 10 20 30 40 -eflag 1 -s 0.3 0.3 0.3 ");
		System.out.println("\n Above example is with selection on three loci and 4 randomly sampled extant units; ");
		System.out.println("\n As epistatic flag (eflag) is SET in the above example you need to provide four epistatic interaction coefficients less than s;");
		System.exit(1);
	}
	
	public static void creteDot(String outFileName) throws FileNotFoundException, UnsupportedEncodingException
	{
		HashMap<Double,HashSet<Integer>> nodepertime = new HashMap<Double, HashSet<Integer>>();
		SortedSet<Double> times = new TreeSet<Double>(Collections.reverseOrder());
		for (Integer v:graphEdges.keySet())
		{
			Edge e=graphEdges.get(v);
			if (nodepertime.containsKey(e.getTime()))
				nodepertime.get(e.getTime()).add(e.getId_fath());
			else
			{
				HashSet<Integer> nodes= new HashSet<Integer>();
				nodes.add(e.getId_fath());
				nodepertime.put(e.getTime(), nodes);
				times.add(e.getTime());
			}
		}


		PrintWriter writer = new PrintWriter(outFileName, "UTF-8");
		//write header
		writer.println("digraph {"); 
		writer.println("{ node [shape=plaintext, fontsize=16];");
		for (Double t:times)
		{
			writer.print(t);
			if (t!=0)
				writer.print("->");
		}
		writer.println(" }");

		for (Double t:times)
		{
			writer.println("{rank=same; "+t+";");
			for (Integer n: nodepertime.get(t))
				if (n!=-1)
					writer.print(n+";");
			writer.println("}");
		}
		writer.println("{rank=same; 0");

		for (int  n=0; n<getExtantUnits()+getExtantUnitsunderSelection(); n++)
			writer.print(n+";");
		writer.println("}");
		//write edges
		for (Integer v:graphEdges.keySet())
		{
			if(graphEdges.get(v).getId_fath()>-1)
			{
				writer.print(graphEdges.get(v).getId_fath()+" -> "+graphEdges.get(v).getId_son());
				if(graphEdges.get(v).isUnderSelection())
					writer.println("[color=red,penwidth=3.0];");
				else
					writer.println(";");
			}
		}
		//end of the files
		writer.println("}");
		writer.close();
	}

	public static void shuffleArray(int[] a) {
		int n = a.length;
		Random random = new Random();
		random.nextInt();
		for (int i = 0; i < n; i++) {
			int change = i + random.nextInt(n - i);
			swap(a, i, change);
		}
	}

	private static void swap(int[] a, int i, int change) {
		int helper = a[i];
		a[i] = a[change];
		a[change] = helper;
	}


	@SuppressWarnings("resource")
	public static Map<String, Double> ComputeFitGroup(double[] AllEpiFit, int k_way, Map<String,Double> Delta,
			int epiflag){

		List<String> deltakey = new ArrayList<String>(Delta.keySet()); 

		//Work with the indices instead of the actual fitness values
		int [] epiidxarr = new int[AllEpiFit.length];
		for (int i = 0; i < AllEpiFit.length; i++)
			epiidxarr[i] = i;

		//Storage handlers of the combination indices for epsistatic effect
		ArrayList<ArrayList<Integer>> IdxComb = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Double>> EpiFitList = new ArrayList<ArrayList<Double>>();

		//Loop governing and generating the combinations 
		for (int j = 1; j <= k_way; j++) {
			ArrayList<ArrayList<Integer>> storecombIdx = new ArrayList<ArrayList<Integer>>(); 
			storecombIdx = findCombinationIdx(AllEpiFit, epiidxarr, AllEpiFit.length, j);	
			for (int i = 0; i < storecombIdx.size(); i++)
				IdxComb.add(storecombIdx.get(i));
		}

		ArrayList<String> Combinations = new ArrayList<String>(); 
		//Obtain the actual epistatic sets from the indices: Relates to Laxmi's tableau
		for (int i = 0; i < IdxComb.size(); i++) {
			ArrayList<Integer> tmpcomb = new ArrayList<Integer>(); 
			ArrayList<Double> tmpvals = new ArrayList<Double>();
			tmpcomb = IdxComb.get(i);
			String combstring = tmpcomb.toString();
			Combinations.add(combstring);
			for(int j = 0; j < tmpcomb.size(); j++) 
				tmpvals.add(AllEpiFit[tmpcomb.get(j)]);
			EpiFitList.add(tmpvals);

		}
		//Print Handlers
		for(int i = 0; i < Combinations.size(); i++)
			System.out.print(Combinations.get(i)+" ---- ");
		System.out.println("");
		//Print Handlers
		for (int i = 0; i < EpiFitList.size(); i++) {
			ArrayList<Double> tmpval = new ArrayList<Double>(); 
			tmpval = EpiFitList.get(i); 
		}

		Map<String,Double> FitGroups = new HashMap<String, Double>(); 
		int numcombs = (int)(Math.pow(2, k_way)-1); 
		for(int ij = 0; ij < EpiFitList.size(); ij++) {
			ArrayList<Double> EachEpi = new ArrayList<Double>(); 
			EachEpi = EpiFitList.get(ij);
			int EpiNum = EachEpi.size();

			/*if(EpiNum > 1 && epiflag == 0) {
				System.out.println("key is "+Combinations.get(ij));
				int [] idxarr = new int[EpiNum]; 
				for (int i = 0; i < EpiNum; i++)
					idxarr[i] = i; 

				double [] eBase = new double[EpiNum+1];
				eBase[0] = 1; 
				double numcombs = 0;
				for (int i = 1; i <= EpiNum; i++)
				{
					numcombs = numcombs + BinomialCoefficient.binomialCoeff(EpiNum,i);
					eBase[i] = findCombinationSum(EachEpi, idxarr, EpiNum,i);
				}

				//double ValtoPut = FunS(sstar+Delta.get(Combinations.get(ij)));
				double ValtoPut = FunS(Delta.get(Combinations.get(ij)));
				FitGroups.put(Combinations.get(ij), ValtoPut);
			}*/
			
			if (EpiNum == 1){
				for (int i = 0; i < EachEpi.size(); i++) 
					FitGroups.put(Combinations.get(ij), FunS(EachEpi.get(i)));
			}
			else if (epiflag == 0) {
				double prod = 1.0, ValtoPut = 0.0;
				double sum = 0.0; 
				//System.out.print("EACHEPI?");
				for(int i = 0; i < EachEpi.size(); i++) {
					//System.out.println(" "+EachEpi.get(i));
					//sum += EachEpi.get(i);
					prod *= FunS(EachEpi.get(i));
				}
				FitGroups.put(Combinations.get(ij), prod);
			}
			else if (epiflag == 1) {
				System.out.println("Epistasis flag is ON");
				System.out.println("--------------------");
				Scanner sc = new Scanner(System.in);
				System.out.print("Enter the interaction value for "+Combinations.get(ij)+": ");
				FitGroups.put(Combinations.get(ij), FunS(Double.parseDouble(sc.nextLine())));
			}	
		}

		 FitGroups.put("[99]",1.0);
		//FitGroups.put("[99]",0.5);
		return FitGroups;
	}

	public static Map<String, Double> ComputeEffPop(Map<String,Double> FitGroups, int epiflag){
		Map<String, Double> EffPop = new HashMap<String,Double>();
		
		for (Map.Entry<String, Double> entry: FitGroups.entrySet()) {
			double singsum = 0.0, dubsum = 0.0, sum3 = 0.0;
			String tmpstr = StrTrun(entry.getKey());
			EffPop.put(entry.getKey(),entry.getValue());
			if(tmpstr.length() == 1) {		
				for (Map.Entry<String, Double> entry2: FitGroups.entrySet()) {
					String newstr = StrTrun(entry2.getKey());
					if (newstr.length() == 2 && newstr.contains(tmpstr)) 
						dubsum += entry2.getValue();
					if(newstr.length() == 3 && newstr.contains(tmpstr))
						sum3 += entry2.getValue();
				}
				EffPop.put(entry.getKey(), (entry.getValue() - dubsum + sum3));
				//System.out.println("Sum is: "+(entry.getValue() - dubsum + sum3));
			}

			if(tmpstr.length() == 2) {
				sum3 = 0.0;
				if (!(tmpstr.equals("99"))) {
					for (Map.Entry<String, Double> entry2: FitGroups.entrySet()) {
						String newstr = StrTrun(entry2.getKey());
						if(newstr.length() > 2 && newstr.contains(tmpstr))
							sum3 = entry2.getValue();
					}
					EffPop.put(entry.getKey(), (entry.getValue()-sum3));
				}	
				else {
					singsum = 0.0; 
					dubsum = 0.0; 
					sum3 = 0.0; 
					
						for (Map.Entry<String, Double> entry2: FitGroups.entrySet()) {
							if (StrTrun(entry2.getKey()).length() == 1) 
								singsum += entry2.getValue(); 
							if (StrTrun(entry2.getKey()).length() == 2 && !(tmpstr.equals(StrTrun(entry2.getKey()))))
								dubsum += entry2.getValue(); 
							if(StrTrun(entry2.getKey()).length() == 3)
								sum3 += entry2.getValue();
						}
					/* 
					if(epiflag == 0) {
						double prod = 1;
						for (Map.Entry<String, Double> entry2: FitGroups.entrySet()) 
							if (StrTrun(entry2.getKey()).length() == 1) 
								prod *= entry2.getValue(); 
	
						dubsum = prod; 
					} */
					System.out.println("Singsum is "+singsum+" dubsum is "+dubsum+" sum3 is "+sum3);
					EffPop.put(entry.getKey(), entry.getValue() - singsum + dubsum - sum3);
					//EffPop.put(entry.getKey(), entry.getValue());
				}
			}

			if(tmpstr.length() == 3)
				EffPop.put(entry.getKey(), entry.getValue()); 	
		}
		 
		return EffPop;	
	}

	public static void SetEdge(Node coal_node, ArrayList<Interval> union, String ipast,  Map<String, ArrayList<Integer>> edgeContain)
	{
		Edge new_edge;
		new_edge = new Edge(graphEdges.size(), coal_node.getId(), -1, union, GetEffPop(ipast));
		new_edge.computeLength();
		new_edge.computeRate(GetEffPop(ipast), getRecomb(), getG());
		edgeContain.get(ipast).add(new_edge.getIdEDGE());
		graphEdges.put(new_edge.getIdEDGE(), new_edge);
		activeEdges.add(new_edge);
		MapActiveEdge.get(ipast).add(new_edge);
		System.out.println("For "+ipast+" map active edge size is "+MapActiveEdge.get(ipast).size());
	}
	
	//Check if has childfriend so [0] -> [0,1], [0,2]
	public static Boolean hasSubFriend(String ipast)
	{
		String tmp = ipast.replace("]",",");
		for(Map.Entry<String, List<Edge>> entry: MapActiveEdge.entrySet()) {
		
			if (entry.getKey().contains(tmp) && entry.getValue().size()>0) //if [0, is contained in [0,1 and [0,1 is not empty
				return true;		
		}

		return false;
	}
	
	public static void SetEdge(Node coal_node, ArrayList<Interval> union, String ipast, String actualKey, Map<String, ArrayList<Integer>> edgeContain, Map<String, ArrayList<Edge>> edgeWait,Map<Double, Edge> edgeTimeWait, Double timew, Double timenow)
	{
		Edge new_edge;
		new_edge = new Edge(graphEdges.size(), coal_node.getId(), -1, union, GetEffPop(ipast));
		new_edge.computeLength();
		new_edge.computeRate(GetEffPop(ipast), getRecomb(), getG());
		
		activeEdges.add(new_edge);
		graphEdges.put(new_edge.getIdEDGE(), new_edge);
		
		
		//if the map is empty and has friend do not move
		if (hasSubFriend(ipast) & MapActiveEdge.get(ipast).isEmpty())
		{
			edgeContain.get(actualKey).add(new_edge.getIdEDGE());
			MapActiveEdge.get(actualKey).add(new_edge);
			return;
		}
		
		if (timenow<timew)
		{
			edgeContain.get(ipast).add(new_edge.getIdEDGE());
			MapActiveEdge.get(ipast).add(new_edge);
		}
		else
		{
		
			if (!edgeWait.containsKey(new_edge))
				edgeWait.put(ipast, new ArrayList<Edge>());
			edgeWait.get(ipast).add(new_edge);		
			edgeTimeWait.put(timenow, new_edge);
			
		}
		
		System.out.println("For "+ipast+" map active edge size is "+MapActiveEdge.get(ipast).size());
	
	}
	
	

	public static void UpdateEffPop(double[] AllEpiFit, List<Double> Epifit, String ipast, int k_way, Map<String,Double> Delta,
			Map<String, Double> FitGroup, Map<String, Double> EffPop, List<Integer> DelFitIdx, int argN, String NS, int epiflag) {
		System.out.println("");
		System.out.println("------- Before Fitlist --------");
		for (int ij = 0; ij < Epifit.size(); ij++)
			System.out.print(" "+Epifit.get(ij));
		System.out.println("");
		int idxtorem = Integer.parseInt(StrTrun(ipast));
		for(int i = 0;i < DelFitIdx.size(); i++) {
			if(DelFitIdx.get(i) < idxtorem) {
				idxtorem--;
				break;
			}
		}
		Epifit.remove(idxtorem);

		System.out.println("------- After Fitlist --------");
		for (int ij = 0; ij < Epifit.size(); ij++)
			System.out.print(" "+Epifit.get(ij));
		System.out.println("");
		//edgeContain.remove(ipast);	
		//MapActiveEdge.remove(ipast);

		double[] revepifit = new double[Epifit.size()];

		for(int i=0; i<Epifit.size(); i++)
			revepifit[i] = Epifit.get(i);
		System.out.println("---- Before Fitness Values ----");
		for(Map.Entry<String, Double> entry: FitGroup.entrySet()) {
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}
		//String newkey = "";
		FitGroup = ComputeFitGroup(revepifit, k_way, Delta,epiflag);
		System.out.println("---- After Fitness Values ----");
		Map<String, Double> tmpfitgroup = new HashMap<String, Double>(); 
		for(Map.Entry<String, Double> entry: FitGroup.entrySet()) {
			/*	String idstr = StrTrun(entry.getKey());
			ArrayList<Integer> newidlist = new ArrayList<Integer>(); 
			if (!entry.getKey().equals(NS)) {
				for(int jk = 0; jk < idstr.length(); jk++) {
					int idx = Character.getNumericValue(idstr.charAt(jk));
					double fitval = revepifit[idx];
					for(int kl = 0; kl < AllEpiFit.length; kl++) {
						if (AllEpiFit[kl] == fitval)
							newidlist.add(kl);
						else
							continue;
					}			
				}
				newkey = newidlist.toString(); 
				tmpfitgroup.put(newkey, entry.getValue());
			}
			else {
				newkey = NS;
				tmpfitgroup.put(newkey, entry.getValue());
			}*/
			tmpfitgroup.put(entry.getKey(), entry.getValue());
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}

		System.out.println("---- Before EffPop Values ----");
		for(Map.Entry<String, Double> entry: EffPop.entrySet()) {
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}
		System.out.println("------------------------");
		EffPop = ComputeEffPop(tmpfitgroup,epiflag);
		for(Map.Entry<String, Double> entry: EffPop.entrySet()) 
			SetEffPop(entry.getKey(), (int) Math.round(argN*entry.getValue()));
		System.out.println("---- After EffPop Values ----");
		for(Map.Entry<String, Double> entry: EffPop.entrySet()) {
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}
		System.out.println("------------------------");
	}

	public static void MapEdgeRemove(String ipast,int id_first_edge) {
		List<Edge> tmpActiveEdge = MapActiveEdge.get(ipast);
		int torem = 0;
		for (int jj = 0; jj < tmpActiveEdge.size(); jj++)
			if (id_first_edge == tmpActiveEdge.get(jj).getIdEDGE())
				torem = jj;
			else 
				continue;
		System.out.println("To remove from mapactiveedge: "+torem);
		System.out.println("Size is : "+tmpActiveEdge.size());
		MapActiveEdge.get(ipast).remove(torem);
	}

	public static void SelectGroups(Map<String,ArrayList<Integer>> IdxGroups, Map<String, ArrayList<Integer>> FinSel, 
			List<String> keys, ArrayList<Integer> immortal, String NS){

		Map<String,ArrayList<Integer>> SelGroups = new HashMap<String,ArrayList<Integer>>();
		for(int jk = 0; jk < keys.size(); jk++) {
			String tmpstr = StrTrun(keys.get(jk));


			if(tmpstr.length() == 1) {

				SelGroups.put(keys.get(jk), IdxGroups.get(keys.get(jk)));
			}
			else if(tmpstr.length() > 1 && !(keys.get(jk).equals(NS))) {
				ArrayList<Integer> tmpid = new ArrayList<Integer>(); 
				for (int i = 0; i < tmpstr.length(); i++) {
					String instr = "["+Character.toString(tmpstr.charAt(i))+"]";

					if (i == 0)
						tmpid.addAll(IdxGroups.get(instr));
					else 
						tmpid.retainAll(IdxGroups.get(instr));
				}
				SelGroups.put(keys.get(jk), tmpid);
			}
			else 
				continue;
		}
		for (Map.Entry<String, ArrayList<Integer>> entry: SelGroups.entrySet()) {
			ArrayList<Integer> tmpid = entry.getValue();
			if(StrTrun(entry.getKey()).length() < 3) {
				ArrayList<Integer> remtmp = new ArrayList<Integer>(); 
				for (Map.Entry<String, ArrayList<Integer>> entry2: SelGroups.entrySet()) {
					if (!(entry.getKey().equals(entry2.getKey()))){
						if(StrTrun(entry2.getKey()).contains(StrTrun(entry.getKey()))) {
							remtmp.addAll(entry2.getValue());
						}
					}
				}	
				Set<Integer> tmpset = new HashSet<Integer>(remtmp);
				ArrayList<Integer> revtmp = new ArrayList<Integer>(tmpset);
				tmpid.removeAll(revtmp);
			}
			FinSel.put(entry.getKey(), tmpid);		
		}

		//Check what is left and assign to 99	
		ArrayList<Integer> remim = new ArrayList<Integer>(); 
		for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
			remim.addAll(entry.getValue());
		}
		immortal.removeAll(remim);

		FinSel.put(NS, immortal);

	}

	/*************************** MAIN Generate the ARG randomly *********************************
	 * Main function implements the ARG sampling algorithm
	 * @param args[0]  integer representing the population size
	 * @param args[1]  integer representing the sample size
	 * @param args[2]  double representing the recombination rate in cM/Mb/gen (x 10^(-8))
	 * @param args[3]  double representing the mutation rate in  mut/bp/gen x 10^(-8)
	 * @param args[4]  integer representing the segment length in Kb
	 * @param args[5]  whole path of the directory where the output files will be stored
	 * @param args[6]  name for the output files
	 * @param args[7] 
	 * @param args[8]  
	 * @param args[9]  
	 * @param args[10] 
	 */

	public static void main(String[] args) {


		int argN;
		int argm;
		double argRecomb;
		double argMut;
		int argG;
		String argPath;
		String argFileName;
		int epiflag; 
		String fits = "";
		
		int runs=1; 
        ArrayList<Integer> mvect = new ArrayList<Integer>();
        ArrayList<Integer> gvect =  new ArrayList<Integer>();
        Map<String, List<String>> params = new HashMap<>();
       
		Runtime runtime = Runtime.getRuntime();
		long start = System.currentTimeMillis();

		List<String> options = null;
		for (int i = 0; i < args.length; i++) {
		    final String cmds = args[i];

		    if (cmds.charAt(0) == '-') {
		        if (cmds.length() < 1) {
		            System.err.println("Error at argument " + cmds);
		            return;
		        }

		        options = new ArrayList<>();
		        params.put(cmds.substring(1), options);
		    }
		    else if (options != null) {
		        options.add(cmds);
		    }
		    else {
		        System.err.println("Illegal parameter usage");
		        return;
		    }
		}
		 
        if(params.containsKey("N") && params.containsKey("g") && params.containsKey("mu")
        		&& params.containsKey("r") && params.containsKey("iter") 
        		&& params.containsKey("eflag") && params.containsKey("m")
        		&& params.containsKey("s")) {
        	System.out.println("All parameters OK.");
        }
        else {
			errorInsertingParameters();
		}
        

		//System.out.println(params.get("N").get(0));
        argN = Integer.parseInt(params.get("N").get(0));
		argRecomb = Double.parseDouble(params.get("r").get(0));
		argMut = Double.parseDouble(params.get("mu").get(0));
		runs = Integer.parseInt(params.get("iter").get(0));
		epiflag = Integer.parseInt(params.get("eflag").get(0));
		
		for (int i = 0; i < params.get("g").size(); i++ ) {
			gvect.add(Integer.parseInt(params.get("g").get(i)));
		}
		for (int i = 0; i < params.get("m").size(); i++ ) {
			mvect.add(Integer.parseInt(params.get("m").get(i)));
		}
		double[] OldEpiFit = new double[params.get("s").size()];
		
		for (int i = 0; i < params.get("s").size(); i++ ) {
			OldEpiFit[i] = Double.parseDouble(params.get("s").get(i));
			fits = fits+params.get("s").get(i);
		}
		int k_way = OldEpiFit.length; //1;	//how many maximum interacting SNPs 
		int numcombs = (int)Math.pow(2,k_way)-k_way-1;
		Map <String, Double> EpiVal = new HashMap<String,Double>(); 	
		
		/*if(epiflag == 1) {
			System.out.println("Epistasis flag is ON");
			System.out.println("--------------------");
			Scanner sc = new Scanner(System.in);
			 int ct1 = 0, ct2 = 1; 
			for (int i = 0; i < numcombs; i++) {
				int ct=0;
				while (ct < k_way) {
					System.out.print("Interaction value for ["+Integer.toString(i)+","+Integer.toString(ct+1)+"] :");
					EpiVal.put("["+Integer.toString(i)+","+Integer.toString(ct+1)+"]", Double.parseDouble(sc.nextLine()));
					ct++;
				}
				numcombs = numcombs - ct; 
			}
			
			System.out.println("---------------");
			
			System.out.println("---- Epi Fitness values ----");
			for(Map.Entry<String, Double> entry: EpiVal.entrySet()) {
				System.out.println(entry.getKey() + " : " + entry.getValue());
			}
			System.out.println("------------------------");
			
		}*/
		
		double [] AllEpiFit = OldEpiFit;	//List of the s1, s2, s3.... values
		//double[] AllEpiFit = new double[] {};
		List<Double> Epifit = new ArrayList<Double>();
		double[] FunFit = new double[AllEpiFit.length];			//Stores the f(s) values

		for (int i = 0; i < AllEpiFit.length; i++) {
			Epifit.add(AllEpiFit[i]);
			FunFit[i] = FunS(AllEpiFit[i]);
		}

		String NS = "[99]";		//Immortal key! 

		Map <String, Double> FitGroup = new HashMap<String,Double>(); 	
		Map <String, Double> EffPop = new HashMap<String,Double>(); 	//This will store f' to calculate N*f'

		//int EpiFlag = 0; //Flag to indicate whether epistasis is ON? We don't need this anymore ?
		

		// This map indicates for which epistasis will be in effect with the delta associated 
		Map <String, Double> Delta = new HashMap<String,Double>();

		
		//Calculating f' for each key
		
		FitGroup = ComputeFitGroup(AllEpiFit, k_way, Delta,epiflag); 
		
		System.out.println("---- Fitness Values ----");
		for(Map.Entry<String, Double> entry: FitGroup.entrySet()) {
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}
		System.out.println("------------------------");
		
	
		
		System.out.println("---------COMPUTING EFFPOP---------------");
		EffPop = ComputeEffPop(FitGroup,epiflag);

		System.out.println("");
		System.out.println("------------------------");

		List<String> keys = new ArrayList<String>(FitGroup.keySet());		//keep a record of the keys being used
		//We will need this to assign keys to argm
		//System.exit(1);
		for(int im=0;im<mvect.size(); im++)
		{
			//for(int indf=0;indf<fvect.length; indf++)
			//{
			for(int indg=0; indg<gvect.size(); indg++)
			{
				//***** Number of run *****
				setTotalRuns(runs);

				argN = 100;
				N = argN;
				argm = mvect.get(im);
				//argRecomb = 1;
				//argMut = 1.5;
				argG = gvect.get(indg);
				argPath = "";
				argFileName = "prova";

				//compute number of extant unit under selection;
				double fns = 0.2 + (int)(Math.random() * ((0.8 - 0.2)));
				fns = 0.2;
				//argmUnderSection = (int) Math.round(fns*argm);    //Math.round(ratioUnderselection*argm);


				String filestast=argPath+"N"+String.valueOf(argN)+"_m"+String.valueOf(argm)+"_m"+ 
						String.valueOf(fns)+"_g"+String.valueOf(argG)+"_s"+fits+"_r"+ 
						String.valueOf(argRecomb)+".txt";


				//***************************** PARAMETER INITIALIZATION  ***************************************

				//***** Path directory output *****
				//Filename Format N#_g#_mu#_r#_m#
				//filename = "N"+getN()+"_g"+getG()+"_mu"+getMu()+"_ro"+getRecomb()+"_m"+getExtantUnits();
				filename = filestast;

				//***** Recombination Rate ****
				setRecomb(argRecomb*Math.pow(10, -8));

				//**** SNP Mutation Rate ****
				mu = argMut*Math.pow(10, -8);

				//**** Segment length g ***
				setG(argG*Math.pow(10, 3));

				//**** whole path of directory for outputfile ***
				//String wholePath = "output_SimRa/";
				String wholePath = argPath;

				//***************************** CREATION OF DIRECTORY FOR OUTPUT FILES  ***********************************
				File dir;
				dir = new File(""+wholePath);
				dir.mkdir();		
				BufferedWriter fwrite = null;
				File file = new File(""+wholePath+filename+"_STATS.txt");
				try {
					fwrite = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file, true), "utf-8"));

				} catch (UnsupportedEncodingException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}


				//***************************** ARG RANDOM GENERATION ALGORITHM   ***********************************

				Random rand = new Random();

				for(int numRun = 1; numRun <= getTotalRuns(); numRun++) {



					BufferedWriter fwritepick = null;
					String tnsName = "stats/selected1_m"+String.valueOf(argm);
					for (int i = 0; i < AllEpiFit.length; i++) 
						tnsName = tnsName+"_s"+String.valueOf(i+1)+"="+String.valueOf(AllEpiFit[i]);
					File filepick = new File(tnsName+"_r="+String.valueOf(argRecomb)+".txt");
					try {
						fwritepick = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filepick, true), "utf-8"));

					} catch (UnsupportedEncodingException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					} catch (FileNotFoundException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}


					//********* INITIALIZATION OF STATIC VARIABLES **************/
					setGeneration(0);
					setRecombinationsNumber(0);
					setMutationsNumber(0);
					setCoalescentNumber(0);
					numberRECandCoal = 0;
					nodeSet = new HashMap<Integer,Node>();
					mutationSet = new HashMap<Double,Mutation>();
					SNPpositionsList = new TreeSet<Double>();
					graphEdges = new HashMap<Integer,Edge>();
					activeEdges = new ArrayList<Edge>();
					splitPoints = new TreeSet<Double>();
					splitPoints.add(0.0);
					splitPoints.add(1.0);

					//SelN = new HashMap<String,Integer>(); 

					//********* START
					ArrayList<String> fitMidx = new ArrayList<String>();		//this will contain the key for each argm
					int[] idxm = new int[argm];
					List<Edge> MutMap = new ArrayList<Edge>();		

					Map<String, ArrayList<Integer>> IdxGroups = new HashMap<String, ArrayList<Integer>>();
					ArrayList<Integer> immortal = new ArrayList<Integer>(); 

					for (int i = 0; i < argm; i++) {
						idxm[i] = i; 
						immortal.add(i);
					}
					Map<String,Integer> selsize = new HashMap<String,Integer>(); 
					for (int i = 0; i < FunFit.length; i++) {
						ArrayList<Integer> tmpid = new ArrayList<Integer>();
						tmpid.add(i);
						selsize.put(tmpid.toString(),(int)Math.round(FunFit[i]*argm));
					}

					//Assignment of each single keys such as s1, s2, s3, etc. according to their respective f(s)
					if(argRecomb > 0) {

						for(int i = 0; i < FunFit.length; i++) {
							shuffleArray(idxm);
							ArrayList<Integer> tmpid = new ArrayList<Integer>();
							ArrayList<Integer> idx = new ArrayList<Integer>(); 

							int  thres = rand.nextInt(argm) + 1;
							//int thres = (int)Math.round(FunFit[i]*argm);
							
							for (int j = 0; j < thres; j++) {
								tmpid.add(idxm[j]);
							}
							idx.add(i);
							IdxGroups.put(idx.toString(), tmpid);
						}
					}
					else {
						Map<Integer, Integer> sizeMap = new TreeMap<Integer, Integer>();

						for(int i = 0; i < FunFit.length; i++) {
							sizeMap.put(i,(int)Math.round(FunFit[i]*argm));//* FunFit[i]*argm  0.25  0.3 change - > 0.5
						}
						
						
						LinkedHashMap<Integer, Integer> reverseSortedMap = new LinkedHashMap<Integer, Integer>();

						//Use Comparator.reverseOrder() for reverse ordering
						sizeMap.entrySet().stream().sorted(Map.Entry.comparingByValue(Comparator.reverseOrder())).forEachOrdered(x -> reverseSortedMap.put(x.getKey(), x.getValue()));

						Map<String, Integer> availablesizeMap = new TreeMap<String, Integer>();
						availablesizeMap.put(NS,immortal.size());

						shuffleArray(idxm);
						ArrayList<Integer> immcopy = immortal;
						Boolean first = true;
						for (Map.Entry<Integer, Integer> entry : reverseSortedMap.entrySet()) {
							ArrayList<Integer> tmpid = new ArrayList<Integer>();
							System.out.println("Key is " + entry.getKey()+" and value is "+entry.getValue());
							if (first) {
								String KeySel = "["+Integer.toString(entry.getKey())+"]";
								Collections.shuffle(immcopy);
								for (int j = 0; j < entry.getValue(); j++) 
									tmpid.add(immcopy.get(j));
								immcopy.removeAll(tmpid);
								IdxGroups.put(KeySel, tmpid);
								availablesizeMap.put(NS, immcopy.size());
								IdxGroups.put(NS, immcopy);
								availablesizeMap.put(KeySel, entry.getValue());
								first = false;
							}
							else
							{
								List<String> Tkeys = new ArrayList<String>(availablesizeMap.keySet());

								Collections.shuffle(Tkeys);
								for (String it: Tkeys)
								{								
									if (entry.getValue() <= availablesizeMap.get(it))
									{
										String KeySel;
										if (it.equals(NS))
										{
											KeySel = "["+Integer.toString(entry.getKey())+"]";
											for (int j = 0; j < entry.getValue(); j++) 
												tmpid.add(immcopy.get(j));
											immcopy.removeAll(tmpid);
										}
										else
										{
											ArrayList<Integer> tlist = IdxGroups.get(it);
											for (int j = 0; j < entry.getValue(); j++) 
												tmpid.add(tlist.get(j));
												
											for (int j = entry.getValue()-1; j>=0; j--) 
												tlist.remove(j);
											IdxGroups.put(it, tlist);
											List<String> items = new ArrayList<String>(Arrays.asList(it.replace("[", "").replaceAll("]", "").split(", ")));
											items.add(Integer.toString(entry.getKey()));
											Collections.sort(items);

											KeySel = "[";
											for (String i:items)
												KeySel =KeySel + i + ", ";

											KeySel = KeySel.substring(0, KeySel.length()-2) + "]";	
										}
										System.out.println("Key is " + KeySel +" and value is "+entry.getValue());
										availablesizeMap.put(it, availablesizeMap.get(it)-entry.getValue());
										availablesizeMap.put(KeySel, entry.getValue());
										IdxGroups.put(KeySel, tmpid);
										break;
									}

								}
							}

						}

					}


					for(Map.Entry<String, ArrayList<Integer>> entry: IdxGroups.entrySet()) {
						System.out.print(" || " + entry.getKey()+ " : ");
						for (int i = 0; i < entry.getValue().size(); i++)
							System.out.print(entry.getValue().get(i)+", ");
					}
					System.out.println("");	

					Map<String, ArrayList<Integer>> FinSel;
					if(argRecomb > 0) {

						//Storage of the selection groups containing keys as strings and the indices assigned as integers. 
						//Note: This will contain redundancies and we have to remove them to get the final selection Map. 
						//Assignment of final selection map
						FinSel = new HashMap<String, ArrayList<Integer>>();
						SelectGroups(IdxGroups, FinSel, keys, immortal, NS);
						for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
							System.out.print("|| "+entry.getKey()+ " : ");
							for (int i = 0; i < entry.getValue().size(); i++)
								System.out.print(entry.getValue().get(i)+", ");
						}
						System.out.println("");	
						for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
							List<Integer> tmplist = entry.getValue();
							for(Map.Entry<String, ArrayList<Integer>> entry2: FinSel.entrySet()) {
								for(int ik = 0; ik < tmplist.size(); ik++) {
									if(!entry2.equals(entry) && entry2.getValue().contains(tmplist.get(ik))) {
										if(StrTrun(entry2.getKey()).length() < StrTrun(entry.getKey()).length())
											entry2.getValue().remove(Integer.valueOf(tmplist.get(ik)));
										else
											entry.getValue().remove(Integer.valueOf(tmplist.get(ik)));
									}
								}
							}
						}

					}
					else
					{
						FinSel = IdxGroups;

					}

					for (int ij = 0; ij < argm; ij++) {
						for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
							if(entry.getValue().contains(ij))
								fitMidx.add(entry.getKey());
						}
					}


					for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
						System.out.print("|| "+entry.getKey()+ " : ");
						for (int i = 0; i < entry.getValue().size(); i++)
							System.out.print(entry.getValue().get(i)+", ");
					}

					System.out.println("");	
					int immnum = FinSel.get(NS).size(); 

					setExtantUnits(immnum); //m = input value - m'
					setExtantUnitsunderSelection(argm - immnum);
					
					System.out.println("---- Eff Pop ----");
					for(Map.Entry<String, Double> entry: EffPop.entrySet()) {
						System.out.println(entry.getKey() + " : " + entry.getValue());
						SetEffPop(entry.getKey(), (int) Math.round(argN*entry.getValue()));
					}


					for(Map.Entry<String, Double> entry: EffPop.entrySet()) {
						setGeneration(entry.getKey(), 0);
					}

					Map<String, Double> SegPos = new HashMap<String, Double>();
					
					//Assign segments positions to keys for recombination
					double[] posset= new double[keys.size()]; 

					for (int i = 0; i < keys.size(); i++) {
						if(i == 0) {
							SegPos.put(keys.get(i), 0.0);
							posset[i] = 0;
						}
						else if(i == keys.size()-1) {
							SegPos.put(keys.get(i), 1.0);
							posset[i] = 1; 
						}
						else {
							double val = Math.random();
							SegPos.put(keys.get(i), val);
							posset[i] = val;
						}
					}
					Arrays.sort(posset);

					System.out.println("------------Segments-------------");
					for(Map.Entry<String, Double> entry: SegPos.entrySet()) {
						System.out.print("|| "+entry.getKey()+ " : "+entry.getValue());
					}
					System.out.println("");


					//********* END
					//edgeContainer for each key
					Map <String, ArrayList<Integer>> edgeContain = new HashMap<String, ArrayList<Integer>>(); 
					//ensuring containers are cleared after each run! 
					MapActiveEdge.clear();
					Epifit.clear();
					MutMap.clear();
					List<Integer> DelFitIdx = new ArrayList<Integer>();
					AllEpiFit = OldEpiFit;
					for (int ik = 0; ik < AllEpiFit.length; ik++) {
						Epifit.add(AllEpiFit[ik]);
					}

					//********* INITIALIZATION OF LEAVE NODES ****************/

					for (int i = 0; i < argm; i++) {

						Node n = new Node(i,false,getGeneration()); 
						ArrayList<Interval> segments = new ArrayList<Interval>();
						segments.add(new Interval(0,1));
						n.setSegments(segments);
						nodeSet.put(i, n);
						ArrayList<Integer> tmparr = new ArrayList<Integer>(); 
						List<Edge> tmpActiveEdge = new ArrayList<Edge>();
						Edge e = null;

						e = new Edge(i,n.getId(), -1, segments, true, GetEffPop(fitMidx.get(i)));
						e.computeLength();
						e.computeDensity();
						e.computeRate(GetEffPop(fitMidx.get(i)), getRecomb(), getG());

						if (!edgeContain.containsKey(fitMidx.get(i))){
							tmparr.add(i);
							edgeContain.put(fitMidx.get(i), tmparr); 
						}
						else {
							tmparr = edgeContain.get(fitMidx.get(i)); 
							tmparr.add(i);
							edgeContain.put(fitMidx.get(i), tmparr); 
						}

						if (GetMapActiveEdge(fitMidx.get(i)) == null) {
							tmpActiveEdge.add(e);
							SetMapActiveEdge(fitMidx.get(i),tmpActiveEdge);	
						}
						else {
							tmpActiveEdge = GetMapActiveEdge(fitMidx.get(i)); 
							tmpActiveEdge.add(e);
							SetMapActiveEdge(fitMidx.get(i),tmpActiveEdge);
						}
						//Add the new lineage to the list of active lineages
						activeEdges.add(e);
						graphEdges.put(e.getIdEDGE(), e);
					}
					//Populate edgeContainer and MapActiveEdges
					for (int kj = 0; kj < keys.size(); kj++) {
						ArrayList<Integer> tmpid = new ArrayList<Integer>();
						ArrayList<Edge> tmpedge = new ArrayList<Edge>();
						if (edgeContain.containsKey(keys.get(kj)))
							continue; 
						else {
							System.out.println("Doesn't contain: "+keys.get(kj));
							edgeContain.put(keys.get(kj), tmpid);
						}
						if (MapActiveEdge.containsKey(keys.get(kj)))
							continue; 
						else
							MapActiveEdge.put(keys.get(kj),tmpedge);
					}

				
					setL(activeEdges.size());

					for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) 
						SetL(entry.getKey(),entry.getValue().size());

					//*************** Construction of the ARG backwards starting from the leaves ****************/
					List<String> tmpkeys = new ArrayList<String>(edgeContain.keySet());	
					while(activeEdges.size() !=1){
						System.out.println("");
						
						System.out.println("----  edges  ----");

						for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) {	
							System.out.println("Key: "+entry.getKey());
							for (int kl = 0; kl < entry.getValue().size(); kl++)
								System.out.print(" "+entry.getValue().get(kl));
							System.out.println("");
						}
						System.out.println("-------------------");
						System.out.println("There are "+activeEdges.size()+" active edges now!!");

						double tpast = 0.0;
						String ipast = "", toRemKey = ""; 
						int flag = 0; 
						double pastlambda=0.0, pastbincoef=0.0; 
						double binCoef_ =0.0, SumRates_ =0.0, lambda_=0.0, t_ =0.0; 
						List<Edge> tmpActiveEdge = new ArrayList<Edge>();  
						List<String> singlekeys = new ArrayList<String>(); 
						List<Edge> saveEdge = new ArrayList<Edge>();
						ArrayList<Integer> saveEC = new ArrayList<Integer>(); 

						Map <String, ArrayList<Double>> timeContain = new HashMap<String, ArrayList<Double>>();
						for(Map.Entry<String, List<Edge>> entry: MapActiveEdge.entrySet()) {

							if(edgeContain.get(entry.getKey()).size()==0)
								continue;

							String Key = entry.getKey(); 
							System.out.println("Key is: "+Key);
							tmpActiveEdge = entry.getValue();
							System.out.println("the map active edges size is : "+tmpActiveEdge.size());
							System.out.println("Edge Contained here is: "+edgeContain.get(entry.getKey()).size());
							//Do this for any key which has at least 1 edge, as it can be picked up for recombination. If index is 0 for 
							//this key then we force it to recombine.
							if(edgeContain.get(entry.getKey()).size() > 1) {
								binCoef_ = BinomialCoefficient.binomialCoeff(GetL(Key), 2);
								SumRates_ = ComputeSumRates(tmpActiveEdge);
								lambda_ = binCoef_ + SumRates_; 						
								ArrayList<Double> tmList = new ArrayList<Double>();
								tmList.add(lambda_);
								tmList.add(binCoef_);
								t_ = GetEffPop(Key)*computeNextTimeEvent(lambda_);
								tmList.add(t_);
								timeContain.put(Key, tmList);

								//t_ = computeNextTimeEvent(lambda_);
								System.out.println("Time generated is : "+t_+" bincoef is "+binCoef_+" Sum Rates is: "+SumRates_);
								if (tpast == 0 || t_ > tpast) { //Take overall max
									tpast = t_;  
									ipast = Key; 
									pastlambda = lambda_; 
									pastbincoef = binCoef_;
								}
								

								if (tpast<0)
									System.out.println("stop");
							}
							//If there is a key such as s1 or s2, which has at most one edge in the container 
							//Then check if there are lineages which contain that key and is still active. 
							//If they are not active then assign this edge to Neutral, 99 case and delete the key. 
							else if(StrTrun(Key).length() == 1 && edgeContain.get(Key).size() < 2) {
								System.out.println("Key is : "+entry.getKey() + " containing "+edgeContain.get(entry.getKey()).size()+ " edges");
								String tmpstr = StrTrun(Key);
								int sinflag = 0; 
								for (int ij = 0; ij < tmpkeys.size(); ij++) {
									if(StrTrun(tmpkeys.get(ij)).contains(tmpstr) && StrTrun(tmpkeys.get(ij)).length() > 1
											&& edgeContain.get(tmpkeys.get(ij)).size() > 1){
										System.out.println("it has its friend: "+tmpkeys.get(ij));
										sinflag = 1; 
										break;
									}
								}
								//System.out.println("sinflag is : "+sinflag+" ");
								if(sinflag == 0 && !edgeContain.get(Key).isEmpty()) {
									ArrayList<Integer> tmpid = edgeContain.get(entry.getKey());
									//edgeContain.get(NS).addAll(tmpid);
									saveEC = tmpid; 
									singlekeys.add(Key);
									toRemKey = Key; 
									List<Edge> tmped = MapActiveEdge.get(Key);
									saveEdge = tmped; 
									for( int mk = 0; mk < tmped.size(); mk++) {
										MutMap.add(tmped.get(mk));
									}
									
									tmpkeys.remove(Key);
								}
							}
							else {
								if(StrTrun(Key).length() > 1 && !Key.equals(NS) && edgeContain.get(Key).size() < 2) {
									ArrayList<Integer> tmpid = edgeContain.get(entry.getKey());
									saveEC = tmpid; 
									saveEdge = tmpActiveEdge; 
									toRemKey = Key;
									
									tmpkeys.remove(Key);
								}
	
								if(ipast.equals(NS) && !edgeContain.get(NS).isEmpty()) {
									if (GetL(NS)==2)
										System.out.println("stop");
									binCoef_ = BinomialCoefficient.binomialCoeff(GetL(NS), 2);
									SumRates_ = ComputeSumRates(MapActiveEdge.get(NS));
									lambda_ = binCoef_ + SumRates_; 
									tpast = GetEffPop(NS)*computeNextTimeEvent(lambda_);
									//tpast = computeNextTimeEvent(lambda_);
									pastbincoef = binCoef_;
									pastlambda = lambda_; 
									ipast = NS;
									ArrayList<Double> tmList = new ArrayList<Double>();
									tmList.add(lambda_);
									tmList.add(binCoef_);
									tmList.add(tpast);
									
									if (tpast<0)
										System.out.println("stop");

									timeContain.put(ipast, tmList);
								}
								else 
									continue;
								
							}

						}


						if(!toRemKey.isEmpty()) {
							System.out.println("Size of Save Edge: "+saveEdge.size());
							MapActiveEdge.get(NS).addAll(saveEdge);
							edgeContain.get(NS).addAll(saveEC);
							for(int kk = 0; kk < saveEdge.size(); kk++) {
								saveEdge.remove(saveEdge.get(kk));
							}
							MapActiveEdge.remove(toRemKey);
							MapActiveEdge.put(toRemKey, saveEdge);
							for(int kk = 0; kk < saveEC.size(); kk++) {
								saveEC.remove(saveEC.get(kk));
							}
							edgeContain.remove(toRemKey);
							edgeContain.put(toRemKey, saveEC);
						}

						if (tpast == 0.0 && pastlambda == 0.0 && pastbincoef == 0) {
							continue;
						}
						System.out.println("Time is: "+tpast+" ipast is "+ipast+" lambda is "+pastlambda+" bincoef is "+pastbincoef);
						System.out.println(" size of container: "+edgeContain.get(ipast).size());
						
						//5) Update the current generation
						setGeneration(computeNextGeneration( timeContain.get(ipast).get(2))); //ipast,

						
						try {
							fwritepick.write(getGeneration()+"\t"+ipast+"\t"+edgeContain.get(ipast).size()+"\t"+numRun+"\n");
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						//}*/
						

						for(Map.Entry<String, List<Edge>> entryn: MapActiveEdge.entrySet()) {


							System.out.println("Key: "+entryn.getKey()+" size: "+edgeContain.get(entryn.getKey()).size());

							if(!timeContain.containsKey(entryn.getKey()))
								continue;

							ipast = entryn.getKey();
							System.out.println(ipast);
							System.out.println(timeContain);
							pastlambda = timeContain.get(ipast).get(0);
							pastbincoef = timeContain.get(ipast).get(1);

							int index_=-1; 
							//System.out.println("Lambda is: "+lambda_+" Bin Coef is "+binCoef_+"  ipast is "+ipast);
							//ARG with recombination nodes

							if(getRecomb() > 0)
								index_ = computeIndexNextEvent(pastlambda, pastbincoef,edgeContain.get(ipast));
							else 
								index_ = 0; 

							if(edgeContain.get(ipast).size() < 2)
								if (!singlekeys.isEmpty()) {
									for (int ij = 0; ij < singlekeys.size(); ij++) {
										System.out.println("KEY to REM is "+singlekeys.get(ij)+" and edge size is: "+MapActiveEdge.get(singlekeys.get(ij)).size());
										MapActiveEdge.get(NS).addAll(MapActiveEdge.get(singlekeys.get(ij)));
										for (int jj = 0; jj < MapActiveEdge.get(singlekeys.get(ij)).size(); jj++)
											MapActiveEdge.get(singlekeys.get(ij)).remove(jj);
										MapActiveEdge.remove(singlekeys.get(ij));
										MapActiveEdge.put(singlekeys.get(ij),MapActiveEdge.get(singlekeys.get(ij)));		
									}
								}

							if (edgeContain.get(ipast).size() < 2 && getRecomb() > 0)
								index_ = 1; 

							int index = index_;
							Random random = new Random();
							//String randKey ="";
							System.out.println("-------------------------------");
							System.out.println("Index is "+index);
							System.out.println("-------------------------------");
							//***** CASE 1: THE NEXT EVENT IS COALESCENT *****
							if(index == 0) {

								/*for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) {
									System.out.println("Key: "+entry.getKey()+" and size is "+entry.getValue().size());
								}*/
								coalescentNumber = coalescentNumber+1;

								//Create a new Coalescent node
								Node coal_node = new Node(nodeSet.size(), false, getGeneration());
								nodeSet.put(coal_node.getId(), coal_node);

								//Pick randomly the first active edge from the set of active edges
								int j; int id_first_edge;

								int pos_j = random.nextInt(edgeContain.get(ipast).size());
								j = activeEdges.indexOf(graphEdges.get(edgeContain.get(ipast).get(pos_j)));
								edgeContain.get(ipast).remove(pos_j);
								System.out.println("POS J IS "+pos_j);

								id_first_edge = activeEdges.get(j).getIdEDGE();


								//Set the father of the active edge selected
								graphEdges.get(id_first_edge).setId_fath(coal_node.getId());
								//Set the time passed for this event (computed by the exponential distribution)
								graphEdges.get(id_first_edge).setTime();
								//System.out.println("VALUE: "+edgeContain.get(randKey).size()+" and "+randKey);

								//to understand how many consecutive coalescent and recombination we have (with the same two lineages)
								int son1 = graphEdges.get(id_first_edge).getId_son();

								//Get the arraylist of Intervals/Segments of the left edge
								ArrayList<Interval> left_segments = activeEdges.get(j).getSegments();

								//Remove the lineage from the set of the active ones
								activeEdges.remove(j);
								MapEdgeRemove(ipast, id_first_edge);

								//Pick randomly the second active edge from the set of active edges

								pos_j = random.nextInt(edgeContain.get(ipast).size());
								j = activeEdges.indexOf(graphEdges.get(edgeContain.get(ipast).get(pos_j)));
								edgeContain.get(ipast).remove(pos_j);	
								System.out.println("POS J IS "+pos_j);
								int id_second_edge = activeEdges.get(j).getIdEDGE();
								int son2 = graphEdges.get(id_second_edge).getId_son();
								if(son1 == son2)
									numberRECandCoal = numberRECandCoal+1;

								//Set the father of the selected active lineage as coalescent node
								graphEdges.get(id_second_edge).setId_fath(coal_node.getId());
								//Set the time passed for this event (computed by the exponential distribution)
								graphEdges.get(id_second_edge).setTime();

								//Get the arraylist of Intervals/Segments of the right edge
								ArrayList<Interval> right_segments = activeEdges.get(j).getSegments();

								//Remove the lineage from the set of active ones
								activeEdges.remove(j);
								MapEdgeRemove(ipast, id_second_edge);
								//Compute the union of the left and the right intervals
								MergeIntervals merge = new MergeIntervals();
								ArrayList<Interval> union = merge.ConcatTwoIntervalList(left_segments, right_segments);
								union = merge.merge(union);

								//Memorize the new set of segments in the coalescent node and the two sons
								coal_node.setSegments(union);
								coal_node.setIDsonsx(son1);
								coal_node.setIDsondx(son2);

								//Create a new active edge outgoing from the new coalescent node
								//Edge new_edge = new Edge(graphEdges.size(), coal_node.getId(), -1, union);
								Edge new_edge;

								if(edgeContain.get(ipast).isEmpty() && StrTrun(ipast).length() > 1 && !StrTrun(ipast).equals(StrTrun(NS))) { 
									String tmpstr1 = "", tmpstr2="";
									String tmpstr = StrTrun(ipast);
									if (tmpstr.length() == 2) {
										tmpstr1 = "["+Character.toString(tmpstr.charAt(0))+"]";
										tmpstr2 = "["+Character.toString(tmpstr.charAt(1))+"]";
									}
									if(tmpstr.length() == 3) {
										double r = random.nextFloat();
										if (r < 0.33) {
											tmpstr1 = "["+Character.toString(tmpstr.charAt(0))+"]";
											tmpstr2 = "["+Character.toString(tmpstr.charAt(1))+", "+Character.toString(tmpstr.charAt(2))+"]";
										}
										else if (r > 0.33 && r < 0.66){
											tmpstr2 = "["+Character.toString(tmpstr.charAt(2))+"]";
											tmpstr1 = "["+Character.toString(tmpstr.charAt(0))+", "+Character.toString(tmpstr.charAt(1))+"]";
										}
										else {
											tmpstr2 = "["+Character.toString(tmpstr.charAt(1))+"]";
											tmpstr1 = "["+Character.toString(tmpstr.charAt(0))+", "+Character.toString(tmpstr.charAt(2))+"]";
										}
									}
									System.out.println("String 1 is: "+tmpstr1+" and String 2 is: "+tmpstr2);

							
									if (argRecomb > 0) {
										recombinationsNumber = recombinationsNumber+1;

										SetEdge(coal_node,union,tmpstr1,edgeContain);
										SetEdge(coal_node,union,tmpstr2,edgeContain);
									}
									else {
										SetEdge(coal_node,union,tmpstr1,edgeContain);
										new_edge = new Edge(graphEdges.size(), coal_node.getId(), -1, union, GetEffPop(tmpstr2));
										new_edge.computeLength();
										new_edge.computeRate(GetEffPop(tmpstr2), getRecomb(), getG());
										edgeContain.get(tmpstr1).add(new_edge.getIdEDGE());
										graphEdges.put(new_edge.getIdEDGE(), new_edge);
										activeEdges.add(new_edge);
										MapActiveEdge.get(tmpstr1).add(new_edge);
									}
									tmpkeys.remove(ipast);
								}
								else if(edgeContain.get(ipast).isEmpty() && StrTrun(ipast).length() == 1) {
									System.out.println("I entered here with "+ipast);
									String tmpstr = StrTrun(ipast);
									for (int ij = 0; ij < tmpkeys.size(); ij++) {
										if(StrTrun(tmpkeys.get(ij)).contains(tmpstr) && StrTrun(tmpkeys.get(ij)).length() > 1
												&& edgeContain.get(tmpkeys.get(ij)).size() > 0){
											System.out.println("it has its friend: "+tmpkeys.get(ij));
											SetEdge(coal_node,union,ipast,edgeContain);
											System.out.println("after size is: "+edgeContain.get(ipast).size());
											flag = 1; 
											break;
										}	
									}
									System.out.println("flag is set to : "+flag);
									if(flag == 0) {
										new_edge = new Edge(graphEdges.size(), coal_node.getId(), -1, union, GetEffPop(NS));
										new_edge.computeLength();
										new_edge.computeRate(GetEffPop(NS), getRecomb(), getG());
										edgeContain.get(NS).add(new_edge.getIdEDGE());
										graphEdges.put(new_edge.getIdEDGE(), new_edge);
										activeEdges.add(new_edge);
										MapActiveEdge.get(NS).add(new_edge);
										MutMap.add(new_edge);
										System.out.println("after size is: "+edgeContain.get(NS).size());
										tmpkeys.remove(ipast);
									}

								}
								else { 
									SetEdge(coal_node,union,ipast,edgeContain);
								}

							}

							//***** CASE 2: THE NEXT EVENT IS RECOMBINATION *****
							else{

								//Increment the number of recombinations events
								recombinationsNumber = recombinationsNumber+1;

								//Create a new recombination node
								Node recomb_node = new Node(nodeSet.size(), true, getGeneration());
								nodeSet.put(recomb_node.getId(), recomb_node);

								//get the ID of the active edge selected
								int id_edge;

								id_edge = edgeContain.get(ipast).get(index-1); 
								edgeContain.get(ipast).remove(index-1);

								//Set the father of the active edge as the ID of the new recomb node 
								graphEdges.get(id_edge).setId_fath(recomb_node.getId());

								//Set the time computed by exponential distribution
								graphEdges.get(id_edge).setTime();

								//Set the list of solid intervals
								int pos_index= activeEdges.indexOf(graphEdges.get(id_edge));
								recomb_node.setSegments(activeEdges.get(pos_index).getSegments());

								//Remove the picked edge from the active ones
								activeEdges.remove(pos_index);
								MapEdgeRemove(ipast, id_edge);
								//Create two new active edges that have as son the new recombination node
								ArrayList<ArrayList<Interval>> splitted = SplittingIntervals.split(recomb_node.getSegments(), recomb_node.getId());
						
								int recflag = 0; 
								if(!(ipast.equals(NS)))
									for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) {
										if(!entry.getKey().equals(NS)) {
											if (entry.getValue().size() > 0) {
												recflag = 1;
												break;
											}	
										}
									}
								String leftkey = "", rightkey = ""; 
								if (recflag == 1) {
									double pos = SplittingIntervals.getRec_pos();
									//System.out.println("pos is : "+pos);
									for(int i=0; i < posset.length-1; i++) {	
										if(pos >= posset[i] && pos < posset[i+1]) {
											for(Map.Entry<String, Double> entry: SegPos.entrySet()) {
												if (entry.getValue() == posset[i])
													leftkey = entry.getKey();
												if(entry.getValue() == posset[i+1])
													rightkey = entry.getKey();

											}
										}
									}
									System.out.println("leftkey is "+leftkey+" and rightkey is "+rightkey);
									if(!(tmpkeys.contains(leftkey)))
										leftkey = NS;
									if(!(tmpkeys.contains(rightkey)))
										rightkey = NS;
								}
								else {
									leftkey = NS;
									rightkey = NS; 
								}

								System.out.println("leftkey is "+leftkey+" and rightkey is "+rightkey);
								//Edge left_edge;
								
								SetEdge(recomb_node,splitted.get(0),leftkey,edgeContain);
								
								SetEdge(recomb_node,splitted.get(0),rightkey,edgeContain);
								
								int remflag = 0; 	
								String tmpstr = StrTrun(ipast);
								for (int ij = 0; ij < tmpkeys.size(); ij++) {
									if(StrTrun(tmpkeys.get(ij)).contains(tmpstr) && StrTrun(tmpkeys.get(ij)).length() > 1
											&& edgeContain.get(tmpkeys.get(ij)).size() > 0){
										System.out.println("it has its friend: "+tmpkeys.get(ij));
										remflag = 1; 
									}
								}
								if(remflag == 0 && edgeContain.get(ipast).size() < 1) {
									tmpkeys.remove(ipast);
								}
							}

							//update L
							setL(activeEdges.size());
							for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) 
								SetL(entry.getKey(),entry.getValue().size());

						}
					} 

					if(activeEdges.size()==1)
						System.out.println("FOUND GMRCA");
					else
						System.out.println("Last generation "+getGeneration()+ "Living Lineages "+getL());
					System.out.println("");
					System.out.println("Size of MutMap is: "+MutMap.size());

					//}
					//***************************** CHECKING SOME PROPERTIES   **********************************
					
					double computed_r = getRecombinationsNumber()/(GetEffPop(NS)*getG()*getGeneration());
					System.out.println(getRecomb()+ " ~ " +computed_r);

					

					//Matrix graph useful for visit the graph and decorate it with mutations SNP and STRs
					DiGraph graph = new DiGraph();
					graph.createDigraph();




					//***************************** DECORATING ARG WITH SNP MUTATIONS   ***********************************



					int totalMutations = DecorateMutationsSNP.computeTotalMutationsNumber();
					System.out.println("Total Number of mutations Y ---> "+totalMutations);
					//System.out.println("*** PROPERTY (26):  mutation rate ~ (Y/ gNT ) ***");
					double computed_mu = totalMutations/(getG()*GetEffPop(NS)*getGeneration());
					System.out.println(getMu()+ " ~ " +computed_mu);
					//setScaledGeneration(getGeneration()*GetEffPop(NS));
					setScaledGeneration(getGeneration());
					///System.exit(1);
					/*
			

			//***************************** DECORATING ARG WITH STRs MUTATIONS  ***********************************

					/*** stats ***/
					try {
						PrintWriter pw = new PrintWriter(new FileOutputStream(
								new File(filestast), 
								true /* append = true */));

						pw.println(getGeneration()+"\t"+getRecombinationsNumber()+"\t"+totalMutations);
						pw.close();
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} 


					//***************************** CREATING FILE STRUCTURES ***********************************
					
					String txtFileL = filename.concat("_L"+".txt");
					CreatingFilesForStructure.createLtxt(wholePath, txtFileL, getGraphEdges(), getNodeSet());

					String txtFileS = filename.concat("_S"+".txt");
					CreatingFilesForStructure.createStxt(wholePath, txtFileS, getGraphEdges(), getNodeSet());


					try {
						fwritepick.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

				} //end numRun

				try {fwrite.close();} catch (Exception ex) {}

				long end = System.currentTimeMillis();
				NumberFormat formatter = new DecimalFormat("#0.00000");
				System.out.print(formatter.format((end - start) / 1000d));
				int mb = 1024;
				runtime = Runtime.getRuntime();
				System.out.println("\t"+(runtime.totalMemory() - runtime.freeMemory()) / mb);
			}
		}
		//	}

	} //end main procedure


	public static void main_wait(String[]args){
		int argN;
		int argm;
		double argRecomb;
		double argMut;
		int argG;
		String argPath;
		String argFileName;
		int epiflag = 1; 
		/**    

		if(args.length < 5) {
			//fare eccezione
			errorInsertingParameters();
		}
		 */
		Runtime runtime = Runtime.getRuntime();
		long start = System.currentTimeMillis();

		/*argN = Integer.parseInt(args[0]);
		argm = Integer.parseInt(args[1]);
		argRecomb = Double.parseDouble(args[2]);
		argMut = Double.parseDouble(args[3]);
		argG = Integer.parseInt(args[4]);
		argPath = args[5].toString();
		argFileName = args[6].toString();*/

		int [] mvect=new int[] {250}; // {20, 50, 80, 120};
		//	int [] mvect=new int[] {50};
		int [] gvect=new int[] {25};
		//double [] fvect=new double[]{0.3}; //, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
		double[] OldEpiFit = new double[] {-0.9};
		double [] AllEpiFit = OldEpiFit;	//List of the s1, s2, s3.... values
		//double[] AllEpiFit = new double[] {};
		List<Double> Epifit = new ArrayList<Double>();
		double[] FunFit = new double[AllEpiFit.length];			//Stores the f(s) values

		for (int i = 0; i < AllEpiFit.length; i++) {
			Epifit.add(AllEpiFit[i]);
			FunFit[i] = FunS(AllEpiFit[i]);
		}

		String NS = "[99]";		//Immortal key! 

		Map <String, Double> FitGroup = new HashMap<String,Double>(); 	
		Map <String, Double> EffPop = new HashMap<String,Double>(); 	//This will store f' to calculate N*f'

		//int EpiFlag = 0; //Flag to indicate whether epistasis is ON? We don't need this anymore ?
		int k_way = 1;	//how many maximum interacting SNPs 

		// This map indicates for which epistasis will be in effect with the delta associated 
		Map <String, Double> Delta = new HashMap<String,Double>();
		//Delta.put("[0, 1]",2.0);
		//Delta.put("[0]", 0.12);
		//Delta.put("[1]", 0.12);*/
		//Delta.put("[0, 1, 2]", 0.05);

		//Calculating f' for each key
		FitGroup = ComputeFitGroup(AllEpiFit, k_way, Delta, epiflag); 
		
		FitGroup.put(NS, 1.0/(2.0+OldEpiFit[0])); //HACK!!!
		
		
		EffPop = ComputeEffPop(FitGroup,epiflag);

		System.out.println("---- Fitness Values ----");
		for(Map.Entry<String, Double> entry: FitGroup.entrySet()) {
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}
		System.out.println("------------------------");
		List<String> keys = new ArrayList<String>(FitGroup.keySet());		//keep a record of the keys being used

		for(int im=0;im<mvect.length; im++)
		{
			for(int indg=0; indg<gvect.length; indg++)
			{
				int runs=100;
				setTotalRuns(runs);
				argN = 1000;
				N = argN;
				argm = mvect[im];
				argRecomb = 0;
				argMut = 1.5;
				argG = gvect[indg];
				argPath = "stats/";
				argFileName = "prova";
				String filestast="stats/N"+String.valueOf(argN)+"_m"+String.valueOf(argm)+"_m1_g"+String.valueOf(argG)+"_s"+"_r"+String.valueOf(argRecomb)+".stats";

				//***************************** PARAMETER INITIALIZATION  ***************************************
				//Filename Format N#_g#_mu#_r#_m#
				//filename = "N"+getN()+"_g"+getG()+"_mu"+getMu()+"_ro"+getRecomb()+"_m"+getExtantUnits();
				filename = filestast;

				//***** Recombination Rate ****
				setRecomb(argRecomb*Math.pow(10, -8));

				//**** SNP Mutation Rate ****
				mu = argMut*Math.pow(10, -8);

				//**** Segment length g ***
				setG(argG*Math.pow(10, 3));

				//**** whole path of directory for outputfile ***
				//String wholePath = "output_SimRa/";
				String wholePath = argPath;

				//***************************** CREATION OF DIRECTORY FOR OUTPUT FILES  ***********************************
				File dir;
				dir = new File(""+wholePath);
				dir.mkdir();		
				BufferedWriter fwrite = null;
				File file = new File(""+wholePath+filename+"_STATS.txt");
				try {
					fwrite = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"));
				} catch (UnsupportedEncodingException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

				//***************************** ARG RANDOM GENERATION ALGORITHM   ***********************************

				Random rand = new Random();

				for(int numRun = 1; numRun <= getTotalRuns(); numRun++) {
					//********* INITIALIZATION OF STATIC VARIABLES **************/
					setGeneration(0); // This may go away
					setRecombinationsNumber(0);
					setMutationsNumber(0);
					setCoalescentNumber(0);
					numberRECandCoal = 0;
					nodeSet = new HashMap<Integer,Node>();
					mutationSet = new HashMap<Double,Mutation>();
					SNPpositionsList = new TreeSet<Double>();
					graphEdges = new HashMap<Integer,Edge>();
					activeEdges = new ArrayList<Edge>();
					splitPoints = new TreeSet<Double>();
					splitPoints.add(0.0);
					splitPoints.add(1.0);

					//********* START
					ArrayList<String> fitMidx = new ArrayList<String>();		//this will contain the key for each argm
					int[] idxm = new int[argm];
					List<Edge> MutMap = new ArrayList<Edge>();

					Map<String, ArrayList<Integer>> IdxGroups = new HashMap<String, ArrayList<Integer>>();
					ArrayList<Integer> immortal = new ArrayList<Integer>(); 

					for (int i = 0; i < argm; i++) {
						idxm[i] = i; 
						immortal.add(i);
					}
					Map<String,Integer> selsize = new HashMap<String,Integer>(); 
					for (int i = 0; i < FunFit.length; i++) {
						ArrayList<Integer> tmpid = new ArrayList<Integer>();
						tmpid.add(i);
						selsize.put(tmpid.toString(),(int)Math.round(FunFit[i]*argm));  //
					}

					/** Assignment of each single keys such as s1, s2, s3, etc. according to their respective f(s) **/
					if(argRecomb > 0) {

						for(int i = 0; i < FunFit.length; i++) {
							shuffleArray(idxm);
							ArrayList<Integer> tmpid = new ArrayList<Integer>();
							ArrayList<Integer> idx = new ArrayList<Integer>(); 

							//int  thres = rand.nextInt(argm) + 1;
							int thres = (int)Math.round(FunFit[i]*argm);
							for (int j = 0; j < thres; j++) {
								tmpid.add(idxm[j]);
							}
							idx.add(i);
							IdxGroups.put(idx.toString(), tmpid);
						}
					}
					else {
						Map<Integer, Integer> sizeMap = new TreeMap<Integer, Integer>();
						for(int i = 0; i < FunFit.length; i++) 
							sizeMap.put(i,(int)Math.round(FunFit[i]*argm)); //
						
						
						LinkedHashMap<Integer, Integer> reverseSortedMap = new LinkedHashMap<Integer, Integer>();
						//Use Comparator.reverseOrder() for reverse ordering
						sizeMap.entrySet().stream().sorted(Map.Entry.comparingByValue(Comparator.reverseOrder())).forEachOrdered(x -> reverseSortedMap.put(x.getKey(), x.getValue()));

						Map<String, Integer> availablesizeMap = new TreeMap<String, Integer>();
						availablesizeMap.put("[99]",immortal.size());
						

						shuffleArray(idxm);
						ArrayList<Integer> immcopy = immortal;
						Boolean first = true;
						for (Map.Entry<Integer, Integer> entry : reverseSortedMap.entrySet()) {
							ArrayList<Integer> tmpid = new ArrayList<Integer>();

							if (first) {
								String KeySel = "["+Integer.toString(entry.getKey())+"]";
								Collections.shuffle(immcopy);
								for (int j = 0; j < entry.getValue(); j++) 
									tmpid.add(immcopy.get(j));
								immcopy.removeAll(tmpid);
								IdxGroups.put(KeySel, tmpid);
								availablesizeMap.put("[99]", immcopy.size());
								IdxGroups.put("[99]", immcopy);
								availablesizeMap.put(KeySel, entry.getValue());
								first = false;
							}
							else
							{
								List<String> Tkeys = new ArrayList<String>(availablesizeMap.keySet());

								Collections.shuffle(Tkeys);
								for (String it: Tkeys)
								{
									if (entry.getValue() <= availablesizeMap.get(it))
									{
										String KeySel;
										if (it.equals("[99]"))
										{
											KeySel = "["+Integer.toString(entry.getKey())+"]";
											for (int j = 0; j < entry.getValue(); j++) 
												tmpid.add(immcopy.get(j));
											immcopy.removeAll(tmpid);
										}
										else
										{
											ArrayList<Integer> tlist = IdxGroups.get(it);
											for (int j = 0; j < entry.getValue(); j++) 
												tmpid.add(tlist.get(j));
											for (int j = entry.getValue()-1; j>=0; j--) 
												tlist.remove(j);
											IdxGroups.put(it, tlist);
											List<String> items = new ArrayList<String>(Arrays.asList(it.replace("[", "").replaceAll("]", "").split(", ")));
											items.add(Integer.toString(entry.getKey()));
											Collections.sort(items);

											KeySel = "[";
											for (String i:items)
												KeySel =KeySel + i + ", ";

											KeySel = KeySel.substring(0, KeySel.length()-2) + "]";	
										}

										availablesizeMap.put(it, availablesizeMap.get(it)-entry.getValue());
										availablesizeMap.put(KeySel, entry.getValue());
										IdxGroups.put(KeySel, tmpid);
										break;
									}
								}
							}
						}
						
						
						//FLIP ****
						 ArrayList<Integer> tmp = IdxGroups.get(NS);
						 IdxGroups.put(NS, IdxGroups.get("[0]"));
						 IdxGroups.put("[0]",tmp);
						
						//FLIP ****
						
					}
					/** END: Assignment of each single keys such as s1, s2, s3, etc. according to their respective f(s) **/

					for(Map.Entry<String, ArrayList<Integer>> entry: IdxGroups.entrySet()) {
						System.out.print("|| "+entry.getKey()+ " : ");
						for (int i = 0; i < entry.getValue().size(); i++)
							System.out.print(entry.getValue().get(i)+", ");
					}
					System.out.println("");	

					Map<String, ArrayList<Integer>> FinSel;
					if(argRecomb > 0) {

						//Storage of the selection groups containing keys as strings and the indices assigned as integers. 
						//Note: This will contain redundancies and we have to remove them to get the final selection Map. 
						//Assignment of final selection map
						FinSel = new HashMap<String, ArrayList<Integer>>();
						SelectGroups(IdxGroups, FinSel, keys, immortal, NS);
						for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
							System.out.print("|| "+entry.getKey()+ " : ");
							for (int i = 0; i < entry.getValue().size(); i++)
								System.out.print(entry.getValue().get(i)+", ");
						}
						System.out.println("");	
						for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
							List<Integer> tmplist = entry.getValue();
							for(Map.Entry<String, ArrayList<Integer>> entry2: FinSel.entrySet()) {
								for(int ik = 0; ik < tmplist.size(); ik++) {
									if(!entry2.equals(entry) && entry2.getValue().contains(tmplist.get(ik))) {
										if(StrTrun(entry2.getKey()).length() < StrTrun(entry.getKey()).length())
											entry2.getValue().remove(Integer.valueOf(tmplist.get(ik)));
										else
											entry.getValue().remove(Integer.valueOf(tmplist.get(ik)));
									}
								}
							}
						}
					}
					else
					{
						FinSel = IdxGroups;
					}

					for (int ij = 0; ij < argm; ij++) {
						for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
							if(entry.getValue().contains(ij))
								fitMidx.add(entry.getKey());
						}
					}

					for(Map.Entry<String, ArrayList<Integer>> entry: FinSel.entrySet()) {
						System.out.print("|| "+entry.getKey()+ " : ");
						for (int i = 0; i < entry.getValue().size(); i++)
							System.out.print(entry.getValue().get(i)+", ");
					}

					int immnum = FinSel.get(NS).size(); 

					setExtantUnits(immnum); //m = input value - m'
					setExtantUnitsunderSelection(argm - immnum);

					System.out.println("---- Eff Pop ----");
					for(Map.Entry<String, Double> entry: EffPop.entrySet()) {
						System.out.println(entry.getKey() + " : " + entry.getValue());
						SetEffPop(entry.getKey(), (int) Math.round(argN*entry.getValue()));
					}

					//For each group we inilizatialize the time
					for(Map.Entry<String, Double> entry: EffPop.entrySet()) {
						setGeneration(entry.getKey(), 0);
					}

					Map<String, Double> SegPos = new HashMap<String, Double>();
					//Assign segments positions to keys for recombination
					double[] posset= new double[keys.size()]; 
					for (int i = 0; i < keys.size(); i++) {
						if(i == 0) {
							SegPos.put(keys.get(i), 0.0);
							posset[i] = 0;
						}
						else if(i == keys.size()-1) {
							SegPos.put(keys.get(i), 1.0);
							posset[i] = 1; 
						}
						else {
							double val = Math.random();
							SegPos.put(keys.get(i), val);
							posset[i] = val;
						}
					}
					Arrays.sort(posset);

					System.out.println("------------Segments-------------");
					for(Map.Entry<String, Double> entry: SegPos.entrySet()) {
						System.out.print("|| "+entry.getKey()+ " : "+entry.getValue());
					}
					System.out.println("");

					Map <String, ArrayList<Integer>> edgeContain = new HashMap<String, ArrayList<Integer>>(); //edgeContainer for each key
					Map <String, ArrayList<Edge>> edgeWaitContain = new HashMap<String, ArrayList<Edge>>(); //edgeContainer for each key
					Map <Double, Edge> edgeTimeWait = new HashMap <Double, Edge>(); 

					MapActiveEdge.clear();
					Epifit.clear();
					MutMap.clear();
					List<Integer> DelFitIdx = new ArrayList<Integer>();
					AllEpiFit = OldEpiFit;
					for (int ik = 0; ik < AllEpiFit.length; ik++) {
						Epifit.add(AllEpiFit[ik]);
					}


					//********* INITIALIZATION OF LEAVE NODES ****************/

					for (int i = 0; i < argm; i++) {

						Node n = new Node(i,false,getGeneration()); 
						ArrayList<Interval> segments = new ArrayList<Interval>();
						segments.add(new Interval(0,1));
						n.setSegments(segments);
						nodeSet.put(i, n);
						ArrayList<Integer> tmparr; 
						List<Edge> tmpActiveEdge;
						Edge e = new Edge(i,n.getId(), -1, segments, true, GetEffPop(fitMidx.get(i)));
						e.computeLength();
						e.computeDensity();
						e.computeRate(GetEffPop(fitMidx.get(i)), getRecomb(), getG());

						if (!edgeContain.containsKey(fitMidx.get(i))){
							tmparr = new ArrayList<Integer>(); 
						}
						else {
							tmparr = edgeContain.get(fitMidx.get(i)); 	
						}
						tmparr.add(i);
						edgeContain.put(fitMidx.get(i), tmparr); 

						if (GetMapActiveEdge(fitMidx.get(i)) == null) {
							 tmpActiveEdge = new ArrayList<Edge>();
						}
						else {
							tmpActiveEdge = GetMapActiveEdge(fitMidx.get(i)); 
						}
						
						tmpActiveEdge.add(e);
						SetMapActiveEdge(fitMidx.get(i),tmpActiveEdge);	
						
						//Add the new lineage to the list of active lineages
						activeEdges.add(e);
						graphEdges.put(e.getIdEDGE(), e);
					}
					//Populate edgeContainer and MapActiveEdges
					for (int kj = 0; kj < keys.size(); kj++) {
						if (!edgeContain.containsKey(keys.get(kj)))
						 {
							System.out.println("Doesn't contain: "+keys.get(kj));
							edgeContain.put(keys.get(kj),  new ArrayList<Integer>());
						}
						if (!MapActiveEdge.containsKey(keys.get(kj)))
							MapActiveEdge.put(keys.get(kj), new ArrayList<Edge>());
					}

					setL(activeEdges.size());

					for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) 
						SetL(entry.getKey(),entry.getValue().size());

					//*************** Construction of the ARG backwards starting from the leaves ****************/
					while(activeEdges.size() !=1){
						
						
						if ((activeEdges.size() % 10)==0)
						{
							
							for(Map.Entry<String, List<Edge>> entry: MapActiveEdge.entrySet()) {
								if (entry.getValue().size() != edgeContain.get(entry.getKey()).size())
									System.out.println(entry.getKey()+" "+entry.getValue().size()+" "+edgeContain.get(entry.getKey()).size());
							}
							
							System.out.println("stop");
						}
						System.out.println("----  edges  ----");

						for(Map.Entry<String, ArrayList<Integer>> entry: edgeContain.entrySet()) {	
							System.out.println("Key: "+entry.getKey());
							for (int kl = 0; kl < entry.getValue().size(); kl++)
								System.out.print(" "+entry.getValue().get(kl));
							System.out.println("");
						}
						System.out.println("-------------------");
						System.out.println("There are "+activeEdges.size()+" active edges now!!");

						//Each key run at every iteration
						for(Map.Entry<String, List<Edge>> entry: MapActiveEdge.entrySet()) {
							
							//if (entry.getKey().equals(NS) & (edgeContain.get("[0]").size()>2))
							//		continue;
							
							if((edgeContain.get(entry.getKey()).size()==0) & !edgeWaitContain.containsKey(entry.getKey()))
								continue;
							else
								if((edgeContain.get(entry.getKey()).size()==0) & edgeWaitContain.containsKey(entry.getKey()))
								{
									//move all the earliest element in wait
									System.out.println("CASE ARE WAITING");
								}
								
							
							

							String Key = entry.getKey(); 
							List<Edge> tmpActiveEdge = entry.getValue();
							
							if((edgeContain.get(Key).size()==1) & !edgeWaitContain.containsKey(Key) & Key.equals(NS))
								continue;
								

							double binCoef_ = BinomialCoefficient.binomialCoeff(GetL(Key), 2);
							double SumRates_ = ComputeSumRates(tmpActiveEdge);
							double lambda_ = binCoef_ + SumRates_; 						
							double t_ = (GetEffPop(Key)/((double)argN))*computeNextTimeEvent(lambda_);
							
							setGeneration(Key, computeNextGeneration(Key, t_));
							
							try {
								PrintWriter pw = new PrintWriter(new FileOutputStream(
										new File("stats/time.txt"), 
										true /* append = true */));
								
								pw.println(activeEdges.size()+"\t"+entry.getKey()+"\t"+getcurrGeneration(Key)+"\t"+edgeContain.get(Key).size());
								pw.close();
							} catch (FileNotFoundException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
							

							//check if the now we can take someone out of the wait list
							Boolean delKey = false;
							ArrayList<Double> tdel  = new ArrayList<Double>();
							for(Map.Entry<Double, Edge> ite: edgeTimeWait.entrySet()) 
							{
								if(getcurrGeneration(Key)>=ite.getKey())
								{
									if (edgeWaitContain.get(Key).contains(ite.getValue()))
									{
										edgeContain.get(Key).add(ite.getValue().getIdEDGE());
										MapActiveEdge.get(Key).add(ite.getValue());
										edgeWaitContain.get(Key).remove(ite.getValue());
										tdel.add(ite.getKey());
										if (edgeWaitContain.get(Key).isEmpty())
											delKey = true;
												
									}
								}
							}
							
							for (Double d: tdel)
								edgeTimeWait.remove(d);
							
							if (delKey)
								edgeWaitContain.remove(Key);
							
							
							if((edgeContain.get(Key).size()==1) & !edgeWaitContain.containsKey(Key) & !Key.equals(NS))
							{
								String tmpstrt =NS;
								if (Key.replace("[","").replace("]", "").length()>1)
								{	//put in wait to the father
									tmpstrt = "[";
									
									String tmpstr = StrTrun(Key);
									List<Character> list = new ArrayList<Character>();
									for(char c : tmpstr.toCharArray())
										list.add(c);
									list.sort(Comparator.naturalOrder());
									for (int iti = 0; iti<list.size()-2; iti++)
										tmpstrt = tmpstrt+list.get(iti)+", ";
									tmpstrt = tmpstrt + list.get(list.size()-2)+"]";
								}
								
								//remove edge from current set and move in wait of the right one
								
								//TODO: Chance that it has to wait!
								if (getcurrGeneration(Key)>=getcurrGeneration(tmpstrt))
								{
									if (!edgeWaitContain.containsKey(MapActiveEdge.get(Key).get(0)))
										edgeWaitContain.put(Key, new ArrayList<Edge>());
									edgeWaitContain.get(Key).add(MapActiveEdge.get(Key).get(0));		
									edgeTimeWait.put(getcurrGeneration(Key), MapActiveEdge.get(Key).get(0));
									continue;
								}
								Integer id_edge = edgeContain.get(Key).get(0);
								edgeContain.get(Key).remove(id_edge);
								edgeContain.get(tmpstrt).add(id_edge);
								
								List<Edge> tmpActiveEdgeT = MapActiveEdge.get(Key);
								
								MapActiveEdge.get(tmpstrt).add(MapActiveEdge.get(Key).get(0));
								MapActiveEdge.get(Key).remove(0);
								
								continue;
							}
							
							if(edgeContain.get(entry.getKey()).size()==1)
								continue;

							//end check
							//check which event
							int index=0; 
							Random random = new Random();

							if(getRecomb() > 0)
								index = computeIndexNextEvent(lambda_, binCoef_,edgeContain.get(Key));


							if(index == 0) {
								//***** CASE 1: THE NEXT EVENT IS COALESCENT *****
								coalescentNumber = coalescentNumber+1;

								//Create a new Coalescent node
								Node coal_node = new Node(nodeSet.size(), false, getcurrGeneration(Key));
								nodeSet.put(coal_node.getId(), coal_node);

								//Pick randomly the first active edge from the set of active edges
								int j; int id_first_edge;

								int pos_j = random.nextInt(edgeContain.get(Key).size());
								j = activeEdges.indexOf(graphEdges.get(edgeContain.get(Key).get(pos_j)));
								edgeContain.get(Key).remove(pos_j);
								System.out.println("POS J IS "+pos_j);

								id_first_edge = activeEdges.get(j).getIdEDGE();

								//Set the father of the active edge selected
								graphEdges.get(id_first_edge).setId_fath(coal_node.getId());
								//Set the time passed for this event (computed by the exponential distribution)
								graphEdges.get(id_first_edge).setTime();

								int son1 = graphEdges.get(id_first_edge).getId_son();

								//Get the arraylist of Intervals/Segments of the left edge
								ArrayList<Interval> left_segments = activeEdges.get(j).getSegments();

								//Remove the lineage from the set of the active ones
								activeEdges.remove(j);
								MapEdgeRemove(Key, id_first_edge);

								//Pick randomly the second active edge from the set of active edges
								pos_j = random.nextInt(edgeContain.get(Key).size());
								j = activeEdges.indexOf(graphEdges.get(edgeContain.get(Key).get(pos_j)));
								edgeContain.get(Key).remove(pos_j);	
								System.out.println("POS J IS "+pos_j);
								int id_second_edge = activeEdges.get(j).getIdEDGE();
								int son2 = graphEdges.get(id_second_edge).getId_son();
								//Set the father of the selected active lineage as coalescent node
								graphEdges.get(id_second_edge).setId_fath(coal_node.getId());
								//Set the time passed for this event (computed by the exponential distribution)
								graphEdges.get(id_second_edge).setTime();

								//Get the arraylist of Intervals/Segments of the right edge
								ArrayList<Interval> right_segments = activeEdges.get(j).getSegments();

								//Remove the lineage from the set of active ones
								activeEdges.remove(j);
								MapEdgeRemove(Key, id_second_edge);
								//Compute the union of the left and the right intervals
								MergeIntervals merge = new MergeIntervals();
								ArrayList<Interval> union = merge.ConcatTwoIntervalList(left_segments, right_segments);
								union = merge.merge(union);

								//Memorize the new set of segments in the coalescent node and the two sons
								coal_node.setSegments(union);
								coal_node.setIDsonsx(son1);
								coal_node.setIDsondx(son2);

								//Create a new active edge outgoing from the new coalescent node
								//Edge new_edge;
								//check status
								if(edgeContain.get(Key).isEmpty() && StrTrun(Key).length() > 1 && !StrTrun(Key).equals(StrTrun(NS))) { 
									//move edge to the right bucket or put in wait list
									//suca
									String tmpstr = StrTrun(Key);
									List<Character> list = new ArrayList<Character>();
									for(char c : tmpstr.toCharArray())
										list.add(c);
									
									Collections.shuffle(list);
									String tmpstr1 = "["+Character.toString(list.get(0))+"]";
									String tmpstr2 = "[";
									List<Character> list1 = new ArrayList<Character>();
									for (int iti = 1; iti<list.size(); iti++)
										list1.add(list.get(iti));
									list1.sort(Comparator.naturalOrder());
									for (int iti = 1; iti<list1.size()-1; iti++)
										tmpstr2 = tmpstr2+list1.get(iti)+", ";
									tmpstr2 = tmpstr2 + list.get(list.size()-1)+"]";

									if (argRecomb > 0) {
										recombinationsNumber = recombinationsNumber+1;
										SetEdge(coal_node,union,tmpstr1,Key, edgeContain,edgeWaitContain, edgeTimeWait, getcurrGeneration(tmpstr1),getcurrGeneration(Key));
										SetEdge(coal_node,union,tmpstr2,Key, edgeContain,edgeWaitContain, edgeTimeWait, getcurrGeneration(tmpstr2),getcurrGeneration(Key));
									}
									else
									{
										//goes to the big one
										String tmpstrt = "[";
										list.sort(Comparator.naturalOrder());
										for (int iti = 0; iti<list.size()-2; iti++)
											tmpstrt = tmpstrt+list.get(iti)+", ";
										tmpstrt = tmpstrt + list.get(list.size()-2)+"]";
										
										if (edgeContain.get(Key).isEmpty() && !edgeWaitContain.containsKey(Key))
											tmpstrt = NS;//TODO: WORKS for 2 SNPs, it should actually recorsively check the father as well. 
										
										SetEdge(coal_node,union,tmpstrt,Key, edgeContain,edgeWaitContain, edgeTimeWait, getcurrGeneration(tmpstrt),getcurrGeneration(Key));
									}
								}
								else if(edgeContain.get(Key).isEmpty() && StrTrun(Key).length() == 1){
									SetEdge(coal_node,union,NS,Key, edgeContain,edgeWaitContain, edgeTimeWait, getcurrGeneration(NS),getcurrGeneration(Key));
										
								} else
								{
									SetEdge(coal_node,union,Key,edgeContain);
								}
								

							}
							else{
								//***** CASE 2: THE NEXT EVENT IS RECOMBINATION *****

							}


						}
					}
					
				
					/*** stats ***/
					try {
						PrintWriter pw = new PrintWriter(new FileOutputStream(
								new File(filestast), 
								true /* append = true */));
						
						int totalMutations = 0; 

						pw.println(getcurrGeneration(NS)+"\t"+getRecombinationsNumber()+"\t"+totalMutations);
						pw.close();
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} 






				} //END RUN

			} //END FOR Gvect
		}//END FOR m-vect

	} //END MAIN

}
