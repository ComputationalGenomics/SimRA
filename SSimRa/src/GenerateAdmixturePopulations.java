import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
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
import org.apache.commons.cli.*;




/**
 * Main class of SimRa.
 * Generates the ARG randomly and then decorates it with SNP mutations and (optionally) STRs mutations.
 * 
 * @author Anna Paola Carrieri
 * @author Filippo Utro
 * @author Aritra Bose
 * 
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


		double  ValtoPut = 0.0;
		double sum = 0.0; 
		for(String  i: items) {
			if (!i.equals("99"))
				sum += (double)SelN.get("["+i+"]")/(double)(2*N);
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
		System.out.println("Random number is: "+r1);
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
			//System.out.println("eval is: "+sum);
			// System.out.println("");
			return;	
		}

		for(int i = start; i <= end && end-i+1 >= k-idx; i++) {
			//System.out.print(i + " || " + idx + " --- ");
			tmp[idx] = arr[i]; 

			findCombination(fitvect, arr, tmp, i+1, end, idx+1, k, order);

			//while(arr[i] == arr[i+1])
			//i++; 
		}

	}
	public static void findcombidx(double fitvect[], int arr[], int tmp[], int start, int end, int idx, int k, ArrayList<ArrayList<Integer>> order) {
		ArrayList<Integer> tmporder = new ArrayList<Integer>();
		if(idx == k) {
			// Current combination is ready to be printed, print it
			for (int j=0; j<k; j++)  
				tmporder.add(tmp[j]);
			//System.out.println("eval is: "+sum);
			// System.out.println("");
			order.add(tmporder);
			return;	
		}
		tmporder.clear();
		for(int i = start; i <= end && end-i+1 >= k-idx; i++) {
			//System.out.print(i + " || " + idx + " --- ");
			tmp[idx] = arr[i]; 

			findcombidx(fitvect, arr, tmp, i+1, end, idx+1, k, order);

			//while(arr[i] == arr[i+1])
			//i++; 
		}

	}
	public static double findCombinationSum(ArrayList<Double> fitvect, int arr[], int n, int k) {
		int[] tmp = new int[k]; 
		double tmpsum = 0.0;
		ArrayList<Double> order = new ArrayList<Double>();
		Arrays.sort(arr);

		findCombination(fitvect,arr,tmp,0,n-1,0,k,order);

		for (int i = 0; i < order.size(); i++) {
			//System.out.print(order.get(i)+" ");
			tmpsum += order.get(i); 
		}
		//System.out.println("");
		//System.out.println("eval is: "+tmpsum);
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
		System.out.println("Numer is: "+numer+" and Denom is: "+denom);
		double sstar = (2*numer)/denom; 

		return sstar; 
	}

	public static double FunS(double val) {
		return ((1+val)/(2+val));
	}
	public static void errorInsertingParameters(){
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


	public static Map<String, Double> ComputeFitGroup(double[] AllEpiFit, int k_way, Map<String,Double> Delta){

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
		/*	for(int i = 0; i < Combinations.size(); i++)
			System.out.print(Combinations.get(i)+" ---- ");
		System.out.println("");*/
		//Print Handlers
		for (int i = 0; i < EpiFitList.size(); i++) {
			ArrayList<Double> tmpval = new ArrayList<Double>(); 
			tmpval = EpiFitList.get(i); 
			for (int j = 0; j < tmpval.size(); j++)
				System.out.print(tmpval.get(j)+" ");
			System.out.println("");
		}

		Map<String,Double> FitGroups = new HashMap<String, Double>(); 

		for(int ij = 0; ij < EpiFitList.size(); ij++) {
			ArrayList<Double> EachEpi = new ArrayList<Double>(); 
			EachEpi = EpiFitList.get(ij);
			int EpiNum = EachEpi.size();

			if(EpiNum > 1 && deltakey.contains(Combinations.get(ij))) {
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

				System.out.println("Base Vector");
				for (int i = 0; i < eBase.length; i++)
					System.out.print(eBase[i] + " ");
				System.out.println("");

				double sstar = ComputeSstar(eBase);
				//double ValtoPut = FunS(sstar+Delta.get(Combinations.get(ij)));
				double ValtoPut = FunS(Delta.get(Combinations.get(ij)));
				FitGroups.put(Combinations.get(ij), ValtoPut);
			}
			else if (EpiNum == 1){
				for (int i = 0; i < EachEpi.size(); i++) 
					FitGroups.put(Combinations.get(ij), FunS(EachEpi.get(i)));
			}
			else {
				double prod = 1.0, ValtoPut = 0.0;
				double sum = 0.0; 
				for(int i = 0; i < EachEpi.size(); i++) {
					sum += EachEpi.get(i);
					//prod *= FunS(val);
				}
				ArrayList<Double> tmp = new ArrayList<Double>();
				double newsum = 0.0;
				for(int i = 0; i < EachEpi.size(); i++) {
					newsum += EachEpi.get(i)/sum;
					tmp.add(newsum);
				}
				int flag = 0;
				double r1 = Math.random();
				System.out.println("Random number is: "+r1);
				for(int i=0; i < EachEpi.size(); i++) {

					if(r1 < tmp.get(i)) {
						flag = 1; 
						ValtoPut = FunS(EachEpi.get(i));
						break;
					}
				}
				if (flag == 0) {
					for (int i = 0; i < EachEpi.size(); i++) {
						if (FunS(EachEpi.get(i)) > ValtoPut)
							ValtoPut = FunS(EachEpi.get(i));
					}
				}
				
				/*String[] tmpitems = Combinations.get(ij).replace("[", "").replace("]", "").replace(" ", "").split(",");
				//System.out.println(EachEpi.get(Integer.valueOf(tmpitems[0])));
				ValtoPut =  FunS(EachEpi.get(Integer.valueOf(tmpitems[0])));
				for (int i =1; i<tmpitems.length; i++)
					ValtoPut = ValtoPut * FunS(EachEpi.get(Integer.valueOf(tmpitems[i])));
				 */
				FitGroups.put(Combinations.get(ij), ValtoPut);
			}
		}

		// FitGroups.put("[99]",1.0);
		FitGroups.put("[99]",0.5);
		return FitGroups;
	}

	public static Map<String, Double> ComputeEffPop(Map<String,Double> FitGroups){
		Map<String, Double> EffPop = new HashMap<String,Double>();

		for (Map.Entry<String, Double> entry: FitGroups.entrySet()) {
			//double singsum = 0.0, dubsum = 0.0, sum3 = 0.0;
			//String tmpstr = StrTrun(entry.getKey());
			EffPop.put(entry.getKey(),entry.getValue());
		}
		//System.out.println("key is: "+tmpstr);
		/*
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
					System.out.println("Singsum is "+singsum+" dubsum is "+dubsum+" sum3 is "+sum3);
					EffPop.put(entry.getKey(), entry.getValue() - singsum + dubsum - sum3);
					EffPop.put(entry.getKey(), entry.getValue());
				}
			}

			if(tmpstr.length() == 3)
				EffPop.put(entry.getKey(), entry.getValue()); 	
		}
		 */
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
			Map<String, Double> FitGroup, Map<String, Double> EffPop, List<Integer> DelFitIdx, int argN, String NS) {
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
		FitGroup = ComputeFitGroup(revepifit, k_way, Delta);
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
		EffPop = ComputeEffPop(tmpfitgroup);
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

	public static void main(String[] args) throws Exception {


		int argN;
		int argm;
		double argRecomb;
		double argMut;
		int argG;
		String argPath;
		String argFileName;

		int runs=1; 
        ArrayList<Integer> mvect = new ArrayList<Integer>();
        ArrayList<Integer> gvect =  new ArrayList<Integer>();
        Map<String, List<String>> params = new HashMap<>();
		
        if(args.length < 5) {
			errorInsertingParameters();
		}
		 
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
		
		//System.out.println(params.get("N").get(0));
		argN = Integer.parseInt(params.get("N").get(0));
		argRecomb = Double.parseDouble(params.get("r").get(0));
		argMut = Double.parseDouble(params.get("mu").get(0));
		runs = Integer.parseInt(params.get("iter").get(0));
		
		for (int i = 0; i < params.get("g").size(); i++ ) {
			gvect.add(Integer.parseInt(params.get("g").get(i)));
		}
		for (int i = 0; i < params.get("m").size(); i++ ) {
			mvect.add(Integer.parseInt(params.get("m").get(i)));
		}
		double[] OldEpiFit = new double[params.get("s").size()];
		
		for (int i = 0; i < params.get("s").size(); i++ ) {
			OldEpiFit[i] = Double.parseDouble(params.get("s").get(i));
		}
		
		/*Options options = new Options();

        Option i1 = new Option("N", "N", true, "Number of Individuals");
        i1.setRequired(true);
        options.addOption(i1);

        Option i2 = new Option("r", "r", true, "Recombination Rate");
        i2.setRequired(true);
        options.addOption(i2);
        
        Option i3 = new Option("mu", "mu", true, "Mutation Rate");
        i3.setRequired(true);
        options.addOption(i3);
        
        Option i4 = new Option("g", false, "Segment Length");
        i4.setRequired(true);
        options.addOption(i4);
        
        Option i5 = new Option("m", false, "# of extant populations");
        i5.setRequired(true);
        options.addOption(i5);
        
        Option i6 = new Option("iter", "iter", true, "# of runs");
        i6.setRequired(true);
        options.addOption(i6);
        
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = parser.parse(options, args);

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("utility-name", options);

            System.exit(1);
        }

        argN = Integer.parseInt(cmd.getOptionValue("N"));
        argRecomb = Double.parseDouble(cmd.getOptionValue("r"));
        argMut =  Double.parseDouble(cmd.getOptionValue("mu"));
        runs = Integer.parseInt(cmd.getOptionValue("iter")); 
        
        if(cmd.hasOption("g")) {
        	System.out.println(args.length);
        	for (int i = 0; i < args.length; i++)
        		gvect.add(Integer.parseInt(args[i]));
        }
        if(cmd.hasOption("m")) {
        	for (int i = 0; i < args.length; i++)
        		mvect.add(Integer.parseInt(args[i]));
        }
        
        System.out.println(argN);
        System.out.println(argRecomb);
*/
		/*argN = Integer.parseInt(args[0]);
		argm = Integer.parseInt(args[1]);
		argRecomb = Double.parseDouble(args[2]);
		argMut = Double.parseDouble(args[3]);
		argG = Integer.parseInt(args[4]);
		argPath = args[5].toString();
		argFileName = args[6].toString();*/

		//int [] mvect=new int[] {20, 50, 80, 120, 150};
		//int [] gvect=new int[] {25, 75};
		//double [] fvect=new double[]{0.3}; //, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
		//double[] OldEpiFit = new double[] {-0.9};
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
		int k_way = OldEpiFit.length; //1;	//how many maximum interacting SNPs 

		// This map indicates for which epistasis will be in effect with the delta associated 
		Map <String, Double> Delta = new HashMap<String,Double>();
		//Delta.put("[0, 1]",-0.3);
		//Delta.put("[0]", 0.12);
		//Delta.put("[1]", 0.12);*/
		//Delta.put("[0, 1, 2]", 0.05);

		
		//Calculating f' for each key
		FitGroup = ComputeFitGroup(AllEpiFit, k_way, Delta); 
		double fsum = 0.0;
		double fp = 0.0;
		double f3 = 0.0;
		for(Map.Entry<String, Double> entry: FitGroup.entrySet()) {
			if (!entry.getKey().equals(NS))
			{
				String line = entry.getKey();
				int count = line.length() - line.replace(",", "").length();
				if (count == 0)
				{	
					fsum -= entry.getValue();
				}
				if (count == 1)
				{
					String t[] = line.replace("[","").replace("]","").split(", ");
					fp += FitGroup.get("["+t[0]+"]")*FitGroup.get("["+t[1]+"]");
					 if (!Delta.containsKey(line))
						 FitGroup.put(line,FitGroup.get("["+t[0]+"]")*FitGroup.get("["+t[1]+"]"));
				}
				if (count == 2)
				{
					f3 = FitGroup.get("[0]")*FitGroup.get("[1]")*FitGroup.get("[2]");
					if (!Delta.containsKey("[0, 1, 2]")) 
						FitGroup.put("[0, 1, 2]",f3);
				}
			}
		}
		
		//fsum = - (FitGroup.get("[0]")+FitGroup.get("[1]")+FitGroup.get("[2]"));
		//fp = FitGroup.get("[0]")*FitGroup.get("[1]");
		//fp += FitGroup.get("[0]")*FitGroup.get("[2]");
		//fp += FitGroup.get("[1]")*FitGroup.get("[2]");	 
		
		FitGroup.put(NS, 1.0+fsum+fp-f3); // +fp ///(2.0+OldEpiFit[0])
		//FitGroup.put("[0, 1]",FitGroup.get("[0]")*FitGroup.get("[1]"));
		//FitGroup.put("[0, 2]",FitGroup.get("[0]")*FitGroup.get("[2]"));
		//FitGroup.put("[1, 2]",FitGroup.get("[1]")*FitGroup.get("[2]"));
		
		
		
		EffPop = ComputeEffPop(FitGroup);

		System.out.println("---- Fitness Values ----");
		for(Map.Entry<String, Double> entry: FitGroup.entrySet()) {
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}
		System.out.println("------------------------");

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
				//int runs=100;
				setTotalRuns(runs);

				//argN = 1000;
				N = argN;
				argm = mvect.get(im);
				//argRecomb = 1.0;
				//argMut = 1.5;
				argG = gvect.get(indg);
				argPath = "stats/";
				argFileName = "prova";

				//compute number of extant unit under selection;
				

				String filestast=argPath+"N"+String.valueOf(argN)+"_m"+String.valueOf(argm)+"_m1"+String.valueOf(fns)+"_g"+String.valueOf(argG)+"_s"+String.valueOf(OldEpiFit[0])+"_r"+String.valueOf(argRecomb)+".stats";


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
					String tnsName = "stats/selected1_m"+String.valueOf(argm)+"_g"+String.valueOf(argG);
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
					
					BufferedWriter fwritesize = null;
					String tnsizeName = "stats/size_m"+String.valueOf(argm)+"_g"+String.valueOf(argG);
					for (int i = 0; i < AllEpiFit.length; i++) 
						tnsizeName = tnsizeName+"_s"+String.valueOf(i+1)+"="+String.valueOf(AllEpiFit[i]);
					File filesize = new File(tnsizeName+"_r="+String.valueOf(argRecomb)+".txt");
					try {
						fwritesize = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filesize, true), "utf-8"));

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

						for(int i = 0; i < FunFit.length; i++) {
							sizeMap.put(i,(int)Math.round(FunFit[i]*argm));//* FunFit[i]*argm  0.25  0.3 change - > 0.5
							//sizeMap.put(i,rand.nextInt(argm) + 1);
						}
						
						
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

					}


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


					//********* END
					
					
					Map <String, ArrayList<Integer>> edgeContain = new HashMap<String, ArrayList<Integer>>(); //edgeContainer for each key
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
					int iteration=-1;
					while(activeEdges.size() !=1){
						iteration=iteration+1;
						


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
							tmpActiveEdge = entry.getValue();
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

								
								if (tpast == 0 || t_ > tpast) { //Take overall max
									tpast = t_;  
									ipast = Key; 
									pastlambda = lambda_; 
									pastbincoef = binCoef_;
								}
								
							}
							//If there is a key such as s1 or s2, which has at most one edge in the container 
							//Then check if there are lineages which contain that key and is still active. 
							//If they are not active then assign this edge to Neutral, 99 case and delete the key. 
							else if(StrTrun(Key).length() == 1 && edgeContain.get(Key).size() < 2) {
								String tmpstr = StrTrun(Key);
								int sinflag = 0; 
								for (int ij = 0; ij < tmpkeys.size(); ij++) {
									if(StrTrun(tmpkeys.get(ij)).contains(tmpstr) && StrTrun(tmpkeys.get(ij)).length() > 1
											&& edgeContain.get(tmpkeys.get(ij)).size() > 1){
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
								if (StrTrun(Key).length() > 1 && !Key.equals(NS) && edgeContain.get(Key).size() < 2) { //(StrTrun(Key).length() > 1 && !Key.equals(NS) && getRecomb() == 0) {
									ArrayList<Integer> tmpid = edgeContain.get(entry.getKey());
									//edgeContain.get(NS).addAll(tmpid);
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
						
						//5) Update the current generation
						setGeneration(computeNextGeneration( timeContain.get(ipast).get(2))); //ipast,

						
						try {
							fwritepick.write(getGeneration()+"\t"+ipast+"\t"+edgeContain.get(ipast).size()+"\t"+numRun+"\n");
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					
						
						try {
							for (Map.Entry<String, ArrayList<Integer>> it_entry : edgeContain.entrySet())
							{
								fwritesize.write(getGeneration()+"\t"+it_entry.getKey()+"\t"+it_entry.getValue().size()+"\t"+numRun+"\t"+iteration+"\n");
							}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						
						

						for(Map.Entry<String, List<Edge>> entryn: MapActiveEdge.entrySet()) {

							if(!timeContain.containsKey(entryn.getKey()))
								continue;

							ipast = entryn.getKey();
							pastlambda = timeContain.get(ipast).get(0);
							pastbincoef = timeContain.get(ipast).get(1);

							int index_=-1; 
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
							
							//***** CASE 1: THE NEXT EVENT IS COALESCENT *****
							if(index == 0) {

				
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
											SetEdge(coal_node,union,ipast,edgeContain);
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
									if(!(tmpkeys.contains(leftkey)))
										leftkey = NS;
									if(!(tmpkeys.contains(rightkey)))
										rightkey = NS;
								}
								else {
									leftkey = NS;
									rightkey = NS; 
								}

								SetEdge(recomb_node,splitted.get(0),leftkey,edgeContain);
								//}
								/*if (edgeContain.get(rightkey).isEmpty()) {
									//Record ID here for mutation
									SetEdge(recomb_node,splitted.get(0),NS,edgeContain);
								}
								else {*/
								SetEdge(recomb_node,splitted.get(0),rightkey,edgeContain);
								//}

								/*left_edge = new Edge(graphEdges.size(), recomb_node.getId(), -1, splitted.get(0), GetEffPop(leftkey));
									left_edge.computeLength();
									left_edge.computeRate(GetEffPop(leftkey), getRecomb(), getG());
									edgeContain.get(leftkey).add(left_edge.getIdEDGE());
									graphEdges.put(left_edge.getIdEDGE(), left_edge);
									activeEdges.add(left_edge);
								}*/
								/*boolean ls=false, rs=false;
								if((SplittingIntervals.getRec_pos()>=segmentposition) && inselection)
								{
									ls=true;
								}
								else
								{
									if (inselection)
										rs=true;
								}
								Edge left_edge;
								if  (ls)
								 left_edge = new Edge(graphEdges.size(), recomb_node.getId(), -1, splitted.get(0), ls, getNs());
								else
									 left_edge = new Edge(graphEdges.size(), recomb_node.getId(), -1, splitted.get(0), ls, getN());

								left_edge.computeLength();
								if (ls)
								{
									left_edge.computeRate(getNs(), getRecomb(), getG());
									edgesUS.add(left_edge.getIdEDGE());
								}
								else
								{
									left_edge.computeRate(getN(), getRecomb(), getG());
									edgesNS.add(left_edge.getIdEDGE());
								}*/


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
					//System.out.println("Total Number of recombinations W ---> "+getRecombinationsNumber());
					//System.out.println("*** PROPERTY (27):  recombination rate ~ (W/ gNT ) ***");
					double computed_r = getRecombinationsNumber()/(GetEffPop(NS)*getG()*getGeneration());
					System.out.println(getRecomb()+ " ~ " +computed_r);

					//System.out.println("Total number of coalescent events Z ---> "+getCoalescentNumber());
					//System.out.println("*** PROPERTY (28): mW  << Z ***");
					//System.out.println(getExtantUnits()*getRecombinationsNumber()+" << "+getCoalescentNumber());

					//Matrix graph useful for visit the graph and decorate it with mutations SNP and STRs
					DiGraph graph = new DiGraph();
					graph.createDigraph();




					//***************************** DECORATING ARG WITH SNP MUTATIONS   ***********************************

					//DecorateMutationsSNP.decoratingWithMutations(getMu(), graph);
					//String filesSNP="N"+String.valueOf(argN)+"_m"+String.valueOf(argm)+"_g"+String.valueOf(argG)+"_f"+"_fN"+"_r"+String.valueOf(argRecomb);


					//TO UNCOMMENT
					//DecorateMutationsSNP.decoratingWithMapMutations(getMu(), graph, MutMap);
					//String filesSNP="N"+String.valueOf(argN)+"_m"+String.valueOf(argm)+"_g"+String.valueOf(argG)+"_f"+"_fN"+"_r"+String.valueOf(argRecomb);

					//String txtFileSNP = filesSNP.concat("_SNP_"+numRun+".txt");
					//DecorateMutationsSNP.createSNP_TxtFile(wholePath, txtFileSNP);

					//END TO UNCOMMENT


					int totalMutations = DecorateMutationsSNP.computeTotalMutationsNumber();
					System.out.println("Total Number of mutations Y ---> "+totalMutations);
					//System.out.println("*** PROPERTY (26):  mutation rate ~ (Y/ gNT ) ***");
					double computed_mu = totalMutations/(getG()*GetEffPop(NS)*getGeneration());
					System.out.println(getMu()+ " ~ " +computed_mu);
					//setScaledGeneration(getGeneration()*GetEffPop(NS));
					setScaledGeneration(getGeneration());
					///System.exit(1);
					/*
			try {
				// PRINTING FILE STATS
				fwrite.write("Recombination rate \t Mutation Rate \t TimeInGenerations \t Time of GMRCA \t #RecombinationNodes \t #CoalescenceNodes \t #SNP Mutations\n");
				fwrite.write(computed_r+"\t"+computed_mu+"\t"+getScaledGeneration()+"\t"+getGeneration()+"\t"+getRecombinationsNumber()+"\t"+getCoalescentNumber()+"\t"+totalMutations+"\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}


			//***************************** DECORATING ARG WITH STRs MUTATIONS  ***********************************

			if(args.length > 7) {

				//Check if the 6th parameter is the string "-STR"
				if(args[7].toString().equalsIgnoreCase("-STR")) {
					//Check if the the number of parameter after "STR" option is correct
					if(args.length == 11) {
						int numSTRs = Integer.parseInt(args[8]);
						int initialeState = Integer.parseInt(args[9]);
						double mutRateSTRs = Double.parseDouble(args[10]);
						DecorateSTRs.inizializeStrs(numSTRs,initialeState,mutRateSTRs*Math.pow(10, -3));
						DecorateSTRs.decoratingDeltaStrsEdges();
						DecorateSTRs.updateStrs();
						String txtFileSTR = filename.concat("_STRs.txt");
						DecorateSTRs.createSTR_TxtFile(wholePath, txtFileSTR, getNodeSet());
					}
					else
						errorInsertingParameters();
				}
				else
					errorInsertingParameters();	
			}	
					 */
					//****** CREATING DOT FILE *****/
					/*	try {
							//creteDot(numRun+".dot");
						} catch (FileNotFoundException e) {
							// TODO Auto-generated catch block
							//e.printStackTrace();
						} catch (UnsupportedEncodingException e) {
							// TODO Auto-generated catch block
							//e.printStackTrace();
						} 
					 */
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
					//Create L file *** FILE L ***
					//String txtFileL = filename.concat("_L_"+numRun+".txt");
					String txtFileL = filename.concat("_L"+".txt");
					CreatingFilesForStructure.createLtxt(wholePath, txtFileL, getGraphEdges(), getNodeSet());

					//Create S file *** FILE S ***
					//String txtFileS = filename.concat("_S_"+numRun+".txt");
					String txtFileS = filename.concat("_S"+".txt");
					CreatingFilesForStructure.createStxt(wholePath, txtFileS, getGraphEdges(), getNodeSet());


					try {
						fwritepick.close();
						fwritesize.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

				} //end numRun

				try {fwrite.close();} catch (Exception ex) {}

				long end = System.currentTimeMillis();
				NumberFormat formatter1 = new DecimalFormat("#0.00000");
				System.out.print(formatter1.format((end - start) / 1000d));
				int mb = 1024;
				runtime = Runtime.getRuntime();
				System.out.println("\t"+(runtime.totalMemory() - runtime.freeMemory()) / mb);
			}
		}
		//	}

	} //end main procedure
}
