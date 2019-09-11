import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.TreeSet;


public class DiGraph {
	
	private int[][] adjMatrix;
	public  int vertices;
	
	
	public int[][] getAdjMatrix() {
		return adjMatrix;
	}

	public void setAdjMatrix(int[][] adjMatrix) {
		this.adjMatrix = adjMatrix;
	}

public void createDigraph(){
	
	vertices = GenerateARG.getNodeSet().size();
	adjMatrix = new int[vertices][vertices];
	
    for (int i = 0; i < vertices; i++ ) {
      for (int j = 0; j < vertices; j++ ) {
	adjMatrix[i][j] = 0;
      }
    }
    
    Iterator<Integer> it_edges = GenerateARG.getGraphEdges().keySet().iterator();  
	
	
	while (it_edges.hasNext()) {  
		
		Integer keyE = it_edges.next();
		if(GenerateARG.getGraphEdges().get(keyE).getId_fath() != -1) {
			
			int p = GenerateARG.getGraphEdges().get(keyE).getId_fath();
			int f = GenerateARG.getGraphEdges().get(keyE).getId_son();
			
			if(adjMatrix[p][f]==0) {
				adjMatrix[p][f] = 1;
			}
			else if(adjMatrix[p][f]==1) {
				adjMatrix[p][f]=2;
			}
			
		}
	}
}


public void printMatrixGraph() {
	
	int i, j;
	System.out.println("Number of vertices :" +adjMatrix.length);
	
	for (i = 0; i < adjMatrix.length; i++) {
		  
	    for (j = 0; j < adjMatrix[0].length; j++) {
	    	if(adjMatrix[i][j]==1)
	    		System.out.println(i+"-->"+j);
	    	if(adjMatrix[i][j]==2) {
	    		System.out.println(i+"-->"+j);
	    		System.out.println(i+"-->"+j);
	    	}
	    }
	    
	  }
}
	
public TreeSet<Integer> computeLeaves(int root) {
	
	TreeSet<Integer> leaves = new TreeSet<Integer>();
	Queue<Integer> q = new LinkedList<Integer>();
	boolean visited[] = new boolean[GenerateARG.getNodeSet().size()];
	
	for(int i = 0; i < visited.length; i++) {
		visited[i] = false;
	}
	
	//If root is a leaf I have to add just that leaf in leaves, no reason to do the visit
	if(root >= 0 && root < GenerateARG.getExtantUnits()) {
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
			for (int u = adjMatrix[w].length-1; u >= 0 ; u--) {
			
			if (adjMatrix[w][u]!=0) { 
				//System.out.println("Visiting edge "+w+"-"+u);
				if(!visited[u]){
		          q.add(new Integer(u));
		          visited[u] = true;
		          //System.out.println("Visiting node :"+u);
		          if(u >= 0 && u < GenerateARG.getExtantUnits()) {
		        	  //u is a leaf node 
		        	  leaves.add(u);
		          }
				 	}
			  	  }	
				}
			}//end while  
	} // end else root is not a leaf
	return leaves;
}

}
