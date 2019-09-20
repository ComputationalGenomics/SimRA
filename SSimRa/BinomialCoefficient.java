
/**
 * The class BinomialCoefficient compute the binomial coefficient (n r)
 * @author Anna Paola Carrieri
 *
 */
public class BinomialCoefficient {

		BinomialCoefficient(){
			
		}
		/**
		  * Function that returns the binomial coefficient (n r)
		  * @param n integer value = number of elements in the set X
		  * @param r integer value = number of elements in each subset of X
		  * @return binomial coefficient = number of distinct k-elements subsets of X
		 */
		public static long binomialCoeff(int n, int k)
		{
		    int[][] C = new int[n+1][k+1];
		    int i, j;
		 
		    // Compute value of Binomial Coefficient in bottom up manner
		    for (i = 0; i <= n; i++)
		    {
		        for (j = 0; j <= min(i, k); j++)
		        {
		            // Base Cases
		            if (j == 0 || j == i)
		                C[i][j] = 1;
		 
		            // Calculate value using previosly stored values
		            else
		                C[i][j] = C[i-1][j-1] + C[i-1][j];
		        }
		    }
		 
		    return C[n][k];
		}
		
		public static int min(int a, int b)
		{
		    return (a<b)? a: b;
		}
	
}
