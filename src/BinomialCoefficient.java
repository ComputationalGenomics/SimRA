/*
  Copyright 2015 IBM Corporation


Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.

You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under 
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF 
ANY KIND, either express or implied. 
See the License for the specific language governing permissions and limitations under the License.


*/

/**
 * The class BinomialCoefficient compute the binomial coefficient (n r)
 * 
 * @author Anna Paola Carrieri
 *
 */
public class BinomialCoefficient {

		BinomialCoefficient(){
			
		}
		/**
		  * Function that returns the binomial coefficient (n r)
		  * @param n  number of elements in the set X
		  * @param k number of elements in each subset of X
		  * @return binomial coefficient, that is the number of distinct k-elements subsets of X
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
		
		/**
		 * The function returns the minimum value between two integer values
		 * @param a integer 
		 * @param b integer
		 * @return the minimum value between the parameters a and b
		 */
		public static int min(int a, int b)
		{
		    return (a<b)? a: b;
		}
	
}
