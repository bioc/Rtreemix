/*
  kmeans.c - kmeans clustering

  Version : 1.1
  Author  : Niko Beerenwinkel
  
  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#include "include/kmeans.h"
#include "include/mtree.h"


int argmin(vector& x, array<vector>& M)
{
  int k_min = -1;
  double d, dist_min = DBL_MAX;

  for (int k=0; k<M.size(); k++)
    {
      if ((d = (x - M[k]).sqr_length()) < dist_min)
	{
	  dist_min = d;
	  k_min = k;
	}
    }
  
  return k_min;
}
    


array<vector> kmeans_init(int K, matrix& X, double min_diff = 1e-10)
{
  int N = X.dim1();
  array<vector> M(K);  // cluster means
  
  // guess initial means
  array<int> sigma = permutation(N);
  int i = 0,  k = 0;
  while ((k < K) & (i < N))
    {
      vector x = X.row(sigma[i]);
      int ok = 1;
      for (int l=0; l<k; l++)
	if ((x - M[l]).sqr_length() < min_diff)
	  {
	    ok = 0;
	    break;
	  }
      if (ok)
	M[k++] = x;
      i++;
    }
      
  if (i >= N)  // failure
    {
      std::cerr << "k-means: Unable to find k = " << K << " sufficiently (min_diff >= " << min_diff << ") different vectors!" << std::endl
		<< "         Try changing k or min_diff." << std::endl;
      exit(1);
    }
  
  return M;
}



double kmeans_iterate(int K, matrix& X, integer_vector& C, array<vector>& M)
{
  int N = X.dim1();
  int L = X.dim2();

  integer_vector C_new(N);  // cluster assignments

  do  // iterate:
    {
      C = C_new;

      // cluster assignment:
      for (int i=0; i<N; i++)
	C_new[i] = argmin(X.row(i), M);

      // junming's code:
      // update means:
      matrix sum (K,L);
      integer_vector count(K);
         	
      for (int i=0; i<N; i++)
      	{
	    int c = C_new[i]; //.to_long();
      	  count[c]++;
      	  sum.row(c) += X.row(i);
      	}
 
      for (int k=0; k<K; k++)
	  M[k] = sum.row(k) * ( 1 / (double) count[k]); //.to_double());
		
    /*
      // niko's code:
      // update means:
      for (int k=0; k<K; k++)
	{
	  vector sum(L);
	  int count = 0;
	  for (int i=0; i<N; i++)
	    {
	      if (C_new[i] == k)
		{
		  count++;
		  sum = sum + X.row(i);
		}
	    }
	  for (int j=0; j<L; j++)  sum[j] /= (double) count;
	  M[k] = sum;
	}
     */	
    } 
  while (C_new != C);

  // calculate within-point scatter W
  double W = 0.0;
  for (int i=0; i<N; i++)
    {
	int c = C[i]; //.to_long();
      W += (X.row(i)-M[c]).sqr_length();
    }
    	
  /* niko's code
    for (int k=0; k<K; k++)
      for (int i=0; i<N; i++)
       if (C[i] == k)
        W += (X.row(i) - M[k]).sqr_length();
  */

  return W;
}



double kmeans(int K, int S, matrix& X, integer_vector& C_opt, array<vector>& M_opt)
{
  int N = X.dim1();  // sample size

  integer_vector C(N);  // cluster assignments
  array<vector> M(K);  // cluster means
  double W = 0.0,  W_opt = DBL_MAX;  // scatter (objective function)

  for (int s=0; s<S; s++)  // try S starting solutions
    {
      M = kmeans_init(K, X);  // initialize randomly: guess means
      W = kmeans_iterate(K, X, C, M);  // iterate
      if (W < W_opt)
	{
	  C_opt = C;
	  M_opt = M;
	  W_opt = W;
	}
    }
  
  return W_opt;
}
