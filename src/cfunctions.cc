/*
  cfunctions.c

  Version : 1.1.04
  Author  : Niko Beerenwinkel

  (c) Max-Planck-Institut fuer Informatik, 2004
*/

#include "include/mtreemix.h"
#include "include/cfunctions.h"


integer_matrix resample(integer_matrix& pattern, matrix& resp, int S, matrix& bresp)
{
  // resample S rows from pattern(,)
  // returns resampled data; corresponding responsibilities are in bresp

  int N = pattern.dim1();  // sample size
  int L = pattern.dim2();  // pattern length
  int K = resp.dim1();

  integer_vector x(L);
  vector gamma(K);  // responsibility of k-th tree 
  integer_matrix R(S,L);

  for (int s=0; s<S; s++)
    {
      int ri = (int) ((double) N * rand() / (RAND_MAX + 1.0));

      // pattern:
      x = pattern.row(ri);
      for (int j=0; j<L; j++)
	R(s,j) = x[j];

      // responsibilities:
      gamma = resp.col(ri);
      for (int k=0; k<K; k++)
        bresp(k,s) = gamma[k]; 
    }

  return R;
}


vector CI(list<double>& L, double alpha)
{
  // calculate (alpha*100)% confidence interval from list of doubles L
  // returns a 2-dim vector whose first entry is the lower limit
  // and whose second entry is the upper limit of the confidence interval

  vector ci(0.0, 1.0);

  int len = L.length();
  if (len >= 2)
    {
      L.sort();  // we can be more efficient than this...

      // compute indices of confidence limits:
      double eps = .5 / (double) len;  // for rounding
      int lower_index = std::max(0, int(alpha*(double)len - eps) - 1);
      int upper_index = std::min(len - 1, int((1.0-alpha)*(double)len + 1.0));
	  
      //ci[0] = L.contents(L.get_item(lower_index));  // lower limit
      //ci[1] = L.contents(L.get_item(upper_index));  // upper limit
      ci[0] = L.contents(lower_index);  // lower limit
      ci[1] = L.contents(upper_index);  // upper limit
    }

  return ci;
}


array< map<edge,double> > mtreemix_bootstrap(array<string>& profile, integer_matrix& pattern, int K, array<graph>& G, array< map<int,node> >& node_no, matrix& resp, int B, double eps, int special_weighing, int uniform_noise, array< map<edge,double> >& lower, array< map<edge,double> >& upper, vector& alpha_lower, vector& alpha_upper, double confidence_level)
{
  // Bootstrap analysis
  // modifies the responsibilities resp(,) ! (see below)

  int N = pattern.dim1();  // sample size
  int L = pattern.dim2();  // pattern length

  vector oneN = ones(N);

  // for running mtreemix_fit1:
  vector                     balpha(1);      // mixture parameter
  array< graph >             bG(1);          // digraph
  array< map<int,node> >     bnode_no(1);    // node of mutation index
  array< map<node,string> >  bevent(1);      // substitution event
  array< map<edge,double> >  bcond_prob(1);  // conditional probabilities
  
  // map nodes onto event indices:
  array< map<node,int> > no_of_node(K);
  for (int k=0; k<K; k++)
    for (int j=0; j<L; j++)
      no_of_node[k][node_no[k][j]] = j;

  // for bootstrapping:
  integer_matrix bpat(N,L);  // one bootstrap sample of size N
  array< map<edge,double> > supp(K);  // bootstrap edge support
  edge e, be;

  for (int k=0; k<K; k++)
    {
      
      list<double> alpha_distr;  // empirical distribution of alpha as estimated by the bootstrap
      array< list<double> > bcond_prob_distr(L);  // empirical distribution of cond_prob as estimated by the bootstrap (edges are indexed by the index of the target node)

      for (int b=0; b<B; b++)
	{
  	  matrix bresp(K,N);
	  bpat = resample(pattern, resp, N, bresp);
          double bN_k = oneN * bresp.row(k);

	  // estimate alpha:
	  alpha_distr.append(bN_k / N);  // == sum(bresp(k,:) / N

          // normalize bresp to obtain distribution over samples (instead of over the tree components)
          for (int i=0; i<N; i++)
	    bresp(k,i) /= bN_k; 
 
	  if (k == 0)
            {
              mtreemix_fit0(profile, bpat, balpha, bG, bnode_no, bevent, bcond_prob, bresp.row(k), uniform_noise, special_weighing);
            }
          else  // k >= 1
	    mtreemix_fit1(profile, bpat, balpha, bG, bnode_no, bevent, bcond_prob, bresp.row(k), eps, special_weighing);
	  
	  for (int j1=0; j1<L; j1++)
	    for (int j2=0; j2<L; j2++)
		if ((be = edge_between(bnode_no[0][j1], bnode_no[0][j2])) != NULL) //nil)  // edge in bootstrap tree
		    if ((e = edge_between(node_no[k][j1], node_no[k][j2])) != NULL) //nil)  // edge in k-th tree of model
            	  {
		    // count edge:
		    // supp[k][e] += (1.0 / (double) B);  // relative counts
		    supp[k][e] += 1.0;  // absolute counts

		    // record cond_prob estimate:
		    bcond_prob_distr[j2].append(bcond_prob[0][be]);
		  }
	}
 
      // calculate confidence intervals:
      vector ci = CI(alpha_distr, confidence_level);
      alpha_lower[k] = ci[0];
      alpha_upper[k] = ci[1];

      forall_edges(e, G[k])
	{
	  list<double> D = bcond_prob_distr[no_of_node[k][target(e)]];  // distrinution of p(e) in G[k]
          ci = CI(D, confidence_level);
	  lower[k][e] = ci[0];
	  upper[k][e] = ci[1];
	}

    }

  bG.clear();  
  return supp;
}

void mtreemix_wait(int L, vector& alpha, array<graph>& G, array< map<edge,double> >& lambda, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob, int n, int sampling_mode, double sampling_param, integer_matrix& pattern, vector& wtime, vector& stime)
{
  // simulate patterns and their waiting times
  
  int K = alpha.dim();

  integer_vector pat(L);

  // map nodes onto event indices:
  array< map<node,int> > no_of_node(K);
  for (int k=0; k<K; k++)
    for (int j=0; j<L; j++)
      no_of_node[k][node_no[k][j]] = j;

  int i = 0;
  while (i < n)
    {
      // draw branching number k:
      int k = discrand(alpha);
      
      // draw/set sampling time:
      if (sampling_mode == EXPONENTIAL)
	stime[i] = expcdf(sampling_param);
      if (sampling_mode == CONSTANT)
	stime[i] = sampling_param;
      
      // wait on branching k until sampling time stime[i]:
      wtime[i] = mtree_wait(G[k], node_no[k][0], cond_prob[k], lambda[k], no_of_node[k], stime[i], pat);
      
      for (int j=0; j<L; j++)  // store pattern
	pattern(i,j) = pat[j];
      
      i++;
    }
  
}

array< list<double> > mtreemix_time(int L, graph& G_k, map<edge,double>& lambda_k, map<int,node>& node_no_k, map<node,int>& no_of_node_k, map<edge,double>& cond_prob_k, int n)
{
  // simulate waiting process on a given tree component G_k
  
  int N = pow2(L-1);

  array< list<double> > wtime(N);  // lists of waiting times for all patterns

  int i = 0;
  while (i < n)
    {      
      mtree_time(G_k, node_no_k[0], cond_prob_k, lambda_k, no_of_node_k, wtime); // wait on branching k

      i++;
    }
  
  return wtime;
}

array< map<edge,double> > rescale_cond_prob(array< graph >& G, array< map<edge,double> >& cond_prob, int sampling_mode, double sampling_param, double output_param)
{
  int K = cond_prob.size();
  edge e;

  array< map<edge,double> > cond_prob_prime(K);

  switch(sampling_mode)
    {
    case CONSTANT: //  t_0  -->  t
      // p(t) = 1 - (1 - p)^(t/t_0)
      for (int k=0; k<K; k++)
	forall_edges(e, G[k])
	  cond_prob_prime[k][e] = 1.0 - pow(1.0 - cond_prob[k][e], output_param / sampling_param); 
      break;
      
    case EXPONENTIAL:  //  lambda_0  -->  lambda
      // p(lambda) = 1 / (1 + (1/p - 1) * lambda/lambda_0)
      for (int k=0; k<K; k++)
	forall_edges(e, G[k])
	  cond_prob_prime[k][e] = 1.0 / (1.0 + ((1.0 / cond_prob[k][e]) - 1.0) * (output_param / sampling_param)); 
      break;
      
    default :
      std::cerr << "Unknown sampling_mode -- " << sampling_mode << std::endl;
      exit(1);
    }
  
  return cond_prob_prime;
}



