/*
  mtreemix.c

  Version : 1.2.00
  Author  : Niko Beerenwinkel
  
  (c) Max-Planck-Institut fuer Informatik, 2004-2005
*/


#include "include/mtree.h"
#include "include/mtreemix.h"
#include "include/Rtreemix_patch.h"
//R includes
#include <R.h>


vector event_freq(integer_matrix& pat)
{
  // compute frequencies of elementary events

  int L = pat.dim2();  // pattern length
  vector freq(L);

  for (int j=0; j<L; j++)
    {
      integer_vector col = pat.col(j);
      freq[j] = nonnegmean(col);  // == rel. frequency
      if (freq[j] < 0.0)  // no data in this column!
	{
	  std::cerr << "No data in column " << j << " !" << std::endl;
	  _Rtreemix_exit(1);
	}
    }

  return freq;
}



void guess_missing_data(integer_matrix& pat, integer_matrix& pat_hat, vector& freq)
{
  int N = pat.dim1();  // sample size
  int L = pat.dim2();  // pattern length
  
  for (int j=0; j<L; j++)
    {
      if (freq[j] < 0.0)  // no data in this column!
	for (int i=0; i<N; i++)
	    pat_hat(i,j) = (myrand() > 0.5) ? 1 : 0;  // random choice
      else
	for (int i=0; i<N; i++)
	  if (pat(i,j) < 0)
	    pat_hat(i,j) = integer(freq[j] + 0.5);  // substitute most common value
	  else
	    pat_hat(i,j) = pat(i,j);
    }
}



void guess_resp(integer_matrix& pattern, int K, int S, matrix& resp)
{
  // guess initial responsibilities via k-means clustering

  int N = pattern.dim1();
  int L = pattern.dim2();

  if (K == 1)
    {
      for (int k=0; k<K; k++)
	for (int i=0; i<N; i++)
	  resp(k,i) = 1.0;
    }
  
  else {  // run (K-1)-means clustering
    
    // make data real:
    matrix X(N,L); 
    for (int i=0; i<N; i++)
      for (int j=0; j<L; j++)
	  X(i,j) = (pattern(i,j) < 0) ? 0.5 : (double) pattern(i,j);  // 0.5 if missing

    integer_vector C(N);  // cluster assignments
    if (K == 1)
      for (int i=0; i<N; i++)
	C[i] = 0;  // no cluster
    else // (K-1)-means:
      {
	array<vector> M(K - 1);  // cluster means
	kmeans(K - 1, S, X, C, M);
      }
    
    // responsibilities:
    double dK = (double) K;
    for (int i=0; i<N; i++)    
      resp(0,i) = 0.01;

    for (int k=1; k<K; k++)
      for (int i=0; i<N; i++)
		resp(k,i) = 0.99 / (dK + 3.0);
    
    for (int k=1; k<K; k++)
      for (int i=0; i<N; i++)
		if (C[i] + 1 == k) 
			resp(k,i) = (3.0 * 0.99) / (3.0 + dK);
  }
  
}



void mtreemix_DOT(array< graph >& G, array< map<node,string> >& mut, array< map<edge,double> >& edge_weight, array< map<int,node> >& node_no, vector& alpha, char *filestem, int edgeweighttype)
{
  // produces output in dot format (graphviz package)
  // see http://www.research.att.com/sw/tools/graphviz/

  node v;  edge e;

  char filename[255];
  sprintf(filename, "%s.dot", filestem);
  std::ofstream dotout(filename);
  
  dotout << "digraph \"" << filestem << "\" {" << std::endl << std::endl;

  // graph attributes:
  dotout << "\tsize = \"7.5, 10\";" << std::endl;
  dotout << "\tratio = \"auto\";" << std::endl;
  dotout << "\tnodesep = \"0.1\";" << std::endl;
  dotout << "\tranksep = \"0.3 equally\";" << std::endl << std::endl;
  
  // node and edge attributes:
  dotout << "\tnode [shape=\"plaintext\", height=\"0.3\", fontsize=\"12\", style=\"filled\", fillcolor=\"white\"];" << std::endl;
  dotout << "\tedge [fontsize=\"10\"];" << std::endl << std::endl;

  for (int k=G.low(); k<=G.high(); k++)
    {
      // cluster name:
      dotout << "\tsubgraph cluster" << k << " {" << std::endl << std::endl;
      // cluster label:
      dotout.precision(2);
      dotout << std::showpoint;
      dotout << "\t\tlabel=\"" <<  alpha[k] << "\";" << std::endl;
      // cluster attributes:
      //dotout << "\t\tstyle=filled;" << std::endl;
      //dotout << "\t\tcolor=lightgrey;" << std::endl << std::endl;
      
      list<edge> E = G[k].all_edges();

      // nodes:
      v = node_no[k][0];
      dotout << "\t\t \"" << v << "\" [label=\"" << mut[k][v] << "\"];" << std::endl;
      forall(e, E)
	{
	  v = target(e);
	  dotout << "\t\t \"" << v << "\" [label=\"" << mut[k][v] << "\"];" << std::endl;
	}
      dotout << std::endl;
      
      // edges:
      forall(e, E)
	{
	  node s = source(e);
	  node t = target(e);
	  dotout << "\t\t \"" << s << "\" -> \"" << t << "\"";
	  if (edgeweighttype == PROBABILITY){
	      dotout.precision(2);
	      dotout << std::showpoint;
	      dotout << " [label=\"" << edge_weight[k][e] << "\"];" << std::endl;
	  }
	  if (edgeweighttype == WAITING_TIME)
	      dotout << " [label=\"" << int(1/edge_weight[k][e]) << "\"];" << std::endl;  // in days
	  if (edgeweighttype == WAITING_TIME_WEEKS)
	      dotout << " [label=\"" << int((1/edge_weight[k][e])/7.0) << "\"];" << std::endl;
	  if (edgeweighttype == BOOTSTRAP_COUNT)
	      dotout << " [label=\"" << int(edge_weight[k][e]) << "\"];" << std::endl;
	}
      dotout << "\t}" << std::endl << std::endl;
    }
  
  dotout << "}" << std::endl;
  
  dotout.close();
}



void mtreemix_DOT(array< graph >& G, array< map<node,string> >& mut, array< map<edge,string> >& edge_label, array< map<int,node> >& node_no, array<string>& tree_label, int uniform_noise, char *filestem)
{
  // produces output in dot format (graphviz package)
  // see http://www.research.att.com/sw/tools/graphviz/
  // this version takes general edge labels (preferred over the previous approach)

  node v;  edge e;

  char filename[255];
  sprintf(filename, "%s.dot", filestem);
  std::ofstream dotout(filename);
  
  dotout << "digraph \"" << filestem << "\" {" << std::endl << std::endl;

  // graph attributes:
  dotout << "\tsize = \"7.5, 10\";" << std::endl;
  dotout << "\tratio = \"auto\";" << std::endl;
  dotout << "\tnodesep = \"0.1\";" << std::endl;
  dotout << "\tranksep = \"0.3 equally\";" << std::endl << std::endl;
  
  // node and edge attributes:
  dotout << "\tnode [shape=\"plaintext\", height=\"0.3\", fontsize=\"12\", style=\"filled\", fillcolor=\"white\"];" << std::endl;
  dotout << "\tedge [fontsize=\"8\"];" << std::endl << std::endl;

  int k = 0;  // noise cluster:
  // cluster name:
  dotout << "\tsubgraph cluster" << k << " {" << std::endl << std::endl;
  // cluster label:
  dotout << "\t\tlabel=\"" << tree_label[k] << "\";" << std::endl;
  // cluster attributes:
  //dotout << "\t\tstyle=filled;" << std::endl;
  //dotout << "\t\tcolor=lightgrey;" << std::endl << std::endl;
      
  list<edge> E = G[k].all_edges();

  // nodes:
  v = node_no[k][0];
  dotout << "\t\t \"" << v << "\" [label=\"" << mut[k][v] << "\"];" << std::endl;
  forall(e, E)
    {
      v = target(e);
      dotout << "\t\t \"" << v << "\" [label=\"" << mut[k][v] << "\"];" << std::endl;
    }
  dotout << std::endl;
      
  // edges:
  while(E.size() >= 1)
    {
      e = E.pop();
      node s = source(e);
      node t = target(e);
      dotout << "\t\t \"" << s << "\" -> \"" << t << "\"";
      if ((! uniform_noise) || (E.size() == 0))
        dotout << " [label=\"" << edge_label[k][e] << "\"];";
      dotout << std::endl;
    }
  dotout << "\t}" << std::endl << std::endl;

  // k >= 1 :
  for (k=G.low()+1; k<=G.high(); k++)
    {
      // cluster name:
      dotout << "\tsubgraph cluster" << k << " {" << std::endl << std::endl;
      // cluster label:
      dotout << "\t\tlabel=\"" << tree_label[k] << "\";" << std::endl;
      // cluster attributes:
      //dotout << "\t\tstyle=filled;" << std::endl;
      //dotout << "\t\tcolor=lightgrey;" << std::endl << std::endl;
      
      list<edge> E = G[k].all_edges();

      // nodes:
      v = node_no[k][0];
      dotout << "\t\t \"" << v << "\" [label=\"" << mut[k][v] << "\"];" << std::endl;
      forall(e, E)
	{
	  v = target(e);
	  dotout << "\t\t \"" << v << "\" [label=\"" << mut[k][v] << "\"];" << std::endl;
	}
      dotout << std::endl;
      
      // edges:
      forall(e, E)
	{
	  node s = source(e);
	  node t = target(e);
	  dotout << "\t\t \"" << s << "\" -> \"" << t << "\"";
	  dotout << " [label=\"" << edge_label[k][e] << "\"];" << std::endl;
	}
      dotout << "\t}" << std::endl << std::endl;
    }
  
  dotout << "}" << std::endl;
  
  dotout.close();
}



list<int> missing_indices(integer_vector pat)
{
  list<int> index;

  for (int j=0; j<pat.dim(); j++)
    if (pat[j] == -1)
      index.append(j);

  return index;
}



double mtreemix_EM(array<string>& profile, integer_matrix& pattern, int K, int M, vector& alpha, array< graph >& G, array< map<int,node> >& node_no, array< map<node,string> >& event, array< map<edge,double> >& cond_prob, integer_matrix& pat_hat, matrix& resp, matrix& wlike, int uniform_noise, double eps, int special_weighing)
{
  // EM algorithm for fitting the parameters 
  // of a K-branchings mixture model
  
  int N = pattern.dim1();  // sample size
  int L = pattern.dim2();  // pattern length
  
  array< matrix >              P(K);     // weighted pair probabilities
  array< edge_array<double> >  dist(K);  // weight functional
  array< list<edge> >          B(K);     // branching

  vector oneN = ones(N);
  vector oneK = ones(K);

  vector freq = event_freq(pattern);  // marginal probabilities of events

  // initial guess:

  // method a) directly guess missing data:
  //guess_missing_data(pattern, pat_hat, freq);

  // method b) modified k-means:
  for (int i=0; i<N; i++)
    for (int j=0; j<L; j++)
      pat_hat(i,j) = pattern(i,j);  // copy
  guess_resp(pat_hat, K, M, resp);
  
  // iterate:
  double logL, logL_new = -DBL_MAX;  // log-likelihood
  int t = 0;  // counter
  do
    {
      t++;
      logL = logL_new;
      
      // "M-like" STEP
      // -------------
      for (int k=0; k<K; k++)
	{
	  // M1) estimate resp-weighted pair probabilities P[k](j1,j2)
	  P[k] = pair_probs(pat_hat, resp.row(k));

	  // M2) estimate branching model B[k] from P[k]
	  mgraph_init(profile, G[k], event[k], dist[k], node_no[k]);

	  if (k == 0)  // noise cluster: star topology
	    {  
	      B[0] = STAR(node_no[0][0]);
	      UNCOVER_BRANCHING(G[0], B[0]);
	      mgraph_weigh(P[0], profile, G[0], dist[0], cond_prob[0], node_no[0], -1.0, special_weighing);
	      
	      if (uniform_noise) {  // set cond_prob to average for all edges:
		int L = profile.size();
		double avg = 0.0;
		edge e;
		forall(e, B[0])
		  avg += cond_prob[0][e];
		avg /= L;  // now, avg == mean_{e in B[0]} (cond_prob[0][e])
		forall(e, B[0])
		  cond_prob[0][e] = avg;
	      }
	    }
	  
	  else  // maximum weight branching topology
	    {
	      mgraph_weigh(P[k], profile, G[k], dist[k], cond_prob[k], node_no[k], eps, special_weighing);
	      B[k] = MAX_WEIGHT_BRANCHING(G[k], event[k], dist[k]);
	      UNCOVER_BRANCHING(G[k], B[k]);
	    }
	  
	  // M3) estimate mixture parameters alpha[k]
	  alpha[k] = (oneN * resp.row(k)) / N;  // == sum(resp(k,:) / N
	}
      
      // E STEP
      // ------
      
      // E1) estimate missing data
      for (int i=0; i<N; i++)
	{
	  // list of indices of missing values in i-th sample:
	  list<int> index_m = missing_indices(pattern.row(i));
	  int L_m = index_m.size(); // number of missing values
	  
	  if (L_m > 0)  // anything missing at all?
	    {
	      double prob, prob_max = 0.0;
	      integer_vector pat_max(L);
	      
	      if (L_m < MISSING_DATA_MAX)  // try all possible patterns:
		{
		  for (int h=0; h<pow2(L_m); h++)
		    {
		      // generate each pattern ...
		      integer_vector pat_m = idx2pat(h, L_m);
		      integer_vector pat = pattern.row(i);
		      int j, j_m = 0;
		      forall(j, index_m)
			pat[j] = pat_m[j_m++];
		      
		      // ... and compute its likelihood:
		      prob = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
		      if (prob > prob_max)  // improvement?
			{
			  prob_max = prob;
			  pat_max = pat;  // remember best one
			}
		    }
		}
	      else  //  vector freq = event_freq(pattern);  // marginal probabilities of events too many possibilities --> heuristic:
		{
		  for (int r=0; r<MD_ML_RUNS; r++)  // repeat
		    {
		      int j;
		      
		      // initially, draw randomly according to freq:
		      integer_vector pat = pattern.row(i);
		      forall(j, index_m)
			pat[j] = (myrand() < freq[j]) ? 1 : 0;
		      prob = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
		      if (prob > prob_max)  // improvement?
			{
			  prob_max = prob;
			  pat_max = pat;
			}

		      // subsequently, improve iteratively:
		      int improve;  // number of improvements
		      do
			{
			  improve = 0; 
			  index_m.permute();  // in random order,
			  forall(j, index_m)  // consider all missing values
			    {
			      pat[j] = (pat[j] + 1) % 2;  // switch at position j
			      prob = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
			      if (prob > prob_max)  // improvement?
				{
				  prob_max = prob;
				  pat_max = pat;  // remeber best one
				  improve++;
				}
			      else
				pat[j] = (pat[j] + 1) % 2;  // undo j-switch
			    }	      
			}
		      while (improve > 0);
		    } 
		}
	      
	      for (int j=0; j<L; j++)  // pat_hat is the ML estimate:
		pat_hat(i,j) = pat_max[j];
	    }
	}
      
      // E2) estimate responsibilities resp(k,i)
      for (int k=0; k<K; k++)
	for (int i=0; i<N; i++)
	  wlike(i,k) = alpha[k] * mtree_like(pat_hat.row(i), G[k], node_no[k], cond_prob[k]);  // "weighted likelihood"
      
      logL_new = 0.0;
      for (int i=0; i<N; i++)
	{
	  double sum = oneK * wlike.row(i);  // sum of weighted likelihoods
	  if (sum <= 0.0)
	    {
	      std::cerr << "EM aborted. Sample no. " << i + 1
			<< " [" << pat_hat.row(i) << "] "
			<< "has likelihood zero!" << std::endl;
	      // mtreemix_save(alpha, G, cond_prob, node_no, "Lzero");  // save model for diagnostics
	      _Rtreemix_exit(1);
	    }
	  
	  for (int k=0; k<K; k++)
	    resp(k,i) = wlike(i,k) / sum;  // responsibility
	  
	  logL_new += log(sum);  // i.i.d. assumption
	}

      //printf("iter = %3d,  logL = %f\n", t, logL_new);  // progress line
    }  while (((logL_new - logL) >= EM_ACCURACY) & (t < EM_MAX_ITER));

  
  
  return logL_new;
}



double mtreemix_E_step(integer_matrix& pattern, int K, vector& alpha, array< graph >& G, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob, integer_matrix& pat_hat, matrix& resp, matrix& wlike)
{
  // E step
  // ------

  int N = pattern.dim1();  // sample size
  int L = pattern.dim2();  // pattern length
  
  vector oneN = ones(N);
  vector oneK = ones(K);
      
  vector freq = event_freq(pattern);  // marginal probabilities of events

  for (int i=0; i<N; i++)
    for (int j=0; j<L; j++)
      pat_hat(i,j) = pattern(i,j);  // copy

  // E1) estimate missing data
      for (int i=0; i<N; i++)
        {
	  // list of indices of missing values in i-th sample:
	  list<int> index_m = missing_indices(pattern.row(i));
	  int L_m = index_m.size(); // number of missing values
	  
	  if (L_m > 0)  // anything missing at all?
	    {
	      double prob, prob_max = 0.0;
	      integer_vector pat_max(L);
	      
	      if (L_m < MISSING_DATA_MAX)  // try all possible patterns:
		{
		  for (int h=0; h<pow2(L_m); h++)
		    {
		      // generate each pattern ...
		      integer_vector pat_m = idx2pat(h, L_m);
		      integer_vector pat = pattern.row(i);
		      int j, j_m = 0;
		      forall(j, index_m)
			pat[j] = pat_m[j_m++];
		      
		      // ... and compute its likelihood:
		      prob = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
		      if (prob > prob_max)  // improvement?
			{
			  prob_max = prob;
			  pat_max = pat;  // remember best one
			}
		    }
		}
	      else  // too many possibilities --> heuristic:
		{
		  for (int r=0; r<MD_ML_RUNS; r++)  // repeat
		    {
		      int j;
		      
		      // initially, draw randomly according to freq:
		      integer_vector pat = pattern.row(i);
		      forall(j, index_m)
			pat[j] = (myrand() < freq[j]) ? 1 : 0;
		      prob = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
		      if (prob > prob_max)  // improvement?
			{
			  prob_max = prob;
			  pat_max = pat;
			}

		      // subsequently, improve iteratively:
		      int improve;  // number of improvements
		      do
			{
			  improve = 0; 
			  index_m.permute();  // in random order,
			  forall(j, index_m)  // consider all missing values
			    {
			      pat[j] = (pat[j] + 1) % 2;  // switch at position j
			      prob = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
			      if (prob > prob_max)  // improvement?
				{
				  prob_max = prob;
				  pat_max = pat;  // remeber best one
				  improve++;
				}
			      else
				pat[j] = (pat[j] + 1) % 2;  // undo j-switch
			    }	      
			}
		      while (improve > 0);
		    } 
		}
	      
	      for (int j=0; j<L; j++)  // pat_hat is the ML estimate:
		pat_hat(i,j) = pat_max[j];
	    }
	}
      
      // E2) estimate responsibilities resp(k,i)
      for (int k=0; k<K; k++)
	for (int i=0; i<N; i++)
	  wlike(i,k) = alpha[k] * mtree_like(pat_hat.row(i), G[k], node_no[k], cond_prob[k]);  // "weighted likelihood"
      
      double logL_new = 0.0;
      for (int i=0; i<N; i++)
	{
	  double sum = oneK * wlike.row(i);  // sum of weighted likelihoods
	  if (sum <= 0.0)
	    {
	      std::cerr << "E-step aborted. Sample no. " << i + 1
			<< " [" << pat_hat.row(i) << "] "
			<< "has likelihood zero!" << std::endl;
	      // mtreemix_save(alpha, G, cond_prob, node_no, "Lzero");  // save model for diagnostics
	      _Rtreemix_exit(1);
	    }
	  
	  for (int k=0; k<K; k++)
	    resp(k,i) = wlike(i,k) / sum;  // responsibility
	  
	  logL_new += log(sum);  // i.i.d. assumption
	}

  return logL_new;
}



double mtreemix_loglike(integer_matrix& pattern, int K, vector& alpha, array<graph>& G, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob)
{
  // Calculate the log-likelihood function  log L(data|model)

  int N = pattern.dim1();  // sample size
  double logL = 0.0;       // log L(data|model)
  matrix wlike(N,K);       // "weighted likelihood" of a sample

  for (int i=0; i<N; i++)
    {
      double sum = 0.0;
      for (int k=0; k<K; k++)
	sum += alpha[k] * mtree_like(pattern.row(i), G[k], node_no[k], cond_prob[k]);

      if (sum <= 0.0)
	{
	  std::cerr << "Warning: The sample: [" << pattern.row(i) << "] has likelihood zero!" << std::endl;
	}
      
      logL += log(sum);
    }
  
  return logL;
}



double mtreemix_prob(integer_vector& pat, int K, vector& alpha, array<graph>& G, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob)
{
  // Calculate the probability (likelihood) of a pattern

  double prob = 0.0;

  for (int k=0; k<K; k++)
    prob += alpha[k] * mtree_like(pat, G[k], node_no[k], cond_prob[k]);

  return prob;
}



void mtreemix_fit0(array<string>& profile, integer_matrix& pattern, vector& alpha, array<graph>& G, array< map<int,node> >& node_no, array< map<node,string> >& event, array< map<edge,double> >& cond_prob, vector& resp, int uniform_noise, int special_weighing)
{
  // Fit independence/noise model (star topology)
  
  alpha[0] = 1.0;

  matrix P = pair_probs(pattern, resp);

  edge_array<double> dist;
  mgraph_init(profile, G[0], event[0], dist, node_no[0]);
  mgraph_weigh(P, profile, G[0], dist, cond_prob[0], node_no[0], -1.0, special_weighing);

  list<edge> B = STAR(node_no[0][0]);
  UNCOVER_BRANCHING(G[0], B);

  if (uniform_noise)  // noise model: set cond_prob to average for all edges
    {
      int L = profile.size();
      double avg = 0.0;
      edge e;
      forall(e, B)
        avg += cond_prob[0][e];
      avg /= L;  // now, avg == mean_{e in B[0]} (cond_prob[0][e])
      forall(e, B)
        cond_prob[0][e] = avg;
    }
}



void mtreemix_fit1(array<string>& profile, integer_matrix& pattern, vector& alpha, array<graph>& G, array< map<int,node> >& node_no, array< map<node,string> >& event, array< map<edge,double> >& cond_prob, vector& resp, double eps, int special_weighing)
{
  // Fit single tree model
  
  alpha[0] = 1.0;

  matrix P = pair_probs(pattern, resp);

  edge_array<double> dist;
  mgraph_init(profile, G[0], event[0], dist, node_no[0]);
  mgraph_weigh(P, profile, G[0], dist, cond_prob[0], node_no[0], eps, special_weighing);

  list<edge> B = MAX_WEIGHT_BRANCHING(G[0], event[0], dist);
  UNCOVER_BRANCHING(G[0], B);

  // keep root node arborescence only:
  list<edge> bfs = mtree_bfs(G[0], node_no[0][0]);
  UNCOVER_BRANCHING(G[0], bfs);
}



double mtreemix_fit(array<string>& profile, integer_matrix& pattern, int K, int M, vector& alpha, array<graph>& G, array< map<int,node> >& node_no, array< map<node,string> >& event, array< map<edge,double> >& cond_prob, integer_matrix& pat_hat, matrix& resp, int uniform_noise, double eps, int special_weighing)
{
  // Fit a K-branchings mixture model to the given data

  int N = pattern.dim1();        // sample size
  matrix wlike(N,K);             // "weighted likelihood"
  
  // run EM:
  double logL = mtreemix_EM(profile, pattern, K, M, alpha, G, node_no, event, cond_prob, pat_hat, resp, wlike, uniform_noise, eps, special_weighing);

  for (int k=0; k<K; k++)  // keep root node arborescences only
    {
      list<edge> bfs = mtree_bfs(G[k], node_no[k][0]);
      UNCOVER_BRANCHING(G[k], bfs);
    }
  return logL;
  
  // [do something with wlike (before returning) ....]
}



void mtreemix_save(vector& alpha, array< graph >& G, array< map<edge,double> >& cond_prob, array< map<int,node> >& node_no, char *filestem)
{
  int L = G[0].number_of_nodes();

  matrix E(L,L);
  
  char filename[255];
  sprintf(filename, "%s.model", filestem);
  std::ofstream mtreemix(filename);
  if (! mtreemix)
    {
      std::cerr << "Can't open output file -- " << filename << std::endl;
      _Rtreemix_exit(1);
    }

  mtreemix << alpha << std::endl;
  for (int k=G.low(); k<=G.high(); k++)
    {
      // build adjacency matrix E:
      for (int i1=0; i1<L; i1++)
	for (int i2=0; i2<L; i2++)
	  {
	    edge e = edge_between(node_no[k][i1], node_no[k][i2]);
	    E(i1,i2) = (e == NULL /*nil*/) ? 0.0 : cond_prob[k][e];
	  }
      
      mtreemix << E;
    }
  
  mtreemix.close();
}



array<string> mtreemix_load(vector& alpha, array< graph >& G, array< map<node,string> >& event, array< map<edge,double> >& cond_prob, array< map<int,node> >& node_no, char *filestem)
{
  matrix E;
  array<string> profile;

  char filename[255];
  sprintf(filename, "%s.model", filestem);
  std::ifstream mtreemix(filename);
  if (! mtreemix)
    {
      std::cerr << "Can't open input file -- " << filename << std::endl;
      _Rtreemix_exit(1);
    }
  
  mtreemix >> alpha >> std::ws;
  int K = alpha.dim();

  G.resize(K);
  event.resize(K);
  cond_prob.resize(K);
  node_no.resize(K);
  
  for (int k=0; k<K; k++)
    {
      // read adjacency matrix E:
      mtreemix >> E;
      
      if (k == 0)  // initialize profile:
	profile = load_profile(filestem, E.dim1());
      
      // add nodes:
      node_no[k].clear();
      event[k].clear();
      for (int i=0; i<E.dim1(); i++)
	{
	  node v = G[k].new_node();
	  node_no[k][i] = v;
	  event[k][v] = profile[i];
	}

      // add edges:
      cond_prob[k].clear();
      for (int i1=0; i1<E.dim1(); i1++)
	for (int i2=0; i2<E.dim2(); i2++)
	  {
	    if (E(i1,i2) > 0)
	      {
		edge e = G[k].new_edge(node_no[k][i1], node_no[k][i2]);
		cond_prob[k][e] = E(i1,i2);
	      }
	  }
    }

  mtreemix.close();

  return profile;
}



integer_matrix mtreemix_draw(int L, vector& alpha, array<graph>& G, array< map<edge,double> >& cond_prob, array< map<int,node> >& node_no, int n, int include_noise)
{
  int K = alpha.dim();

  // map nodes onto event indices:
  array< map<node,int> > no_of_node(K);
  for (int k=0; k<K; k++)
    for (int j=0; j<L; j++)
      no_of_node[k][node_no[k][j]] = j;

  integer_matrix pattern(n,L);
  integer_vector pat(L);

  int i = 0;
  while (i < n)
    {
      // draw branching number k:
      int k = discrand(alpha);

      // draw sample from branching k:
      if (include_noise | (k > 0))
	{
	  pat = mtree_draw(L, G[k], node_no[k][0], cond_prob[k], no_of_node[k]);
	  for (int j=0; j<L; j++)
	    pattern(i,j) = pat[j];
	  i++;
	}
    }

  return pattern;
}		    
  


array< map<edge,double> > waiting_times(array< map<edge,double> >& cond_prob, int sampling_mode, double sampling_param)
{
  // compute waiting times (parameter lambda
  // of the exponential distribution) on each edge

  int K = cond_prob.size();
  array< map<edge,double> > lambda(K);
  edge e;
  
  switch(sampling_mode)
    {
    case CONSTANT: // sampling_param == t_0
      // lambda = (lambda_T * p) / (1 - p)
      for (int k=cond_prob.low(); k<=cond_prob.high(); k++)
	{
	  forall_defined(e, cond_prob[k])  // lambda = -log(1 - p) / t_0
	    lambda[k][e] = -log(1 - cond_prob[k][e]) / sampling_param;
	}
      break;

    case EXPONENTIAL:  // sampling_param == lambda_s
      for (int k=cond_prob.low(); k<=cond_prob.high(); k++)
	{
	  forall_defined(e, cond_prob[k])  {// lambda = (lambda_s * p) / (1 - p)
	    lambda[k][e] = (sampling_param * cond_prob[k][e]) / (1 - cond_prob[k][e]);
	  }
	}
      break;
    
    default :
      std::cerr << "Unknown sampling_mode -- " << sampling_mode << std::endl;
      _Rtreemix_exit(1);
    }
  return lambda;
}





vector mtreemix_distr(int L, vector& alpha, array<graph>& G, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob)
{
  // compute entire distribution induced by the model

  int K = alpha.dim();
  
  // generate all patterns:
  int n = pow2(L-1);
  vector prob(n);  // vector of probabilities of all patterns

  for (int i=0; i<n; i++)
    {
      integer_vector pat = index2pattern(i, L);
      prob[i] = mtreemix_prob(pat, K, alpha, G, node_no, cond_prob);
    }

  return prob;
}




int has_missing_values(integer_matrix& pat)
{
  int N = pat.dim1();  // sample size
  int L = pat.dim2();  // pattern length

  for (int i=0; i<N; i++)
    for (int j=0; j<L; j++)
      if (pat(i,j) == -1)
	return 1;

  return 0;
}

  

//  Author: Junming Yin
//  (c) Max-Planck-Institut fuer Informatik, 2005


void alpha_random(vector& alpha, int& K)
{
  // randomly generate weight of each single tree

  double sum = 0.0; 
  
  for (int i=0;i<K;i++)
    {
      alpha[i] = (double) rand()/ (double) RAND_MAX;
      sum +=alpha[i]; 
    }
      
  //normalize alpha, make sure sum up to one 
  for (int i=0;i<K;i++)
    alpha[i] = alpha[i]/sum;
    
}


void mtreemix_random(int K, int L, array<string>& profile, vector& alpha, array <graph>& G, array< map<node,string> >& event, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob, int& star, int& uniform, double& min, double& max)
{

  //randomly generate a mixture tree
  //K: number of tree components
  //L: number of nodes in each tree
  //star = 1: the first component is a star
  //uniform = 1: the first component is a uniform star
  //min & max: the minimum and the maximum of randomly drawn conditional probability associated with edges
  
  profile.resize(L);
  for (int j=0; j<L; j++)
    profile[j] = tostring("%d", j);
  
  //generate value for alpha
  alpha_random(alpha, K);  
  
  // the string of length (L-2) used to generate the tree
  vector s(L-2);

  if (star == 1)
    {
      //the sequence to generate the star
      for(int i=0;i<L-2;i++)
	s[i] = 0;
      mtree_random(L, profile, G[0], event[0], node_no[0], cond_prob[0], s, 0, min, max);

      if(uniform == 1) {
        // a uniform star
        edge e;
        double weight = ((double) rand()/ (double) RAND_MAX)*(max-min) + min; 
        forall_edges(e, G[0])
	  cond_prob[0][e] = weight;
      }
      //randomly generate other K-1 tree components
      for (int k=1; k<K; k++)
	mtree_random(L, profile, G[k], event[k], node_no[k], cond_prob[k], s, 1, min, max);
    }
  else
    {
      //randomly generate K tree components
      for (int k=0; k<K; k++)
	mtree_random(L, profile, G[k], event[k], node_no[k], cond_prob[k], s, 1, min, max);
    }
}


double mtree_state(map<edge, int>& state, integer_vector& pattern, graph& G, map<int,node>& node_no, map<edge,double>& prob_cond)
{
  // Maps a sample on a mutagenetic tree and computes its likelihood.
  // G has to be a branching!
  // Return the state for each edge
  // state[e] = 1:  prob_cond[e] appears in the calculation of P(pattern)
  // state[e] = -1: (1 - prob_cond[e]) appears in the calculation of P(pattern)
  // state[e] = 0:  prob_cond[e] does not appear in the calculation of P(pattern)

  //initial state of each edge
  edge e;
  forall_edges(e, G)
    state[e] = 0;

  if (pattern[0] < 1)  // initial event has to occur!
    return 0.0;

  int L = pattern.dim();  // pattern length

  node v, w; 
  node_array<int> dist;
  queue<node> Q;

  // put pattern in set:
  node_set Pat(G);
  for (int j=0; j<L; j++)
    if (pattern[j] > 0)
      Pat.insert(node_no[j]);

  // initialize dist[] for BFS
  dist.init(G);
  forall_nodes(w, G)
    if (Pat.member(w))
      dist[w] = -1;
    else
      dist[w] = 0;  // Do not visit events that are not part of pattern!

  int visited_events = 1;
  double like = 1.0;

  // BFS on pattern nodes:
  node root = node_no[0];
  Q.append(root);
  dist[root] = 0;

  while (! Q.empty())
    {
      v = Q.pop();
      forall_out_edges(e, v)
	{
	  w = target(e);
	  if (dist[w] < 0)
	    {
	      Q.append(w);
	      dist[w] = dist[v] + 1;
	      visited_events++;
	      like *= prob_cond[e];
	      state[e] = 1;
	    }
	  else
	    {
	      state[e] = -1;
	      like *= (1 - prob_cond[e]);
	    }
	}
    }
  
  if (visited_events < Pat.size())  // ie. pattern does not fit onto branching
    like = 0.0;      // hence likelihood zero
  
  return like;
}


integer_vector myindex2pattern(int& num_nonzeros, int index, int L)
{
  // pattern \in 2^{1, ..., L} of a given index \in [0 .. 2^(L-1)]
  // Note: We always assume pat[0] == 1 !
  // num_nonzero: number of nonzeros (except the first one)  
  integer_vector pat(L);

  pat[0] = 1; 
  num_nonzeros = 0;

  for (int j=1; j<L; j++)
    {
      int mod2 = index % 2;
      num_nonzeros += (mod2 == 1);
      pat[j] = integer(mod2);
      index = (index - mod2) / 2;
    }
     
  return pat;
}


double power(double m, int n)
{
  // m^n

  double p = 1;

  for (int i=0; i<n; i++)
    p *= m;

  return p;
}

matrix mtreemix_distance(int L, int K1, array< graph >& G1, array< map<int,node> >& node_no1, int K2, array< graph >& G2, array< map<int,node> >& node_no2)
{
  // return the distance matrix of two mixture models
  // dist(k1,k2) = 0 means k1th component of G1 and k2th component of G2 are exactly the same
  // dist(k1,k2) = 1 means k1th component of G1 and k2th component of G2 are totally "different"

  // build the distance matrix
  matrix dist(K1,K2);

  for(int k1=0; k1 < K1; k1++)
    for(int k2=0; k2 < K2; k2++)
      dist(k1,k2) = mtree_distance(L, G1[k1], node_no1[k1], G2[k2], node_no2[k2]);

  return dist;  
}
