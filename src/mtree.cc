/*
  mtree.c

  Version : 1.2.00
  Author  : Niko Beerenwinkel

  (c) Max-Planck-Institut fuer Informatik, 2004-2005
*/



#include "include/mtree.h"


array<string> load_profile(char *filestem, int L)
{
  // load profile (list of events) from file

  array<string> profile;

  // load from file
  char filename[255];
  sprintf(filename, "%s.prf", filestem);
  std::ifstream prf(filename);
  if (prf)
    {
      int j = 0;
      while (prf)
	{
	  //string s = read_line(prf);
	    string s;
	    getline(prf,s);
	  if (s.length() > 0)
	    {
	      profile.resize(j+1);
	      profile[j++] = s;
	    }
	}
      prf.close();
      
      if (j != L)
	{
	  std::cerr << "Number of profile labels does not coincide with number of data columns and/or model dimensions!" << std::endl;
	  exit(1);
	}
    }

  else  // generic labels:
    {
      profile.resize(L);
      for (int j=0; j<L; j++)
	profile[j] = tostring("%d", j);
    }

  return profile;
}


void save_profile(array<string>& profile, char* filestem)
{
  // save profile
  
  int L = profile.size();

  char filename[255];
  sprintf(filename, "%s.prf", filestem);
  std::ofstream prf(filename);
  if (! prf)
    {
      std::cerr << "Can't open output file -- " << filename << std::endl;
      exit(1);
    }
  
  for (int j=0; j<L; j++)
     prf << profile[j] << std::endl;

  prf.close();

}



integer_matrix load_pattern(char *filestem)
{
  integer_matrix pattern;
  
  char filename[255];
  sprintf(filename, "%s.pat", filestem);

  std::ifstream patfile(filename);
  if (! patfile)
    {
      std::cerr << "Can't open input file -- " << filename << std::endl;
      exit(1);
    }
  patfile >> pattern;
  patfile.close();
 
  return pattern;
}



void save_pattern(integer_matrix& pattern, char *filestem)
{
  char filename[255];
  sprintf(filename, "%s.pat", filestem);

  std::ofstream patfile(filename);
  if (! patfile)
    {
      std::cerr << "Can't open output file -- " << filename << std::endl;
      exit(1);
    }
  patfile << pattern;
  patfile.close();
 
}



matrix pair_probs(integer_matrix& pat, vector& resp)
{
  // calculate all pairwise probabilities

  int N = pat.dim1();  // sample size
  int L = pat.dim2();  // pattern length
  
  matrix P(L, L);  // Note: P(0,:), P(:,0) and the diagonal
  // all give the marginal probability of the single event!

  for (int j1=0; j1<L; j1++)
    for (int j2=j1; j2<L; j2++)
      {
	int count = 0;          // number of full data pairs
	double wcount = 0.0;    // resp-weighted count
	double resp_sum = 0.0;  // sum of responsibilities ( == N_k ) 
	
	for (int i=0; i<N; i++)
	  if ((pat(i,j1) >= 0) && (pat(i,j2) >= 0))  // ie. both non-missing
	    {
	      count++;
	      //wcount += resp[i] * (to_double(pat(i,j1)) * to_double(pat(i,j2)));
	      wcount += resp[i] * (double)(pat(i,j1) * pat(i,j2));
	      resp_sum += resp[i];
	    }
	
	if (count == 0)  // no data in the column pair (j1,j2)
	  {
	    std::cerr << "Warning: No data in column pair (" << j1 << ", " << j2 << ")! Assuming independence." << std::endl;
	    wcount = P(0,j1) * P(0,j2);  // assuming independence
	  }
	
	P(j1,j2) = P(j2,j1) = (wcount / resp_sum) + PSEUDO_COUNT;
      }
  
  return P;
}



void mgraph_init(array<string>& profile, graph& G, map<node,string>& event, edge_array<double>& dist, map<int,node>& node_no)
{

  node v, w;

  G.del_all_nodes();
  G.del_all_edges();
  event.clear();  
  node_no.clear();

  // for each event a node
  for (int i=0; i<profile.size(); i++)
    {
      v = G.new_node();
      node_no[i] = v;
      event[v] = profile[i];
    }


  forall_nodes(v, G)
    forall_nodes(w, G)
      if (v != w && w != node_no[0])
	  G.new_edge(v, w);
      

  dist.init(G); 


}



void mgraph_weigh(matrix& P, array<string>& profile, graph& G, edge_array<double>& dist, map<edge,double>& prob_cond, map<int,node>& node_no, double eps, int special_weighing)
{
  edge e;
  map<edge,double> w;
  double w_min = DBL_MAX;

  prob_cond.clear();

  // weigh edges:
  for (int i=0; i<profile.size(); i++)
    for (int j=1; j<profile.size(); j++) // omit edges into root==node_no[0]
	if ((e = edge_between(node_no[i], node_no[j])) != NULL)//nil)  // e = (i,j)
	{
	  prob_cond[e] = P(j,i) / P(i,i);  // == P(j|i)
	  if ((eps > 0) && (prob_cond[e] < eps))  // delete light edges
	    // use negative eps to avoid this, eg. in the noise component
	    G.del_edge(e);
	  else
	    {
	      if ((special_weighing == 1) && (i == 0))
		w[e] = log(P(j,j));
	      else  // Desper weight functional:
		w[e] = log(P(i,j)) - log(P(i,i) + P(j,j)) - log(P(j,j));
	      w_min = std::min(w_min, w[e]);
	    }
	}
  
  forall_edges(e, G)
    dist[e] = 1.0 + -w_min + w[e];  // shift to positive reals

}



edge edge_between(node& v, node& w)
{
  edge e;
  graph *G = graph_of(v);

  forall_edges(e, (*G))
    {
      if (source(e)==v && target(e)==w)
	return e;
    }

  return NULL;//nil;
}



double mstar_like(int *pattern, int pattern_length, matrix& P)
{
  double like = 1.0;

  for (int k=0; k<pattern_length; k++)
    like *= (pattern[k] == 1) ? P(k,k) : (1 - P(k,k));
  
  return like;
}



list<edge> mtree_bfs(graph& G, node& root)
{
  // BFS search in G starting from root
  // return all visited edges.

  list<edge> L;
  node v, w;  edge e;
  node_array<int> dist;
  queue<node> Q;

  // initialize:
  dist.init(G);
  forall_nodes(w, G)
    dist[w] = -1;

  // BFS:
  Q.append(root);
  dist[root] = 0;

  while (! Q.empty())
    {
      v = Q.pop();
      forall_out_edges(e, v)
	{
	  L.append(e);
	  w = target(e);
	  if (dist[w] < 0)
	    {
	      Q.append(w);
	      dist[w] = dist[v] + 1;
	    }
	}
    }
  
  return L;
}


double mtree_like(integer_vector& pattern, graph& G, map<int,node>& node_no, map<edge,double>& prob_cond)
{
  // Maps a sample on a mutagenetic tree and computes its likelihood.
  //
  // G has to be a branching!

  if (pattern[0] < 1)  // initial event has to occur!
    return 0.0;

  int L = pattern.dim();  // pattern length

  node v, w;  edge e;
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
	    }
	  else
	    like *= (1 - prob_cond[e]);
	}
    }
  
  if (visited_events < Pat.size())  // ie. pattern does not fit onto branching
    like = 0.0;      // hence likelihood zero
  
  return like;
}



array<int> permutation(int n)
{
  array<int> sigma(n);
  
  for (int k=0; k<n; k++)
    sigma[k] = k;
  sigma.permute();  // draw a permutation sigma \in S(n)

  return sigma;
}


list<edge> STAR(node& root)
{
  // return edges of star topology with given root
  
  edge e;
  list<edge> star;
  
  
  forall_out_edges(e, root)
      star.append(e);

  
  return star;
}



double myrand()
{
    // srand( (unsigned)time( NULL ) );
    return (double) rand() / (double) RAND_MAX;
}


int discrand(vector& P)
{
  int K = P.dim();

  // draw a random integer from {0, ..., K-1}
  // according to the discrete probability distribution P

  double s = myrand();
  int k = 0;
  double cP = P[0];  // cumulative P
  
  while ((cP < s) & (k < K-1))
    cP += P[++k];

  return k;
}


double expcdf(double lambda)
{
  // exponential CDF
  
  double Log;
  double u = myrand();

  while (! ((Log = -log(1 - u)) <= DBL_MAX))  // avoid log(zero)
    u = myrand();
  
  return Log / lambda;
}



integer_vector mtree_draw(int L, graph& G, node& root, map<edge,double>& cond_prob, map<node,int>& no_of_node)
{
  // draw a sample from branching G according to cond_prob

  integer_vector pat(L);

  // data structures:
  node v, w;  
  edge e;
  node_array<int> dist;
  queue<node> Q;

  // initialize:
  dist.init(G);
  forall_nodes(w, G)
    dist[w] = -1;

  Q.append(root);
  dist[root] = 0;
  pat[0] = 1;  // initial event

  // BFS:
  while (! Q.empty())
    {
      v = Q.pop();
      forall_out_edges(e, v)
	{
	  w = target(e);
	  if (dist[w] < 0)
	    {
	      double s = myrand();
	      if (s < cond_prob[e])  // i.e. draw e randomly
		{
		  Q.append(w);
		  pat[no_of_node[w]] = 1;
		}
	      dist[w] = dist[v] + 1;
	    }
	}
    }

  return pat;
}



double mtree_wait(graph& G, node& root, replaceleda::map<edge,double>& cond_prob, map<edge,double>& lambda, map<node,int>& no_of_node, double stime, integer_vector& cpattern)
{
  // draw a sample from branching G according to 
  // exponential(lambda) waiting times on E(G)

  int L = G.number_of_nodes();
  edge e, e1;
  
  // draw random waiting times
  map<edge,double> time;
  forall_defined(e, lambda)
    time[e] = expcdf(lambda[e]);

  forall_defined(e, cond_prob) 
    if(cond_prob[e] == 1.0)
      time[e] = 0.0;
  
  // initialize data structures
  edge cedge;  // current edge
  // current pattern (null event):
  cpattern[0] = 1;
  for (int j=1; j<L; j++)
    cpattern[j] = 0;
  double ctime = 0.0; // current time
  double wtime = 0.0; // waiting time of the pattern

  // initialize Q
  p_queue<double,edge> Q;  // edge queue, waiting time prioritized
  forall_out_edges(e, root) 
    Q.insert(time[e], e);
  
  
  // traverse (wait along the branching G):
  while ((! Q.empty()) && ctime < stime)  // current time < sampling time ?
    {
      pq_item min_it = Q.find_min();  // edge with minimal waiting time
      //pq_elem<double,edge> min_it = Q.find_min();
      cedge = Q.inf(min_it);
      Q.del_item(min_it);  // remove from Q
      
      forall_out_edges(e, target(cedge))  // insert child edges of cedge
	{
          if(cond_prob[e] == 1.0) {
            cpattern[no_of_node[target(e)]] = 1;  // update cpattern
            time[e] += time[cedge];
       
            forall_out_edges(e1, target(e)) {
              time[e1] += time[e];
              Q.insert(time[e1], e1);
            }
          } else {
	    time[e] += time[cedge];  // waiting times of subsequent edges add up
	    Q.insert(time[e], e);
          }
	}
      
      ctime = std::max(ctime, time[cedge]);  // update ctime
      if (ctime <= stime)
	{
	  cpattern[no_of_node[target(cedge)]] = 1;  // update cpattern
	  wtime = ctime;
	}
    }
  return wtime;
}


void mtree_time(graph& G, node& root, replaceleda::map<edge,double>& cond_prob, map<edge,double>& lambda, map<node,int>& no_of_node, array< list<double> >& wtime)
{
  // draw a sample from branching G according to 
  // exponential(lambda) waiting times on E(G)
  // store patterns and their waiting times in wtime

  int L = G.number_of_nodes();
  edge e, e1;
  
  // draw random waiting times on edges
  map<edge,double> time;
  forall_defined(e, lambda)
    time[e] = expcdf(lambda[e]);

  forall_defined(e, cond_prob) 
    if(cond_prob[e] == 1.0)
      time[e] = 0.0;
    
  // initialize data structures
  edge cedge;                  // current edge
  integer_vector cpattern(L);  // current pattern
  cpattern[0] = 1;             // null event
  double ctime = 0.0;          // current time
  wtime[pattern2index(cpattern)].append(ctime);  // record initial null pattern

  // initialize Q
  p_queue<double,edge> Q;  // edge queue, waiting time prioritized
  forall_out_edges(e, root)
    Q.insert(time[e], e);
  
  // traverse (wait along the branching G):
  while (! Q.empty())
    {
      pq_item min_it = Q.find_min();  // edge with minimal waiting time
      cedge = Q.inf(min_it);
      Q.del_item(min_it);  // remove from Q
      
      forall_out_edges(e, target(cedge))  // insert child edges of cedge
	{
          if(cond_prob[e] == 1.0) {
            cpattern[no_of_node[target(e)]] = 1;  // update cpattern
            time[e] += time[cedge];
       
            forall_out_edges(e1, target(e)) {
              time[e1] += time[e];
              Q.insert(time[e1], e1);
            }
          } else {
	    time[e] += time[cedge];  // waiting times of subsequent edges add up
	    Q.insert(time[e], e);
          }
	}
      
      cpattern[no_of_node[target(cedge)]] = 1;  // update cpattern
      ctime = std::max(ctime, time[cedge]);          // update ctime
      wtime[pattern2index(cpattern)].append(ctime);  // record waiting time
    }

}



int pow2(int n)
{
  // 2^n

  int p = 1;

  for (int i=0; i<n; i++)
    p *= 2;

  return p;
}



double log2(double x)
{
	return log((double)x)/log((double)2);
}


vector ones(int N)
{
  vector one(N);

  for (int i=0; i<N; i++)
    one[i] = 1.0;

  return one;
}



int pattern2index(integer_vector& pat)
{
  // index \in [0 .. 2^(L-1)] of a given pattern \in 2^{1, ..., L}

  int L = pat.dim();

  int index = 0;
  for (int j=1; j<L; j++)  // we always assume pat[0] == 1 !
    index += (pat[j] == 1) * pow2(j-1);

  return index;
}


integer_vector index2pattern(int index, int L)
{
  // pattern \in 2^{1, ..., L} of a given index \in [0 .. 2^(L-1)]
  // Note: We always assume pat[0] == 1 !
  
  integer_vector pat(L);

  pat[0] = 1; 
  for (int j=1; j<L; j++)
    {
      int mod2 = index % 2;
      pat[j] = integer(mod2);
      index = (index - mod2) / 2;
    }
     
  return pat;
}



int pat2idx(integer_vector& pat)  // not yet tested!!!
{
  // index \in [0 .. 2^(L-1)] of a given pattern \in 2^{1, ..., L}
  // no constraints here!

  int L = pat.dim();

  int index = 0;
  for (int j=0; j<L; j++)  // we always assume pat[0] == 1 !
    index += (pat[j] == 1) * pow2(j-1);

  return index;
}



integer_vector idx2pat(int index, int L)
{
  // pattern \in 2^{1, ..., L} of a given index \in [0 .. 2^(L-1)]
  // no constraints here!

  integer_vector pat(L);

  for (int j=0; j<L; j++)
    {
      int mod2 = index % 2;
      pat[j] = integer(mod2);
      index = (index - mod2) / 2;
    }
     
  return pat;
}



double nonnegmean(vector& X)
{
  // mean of non-negative entries of X

  vector is_nonneg = ones(X.dim());  // indicates non-negative entries in X
  int N = 0;                         // number of non-negative entries in X
  
  for (int i=0; i<X.dim(); i++)
    {
      if (X[i] >= 0.0)
	N++;
      else
	is_nonneg[i] = 0.0;
    }

  return (is_nonneg * X) / (double) N;
}


double nonnegmean(list<double>& L)
{
  // mean of non-negative entries of L

  int N = 0;         // number of non-negative entries in X
  double sum = 0.0;  // sum of non-negative entries in X
  
  double x;
  forall(x, L)
    {
      if (x >= 0.0)
	{
	  sum += x;
	  N++;
	}
    }

  return sum / (double) N;
}


double nonnegmean(integer_vector& X)
{
  // mean of non-negative entries of X

  vector XX(X.dim()); // double X
  vector is_nonneg = ones(X.dim());  // indicates non-negative entries in X
  int N = 0;                         // number of non-negative entries in X
  
  for (int i=0; i<X.dim(); i++)
    {
      if (X[i] >= 0)
	{
	    XX[i] = (double) X[i];
	  N++;
	}
      else
	is_nonneg[i] = 0.0;
    }

  return (N == 0) ? -1.0 : (is_nonneg * XX) / (double) N;
}



//  Author: Junming Yin
//  (c) Max-Planck-Institut fuer Informatik, 2005


void mtree_directed(graph& G, node& root, map<edge,double>& cond_prob, double min, double max)
{
  // traverse the undirected rooted tree with BFS and delete all the reverse edges

  node v, w;  edge e, rev_e;
  node_array<int> dist;  //mapping from node to int 
  queue<node> Q;

  // initialize:
  dist.init(G);
  forall_nodes(w, G)
    dist[w] = -1;

  // BFS:
  Q.append(root);
  dist[root] = 0;

  //  std::cout << G.number_of_edges() << "\n";
  while (! Q.empty()) //Q is not empty   Q.empty() == false
    {
      v = Q.pop();
      forall_out_edges(e, v)  
    	{
	     w = target(e);
	     rev_e = edge_between(w,v);
	     G.del_edge(rev_e);  //delete the reverse edge
	     cond_prob[e] =  ((double) rand()/ (double) RAND_MAX)*(max-min) + min; //randomly draw from min to max
	     if (dist[w] < 0)
	      {
	        Q.append(w);
   	        dist[w] = dist[v] + 1;
	      }
    	}
    }
}


void mtree_random(int L, array<string>& profile, graph& G, map<node,string>& event, map<int,node>& node_no, map<edge,double>& cond_prob, vector& s, int random, double min, double max)
{
  // generate a random directed tree

  // vector s: string of length L-2 use to generate the topology of the tree
  // random = 1: randomly generate sequence s
  // random = 0: use a specified sequence s
  // min & max: the minimum and the maximum of randomly drawn conditional probability

  // the degree of each node
  // deg(j) = #(j appears in s) + 1
  vector deg = ones(L); 

  // the priority quene storing all the nodes with degree 1
  p_queue<int,node> Q;

  node v, v1, v2;
  int i;

  if (random==1) //randomly generate a tree
    {
      for (i=0;i<L-2;i++)
	{
	  s[i] = rand() % L;  //random draw integer from 0 to L-1
	  deg[(uint) s[i]]++; // increase the degree of that node
	}
    }
  else // use the specified vector a to generate the tree
    {
      for (i=0;i<L-2;i++)
	  deg[(uint)s[i]]++; // increase the degree of that node
    }
    
  //intial the graph
  G.del_all_nodes(); G.del_all_edges();
  event.clear();
  node_no.clear();

  // for each event a node
  for(i=0;i<L;i++)
  {
    v = G.new_node();
    node_no[i] = v;
    event[v] = profile[i];
  }  

  //insert deg[i] = 1 to priority queue
  for (i=0; i<L; i++)   
     if(deg[i]==1)
       Q.insert(i,node_no[i]);

  //build the graph
  for (i=0;i<L-2;i++)
  {
   //retrieve the node with minimal number in priority queue

   pq_item min_it = Q.find_min();  // node with minimal number
   v = Q.inf(min_it);
   Q.del_item(min_it);  // remove from Q

   if (--deg[(uint)s[i]]==1) // decrease the degree
       Q.insert((uint)s[i],node_no[(uint)s[i]]);

   // make a new edge
   G.new_edge(node_no[(uint)s[i]],v);
   G.new_edge(v, node_no[(uint)s[i]]);

  }
  
  //last two items in the Q
  pq_item min_it = Q.find_min();
  v1 = Q.inf(min_it);
  Q.del_item(min_it);

  min_it = Q.find_min();
  v2 = Q.inf(min_it);
  Q.del_item(min_it);

  // make a new edge
  G.new_edge(v1,v2);
  G.new_edge(v2,v1);

  //make the tree directed by BFS
  mtree_directed(G, node_no[0], cond_prob, min, max);
}


double infinity_norm(int L, integer_matrix m)
{
  // caculate the infinity norm (maximum absolute row norm) of a matrix
  // || A ||_{\inf} = max_i \sum_j |A_{ij}|

  double max = 0;

  for (int i=0; i < L; i++)
    {
      double sum_row = 0.0;
      for(int j=0; j < L; j++)
	  sum_row = sum_row + (double)abs(m(i,j));

      if(sum_row > max)
      	max = sum_row;
    }
      
  return max;
}


double mtree_distance(int L, graph& G1, map<int,node>& node_no1, graph& G2, map<int,node>& node_no2)
{
  // return the distance between two tree components
  // dist = 0 means G1 and G2 are exactly the same
  // dist = 1 means G2 and G2 are totally different

  double dist;
  integer_matrix adj1(L,L), adj2(L,L); // the adjacency matrix
  edge e;

  map<node,int> no_of_node1, no_of_node2;

  for (int j=0; j<L; j++)
      no_of_node1[node_no1[j]] = j;
  for (int j=0; j<L; j++)
      no_of_node2[node_no2[j]] = j;

  // build the adjacency matrix:    
  forall_edges(e, G1)
    { adj1(no_of_node1[source(e)],no_of_node1[target(e)]) = 1; }

  forall_edges(e, G2)
    { adj2(no_of_node2[source(e)],no_of_node2[target(e)]) = 1; }

  // caculate the distance
  //          || adj1 - adj2 ||_{\inf}
  // dist = --------------------------
  //                   L - 1

  integer_matrix  diff_adj = adj1 - adj2;
  double norm = infinity_norm(L, diff_adj);		
  dist = norm/(L-1); // normalized distance value

  return dist;
}
