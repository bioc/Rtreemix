/*
  mtree.h

  Version : 1.1
  Author  : Niko Beerenwinkel

  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#ifndef _MTREE_H
#define _MTREE_H

#include "max_weight_branch.h"

/*
  #include <LEDA/array.h>
  #include <LEDA/vector.h>
  #include <LEDA/matrix.h>
  #include <LEDA/integer.h>
  #include <LEDA/integer_vector.h>
  #include <LEDA/integer_matrix.h>
  #include <LEDA/list.h>
  #include <LEDA/map.h>
  #include <LEDA/p_queue.h>
  #include <LEDA/queue.h>
  #include <LEDA/graph.h>
*/
#include "replaceleda.hh"

// data:

array<replaceleda::string> load_profile(char *filestem, int N);   // profiles
void save_profile(array<replaceleda::string>& profile, char* filestem);

integer_matrix load_pattern(char *filestem);  // patterns
void save_pattern(integer_matrix& pat, char *filestem);


// joint probabilities:

#define PSEUDO_COUNT 1e-5

matrix pair_probs(integer_matrix& pat, replaceleda::vector& resp);


// mutational patterns:

int pattern2index(integer_vector& pat);
int pat2idx(integer_vector& pat);
// returns index of a given pattern

integer_vector index2pattern(int index, int L);
integer_vector idx2pat(int index, int L);
// returns pattern to a given index


// mutagenetic trees:

void mgraph_init(array<replaceleda::string>& profile, graph& G, replaceleda::map<node,replaceleda::string>& event, edge_array<double>& dist, replaceleda::map<int,node>& node_no);

void mgraph_weigh(matrix& P, array<replaceleda::string>& mut, graph& G, edge_array<double>& dist, replaceleda::map<edge,double>& edge_label, replaceleda::map<int,node>& node_no, double eps, int special_weighing);

replaceleda::list<edge> mtree_bfs(graph& G, node& root);
// breadth first search

double mtree_like(integer_vector& pattern, graph& G, replaceleda::map<int,node>& node_no, replaceleda::map<edge,double>& prob_cond);
// Compute the likelihood of a pattern in a given branching

double mstar_like(int *pattern, int pattern_length, matrix& P);
// Likelihood in star

integer_vector mtree_draw(int L, graph& G, node& root, replaceleda::map<edge,double>& cond_prob, replaceleda::map<node,int>& no_of_node);
// draw a sample from G according to cond_prob

double mtree_wait(graph& G, node& root, replaceleda::map<edge,double>& cond_prob, replaceleda::map<edge,double>& lambda, replaceleda::map<node,int>& no_of_node, double stime, integer_vector& cpattern);
// draw a sample from branching G according to exp(lambda) waiting times on E(G) at sampling time stime

void mtree_time(graph& G, node& root, replaceleda::map<edge,double>& cond_prob, replaceleda::map<edge,double>& lambda, replaceleda::map<node,int>& no_of_node, array< replaceleda::list<double> >& wtime);
// estimate pattern waiting times on G


void mtree_directed(graph& G, node& root, replaceleda::map<edge,double>& cond_prob, double min, double max);
void mtree_random(int L, array<replaceleda::string>& profile, graph& G, replaceleda::map<node,replaceleda::string>& event, replaceleda::map<int,node>& node_no, replaceleda::map<edge,double>& cond_prob, replaceleda::vector& s, int random, double min, double max);

double infinity_norm(int L, integer_matrix m);
double mtree_distance(int L, graph& G1, replaceleda::map<int,node>& node_no1, graph& G2, replaceleda::map<int,node>& node_no2);

// utilities:

#define NaN 0.0/0.0;

edge edge_between(node& v, node& w);

array<int> permutation(int n);

replaceleda::list<edge> STAR(node& root);

double myrand();
int discrand(replaceleda::vector& P);

double expcdf(double lambda);  // exponential CDF

int pow2(int n);
double log2(double x);

replaceleda::vector ones(int N);

double nonnegmean(replaceleda::vector& X);  // mean of non-negative entries of X
double nonnegmean(integer_vector& X);
double nonnegmean(replaceleda::list<double>& L);

#endif /* _MTREE_H */
