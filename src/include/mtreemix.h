/*
  mtreemix.h

  Version : 1.1
  Author  : Niko Beerenwinkel

  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#ifndef _MTREEMIX_H
#define _MTREEMIX_H

#include "mtree.h"
#include "kmeans.h"


// K-branchings mixture models:

// EM stopping criteria:
#define EM_ACCURACY 1e-5  // difference in log-likelihood
#define EM_MAX_ITER 1000  // max. number of iterations

#define MISSING_DATA_MAX 10  // max number of missing values per pattern to explore all possibilities. If this number is exceeded, a heuristic is applied

#define MD_ML_RUNS 100 // number of starting solutions for missing data ML heuristic

enum { PROBABILITY,
       WAITING_TIME, 
       WAITING_TIME_WEEKS, 
       BOOTSTRAP_COUNT };  // edge weight types

enum { CONSTANT,
       EXPONENTIAL };  // sampling modes


array<replaceleda::string> mtreemix_load(replaceleda::vector& alpha, array< graph >& G, array< replaceleda::map<node,replaceleda::string> >& event, array< replaceleda::map<edge,double> >& cond_prob, array< replaceleda::map<int,node> >& node_no, char *filestem);
// load model

void mtreemix_save(replaceleda::vector& alpha, array< graph >& G, array< replaceleda::map<edge,double> >& cond_prob, array< replaceleda::map<int,node> >& node_no, char *filestem);
// save model

void mtreemix_DOT(array< graph >& G, array< replaceleda::map<node,replaceleda::string> >& mut, array< replaceleda::map<edge,double> >& edge_weight, array< replaceleda::map<int,node> >& node_no, replaceleda::vector& alpha, char *filestem, int edgeweighttype);
// print model to DOT file, specific double edge labels

void mtreemix_DOT(array< graph >& G, array< replaceleda::map<node,replaceleda::string> >& mut, array< replaceleda::map<edge,replaceleda::string> >& edge_label, array< replaceleda::map<int,node> >& node_no, array<replaceleda::string>& tree_label, int uniform_noise, char *filestem);
// print model to DOT file, general edge and tree labels


double mtreemix_prob(integer_vector& pat, int K, replaceleda::vector& alpha, array<graph>& G, array< replaceleda::map<int,node> >& node_no, array< replaceleda::map<edge,double> >& cond_prob);
// Calculate the probability (likelihood) of a pattern

double mtreemix_loglike(integer_matrix& pattern, int K, replaceleda::vector& alpha, array<graph>& G, array< replaceleda::map<int,node> >& node_no, array< replaceleda::map<edge,double> >& cond_prob);
// Calculate the log-likelihood function  log L(data|model)


void mtreemix_fit0(array<replaceleda::string>& profile, integer_matrix& pattern, replaceleda::vector& alpha, array<graph>& G, array< replaceleda::map<int,node> >& node_no, array< replaceleda::map<node,replaceleda::string> >& event, array< replaceleda::map<edge,double> >& cond_prob, replaceleda::vector& resp, int uniform_noise, int special_weighing);
// Fit independence/noise model (star topology)

void mtreemix_fit1(array<replaceleda::string>& profile, integer_matrix& pattern, replaceleda::vector& alpha_opt, array<graph>& G_opt, array< replaceleda::map<int,node> >& node_no_opt, array< replaceleda::map<node,replaceleda::string> >& event_opt, array< replaceleda::map<edge,double> >& cond_prob_opt, replaceleda::vector& resp, double eps, int special_weighing);
// Fit single tree model
  
double mtreemix_fit(array<replaceleda::string>& profile, integer_matrix& pattern, int K, int M, replaceleda::vector& alpha, array<graph>& G, array< replaceleda::map<int,node> >& node_no, array< replaceleda::map<node,replaceleda::string> >& event, array< replaceleda::map<edge,double> >& cond_prob, integer_matrix& pat_hat, matrix& resp, int uniform_noise, double eps, int special_weighing);
// Fit a K-branchings mixture model to the given data


integer_matrix mtreemix_draw(int L, replaceleda::vector& alpha, array<graph>& G, array< replaceleda::map<edge,double> >& cond_prob, array< replaceleda::map<int,node> >& node_no, int n, int include_noise);
// draw sample from model


array< replaceleda::map<edge,double> > waiting_times(array< replaceleda::map<edge,double> >& cond_prob, int sampling_mode, double sampling_param);
// estimate edge waiting times


replaceleda::vector mtreemix_distr(int L, replaceleda::vector& alpha, array<graph>& G, array< replaceleda::map<int,node> >& node_no, array< replaceleda::map<edge,double> >& cond_prob);
// compute entire distribution induced by the model


int has_missing_values(integer_matrix& pat);

void alpha_random(replaceleda::vector& alpha, int& K);

void mtreemix_random(int K, int L, array<replaceleda::string>& profile, replaceleda::vector& alpha, array <graph>& G, array< replaceleda::map<node,replaceleda::string> >& event, array< replaceleda::map<int,node> >& node_no, array< replaceleda::map<edge,double> >& cond_prob, int& star, int& uniform, double& min, double& max);

double mtree_state(replaceleda::map<edge, int>& state, integer_vector& pattern, graph& G, replaceleda::map<int,node>& node_no, replaceleda::map<edge,double>& prob_cond);

integer_vector myindex2pattern(int& num_nonzeros, int index, int L);

double power(double m, int n);

matrix mtreemix_distance(int L, int K1, array< graph >& G1, array< replaceleda::map<int,node> >& node_no1, int K2, array< graph >& G2, array< replaceleda::map<int,node> >& node_no2);

#endif /* _MTREEMIX_H */
