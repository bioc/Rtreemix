/*
  Cfunctions.h

   Version : 1.1.04
   Author  : Niko Beerenwinkel

  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#ifndef _CFUNCTIONS_H
#define _CFUNCTIONS_H

#include "mtreemix.h"

integer_matrix resample(integer_matrix& pattern, matrix& resp, int S, matrix& bresp);
// resample S rows from the pattern

vector CI(list<double>& L, double alpha);
// calculate (alpha*100)% confidence interval from list of doubles L

array< map<edge,double> > mtreemix_bootstrap(array<string>& profile, integer_matrix& pattern, int K, array<graph>& G, array< map<int,node> >& node_no, matrix& resp, int B, double eps, int special_weighing, int uniform_noise, array< map<edge,double> >& lower, array< map<edge,double> >& upper, vector& alpha_lower, vector& alpha_upper, double confidence_level);
// Bootstrap analysis

void mtreemix_wait(int L, vector& alpha, array<graph>& G, array< map<edge,double> >& lambda, array< map<int,node> >& node_no, array< map<edge,double> >& cond_prob, int n, int sampling_mode, double sampling_param, integer_matrix& pattern, vector& wtime, vector& stime);
// simulate patterns and their waiting times

array< list<double> > mtreemix_time(int L, graph& G_k, map<edge,double>& lambda_k, map<int,node>& node_no_k, map<node,int>& no_of_node_k, map<edge,double>& cond_prob_k, int n);
// estimate pattern waiting times for a given tree component

array< map<edge,double> > rescale_cond_prob(array< graph >& G, array< map<edge,double> >& cond_prob, int sampling_mode, double sampling_param, double output_param);

#endif /* _CFUNCTIONS_H */
