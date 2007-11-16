/*
  max_weight_branch.h

  Version : 1.1
  Author  : Niko Beerenwinkel
  
  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#ifndef _MAX_WEIGHT_BRANCH_H
#define _MAX_WEIGHT_BRANCH_H

/*
  #include <LEDA/graph.h>
  #include <LEDA/graph_gen.h>
  #include <LEDA/graph_misc.h>
  #include <LEDA/node_array.h>
  #include <LEDA/edge_array.h>
  #include <LEDA/node_set.h>
  #include <LEDA/edge_set.h>
  #include <LEDA/graphwin.h>
*/

#include "replaceleda.hh"

#include <stdlib.h>
#include <math.h>

using namespace replaceleda;

replaceleda::list<edge> MAX_WEIGHT_BRANCHING(graph& G0, replaceleda::map<node,replaceleda::string>& node_label, edge_array<double>& weight);
// Compute a maximum weight branching of the weighted digraph (G, weight).
// The list of edges of the branching is returned.

double BRANCHING_WEIGHT(replaceleda::list<edge>& B, edge_array<double>& weight);
// total weight of branching B

void UNCOVER_BRANCHING(graph& G, replaceleda::list<edge>& B);
// delete all non-branching edges (those not in B) from G

void DOT(graph& G, replaceleda::map<node,replaceleda::string>& node_label, replaceleda::map<edge,double>& cond_prob, char *filename);
// output G in dot format


#endif /* _MAX_WEIGHT_BRANCH_H */
