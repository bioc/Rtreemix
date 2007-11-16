/*
  kmeans.h

  Version : 1.1
  Author  : Niko Beerenwinkel

  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#ifndef _KMEANS_H
#define _KMEANS_H

/*
  #include <LEDA/array.h>
  #include <LEDA/vector.h>
  #include <LEDA/matrix.h>
  #include <LEDA/integer_vector.h>
*/
#include "replaceleda.hh"

using namespace replaceleda;


double kmeans(int K, int S, matrix& X, integer_vector& C_opt, array<replaceleda::vector>& M_opt);
/*
  k-means clustering of data in X (columns = features, rows = samples)
  K = number of clusters
  S = number of random starting solutions
  C_opt = class indicator
  M_opt = means
*/


#endif /* _KMEANS_H */
