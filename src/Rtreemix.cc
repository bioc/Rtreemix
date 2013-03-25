// C includes

#include "include/mtreemix.h"

#include "include/Node.hh"

#include "include/cfunctions.h"



#include <iostream>

#include <vector>

#include <sstream>

#include <time.h>



// R includes 

#include <R.h>

#include <Rinternals.h>

#include <Rdefines.h>

#include <Rmath.h>

#include <R_ext/Rdynload.h>





extern "C" {



// A function for retrieving the (i, j)-th element of an integer matrix m   

#define R_INT_MATRIX(m,i,j) (INTEGER(m)[INTEGER(getAttrib(m, R_DimSymbol))[0] * j + i ]);



// List of all functions with their signatures:  

array<string> C_get_profile(SEXP R_events);

integer_matrix C_get_pattern(SEXP R_mat);



SEXP R_int_matrix(integer_matrix C_mat);

SEXP R_all_patterns(SEXP R_no_events);

SEXP R_real_matrix(matrix C_mat);

SEXP R_real_vector(vector v);



SEXP R_scalarString(const char *v);

  

  

SEXP R_fit(SEXP R_profile, SEXP R_pattern, SEXP R_K, SEXP R_M, 

	   SEXP R_uniform_noise, SEXP R_eps, SEXP R_weighing, SEXP R_seed);



SEXP R_fit1(SEXP R_profile, SEXP R_pattern, SEXP R_eps, SEXP R_weighing);



SEXP R_fit0(SEXP R_profile, SEXP R_pattern, SEXP R_uniform_noise, 

	    SEXP R_weighing);



SEXP R_bootstrap(SEXP R_profile, SEXP R_pattern, SEXP R_K, SEXP R_M, 

	   SEXP R_uniform_noise, SEXP R_eps, SEXP R_weighing, SEXP R_seed, SEXP R_B, SEXP R_conf);



int get_index(SEXP listNames, const char *str);



void R_get_graph(SEXP R_alpha, SEXP listG, vector& alpha, array<graph>& G, array< map<node,string> >& event, 

		 array< map<edge,double> >& cond_prob, array< map<int,node> >&  node_no);



SEXP R_draw(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_n, SEXP R_seed);



SEXP R_simulate(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_sampling_mode, 

		SEXP R_sampling_param, SEXP R_n, SEXP R_seed);



SEXP R_time(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_sampling_mode, 

	    SEXP R_sampling_param, SEXP R_n, SEXP R_seed);



SEXP R_likelihood(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_pattern);



SEXP R_random(SEXP R_K, SEXP R_L, SEXP R_star, SEXP R_uniform, SEXP R_min, SEXP R_max, SEXP R_seed);



SEXP R_distr(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_sampling_mode, 

	     SEXP R_sampling_param, SEXP R_out_param);







// Function that converts an R character vector to an array of strings

array<string> C_get_profile(SEXP R_events) {

  // Get the length of the R vector

  int L = length(R_events);

  

  array<string> C_profile(L);

  

  // Protect and coerce the R object 

  PROTECT(R_events = coerceVector(R_events, STRSXP));

  

  // Fill the C array with the corresponding elements of the R vector

  for(int i = 0; i < L; i ++)  

    C_profile[i] = string (CHAR(STRING_ELT(R_events, i)));

  

  

  UNPROTECT(1);

  return(C_profile);

}

  

  

  

// Function for converting R integer matrix to the C structure integer_matrix

// It uses the function R_INT_MATRIX 

integer_matrix C_get_pattern(SEXP R_mat) {

 

  // Get the dimensions of the R matrix

  SEXP Rdim = getAttrib(R_mat, R_DimSymbol);

  int nrows = INTEGER(Rdim)[0]; 

  int ncols = INTEGER(Rdim)[1];

  

  // Coerce it to a vector

  PROTECT(R_mat = coerceVector(R_mat,INTSXP));

 

  // Define a C matrix

  integer_matrix C_mat(nrows, ncols);



  // Fill it with the R matrix elements

  for(int i = 0; i < nrows; i ++) {

    for (int j = 0; j < ncols; j++) {

      C_mat[i][j] = R_INT_MATRIX(R_mat, i, j);

    }

  } 

  

  UNPROTECT(1);

  

  return(C_mat);

}



// Function for converting from the C structure integer_matrix to an R matrix of integers

SEXP R_int_matrix(integer_matrix C_mat) {   

  SEXP R_mat;

  

  // Coerce the SEXP to an R matrix of integers 

  PROTECT(R_mat = allocMatrix(INTSXP, C_mat.dim1(), C_mat.dim2()));

  // Fill the R integer matrix with the elements of the C integer_matrix.

  for(int i = 0; i < C_mat.dim1(); i ++) 

    for (int j = 0; j < C_mat.dim2(); j++) 

      INTEGER(R_mat)[i + j*C_mat.dim1()]= C_mat[i][j];

  

  

  UNPROTECT(1);



  return(R_mat);

}



// Function for creating an R integer matrix containing all possible patterns 

//for a specified number of events

SEXP R_all_patterns(SEXP R_no_events) {

  int L = INTEGER_VALUE(R_no_events);

  integer_matrix pat(pow2(L-1), L);

  // Generate all patterns:

  for (int i=0; i < pow2(L-1); i++)

    pat[i] = index2pattern(i, L);



  return(R_int_matrix(pat));

}



// Function for converting from the C structure matrix to an R real matrix 

SEXP R_real_matrix(matrix C_mat) {   

  SEXP R_mat;

  // Coerce the SEXP to an R matrix of reals

  PROTECT(R_mat = allocMatrix(REALSXP, C_mat.dim1(), C_mat.dim2()));

  // Fill the R real matrix with the elements of the C real matrix

  for(int i = 0; i < C_mat.dim1(); i ++) 

    for (int j = 0; j < C_mat.dim2(); j++) 

      REAL(R_mat)[i + j*C_mat.dim1()]= C_mat[i][j];

  

  

  UNPROTECT(1);



  return(R_mat);

}



// Function for converting C vectors to R vectors 

SEXP R_real_vector(vector v) {

  SEXP R_vec;

  // Coerce the SEXP to a vector

  PROTECT(R_vec = allocVector(REALSXP, v.size()));

  // Fill the R vector with the elements of the C vector

  for(unsigned int i = 0; i < v.size(); i++)

    REAL(R_vec)[i] = v[i];



  UNPROTECT(1);

  return(R_vec);

}



// Helper function for converting C (const char *) to an R character vector

SEXP R_scalarString(const char *v) {

  SEXP ans;

  // Coerce the SEXP in an R character vector

  PROTECT(ans = allocVector(STRSXP, 1));

  // Fill the elements of ans with the corresponding elements from the C (const char *)

  if(v)

    SET_STRING_ELT(ans, 0, mkChar(v));

  UNPROTECT(1);

  return(ans);

}







// Fitting an Rtreemix model with R_K tree components to the given set of patterns

SEXP R_fit(SEXP R_profile, SEXP R_pattern, SEXP R_K, SEXP R_M, 

	   SEXP R_uniform_noise, SEXP R_eps, SEXP R_weighing, SEXP R_seed) {



  // Load the necessary data from the arguments in their corresponding C structures



  int K = INTEGER_VALUE(R_K); // number of tree components

  int M = INTEGER_VALUE(R_M); // number of starting solutions for the k-means clustering

  int uniform_noise = INTEGER_VALUE(R_uniform_noise); // equal (=1) or unequal (=0) edge weights in noise component

  int weighing = INTEGER_VALUE(R_weighing); // use special weights log(Pr(v)) for edges (root, v) (=1) or not (=0)

  double eps = NUMERIC_VALUE(R_eps); // minimum conditional probability eps to include edge



  // Setting the seed for srand

  if(INTEGER_VALUE(R_seed) == -1) 

    srand((unsigned)time(NULL));

  else 

    srand((unsigned)INTEGER_VALUE(R_seed));



  

  integer_matrix pattern = C_get_pattern(R_pattern);

  array<string> profile = C_get_profile(R_profile);



  // Fit mtreemix model:

  vector                     alpha(K);      // mixture parameter

  array< graph >             G(K);          // digraph

  array< map<int,node> >     node_no(K);    // node of mutation index

  array< map<node,string> >  event(K);      // substitution event

  array< map<edge,double> >  cond_prob(K);  // conditional probabilities

  integer_matrix             pat_hat(pattern.dim1(), pattern.dim2());  // estimated full data

  matrix                     resp(K, pattern.dim1());     // responsibilities   



  // Call to the C function for fitting an mtreemix model

  mtreemix_fit(profile, pattern, K, M, alpha, G, node_no, event, cond_prob, pat_hat, 

	       resp, uniform_noise, eps, weighing);

 

  int j, curEl, num_nodes, num_edges;

  int curEd;

  

  // SEXP variables needed for converting all the structures connected with the mdoel from C to R

  SEXP list_graphs, one_graph, klass;

  SEXP newNames, graphN, graphE;

  SEXP curRval, curEdges, curWeights;



  SEXP attrKlass, edObj, edData, edNames, wList, edWeights; 



  SEXP result, list_names, edg_mode;

  

  // The list that will contain the function result

  PROTECT(result = allocVector(VECSXP, 4));

  // Names of the lists that are part of the output 

  PROTECT(list_names = allocVector(STRSXP, 4));

  SET_STRING_ELT(list_names, 0, mkChar("alpha")); // the weight vector of the model

  SET_STRING_ELT(list_names, 1, mkChar("resp")); // the responsibilities

  SET_STRING_ELT(list_names, 2, mkChar("pat.hat")); // the complete sample matrix (if there were some missing data)

  SET_STRING_ELT(list_names, 3, mkChar("graphs.mixture")); // the list of the graphs each for every tree component

  

  setAttrib(result, R_NamesSymbol, list_names);



  SET_VECTOR_ELT(result, 0, R_real_vector(alpha));



  SET_VECTOR_ELT(result, 1, R_real_matrix(resp));



  if (has_missing_values(pattern)) 

    SET_VECTOR_ELT(result, 2, R_int_matrix(pat_hat));

  else

    SET_VECTOR_ELT(result, 2, allocMatrix(REALSXP, 0, 0)); 	

    

  // Build the list of graphs that give the trees composing the mixture model:

  PROTECT(list_graphs = allocVector(VECSXP, K));



  klass = MAKE_CLASS("graphNEL");



  // Needed for the edgeData construction

  attrKlass = MAKE_CLASS("attrData");



  PROTECT(newNames = allocVector(STRSXP, 2));

  SET_STRING_ELT(newNames, 0, mkChar("edges"));

  SET_STRING_ELT(newNames, 1, mkChar("weights"));



  // Build each of the K graphs

  for (int k = 0; k < K; k++) {    

    PROTECT(one_graph = NEW_OBJECT(klass));

    // R version >= 2.7.0 - edgemode slot deprecated
	  PROTECT(edg_mode = allocVector(VECSXP, 1));  
	  
	  setAttrib(edg_mode, R_NamesSymbol, R_scalarString("edgemode"));
	  
    SET_VECTOR_ELT(edg_mode, 0, R_scalarString("directed"));  

    SET_SLOT(one_graph, Rf_install("graphData"), edg_mode);   

    num_nodes = G[k].number_of_nodes();

    num_edges = G[k].number_of_edges();

  

    if (num_nodes == 0) {

      SET_SLOT(one_graph, Rf_install("nodes"), allocVector(STRSXP, 0));

      SET_SLOT(one_graph, Rf_install("edgeL"), allocVector(VECSXP, 0)); 

    } else {

      // Make the node and edges lists of the graph object

      PROTECT(graphN = allocVector(STRSXP, num_nodes));

      PROTECT(graphE = allocVector(VECSXP, num_nodes));



      // Structures for the edgeData construction

      PROTECT(edObj = NEW_OBJECT(attrKlass));

      PROTECT(edData = allocVector(VECSXP, num_edges)); 

      PROTECT(edNames = allocVector(STRSXP, num_edges));

      curEd = 0;



      node n;

      edge e;

      j = 0;

      

      // Build the node list in R

      forall_nodes (n, G[k]) {    

	SET_STRING_ELT(graphN, j, STRING_ELT(R_scalarString(event[k][n].c_str()), 0));

  

	PROTECT(curRval = allocVector(VECSXP, 2));

	setAttrib(curRval, R_NamesSymbol, newNames);

	if(G[k].outdeg(n) == 0) {

	  SET_VECTOR_ELT(curRval, 0, allocVector(INTSXP, 0));

	  SET_VECTOR_ELT(curRval, 1, allocVector(REALSXP, 0));

	} else {  

	  curEl = 0;

	  // The edge list for the node n

	  PROTECT(curEdges = allocVector(INTSXP, G[k].outdeg(n)));

	 

	  // The corresponding weights list for the node n

	  PROTECT(curWeights = allocVector(REALSXP, G[k].outdeg(n)));

	  

          // Construct the edgeData:

	  forall_out_edges(e, n) {

	   SET_STRING_ELT(edNames, curEd, 

			  STRING_ELT(R_scalarString((event[k][n] + "|" 

						     + event[k][G[k].target(e)]).c_str()), 0));   



	   PROTECT(wList = allocVector(VECSXP, 1));

	   PROTECT(edWeights = allocVector(REALSXP, 1));

	   setAttrib(wList, R_NamesSymbol, R_scalarString("weight"));

	   REAL (edWeights)[0] = cond_prob[k][e];

 

	   SET_VECTOR_ELT(wList, 0, edWeights);



	   SET_VECTOR_ELT(edData, curEd, wList);



	   //SET_STRING_ELT(curEdges, curEl, STRING_ELT(R_scalarString(event[k][G[k].target(e)].c_str()), 0));

	   INTEGER (curEdges)[curEl] = (int) (G[k].target(e)->getIndex()) + 1;

	   REAL (curWeights)[curEl] = cond_prob[k][e];

	  	   

	   UNPROTECT(2); //edWeights, wList

	   curEd++;

           curEl ++; 

	  }

	 

	  SET_VECTOR_ELT(curRval, 0, curEdges);

	  SET_VECTOR_ELT(curRval, 1, curWeights);

	  UNPROTECT(2); // curEdges, curWeights

	}

	SET_VECTOR_ELT(graphE, j, curRval);

	UNPROTECT(1); //curRval

	

	j ++;

      }

      setAttrib(graphE, R_NamesSymbol, graphN); 

      setAttrib(edData, R_NamesSymbol, edNames);

      // Install the necessary slots an R graph needs:

      SET_SLOT(edObj, Rf_install("default"), allocVector(VECSXP, 0));

      SET_SLOT(edObj, Rf_install("data"), edData);



      SET_SLOT(one_graph, Rf_install("edgeData"), edObj);

      SET_SLOT(one_graph, Rf_install("edgeL"), graphE);

      SET_SLOT(one_graph, Rf_install("nodes"), graphN);



      

      UNPROTECT(5); //graphN, graphE, edNames, edData, edObj

    }

    

    SET_VECTOR_ELT(list_graphs, k, one_graph);

    

    UNPROTECT(2); // one_graph, edg_mode

  }    



  SET_VECTOR_ELT(result, 3, list_graphs);

  

  UNPROTECT(2); // newNames, list_graphs

  

  UNPROTECT(2); //result, list_names

 

  return(result); //put result

  

}

  

// Fitting a single tree model to the given set of patterns 

SEXP R_fit1(SEXP R_profile, SEXP R_pattern, SEXP R_eps, SEXP R_weighing) {



  // Load the necessary data in their corresponding C structures

  int K = 1; // number of tree components

  int weighing = INTEGER_VALUE(R_weighing); // use special weights log(Pr(v)) for edges (root, v) (=1) or not (=0)

  double eps = NUMERIC_VALUE(R_eps); // minimum conditional probability eps to include edge



  integer_matrix pattern = C_get_pattern(R_pattern);

  array<string> profile = C_get_profile(R_profile);



  // Fit a mixture model:

  vector                     alpha(1);      // mixture parameter

  array< graph >             G(1);          // digraph

  array< map<int,node> >     node_no(1);    // node of mutation index

  array< map<node,string> >  event(1);      // substitution event

  array< map<edge,double> >  cond_prob(1);  // conditional probabilities

  

  vector resp(pattern.dim1());

  for (int i = 0; i < pattern.dim1(); i ++)

    resp[i] = 1.0;

  

  mtreemix_fit1(profile, pattern, alpha, G, node_no, event, cond_prob, 

	      resp, eps, weighing);



  int j, curEl, num_nodes, num_edges;

  int curEd;



  // SEXP variables needed for converting all the structures connected with the mdoel from C to R

  SEXP list_graphs, one_graph, klass;

  SEXP newNames, graphN, graphE;

  SEXP curRval, curEdges, curWeights;



  SEXP attrKlass, edObj, edData, edNames, wList, edWeights; 



  SEXP result, list_names, edg_mode;

  

  //The list that will contain the function result

  PROTECT(result = allocVector(VECSXP, 3));

  // Names of the lists that are part of the output

  PROTECT(list_names = allocVector(STRSXP, 3));

  SET_STRING_ELT(list_names, 0, mkChar("alpha")); // the weight of the tree

  SET_STRING_ELT(list_names, 1, mkChar("resp")); // the responsibility of the tree

  SET_STRING_ELT(list_names, 2, mkChar("graphs.mixture")); // the list of the graphs each for every tree component

  

  setAttrib(result, R_NamesSymbol, list_names);



  SET_VECTOR_ELT(result, 0, R_real_vector(alpha));



  SET_VECTOR_ELT(result, 1, R_real_vector(resp));

  

  // Build the list that contains the graph corresponding to the single tree component:  

  PROTECT(list_graphs = allocVector(VECSXP, K));

  

  klass = MAKE_CLASS("graphNEL");



  // Needed for the edgeData construction

  attrKlass = MAKE_CLASS("attrData");



  PROTECT(newNames = allocVector(STRSXP, 2));

  SET_STRING_ELT(newNames, 0, mkChar("edges"));

  SET_STRING_ELT(newNames, 1, mkChar("weights"));

  

  // Build the graph

  for (int k = 0; k < K; k++) {    

    PROTECT(one_graph = NEW_OBJECT(klass));

    // R version >= 2.7.0 - edgemode slot deprecated
	  PROTECT(edg_mode = allocVector(VECSXP, 1));  
	  
	  setAttrib(edg_mode, R_NamesSymbol, R_scalarString("edgemode"));
	  
    SET_VECTOR_ELT(edg_mode, 0, R_scalarString("directed"));  

    SET_SLOT(one_graph, Rf_install("graphData"), edg_mode);   
    

    num_nodes = G[k].number_of_nodes();

    num_edges = G[k].number_of_edges();

  

    if (num_nodes == 0) {

      SET_SLOT(one_graph, Rf_install("nodes"), allocVector(STRSXP, 0));

      SET_SLOT(one_graph, Rf_install("edgeL"), allocVector(VECSXP, 0)); 

    } else {

      // Make the node and edges lists of the graph object

      PROTECT(graphN = allocVector(STRSXP, num_nodes));

      PROTECT(graphE = allocVector(VECSXP, num_nodes));



      // Structures for the edgeData construction

      PROTECT(edObj = NEW_OBJECT(attrKlass));

      PROTECT(edData = allocVector(VECSXP, num_edges)); 

      PROTECT(edNames = allocVector(STRSXP, num_edges));

      curEd = 0;



      node n;

      edge e;

      j = 0;



      // Build the node list in R

      forall_nodes (n, G[k]) {    

	SET_STRING_ELT(graphN, j, STRING_ELT(R_scalarString(event[k][n].c_str()), 0));

  

	PROTECT(curRval = allocVector(VECSXP, 2));

	setAttrib(curRval, R_NamesSymbol, newNames);

	if(G[k].outdeg(n) == 0) {

	  SET_VECTOR_ELT(curRval, 0, allocVector(INTSXP, 0));

	  SET_VECTOR_ELT(curRval, 1, allocVector(REALSXP, 0));

	} else {  

	  curEl = 0;

	  // The edge list for the node n

	  PROTECT(curEdges = allocVector(INTSXP, G[k].outdeg(n)));

	 

	  // The corresponding weights list for the node n

	  PROTECT(curWeights = allocVector(REALSXP, G[k].outdeg(n)));

	  // Construct the edgeData:

	  forall_out_edges(e, n) {

	   SET_STRING_ELT(edNames, curEd, 

			  STRING_ELT(R_scalarString((event[k][n] + "|" 

						     + event[k][G[k].target(e)]).c_str()), 0));   



	   PROTECT(wList = allocVector(VECSXP, 1));

	   PROTECT(edWeights = allocVector(REALSXP, 1));

	   setAttrib(wList, R_NamesSymbol, R_scalarString("weight"));

	   REAL (edWeights)[0] = cond_prob[k][e];

	   SET_VECTOR_ELT(wList, 0, edWeights);



	   SET_VECTOR_ELT(edData, curEd, wList);





	   INTEGER (curEdges)[curEl] = (int) (G[k].target(e)->getIndex()) + 1;

	   REAL (curWeights)[curEl] = cond_prob[k][e];

	  	   

	   UNPROTECT(2); //edWeights, wList

	   curEd++;

           curEl ++; 

	  }

	 

	  SET_VECTOR_ELT(curRval, 0, curEdges);

	  SET_VECTOR_ELT(curRval, 1, curWeights);

	  UNPROTECT(2); // curEdges, curWeights

	}

	SET_VECTOR_ELT(graphE, j, curRval);

	UNPROTECT(1); //curRval

	

	j ++;

      }

      setAttrib(graphE, R_NamesSymbol, graphN); 

      setAttrib(edData, R_NamesSymbol, edNames);

      // Install the necessary slots an R graph needs:

      SET_SLOT(edObj, Rf_install("default"), allocVector(VECSXP, 0));

      SET_SLOT(edObj, Rf_install("data"), edData);



      SET_SLOT(one_graph, Rf_install("edgeData"), edObj);

      SET_SLOT(one_graph, Rf_install("edgeL"), graphE);

      SET_SLOT(one_graph, Rf_install("nodes"), graphN);



      

      UNPROTECT(5); //graphN, graphE, edNames, edData, edObj

    }

    

    SET_VECTOR_ELT(list_graphs, k, one_graph);

    

    UNPROTECT(2); // one_graph, edg_mode

  }    



  SET_VECTOR_ELT(result, 2, list_graphs);

  

  UNPROTECT(2); // newNames, list_graphs

  

  UNPROTECT(2); //result, list_names



  return(result);

  

}



// Fitting a single star model to the data (independence of events is assumed)

SEXP R_fit0(SEXP R_profile, SEXP R_pattern, SEXP R_uniform_noise, 

	    SEXP R_weighing) {



  // Load the necessary data in their corresponding C structures

  int K = 1; // number of tree components

  int uniform_noise = INTEGER_VALUE(R_uniform_noise); // equal (=1) or unequal (=0) edge weights in noise component

  int weighing = INTEGER_VALUE(R_weighing); // use special weights log(Pr(v)) for edges (root, v) (=1) or not (=0)



  integer_matrix pattern = C_get_pattern(R_pattern);

  array<string> profile = C_get_profile(R_profile);





  // Fit  a mixture model:

  vector                     alpha(1);      // mixture parameter

  array< graph >             G(1);          // digraph

  array< map<int,node> >     node_no(1);    // node of mutation index

  array< map<node,string> >  event(1);      // substitution event

  array< map<edge,double> >  cond_prob(1);  // conditional probabilities

  

  vector resp(pattern.dim1());

  for (int i = 0; i < pattern.dim1(); i ++)

    resp[i] = 1.0;

  

  mtreemix_fit0(profile, pattern, alpha, G, node_no, event, cond_prob, 

	      resp, uniform_noise, weighing);



  int j, curEl, num_nodes, num_edges;

  int curEd;



  // SEXP variables needed for converting all the structures connected with the mdoel from C to R

  SEXP list_graphs, one_graph, klass;

  SEXP newNames, graphN, graphE;

  SEXP curRval, curEdges, curWeights;



  SEXP attrKlass, edObj, edData, edNames, wList, edWeights; 



  SEXP result, list_names, edg_mode;

  

  //The list that will contain the returned data

  PROTECT(result = allocVector(VECSXP, 3));

  // Names of the lists that are part of the output

  PROTECT(list_names = allocVector(STRSXP, 3));

  SET_STRING_ELT(list_names, 0, mkChar("alpha")); // the weight of the star tree

  SET_STRING_ELT(list_names, 1, mkChar("resp")); // the responsibility of the star tree

  SET_STRING_ELT(list_names, 2, mkChar("graphs.mixture")); // the graph corresponding to the single star tree component

  

  setAttrib(result, R_NamesSymbol, list_names);



  SET_VECTOR_ELT(result, 0, R_real_vector(alpha));



  SET_VECTOR_ELT(result, 1, R_real_vector(resp));



  // Build the list that contains the graph corresponding to the single star tree component:    

  PROTECT(list_graphs = allocVector(VECSXP, K));



  klass = MAKE_CLASS("graphNEL");



  // Needed for the edgeData construction

  attrKlass = MAKE_CLASS("attrData");



  PROTECT(newNames = allocVector(STRSXP, 2));

  SET_STRING_ELT(newNames, 0, mkChar("edges"));

  SET_STRING_ELT(newNames, 1, mkChar("weights"));



  // Build the graph

  for (int k = 0; k < K; k++) {    

    PROTECT(one_graph = NEW_OBJECT(klass));

    // R version >= 2.7.0 - edgemode slot deprecated
	  PROTECT(edg_mode = allocVector(VECSXP, 1));  
	  
	  setAttrib(edg_mode, R_NamesSymbol, R_scalarString("edgemode"));
	  
    SET_VECTOR_ELT(edg_mode, 0, R_scalarString("directed"));  

    SET_SLOT(one_graph, Rf_install("graphData"), edg_mode);   

    

    num_nodes = G[k].number_of_nodes();

    num_edges = G[k].number_of_edges();

  

    if (num_nodes == 0) {

      SET_SLOT(one_graph, Rf_install("nodes"), allocVector(STRSXP, 0));

      SET_SLOT(one_graph, Rf_install("edgeL"), allocVector(VECSXP, 0)); 

    } else {

      // Make the node and edges lists of the graph object.

      PROTECT(graphN = allocVector(STRSXP, num_nodes));

      PROTECT(graphE = allocVector(VECSXP, num_nodes));



      // Structures for the edgeData construction.

      PROTECT(edObj = NEW_OBJECT(attrKlass));

      PROTECT(edData = allocVector(VECSXP, num_edges)); 

      PROTECT(edNames = allocVector(STRSXP, num_edges));

      curEd = 0;



      node n;

      edge e;

      j = 0;



      // Build the node list in R

      forall_nodes (n, G[k]) {    

	SET_STRING_ELT(graphN, j, STRING_ELT(R_scalarString(event[k][n].c_str()), 0));

  

	PROTECT(curRval = allocVector(VECSXP, 2));

	setAttrib(curRval, R_NamesSymbol, newNames);

	if(G[k].outdeg(n) == 0) {

	  SET_VECTOR_ELT(curRval, 0, allocVector(INTSXP, 0));

	  SET_VECTOR_ELT(curRval, 1, allocVector(REALSXP, 0));

	} else {  

	  curEl = 0;

	  // The edge list for the node n

	  PROTECT(curEdges = allocVector(INTSXP, G[k].outdeg(n)));

	 

	  // The corresponding weights list for the node n

	  PROTECT(curWeights = allocVector(REALSXP, G[k].outdeg(n)));

	  // Construct the edgeData:

	  forall_out_edges(e, n) {

	   SET_STRING_ELT(edNames, curEd, 

			  STRING_ELT(R_scalarString((event[k][n] + "|" 

						     + event[k][G[k].target(e)]).c_str()), 0));   



	   PROTECT(wList = allocVector(VECSXP, 1));

	   PROTECT(edWeights = allocVector(REALSXP, 1));

	   setAttrib(wList, R_NamesSymbol, R_scalarString("weight"));

	   REAL (edWeights)[0] = cond_prob[k][e];

	   SET_VECTOR_ELT(wList, 0, edWeights);



	   SET_VECTOR_ELT(edData, curEd, wList);





	   INTEGER (curEdges)[curEl] = (int) (G[k].target(e)->getIndex()) + 1;

	   REAL (curWeights)[curEl] = cond_prob[k][e];

	  	   

	   UNPROTECT(2); //edWeights, wList

	   curEd++;

           curEl ++; 

	  }

	 

	  SET_VECTOR_ELT(curRval, 0, curEdges);

	  SET_VECTOR_ELT(curRval, 1, curWeights);

	  UNPROTECT(2); // curEdges, curWeights

	}

	SET_VECTOR_ELT(graphE, j, curRval);

	UNPROTECT(1); //curRval

	

	j ++;

      }

      setAttrib(graphE, R_NamesSymbol, graphN); 

      setAttrib(edData, R_NamesSymbol, edNames);

      // Install the necessary slots an R graph needs:

      SET_SLOT(edObj, Rf_install("default"), allocVector(VECSXP, 0));

      SET_SLOT(edObj, Rf_install("data"), edData);



      SET_SLOT(one_graph, Rf_install("edgeData"), edObj);

      SET_SLOT(one_graph, Rf_install("edgeL"), graphE);

      SET_SLOT(one_graph, Rf_install("nodes"), graphN);



      

      UNPROTECT(5); //graphN, graphE, edNames, edData, edObj

    }

    

    SET_VECTOR_ELT(list_graphs, k, one_graph);

    

    UNPROTECT(2); // one_graph, edg_mode

  }    



  SET_VECTOR_ELT(result, 2, list_graphs);

  

  UNPROTECT(2); // newNames, list_graphs

  

  UNPROTECT(2); //result, list_names



  return(result);

  

}



 



// Fitting Rtreemix model to given set of patterns and analyzing its variance with the bootstrap method

SEXP R_bootstrap(SEXP R_profile, SEXP R_pattern, SEXP R_K, SEXP R_M, 

		 SEXP R_uniform_noise, SEXP R_eps, SEXP R_weighing, SEXP R_seed, SEXP R_B, SEXP R_conf) {



  // Load the necessary data in their corresponding C structures



  int K = INTEGER_VALUE(R_K); // number of tree components

  int M = INTEGER_VALUE(R_M); // number of starting solutions for the k-means clustering

  int uniform_noise = INTEGER_VALUE(R_uniform_noise); // equal (=1) or unequal (=0) edge weights in noise component

  int weighing = INTEGER_VALUE(R_weighing); // use special weights log(Pr(v)) for edges (root, v) (=1) or not (=0)

  double eps = NUMERIC_VALUE(R_eps); // minimum conditional probability eps to include edge

  int B = INTEGER_VALUE(R_B); // number of bootstrap samples

  double confidence = NUMERIC_VALUE(R_conf); // confidence level for intervals

  

  // Setting the seed for srand

  if(INTEGER_VALUE(R_seed) == -1) 

    srand((unsigned)time(NULL));

  else 

    srand((unsigned)INTEGER_VALUE(R_seed));



  

  integer_matrix pattern = C_get_pattern(R_pattern);

  array<string> profile = C_get_profile(R_profile);



  // Fit mtreemix model:

  vector                     alpha(K);      // mixture parameter

  array< graph >             G(K);          // digraph

  array< map<int,node> >     node_no(K);    // node of mutation index

  array< map<node,string> >  event(K);      // substitution event

  array< map<edge,double> >  cond_prob(K);  // conditional probabilities

  integer_matrix             pat_hat(pattern.dim1(), pattern.dim2());  // estimated full data

  matrix                     resp(K, pattern.dim1());     // responsibilities   

  

  mtreemix_fit(profile, pattern, K, M, alpha, G, node_no, event, cond_prob, pat_hat, 

	      resp, uniform_noise, eps, weighing);





  // Bootstrap:

  array< map<edge,double> > lower(K);  // lower bound of CI for cond_prob

  array< map<edge,double> > upper(K);  // upper bound of CI for cond_prob

  vector alpha_lower(K), alpha_upper(K);  // CI for alpha

  array< map<edge,double> > supp = mtreemix_bootstrap(profile, pattern, K, G, node_no, 

						      resp, B, eps, weighing, 

						      uniform_noise, lower, upper, alpha_lower, 

						      alpha_upper, confidence);



  int j, curEl, num_nodes, num_edges;

  int curEd;

  

  // SEXP variables needed for converting all the structures connected with the mdoel from C to R

  SEXP list_graphs, one_graph, klass;

  SEXP newNames, graphN, graphE;

  SEXP curRval, curEdges, curWeights;



  SEXP attrKlass, edObj, edData, edNames, wList, edWeights, edCI, edCINames, edAttr;



  SEXP alphaCI, alphaVec;



  SEXP result, list_names, edg_mode;

  

  //The list that will contain the returned data

  PROTECT(result = allocVector(VECSXP, 5));

  // Names of the lists that are part of the output

  PROTECT(list_names = allocVector(STRSXP, 5));

  SET_STRING_ELT(list_names, 0, mkChar("alpha")); // the weight vector of the model

  SET_STRING_ELT(list_names, 1, mkChar("resp")); // the responsibilities

  SET_STRING_ELT(list_names, 2, mkChar("pat.hat")); // the complete sample matrix (if there were some missing data)

  SET_STRING_ELT(list_names, 3, mkChar("graphs.mixture")); // the list of the graphs each for every tree component

  SET_STRING_ELT(list_names, 4, mkChar("alpha.ci")); // te list of the confidence intervals for the weights

  

  setAttrib(result, R_NamesSymbol, list_names);



  SET_VECTOR_ELT(result, 0, R_real_vector(alpha));



  SET_VECTOR_ELT(result, 1, R_real_matrix(resp));



  if (has_missing_values(pattern)) 

    SET_VECTOR_ELT(result, 2, R_int_matrix(pat_hat));

  else

    SET_VECTOR_ELT(result, 2, allocMatrix(REALSXP, 0, 0)); 

    



  PROTECT(alphaCI = allocVector(VECSXP, K)); // the list holding the confidence intervals for alpha



  PROTECT(list_graphs = allocVector(VECSXP, K)); // the list for the trees



  klass = MAKE_CLASS("graphNEL");



  // Needed for the edgeData construction

  attrKlass = MAKE_CLASS("attrData");



  PROTECT(newNames = allocVector(STRSXP, 2));

  SET_STRING_ELT(newNames, 0, mkChar("edges"));

  SET_STRING_ELT(newNames, 1, mkChar("weights"));



  PROTECT(edAttr = allocVector(STRSXP, 2));

  SET_STRING_ELT(edAttr, 0, mkChar("weight"));

  SET_STRING_ELT(edAttr, 1, mkChar("ci"));



  PROTECT(edCINames = allocVector(STRSXP, 3));

  SET_STRING_ELT(edCINames, 0, mkChar("lower"));

  SET_STRING_ELT(edCINames, 1, mkChar("upper"));

  SET_STRING_ELT(edCINames, 2, mkChar("supp"));



  // Build each of the K graphs with the confidence intervals:

  for (int k = 0; k < K; k++) {

    

    // The confidence intervals for alpha:

    PROTECT(alphaVec = allocVector(REALSXP, 2));

    REAL (alphaVec) [0] = alpha_lower[k];

    REAL (alphaVec) [1] = alpha_upper[k];



    SET_VECTOR_ELT(alphaCI, k, alphaVec);

    

    UNPROTECT(1); // alphaVec



    

    PROTECT(one_graph = NEW_OBJECT(klass)); // the k-th tree

    // R version >= 2.7.0 - edgemode slot deprecated
	  PROTECT(edg_mode = allocVector(VECSXP, 1));  
	  
	  setAttrib(edg_mode, R_NamesSymbol, R_scalarString("edgemode"));
	  
    SET_VECTOR_ELT(edg_mode, 0, R_scalarString("directed"));  

    SET_SLOT(one_graph, Rf_install("graphData"), edg_mode);   

    

    num_nodes = G[k].number_of_nodes();

    num_edges = G[k].number_of_edges();

  

    if (num_nodes == 0) {

      SET_SLOT(one_graph, Rf_install("nodes"), allocVector(STRSXP, 0));

      SET_SLOT(one_graph, Rf_install("edgeL"), allocVector(VECSXP, 0)); 

    } else {

      // Make the node and edges lists of the graph object

      PROTECT(graphN = allocVector(STRSXP, num_nodes));

      PROTECT(graphE = allocVector(VECSXP, num_nodes));



      // Structures for the edgeData construction

      PROTECT(edObj = NEW_OBJECT(attrKlass));

      PROTECT(edData = allocVector(VECSXP, num_edges)); 

      PROTECT(edNames = allocVector(STRSXP, num_edges));

      curEd = 0;



      node n;

      edge e;

      j = 0;



      // Build the node list in R

      forall_nodes (n, G[k]) {    

	SET_STRING_ELT(graphN, j, STRING_ELT(R_scalarString(event[k][n].c_str()), 0));

  

	PROTECT(curRval = allocVector(VECSXP, 2));

	setAttrib(curRval, R_NamesSymbol, newNames);

	if(G[k].outdeg(n) == 0) {

	  SET_VECTOR_ELT(curRval, 0, allocVector(INTSXP, 0));

	  SET_VECTOR_ELT(curRval, 1, allocVector(REALSXP, 0));

	} else {  

	  curEl = 0;

	  // The edge list for the node n

	  PROTECT(curEdges = allocVector(INTSXP, G[k].outdeg(n)));

	 

	  // The corresponding weights list for the node n

	  PROTECT(curWeights = allocVector(REALSXP, G[k].outdeg(n)));

	  

	  forall_out_edges(e, n) {

	   SET_STRING_ELT(edNames, curEd, 

			  STRING_ELT(R_scalarString((event[k][n] + "|" 

						     + event[k][G[k].target(e)]).c_str()), 0));   

	   // Construct the edgeData containing the weight and ci attributes

	   PROTECT(wList = allocVector(VECSXP, 2));

	   setAttrib(wList, R_NamesSymbol, edAttr);

           // Edge weights:

	   PROTECT(edWeights = allocVector(REALSXP, 1));

	   REAL (edWeights)[0] = cond_prob[k][e];

           // Confidence intervals for edges:

	   PROTECT(edCI = allocVector(REALSXP, 3));

	   setAttrib(edCI, R_NamesSymbol, edCINames);

	   REAL (edCI)[0] = lower[k][e];

	   REAL (edCI)[1] = upper[k][e];

	   REAL (edCI)[2] = supp[k][e];

 

	   SET_VECTOR_ELT(wList, 0, edWeights);

	   SET_VECTOR_ELT(wList, 1, edCI);



	   SET_VECTOR_ELT(edData, curEd, wList);





	   INTEGER (curEdges)[curEl] = (int) (G[k].target(e)->getIndex()) + 1;

	   REAL (curWeights)[curEl] = cond_prob[k][e];

	  	   

	   UNPROTECT(3); //edWeights, wList, edCI

	   curEd++;

           curEl ++; 

	  }

	 

	  SET_VECTOR_ELT(curRval, 0, curEdges);

	  SET_VECTOR_ELT(curRval, 1, curWeights);

	  UNPROTECT(2); // curEdges, curWeights

	}

	SET_VECTOR_ELT(graphE, j, curRval);

	UNPROTECT(1); //curRval

	

	j ++;

      }

      setAttrib(graphE, R_NamesSymbol, graphN); 

      setAttrib(edData, R_NamesSymbol, edNames);

      // Install the necessary slots an R graph needs:

      SET_SLOT(edObj, Rf_install("default"), allocVector(VECSXP, 0));

      SET_SLOT(edObj, Rf_install("data"), edData);



      

      SET_SLOT(one_graph, Rf_install("edgeL"), graphE);

      SET_SLOT(one_graph, Rf_install("edgeData"), edObj);

      SET_SLOT(one_graph, Rf_install("nodes"), graphN);



      

      UNPROTECT(5); //graphN, graphE, edNames, edData, edObj

    }

    

    SET_VECTOR_ELT(list_graphs, k, one_graph);

    

    UNPROTECT(2); // one_graph, edg_mode

  }    



  SET_VECTOR_ELT(result, 3, list_graphs);

  SET_VECTOR_ELT(result, 4, alphaCI);

  

  UNPROTECT(4); // newNames, list_graphs, alphaCI, edAttr

  

  UNPROTECT(3); //result, list_names, edCINames



  return(result); 

  

}



// Get the index of the string str in the list of names.

int get_index(SEXP listNames, const char *str) {

  int index = -1;



  for (int i = 0; i < length(listNames); i++) {

    if(strcmp(CHAR(STRING_ELT(listNames, i)), str) == 0) {

      index = i;

      break;

    }

  }

  return index;

}





// Get the C graph structrure from a given R graph

void R_get_graph(SEXP R_alpha, SEXP listG, vector& alpha, array<graph>& G, array< map<node,string> >& event, 

		 array< map<edge,double> >& cond_prob, array< map<int,node> >&  node_no) {

  node v;

  edge e;



  PROTECT(R_alpha = coerceVector(R_alpha, REALSXP));

  PROTECT(listG = coerceVector(listG, VECSXP));

  SEXP one_graph, gN, gEW, ew, names;

  

  G.resize(length(listG));

  event.resize(length(listG));

  cond_prob.resize(length(listG));

  node_no.resize(length(listG));

  

  for(int k = 0; k < length(listG); k++) {

    alpha[k] = REAL (R_alpha)[k]; // get the weights of the tree components



    PROTECT(one_graph = coerceVector(VECTOR_ELT(listG, k), VECSXP));

    

    PROTECT(gN = AS_CHARACTER(VECTOR_ELT(one_graph, 0)));

    

    // Build the node structure in C by using the node list from R:

    node_no[k].clear();

    event[k].clear();   

    for (int l = 0; l < length(gN); l++) {

      v = G[k].new_node();

      event[k][v] = string (CHAR(STRING_ELT(gN, l)));

      node_no[k][l] = v;      

    }

	    

    PROTECT(gEW = coerceVector(VECTOR_ELT(one_graph, 1), VECSXP));

    // Build the edge structure of the graphs in C by using the edge structure in R:

    cond_prob[k].clear();

    for (int i = 0; i < length(gN); i++) {

      PROTECT(ew = coerceVector(VECTOR_ELT(gEW, i), REALSXP));

      

      if(length(ew) != 0) {

	names = AS_CHARACTER(getAttrib(ew, R_NamesSymbol));

	// Edge weights:

	for(int j = 0; j < length(ew); j++) {

	  e = G[k].new_edge(G[k].getNode(i), G[k].getNode(get_index(gN, CHAR(STRING_ELT(names, j)))));

	  cond_prob[k][e] = REAL(ew)[j];

	}

      }

      UNPROTECT(1); // ew

    }

  

    

    UNPROTECT(2); // gN, gEW

    UNPROTECT(1); //one_graph

  }

   

  UNPROTECT(2); //listG, R_alpha

}



// Function for drawing samples from a given mixture model

SEXP R_draw(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_n, SEXP R_seed) {

  // Load the necessary data in their corresponding C structures

  int L = INTEGER_VALUE(R_L); // number of genetic events

  int n = INTEGER_VALUE(R_n); // number of samples to draw



  // Setting the seed for srand

  if(INTEGER_VALUE(R_seed) == -1) 

    srand((unsigned)time(NULL));

  else 

    srand((unsigned)INTEGER_VALUE(R_seed));



  // Load the model in R_trees in the corresponding C structures:

  vector alpha(length(R_trees));      // mixture parameter

  array<graph> G; // digraph

  array< map<node,string> > event; // genetic event

  array< map<edge,double> > cond_prob; // conditional probabilities

  array< map<int,node> > node_no; // node of mutation index



  R_get_graph(R_alpha, R_trees, alpha, G, event, cond_prob, node_no);

  

  // Simulate data:

  integer_matrix patsim = mtreemix_draw(L, alpha, G, cond_prob, node_no, n, 1);



  return(R_int_matrix(patsim));



}



// Function for simulating patterns and their waiting times from a given mixture model

SEXP R_simulate(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_sampling_mode, 

		SEXP R_sampling_param, SEXP R_n, SEXP R_seed){



  // Load the necessary data in their corresponding C structures

  int L = INTEGER_VALUE(R_L); // number of genetic events 

  int  sampling_mode = (INTEGER_VALUE(R_sampling_mode) != 0)?EXPONENTIAL:CONSTANT; // sampling mode for the simulations: exponential or constant

  double sampling_param = NUMERIC_VALUE(R_sampling_param); // sampling parameter that corresponds to the sampling mode

  int n = INTEGER_VALUE(R_n); // number of simulations



  // Setting the seed for srand

  if(INTEGER_VALUE(R_seed) == -1) 

    srand((unsigned)time(NULL));

  else 

    srand((unsigned)INTEGER_VALUE(R_seed));

  

  // Load the model in R_trees in the corresponding C structures:

  vector alpha(length(R_trees));      // mixture parameter

  array<graph> G; // digraph

  array< map<node,string> > event; // genetic event

  array< map<edge,double> > cond_prob; // conditional probabilities

  array< map<int,node> > node_no; // node of mutation index



  R_get_graph(R_alpha, R_trees, alpha, G, event, cond_prob, node_no);

  

  // Estimate exponential waiting times on edges:

  array< map<edge,double> >  lambda;  // exponential parameter lambda_i

  lambda = waiting_times(cond_prob, sampling_mode, sampling_param);

  

  // Simulate data:

  integer_matrix pattern(n, L);  // simulated patterns,

  vector wtime(n);               //    their waiting times,

  vector stime(n);               //    and their sampling times



  mtreemix_wait(L, alpha, G, lambda, node_no, cond_prob, n, sampling_mode, 

		NUMERIC_VALUE(R_sampling_param), pattern, wtime, stime);



  SEXP result, listNames;



  PROTECT(listNames = allocVector(STRSXP, 3));

  SET_STRING_ELT(listNames, 0, mkChar("patterns"));

  SET_STRING_ELT(listNames, 1, mkChar("wtimes"));

  SET_STRING_ELT(listNames, 2, mkChar("stimes"));



  PROTECT(result = allocVector(VECSXP, 3));

  setAttrib(result, R_NamesSymbol, listNames);

  UNPROTECT(1); // listNames

  

  SET_VECTOR_ELT(result, 0, R_int_matrix(pattern));

  SET_VECTOR_ELT(result, 1, R_real_vector(wtime));

  SET_VECTOR_ELT(result, 2, R_real_vector(stime));



  UNPROTECT(1); //result



  return(result);

   

}



// Function for estimating pattern waiting times with respect to a given mixture model

SEXP R_time(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_sampling_mode, 

		SEXP R_sampling_param, SEXP R_n, SEXP R_seed){

  // Load the necessary data in their corresponding C structures

  int L = INTEGER_VALUE(R_L); // number of genetic events 

  int  sampling_mode = (INTEGER_VALUE(R_sampling_mode) != 0)?EXPONENTIAL:CONSTANT; // sampling mode for the simulations: exponential or constant

  double sampling_param = NUMERIC_VALUE(R_sampling_param); // sampling parameter that corresponds to the sampling mode

  int n = INTEGER_VALUE(R_n); // number of simulations 



  // Setting the seed for srand

  if(INTEGER_VALUE(R_seed) == -1) 

    srand((unsigned)time(NULL));

  else 

    srand((unsigned)INTEGER_VALUE(R_seed));

  

  // Load the model in R_trees:

  vector alpha(length(R_trees));      // mixture parameter

  array<graph> G; // digraph

  array< map<node,string> > event; // genetic event

  array< map<edge,double> > cond_prob; // conditional probabilities

  array< map<int,node> > node_no; // node of mutation index



  R_get_graph(R_alpha, R_trees, alpha, G, event, cond_prob, node_no);



  // Map nodes onto event indices:

  array< map<node,int> > no_of_node(length(R_trees));

  for (int k = 0; k < length(R_trees); k++)

    for (int j = 0; j < L; j++)

      no_of_node[k][node_no[k][j]] = j;



  // Estimate exponential waiting times on edges:

  array< map<edge,double> >  lambda;  // exponential parameter lambda_i

  lambda = waiting_times(cond_prob, sampling_mode, sampling_param);



  SEXP R_wtime, R_wtime_comp;

  PROTECT(R_wtime = allocVector(VECSXP, length(R_trees)));

  

  for (int k = 0; k < length(R_trees); k++) { 

    // Estimate pattern waiting times for the branching k:

    array< list<double> > wtime = mtreemix_time(L, G[k], lambda[k], 

			  node_no[k], no_of_node[k], cond_prob[k], n);



    PROTECT(R_wtime_comp = allocVector(REALSXP, pow2(L-1)));

    // Compute the expected waiting times for all possible patterns for L events:

    for (int i=0; i < pow2(L-1); i++)	

	// Expected waiting time:

	REAL (R_wtime_comp) [i] = nonnegmean(wtime[i]);

      

    SET_VECTOR_ELT(R_wtime, k, R_wtime_comp);

    UNPROTECT(1); // R_wtime_comp

    // Clear the array for the new cycle

    wtime.clear();

  }

  // Create a matrix of all possible patterns for L events:

  integer_matrix pat(pow2(L-1), L);

  for (int i=0; i < pow2(L-1); i++)

      // Pattern:

      pat[i] = index2pattern(i, L);



  SEXP result;

  PROTECT(result = allocVector(VECSXP, 2));

  

  SET_VECTOR_ELT(result, 0, R_int_matrix(pat));

  SET_VECTOR_ELT(result, 1, R_wtime);



  UNPROTECT(2); // R_wtime, result



  return(result);

}



// Function for computing the (log-)likelihoods of a given mixture model and a set of patterns

SEXP R_likelihood(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_pattern) {

  // Load the necessary data in their corresponding C structures

  int L = INTEGER_VALUE(R_L); // number of genetic events 

  integer_matrix pattern = C_get_pattern(R_pattern); // the set of patterns

  int N = pattern.dim1(); // number of patterns (samples)



  // Load the model in R_trees in the corresponding C structures:

  vector alpha(length(R_trees));      // mixture parameter

  array<graph> G; // digraph

  array< map<node,string> > event; // genetic event

  array< map<edge,double> > cond_prob; // conditional probabilities

  array< map<int,node> > node_no; // node of mutation index



  R_get_graph(R_alpha, R_trees, alpha, G, event, cond_prob, node_no);



  // Compute log-likelihoods and weighted likelihoods for the given set of patterns:

  vector logL(N);

  matrix wlike(N,length(R_trees));



  for (int i=0; i<N; i++)

    {

      integer_matrix pat(1,L);

      for (int j=0; j<L; j++)

	pat(0,j) = pattern(i,j);

      

      // Log-likelihood of sample i

      logL[i] = mtreemix_loglike(pat, length(R_trees), alpha, G, node_no, cond_prob);

      

      // Weighted likelihoods in the component trees 1, ..., K

      for (int k=0; k<length(R_trees); k++)

	wlike(i,k) = alpha[k] * mtree_like(pat.row(0), G[k], node_no[k], cond_prob[k]);

    }

  

  SEXP result;

  PROTECT(result = allocVector(VECSXP, 2));

  

  SET_VECTOR_ELT(result, 0, R_real_matrix(wlike));

  SET_VECTOR_ELT(result, 1, R_real_vector(logL));



  UNPROTECT(1); //result



  return(result);

}



// Function fro generating a random Rtreemix model with R_K components and R_L events

SEXP R_random(SEXP R_K, SEXP R_L, SEXP R_star, SEXP R_uniform, SEXP R_min, SEXP R_max, SEXP R_seed) {



  int K = INTEGER_VALUE(R_K); // number of tree components

  int L = INTEGER_VALUE(R_L); // number of starting solutions for the k-means clustering

  int star = INTEGER_VALUE(R_star); // (=1) with star component, (=0) without

  int uniform = INTEGER_VALUE(R_uniform); // equal (=1) or unequal (=0) edge weights in noise component

  double min = NUMERIC_VALUE(R_min); // minimum value for the conditional probabilities on the tree edges

  double max = NUMERIC_VALUE(R_max); // maximum value for the conditional probabilities on the tree edges



  // Setting the seed for srand

  if(INTEGER_VALUE(R_seed) == -1) 

    srand((unsigned)time(NULL));

  else 

    srand((unsigned)INTEGER_VALUE(R_seed));



  array<string>              profile;       // names of the genetic events

  vector                     alpha(K);      // mixture parameter

  array< graph >             G(K);          // digraph

  array< map<int,node> >     node_no(K);    // node of mutation index

  array< map<node,string> >  event(K);      // substitution event

  array< map<edge,double> >  cond_prob(K);  // conditional probabilities



  // Generate random mixture tree

  mtreemix_random(K, L, profile, alpha, G, event, node_no, 

		  cond_prob, star, uniform, min, max);





  int j, curEl, num_nodes, num_edges;

  int curEd;

  // SEXP variables needed for converting all the structures connected with the mdoel from C to R  

  SEXP list_graphs, one_graph, klass;

  SEXP newNames, graphN, graphE;

  SEXP curRval, curEdges, curWeights;



  SEXP attrKlass, edObj, edData, edNames, wList, edWeights; 



  SEXP result, list_names, edg_mode;

  

  // The list that will contain the function result

  PROTECT(result = allocVector(VECSXP, 2));

  // Names of the lists that are part of the output

  PROTECT(list_names = allocVector(STRSXP, 2));

  SET_STRING_ELT(list_names, 0, mkChar("alpha")); // the weight vector of the model

  SET_STRING_ELT(list_names, 1, mkChar("graphs.mixture")); // the list of the graphs each for every tree component

  

  setAttrib(result, R_NamesSymbol, list_names);



  SET_VECTOR_ELT(result, 0, R_real_vector(alpha));



  // Build the list of graphs that give the trees composing the mixture model:    

  PROTECT(list_graphs = allocVector(VECSXP, K));



  klass = MAKE_CLASS("graphNEL");



  // Needed for the edgeData construction

  attrKlass = MAKE_CLASS("attrData");



  PROTECT(newNames = allocVector(STRSXP, 2));

  SET_STRING_ELT(newNames, 0, mkChar("edges"));

  SET_STRING_ELT(newNames, 1, mkChar("weights"));



  // Build each of the K graphs

  for (int k = 0; k < K; k++) {    

    PROTECT(one_graph = NEW_OBJECT(klass));

    // R version >= 2.7.0 - edgemode slot deprecated
	  PROTECT(edg_mode = allocVector(VECSXP, 1));  
	  
	  setAttrib(edg_mode, R_NamesSymbol, R_scalarString("edgemode"));
	  
    SET_VECTOR_ELT(edg_mode, 0, R_scalarString("directed"));  

    SET_SLOT(one_graph, Rf_install("graphData"), edg_mode);   

    

    num_nodes = G[k].number_of_nodes();

    num_edges = G[k].number_of_edges();

  

    if (num_nodes == 0) {

      SET_SLOT(one_graph, Rf_install("nodes"), allocVector(STRSXP, 0));

      SET_SLOT(one_graph, Rf_install("edgeL"), allocVector(VECSXP, 0)); 

    } else {

      // Make the node and edges lists of the graph object

      PROTECT(graphN = allocVector(STRSXP, num_nodes));

      PROTECT(graphE = allocVector(VECSXP, num_nodes));



      // Structures for the edgeData construction

      PROTECT(edObj = NEW_OBJECT(attrKlass));

      PROTECT(edData = allocVector(VECSXP, num_edges)); 

      PROTECT(edNames = allocVector(STRSXP, num_edges));

      curEd = 0;



      node n;

      edge e;

      j = 0;



      // Build the node list in R

      forall_nodes (n, G[k]) {    

	SET_STRING_ELT(graphN, j, STRING_ELT(R_scalarString(("E" + event[k][n]).c_str()), 0));

  

	PROTECT(curRval = allocVector(VECSXP, 2));

	setAttrib(curRval, R_NamesSymbol, newNames);

	if(G[k].outdeg(n) == 0) {

	  SET_VECTOR_ELT(curRval, 0, allocVector(INTSXP, 0));

	  SET_VECTOR_ELT(curRval, 1, allocVector(REALSXP, 0));

	} else {  

	  curEl = 0;

	  // The edge list for the node n

	  PROTECT(curEdges = allocVector(INTSXP, G[k].outdeg(n)));

	 

	  // The corresponding weights list for the node n

	  PROTECT(curWeights = allocVector(REALSXP, G[k].outdeg(n)));



	  // Construct the edgeData:

	  forall_out_edges(e, n) {

	   SET_STRING_ELT(edNames, curEd, 

			  STRING_ELT(R_scalarString(("E" + event[k][n] + "|" 

						     + "E" + event[k][G[k].target(e)]).c_str()), 0));   



	   PROTECT(wList = allocVector(VECSXP, 1));

	   PROTECT(edWeights = allocVector(REALSXP, 1));

	   setAttrib(wList, R_NamesSymbol, R_scalarString("weight"));

	   REAL (edWeights)[0] = cond_prob[k][e];

	   SET_VECTOR_ELT(wList, 0, edWeights);



	   SET_VECTOR_ELT(edData, curEd, wList);





	   INTEGER (curEdges)[curEl] = (int) (G[k].target(e)->getIndex()) + 1;

	   REAL (curWeights)[curEl] = cond_prob[k][e];

	  	   

	   UNPROTECT(2); //edWeights, wList

	   curEd++;

           curEl ++; 

	  }

	 

	  SET_VECTOR_ELT(curRval, 0, curEdges);

	  SET_VECTOR_ELT(curRval, 1, curWeights);

	  UNPROTECT(2); // curEdges, curWeights

	}

	SET_VECTOR_ELT(graphE, j, curRval);

	UNPROTECT(1); //curRval

	

	j ++;

      }

      setAttrib(graphE, R_NamesSymbol, graphN); 

      setAttrib(edData, R_NamesSymbol, edNames);

      // Install the necessary slots an R graph needs:

      SET_SLOT(edObj, Rf_install("default"), allocVector(VECSXP, 0));

      SET_SLOT(edObj, Rf_install("data"), edData);



      SET_SLOT(one_graph, Rf_install("edgeData"), edObj);

      SET_SLOT(one_graph, Rf_install("edgeL"), graphE);

      SET_SLOT(one_graph, Rf_install("nodes"), graphN);



      

      UNPROTECT(5); //graphN, graphE, edNames, edData, edObj

    }

    

    SET_VECTOR_ELT(list_graphs, k, one_graph);

    

    UNPROTECT(2); // one_graph, edg_mode

  }    



  SET_VECTOR_ELT(result, 1, list_graphs);

  

  UNPROTECT(2); // newNames, list_graphs

  

  UNPROTECT(2); //result, list_names



  return(result);

}



// Function for computing the entire distribution encoded by a given mixture model

SEXP R_distr(SEXP R_L, SEXP R_trees, SEXP R_alpha, SEXP R_sampling_mode, 

	     SEXP R_sampling_param, SEXP R_out_param) {

  // Load the necessary data in their corresponding C structures

  int L = INTEGER_VALUE(R_L); // number of genetic events 

  int  sampling_mode = INTEGER_VALUE(R_sampling_mode); // specify the mode (exponential or constant) for the sampling times for the observed input and output probabilities

  double sampling_param = NUMERIC_VALUE(R_sampling_param); // sampling parameter that corresponds to the sampling mode for the input probabilities

  double output_param = NUMERIC_VALUE(R_out_param); // sampling parameter that corresponds to the sampling mode for the output probabilities



  // Load the model in R_trees:

  vector                     alpha(length(R_trees));      // mixture parameter

  array< graph >             G;          // digraph

  array< map<int,node> >     node_no;    // node of mutation index

  array< map<node,string> >  event;      // genetic event

  array< map<edge,double> >  cond_prob;  // conditional probabilities



  R_get_graph(R_alpha, R_trees, alpha, G, event, cond_prob, node_no);



  if (sampling_mode != -1) {

    sampling_mode = (sampling_mode != 0)?EXPONENTIAL:CONSTANT;

    cond_prob = rescale_cond_prob(G, cond_prob, sampling_mode, 

				  sampling_param, output_param);

  }



  // Calculate the probabilities:

  vector prob = mtreemix_distr(L, alpha, G, node_no, cond_prob);

  SEXP R_prob;



  PROTECT(R_prob = allocVector(REALSXP, prob.dim()));

  

  integer_matrix pat(prob.dim(), L);

  for (int i=0; i < prob.dim(); i++)

    {

      // pattern:

      pat[i] = index2pattern(i, L);

      // Probability induced with the tree mixture model:

      REAL (R_prob) [i] = prob[i];

    }



  SEXP result;

  PROTECT(result = allocVector(VECSXP, 2));

  

  SET_VECTOR_ELT(result, 0, R_int_matrix(pat));

  SET_VECTOR_ELT(result, 1, R_prob);



  UNPROTECT(2); // R_prob, result





  return(result);

} 

} // extern "C"



