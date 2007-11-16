#include "include/replaceleda.hh"

//#include "include/Graph.hh"
//#include "include/WrapGraph.hh"

//using namespace replaceleda;

replaceleda::list<replaceleda::edge> replaceleda::graphwrapper_all_edges(replaceleda::graph &G){

    return G.all_edges();
}

replaceleda::list<replaceleda::node> replaceleda::graphwrapper_all_nodes(replaceleda::graph &G){

    return G.all_nodes();
}
