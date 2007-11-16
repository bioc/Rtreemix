#ifndef _WRAP_GRAPH_HH
#define _WRAP_GRAPH_HH

#include "Edge.hh"
#include "List.hh"
#include "Node.hh"

namespace replaceleda{

    class graph;
    
    list<edge> graphwrapper_all_edges(graph &G);
    list<node> graphwrapper_all_nodes(graph &G);

}

#endif
