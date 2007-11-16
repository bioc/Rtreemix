#ifndef _REPLACELEDA_HH
#define _REPLACELEDA_HH

typedef unsigned int uint;

#include "Array.hh"
#include "Graph.hh"
#include "List.hh"
#include "Map.hh"
#include "Matrix.hh"
#include "Queue.hh"
#include "Vector.hh"

#include <fstream>
#include <climits>


namespace replaceleda {

    const double DBL_MAX = 1.79769313486231470e+308;
    

    typedef mmatrix<int> integer_matrix;
    typedef mmatrix<double> matrix;
    
    typedef mvector<int> integer_vector;
    typedef mvector<double> vector;

    typedef std::string string;

    typedef int integer;    

//    template<class T>
    bool member(std::set<node> s, node e);
    bool member(std::set<edge> s, edge e);    

//    template <class T>
    std::vector<int> permute(std::vector<int> v);

//    template <class T>
    list<int> permute(list<int> l);

    graph *graph_of(node n);

    graph *graph_of(edge e);

    node source(edge e);

    node target(edge e);

    void mydfs(graph &G, list<edge> &cycle, std::set<node> &unseen, std::set<node> &on_path, node v);
    bool Is_Acyclic(graph &G, list<edge> &cycle);

    void CopyGraph(graph &g1, graph &g2);
    void CopyGraph(GRAPH<node, edge>  &g_target, graph &g_source);

    //split string
    std::vector<replaceleda::string> strsplit(replaceleda::string inputstr, replaceleda::string splitseq);

    //does something like the formated output in C
    replaceleda::string tostring(replaceleda::string inputstr, ...);

    //from Matrix.hh
    //template <class T>
    std::ostream& operator<< (std::ostream &os, mmatrix<double> &m);
    std::ostream& operator<< (std::ostream &os, mmatrix<int> &m);

//    template <class T>
    std::istream& operator>> (std::istream &is, mmatrix<double> &m);
    std::istream& operator>> (std::istream &is, mmatrix<int> &m);
   
//    template<class T>
    mmatrix<double> transpose(mmatrix<double> &m);
    mmatrix<int> transpose(mmatrix<int> &m);

    std::ostream& operator<< (std::ostream &os, std::set<node> s);

    //from Vector.hh    
//    template <class T>
    std::ostream& operator<< (std::ostream &os, mvector<double> v);
    std::ostream& operator<< (std::ostream &os, mvector<int> v);
    
//    template <class T>
    std::istream& operator>> (std::istream &is, mvector<double> &v);    
    std::istream& operator>> (std::istream &is, mvector<int> &v);    

    std::ostream& operator<< (std::ostream &os, graph &g);

    void printGraph (graph &g, edge_array<double>  &weights);

}

//#include "replaceleda.cc"

#endif // _REPLACELEDA_HH
