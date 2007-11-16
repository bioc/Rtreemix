#ifndef _EDGE_ARRAY_HH
#define _EDGE_ARRAY_HH

#include "Edge.hh"
//#include "NodeN.hh"
#include "List.hh"
#include "WrapGraph.hh"
#include <map>

namespace replaceleda{

    template <class T>
    class edge_array {
    public:
	edge_array() : mygraph(0) {}
	virtual ~edge_array(){
	    mygraph = NULL;
	    storage.clear();
	}
	
	T operator[](edge e) const { return storage[e];}
	T& operator[](edge e) {
	    if (storage.find(e) == storage.end())  
		storage.insert(std::make_pair(e, T()));
	    return storage[e];
	}

	void set(edge e, T elem){
	    storage[e] = elem;
	}

	graph *get_graph(){ return mygraph;}


	void init(graph &G){
	    edge e;
	    T elem = T(0);
	    //list<edge> el = G.all_edges();
	    list<edge> el = graphwrapper_all_edges(G);
	    for(uint n = 0; n < el.size(); ++n){
		e = el[n];	   
		storage.insert(std::make_pair(e,elem));
	    }
	}

//	void init (graph &G);
	void clear(){  storage.clear(); }
	void erase(T item){ storage.erase(item);}


	void init(graph &G, T def){
	    edge e;
	    //list<edge> el = G.all_edges();
	    list<edge> el = graphwrapper_all_edges(G);
	    for(uint n = 0; n < el.size(); ++n){
		e = el[n];	   
		storage.insert(std::make_pair(e,def));
	    }
	}

//	void init(graph &G, T def);
    protected:
	std::map<edge, T> storage;
	graph *mygraph;

    };

}


#endif // _Edge_Array_hh
