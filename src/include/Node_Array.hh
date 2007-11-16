#ifndef _NODE_ARRAY_HH
#define _NODE_ARRAY_HH

#include "List.hh"
#include "Node.hh"
//#include "NodeN.hh"
#include "WrapGraph.hh"
#include <map>

namespace replaceleda{
   
    template <class T>
    class node_array {
    public:
	node_array() : mygraph(0) {}
	virtual ~node_array(){
	    mygraph = NULL;
	    storage.clear();
	}
	
	T operator[] (node n) const { return storage[n];}
	T& operator[] (node n) {
	    if(storage.find(n) == storage.end())
		storage.insert(std::make_pair(n, T()));
	    return storage[n];
	}

	graph* get_graph(){ return mygraph;}

	void set(node n, T elem){
	    storage[n] = elem;
	}

	void clear(){ storage.clear(); }


	void init(graph &G){
	    node v;
	    T e = T(0);
	    //list<node> nl = G.all_nodes();
	    list<node> nl = graphwrapper_all_nodes(G);
	    for(uint n = 0; n < nl.size(); ++n){
		v = nl[n];	   
		storage.insert(std::make_pair(v,e));
	    }
	}

	void init(graph &G, T def){
	    node v;	    
	    //list<node> nl = G.all_nodes();
	    list<node> nl = graphwrapper_all_nodes(G);
	    for(uint n = 0; n < nl.size(); ++n){
		v = nl[n];	   
		storage.insert(std::make_pair(v,def));
	    }
	}

    protected:
	std::map<node, T> storage;
	graph *mygraph;
    };

}


#endif //_Node_Array_hh
