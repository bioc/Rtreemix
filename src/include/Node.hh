
#ifndef _NODE_HH
#define _NODE_HH

//#include <list>
#include <set>
#include <iostream>

#include "Edge.hh"
#include "List.hh"
#include "RefCounted.hh"

namespace replaceleda{

    class graph;

    class Node : public RefCounted {
    public:
	Node () {}
	Node(graph *g, uint _idx = 0) :
	    index(_idx), mygraph(g)
	    {
		adj_e.clear();
		in_e.clear();
		out_e.clear();
	    }

	virtual ~Node() {	    

            //time of destruction has come!
	    //delete all out_arcs!
	    while (! out_e.empty()){
		edge e = out_e.front();
		Node* partner = e->getTarget();
		partner->del_edge_adj(e); del_edge_adj(e);
		partner->del_edge_in(e); //del_edge_out(e);
		out_e.remove(e);
		//delete e;
	    }
	    //delete all in_arcs!
	    while (! in_e.empty()){
		edge e = in_e.front();
		Node* partner = e->getSource();
		partner->del_edge_adj(e); del_edge_adj(e);
		partner->del_edge_out(e); //del_edge_in(e);
		in_e.remove(e);
		//delete e;
	    }

	    mygraph = NULL;
	    adj_e.clear();
	    in_e.clear();
	    out_e.clear();

	}

	graph* get_graph(){ return mygraph;}
	
	uint getIndex(){ return index;}

	void add_edge_adj(edge e){add_edge(e,1);}
	void add_edge_in(edge  e){add_edge(e,2);}
	void add_edge_out(edge e){add_edge(e,3);}

	void del_edge_adj(edge e){del_edge(e,1);}
	void del_edge_in(edge e) {del_edge(e,2);}
	void del_edge_out(edge e){del_edge(e,3);}


	list<edge> adj_edges(){return adj_e;}
	list<edge> in_edges() {return in_e; }
	list<edge> out_edges(){return out_e;}
	
	void del_all_adj(){ adj_e.clear(); }
	void del_all_in(){ in_e.clear(); }
	void del_all_out(){ out_e.clear(); }

	void del_all(){ del_all_adj(); del_all_in(); del_all_out();}

    protected:       
	uint index;	
	graph *mygraph;
	list<edge> adj_e, in_e, out_e;

    protected:
	void add_edge(edge e, uint eType){
	    switch(eType){
	    case 1: adj_e.push_back(e); break;
	    case 2: in_e.push_back(e); break;
	    case 3: out_e.push_back(e); break;
	    }
	}

	//here the edges in the partner should be removed!
	void del_edge(edge e, uint eType){	    
	    switch(eType){
	    case 1: adj_e.remove(e); break;
	    case 2: in_e.remove(e); break;
	    case 3: out_e.remove(e); break;
	    }
	}	
    };
    
    struct int_node {
	int i;
	node n;
	
	bool operator<( const int_node &in) const {
	    //return i < in.i;
	    return smallontop(in);
	}
	bool operator>( const int_node &in) const {
	    return i > in.i;
	}
	bool smallontop( const int_node &in) const {
	    return i > in.i;
	}
	bool bigontop( const int_node &in) const {
	    return i < in.i;
	}
    };

    class node_set : public set<node> {
    public: 
	node_set() {}
	node_set(graph &G){}
	virtual ~node_set() {}
    protected:

    };
}




#endif //_Node_hh
