
#ifndef _GRAPH_HH
#define _GRAPH_HH

#include <vector>
#include <map>

#include "Edge.hh"
#include "Edge_Array.hh"
#include "List.hh"
#include "Node.hh"
#include "Node_Array.hh"


//caution just a dummy until I decide for a final layout of class Graph!!!
#define forall_nodes(n,G) for(uint icount = 0; icount < G.number_of_nodes() ,  n = G.getNode(icount);++icount)
#define forall_edges(e,G) for(uint jcount = 0; e = G.getEdge(jcount), jcount < G.number_of_edges();++jcount)
#define forall_out_edges(e,n) for(uint kcount = 0; e = n->out_edges().get_item(kcount) , kcount < n->out_edges().size()/*graph_of(n)->outdeg(n)*/; kcount++)
#define forall_in_edges(e,n) for(uint kcount = 0; e = n->in_edges().get_item(kcount), kcount < graph_of(n)->indeg(n); ++kcount)
#define forall_defined(el,map) for(uint kcount = 0; el = map.getDefined().get_item(kcount), kcount < map.getDefined().size(); ++kcount)
#define forall(e,L) for(uint ecount = 0; e = L.get_item(ecount), ecount < L.size(); ++ecount)

namespace replaceleda{

    //typedef std::set<node> node_set;
    
    class graph {
    public:
	graph() : directed(1), maxnodeindex(0){
	    mynodes.clear();
	    myedges.clear(); 
	}

	virtual ~graph(){
	    clear();
	}

	virtual void clear(){
	    mynodes.clear();
	    myedges.clear();
	    maxnodeindex=0;
	}

	//clean the graph
	//private?
	void del_all_nodes(){
	    //while(! mynodes.empty()){
		//node n = mynodes.front();
		//del_node(n);
	    //}
	    mynodes.clear();
	    maxnodeindex = 0;
	    updateEdgesInGraph();
	}

	void del_all_edges(){
	    while(! myedges.empty()){
		edge e = myedges.front();
		del_edge(e);		
	    }
	    myedges.clear();
	}

	//node operations
	node new_node(){
	    node result;
	    result = new Node(this, maxnodeindex++);
	    mynodes.push_back(result);	    
	    return result;
	}

	void del_node(node n){
	    mynodes.remove(n);
	    n=0;
	    //eges are handeled in destructor of node!
            updateEdgesInGraph();
	}

	node getNode(uint n){ return mynodes.get_item(n);}

	uint outdeg(node n){ return n->out_edges().size();}
	uint indeg(node n ){
	    if (directed)
		return n->in_edges().size();
	    else
		return 0;
	}
	uint degree(node n){ return n->adj_edges().size();}
	
	list<edge> in_edges(node v) { 
	    return v->in_edges();
	}
	list<edge> out_edges(node v){
	    return v->out_edges();
	}
	list<edge> adj_edges(node v){
	    return v->adj_edges();
	}

	list<node> adj_nodes(node v){
	    list<node> result;
	    for(list<edge>::iterator eit = v->adj_edges().begin(); eit != v->adj_edges().end(); ++eit)
		result.push_back(opposite(v, *eit));
	    return result;
	}

	list<node> all_nodes(){return mynodes;}
	node first_node(){ return mynodes.front();}

	//edge operations
	//new edge between two nodes v and w
	edge new_edge(node v, node w){
	    //create edge between v, w: v->w
	    edge result = new Edge(v,w,this);
	    myedges.push_back(result);
	    updateEdgesInNodes(v,w,result);
	    return result;
	}
	
	//new edge between source of e and w
	edge new_edge(edge e, node w){
	    //create edge between v, w: v->w
	    node v = e->getSource();
	    edge result = new Edge(v,w,this);
	    myedges.push_back(result);
	    updateEdgesInNodes(v,w,result);	    
	    return result;
	}

	//new edge between v and target of e
	edge new_edge(node v, edge e){
	    //create edge between v, w: v->w
	    node w = e->getTarget();
	    edge result = new Edge(v,w,this);
	    myedges.push_back(result);
	    updateEdgesInNodes(v,w,result);	    
	    return result;
	}

	virtual void del_edge(edge e){
	    //del entry in nodes
	    node s = e->getSource(), t = e->getTarget();
	    if (s){
		s->del_edge_out(e); s->del_edge_adj(e);
	    }
	    if (t){
		t->del_edge_in(e); t->del_edge_adj(e);
	    }
	    if(! directed){
		if (s){
		s->del_edge_in(e);
		}
		if (t){
		    t->del_edge_out(e);
		}
	    }
	    //del edge
	    myedges.remove(e);
	    //delete e;
	    //e = 0;
	}

	edge first_edge () { return myedges.front();}
	node source (edge e){ return (e->getSource());}
	node target (edge e){return (e->getTarget());}

	edge getEdge(uint e){ return myedges.get_item(e);}

	list<edge> all_edges(){
	    return myedges;
	}

	node opposite(node v, edge e){
	    node result = e->getSource();
	    if( v == result)
		return e->getTarget();
	    return result;
	}
	
	uint number_of_nodes() {return mynodes.size();}
	uint number_of_edges() {return myedges.size();}

    protected:
	bool directed;
	uint maxnodeindex;

	list<node> mynodes;
	list<edge> myedges;


    protected:
	void updateEdgesInNodes(node v, node w, edge e){
	    //let v, w know that they have a new edge
	    v->add_edge_adj(e); w->add_edge_adj(e); 
	    v->add_edge_out(e); w->add_edge_in(e);
	    if(!directed){
		v->add_edge_in(e); w->add_edge_out(e);
	    }
	}

	void updateEdgesInGraph(){
	    myedges.clear();
	    edge ce;
	    //just works with directed graphs... for now
	    for(list<node>::iterator cn = mynodes.begin(); cn != mynodes.end(); ++cn)
		//for(list<edge>::iterator ce = (*cn)->out_edges().begin(); ce != (*cn)->out_edges().end(); ++ce)
		for(uint kcount = 0; ce = (*cn)->out_edges().get_item(kcount), kcount < (*cn)->out_edges().size(); ++kcount)
		    myedges.push_back(ce);
		
	}	


};


    //Parametreized GRAPH
    template <class T1, class T2>
    class GRAPH : public graph {
	//carries additional information for every node (T1) and every edge(T2)
    public:
	GRAPH () : graph(){
	    mynode_info.clear();
	    myedge_info.clear();
	    maxnodeindex = 0;
	}
	virtual ~GRAPH(){
	    mynode_info.clear();
	    myedge_info.clear();
	}

	virtual void clear(){
	    graph::clear();
	    mynode_info.clear();
	    myedge_info.clear();
	}

	T1 inf(node n){ return mynode_info[n];}
	T2 inf(edge e){ return myedge_info[e];}

	T1 operator[](node n) const {return mynode_info[n];}
	T2 operator[](edge e) const {return myedge_info[e];}
	T1& operator[](node n){ return mynode_info[n];}
	T2& operator[](edge e){ return myedge_info[e];}

	void assign(node n, T1 elem){ mynode_info[n] = elem;}
	void assign(edge e, T2 elem){ myedge_info[e] = elem;}
	
	node_array<T1> node_data(){ return mynode_info;}
	edge_array<T2> edge_data(){ return myedge_info;}
	
    protected:
	node_array<T1> mynode_info;
	edge_array<T2> myedge_info;
    };

}

#endif // _Graph_hh
