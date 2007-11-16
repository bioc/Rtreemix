
#ifndef _NODEN_HH
#define _NODEN_HH

#include <string>
#include <set>
#include <iostream>

#include "Set.hh"
#include "List.hh"

namespace replaceleda{

    class graph;

    class Edge;
    typedef Edge* edge;

    class Node;
    class node;

    class Node {
    public:
	Node () : m_refCount(0) {}
	Node(graph *g, uint _idx = 0) :
	    index(_idx), mygraph(g), adj_e(0), in_e(0), out_e(0), m_refCount(0) {}

	virtual ~Node(){
	    //ASSERT(0==m_refCount);
	}
	
	//this this needed?
	//void set_graph(graph *g){ mygraph = g;}
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

	void del_edge(edge e, uint eType){	    
	    switch(eType){
	    case 1: adj_e.remove(e); break;
	    case 2: in_e.remove(e); break;
	    case 3: out_e.remove(e); break;
	    }
	}	
    private:
	uint mutable m_refCount;
	void AddRef() {
	    ++m_refCount;
	}
	void Release() {
	    if (0 == --m_refCount)
		delete this;
	}
	
	friend class node;
    };

    class node{
    public:
	node(Node* pTarget = 0) : m_pTarget(pTarget) {
	    if (m_pTarget)
		m_pTarget->AddRef();
	}
	node(node& v) : m_pTarget(v.m_pTarget) {
	    if (m_pTarget)
		m_pTarget->AddRef();
	}
	~node (){
	    if (m_pTarget)
		m_pTarget->Release();
	}
	node& operator=(node& v) {	   
	    if (v.m_pTarget)
		v.m_pTarget->AddRef();
	    if (m_pTarget)
		m_pTarget->Release();
	    m_pTarget = v.m_pTarget;
	    return *this;
	}
	node& operator=(Node* pTarget) {
	    if(pTarget)
		pTarget->AddRef();
	    if(m_pTarget)
		m_pTarget->AddRef();
	    m_pTarget = pTarget;
	    return *this;
	}
	
	Node& operator*() {
	    return *m_pTarget;
	}
	Node* operator->() {
	    return m_pTarget;
	}

/*	bool operator==(node v){
	    return (m_pTarget==v.operator->());
	}
	bool operator!=(node v){
	    return (!operator==(v));
	}
*/
    private:
	Node*  m_pTarget;
    };
    
    //graph *graph_of(node n);
    //std::ostream& operator<<(std::ostream &os, node v);

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

    class Edge {
    public:
	Edge(node _source = NULL, node _target = NULL, graph *_graph = NULL, double _weight = 0.0) : 
	    weight(_weight), source(_source), target(_target), mygraph(_graph) {}
	virtual ~Edge(){}
	    
	void set_graph(graph *g){ mygraph = g;}
	graph* get_graph(){ return mygraph;}
	
	node getSource(){ return source;}
	node getTarget(){ return target;}
	
	void setWeight(double _w)  { weight = _w;}
	double getWeight(){ return weight;}


	//some overloaded operators
	bool operator==(Edge e){
	    if (mygraph != e.get_graph())/* or 
		target != e.getTarget() or
		source != e.getSource())*/
		return false;
	    return true;
	}
	bool operator!=(Edge e){
	    return (!operator==(e));
	}
	bool operator<(const Edge e) const {
	    return weight < e.weight;
	}
	bool operator>(const Edge e) const {
	    return weight > e.weight;
	}

    public:
	std::string label;
	double weight;
    protected:

	node source, target;
	graph *mygraph;
    };

    const edge nullEdge = NULL;

    //for priority queue, we want small items on top !
    struct double_edge {
	double d;
	edge e;

	bool operator<(const double_edge &de) const{
	    //return d < de.d;
	    return smalltotop(de);
	}
	bool operator>(const double_edge &de) const{
	    return d > de.d;
	}
	bool smalltotop(const double_edge &de) const{
	    return d > de.d;
	}
	bool bigtotop(const double_edge &de) const{
	    return d < de.d;
	}
    };


    class edge_set : public set<edge> {
    public:
	edge_set() {}
	edge_set(graph &G){}
	virtual ~edge_set(){}
    protected:

    };

}

#endif //_NodeN_hh
