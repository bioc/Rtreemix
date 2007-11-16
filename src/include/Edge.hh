
#ifndef _EDGE_HH
#define _EDGE_HH

#include "RefCounted.hh"
#include "Set.hh"
#include <iostream>
#include <string>

namespace replaceleda{

    class graph;
    class Node;
    //class node;
    //typedef Node* node;
    typedef RefCountPtr<Node> node;

    class Edge : public RefCounted {
    public:
	Edge(Node* _source = NULL, Node* _target = NULL, graph *_graph = NULL, double _weight = 0.0) : 
	    weight(_weight), source(_source), target(_target), mygraph(_graph) {}
	virtual ~Edge(){
	    source = 0;
	    target = 0;
	    mygraph = 0;
        }
	    
	void set_graph(graph *g){ mygraph = g;}
	graph* get_graph(){ return mygraph;}
	
	Node* getSource(){ return source;}
	Node* getTarget(){ return target;}
	
	void setWeight(double _w)  { weight = _w;}
	double getWeight(){ return weight;}

	//some overloaded operators
	bool operator==(Edge e){
	    if (mygraph != e.get_graph() ||
		target != e.getTarget() ||
		source != e.getSource())
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

	//node source, target;
	Node *source, *target;
	graph *mygraph;
    };

    typedef RefCountPtr<Edge> edge;

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

#endif //_Edge_hh
