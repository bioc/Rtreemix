
#include "include/replaceleda.hh"

#include <cstdarg>
#include <sstream>

//R includes
#include <R.h>

//using namespace replaceleda;

//template<class T>
bool replaceleda::member(std::set<replaceleda::node> s, replaceleda::node e){
    return (s.find(e) != s.end());
}

bool replaceleda::member(std::set<replaceleda::edge> s, replaceleda::edge e){
    return (s.find(e) != s.end());
}


//template <class T>
std::vector<int> replaceleda::permute(std::vector<int> v){
    std::vector<int> dummy = v, result;
    //srand((unsigned)time(NULL));
    while(dummy.size() > 0){
	int r = int(dummy.size() * (float) rand() / (float) (RAND_MAX + 1.0));
	result.push_back(dummy[r]);
	dummy.erase(dummy.begin()+r);
    }
    return result;
}

//template <class T>
replaceleda::list<int> replaceleda::permute(list<int> l){
    replaceleda::list<int> result = l;
    std::vector<int> dummy;
    while(result.size() > 0){
	dummy.push_back(result.front());
	result.pop_front();
    }
    dummy = permute(dummy);
    while(dummy.size() > 0){	    
	result.push_front(dummy.back());
	dummy.pop_back();
    }
    return result;
}

//algorithms from Graph.hh

replaceleda::graph *replaceleda::graph_of(replaceleda::node n){
    return n->get_graph();
}

replaceleda::graph *replaceleda::graph_of(replaceleda::edge e){
    return e->get_graph();
}


replaceleda::node replaceleda::source(replaceleda::edge e) {
    return e->getSource();
}

replaceleda::node replaceleda::target(replaceleda::edge e) {
    return e->getTarget();
}     

void replaceleda::mydfs(replaceleda::graph& G, replaceleda::list<replaceleda::edge>& cycle, std::set<replaceleda::node>& unseen, std::set<replaceleda::node>& on_path, replaceleda::node v){

}

bool replaceleda::Is_Acyclic(replaceleda::graph &G, replaceleda::list<replaceleda::edge> &cycle){
    
    //use depth first search to detect cycles

    std::set<replaceleda::node> unseen, on_path;
    replaceleda::node v;

    std::vector<node> stack, tb_stack;
    forall_nodes(v,G)
	unseen.insert(v);
    on_path.clear();
    while (! unseen.empty()){
	on_path.clear();
	v = *(unseen.begin());
	stack.push_back(v);
	while (!stack.empty()){	    
	    v = stack.back();
	    stack.pop_back();
	    while ((v == NULL) && !stack.empty()){
		if (v == NULL){
		    on_path.erase(tb_stack.back());
		    tb_stack.pop_back();
		}
		v = stack.back();
		stack.pop_back();
	    }
	    if (v != NULL){
		unseen.erase(v);		
		on_path.insert(v);
		stack.push_back(NULL);
		tb_stack.push_back(v);
		edge e;
		
		forall_out_edges(e,v){
		    node w = target(e);
		    if( member(on_path, w)){
			cycle.push_back(e);
		    } else if (member(unseen, w))
			stack.push_back(w);		
		}
	    }
	}
    }
    return cycle.empty();
}


void replaceleda::CopyGraph(replaceleda::graph &g_target, replaceleda::graph &g_source){

    node v, w;
    edge e;
    node_array<node> mapping;

    g_target.clear();
    //create nodes!
    forall_nodes(v, g_source){
	w = g_target.new_node();
	mapping[v] = w;
    }

    forall_edges(e,g_source){
	v = mapping[source(e)];
	w = mapping[target(e)];
	g_target.new_edge(v,w);
    }

};

void replaceleda::CopyGraph(replaceleda::GRAPH<replaceleda::node, replaceleda::edge>  &g_target, replaceleda::graph &g_source){

    node v, w;
    edge e, ne;
    node_array<node> mapping;

    g_target.clear();
    //create nodes!
    forall_nodes(v, g_source){
	w = g_target.new_node();
	mapping[v] = w;
	g_target[w] = v;
    }

    forall_edges(e,g_source){
	v = mapping[source(e)];
	w = mapping[target(e)];
	ne = g_target.new_edge(v,w);
	g_target[ne] = e;
    }

};


std::vector<std::string> replaceleda::strsplit(replaceleda::string inputstr, std::string splitseq){

    std::vector<std::string> result;
    std::string::size_type pos, opos = 0;
    
    pos = inputstr.find(splitseq);
    
    while (pos != std::string::npos){
        std::string sub = inputstr.substr(opos, pos - opos);
        result.push_back(sub);
        opos = pos + 1;
        pos = inputstr.find(splitseq, opos+1);
    }

    result.push_back(inputstr.substr(opos, pos - opos));

    return result;

}

replaceleda::string replaceleda::tostring(replaceleda::string inputstr, ...){

    std::vector<replaceleda::string> token;
    va_list arguments;

    va_start(arguments, inputstr);

    std::ostringstream result;

    token = strsplit(inputstr, " ");
    for(std::vector<replaceleda::string>::iterator iter = token.begin(); iter != token.end(); ++iter){  
            if (*iter == "%d"){
                int value = va_arg(arguments, int);
                result << " " << value;
            } else if (*iter == "%f"){
                double value = va_arg(arguments, double);
                result << " " << value;             
            } else
                result << " " << *iter;     
    }

    return result.str().substr(1);
}



//from Matrix.hh

//template <class T>
std::ostream& replaceleda::operator<< (std::ostream &os, replaceleda::mmatrix<double> &m){
    os << m.dim1() << " " << m.dim2() << std::endl;
    for(int i=0;i<m.dim1();++i){
	for(int j=0;j<m.dim2();++j)
	    //os << m.elem(i,j) << " ";
	    os << m[i][j] << " ";
	os << std::endl;	    
    }
    return os;
}

std::ostream& replaceleda::operator<< (std::ostream &os, replaceleda::mmatrix<int> &m){
    os << m.dim1() << " " << m.dim2() << std::endl;
    for(int i=0;i<m.dim1();++i){
	for(int j=0;j<m.dim2();++j)
	    //os << m.elem(i,j) << " ";
	    os << m[i][j] << " ";
	os << std::endl;	    
    }
    return os;
}


//template <class T>
std::istream& replaceleda::operator>> (std::istream &is, replaceleda::mmatrix<double> &m){
    int N,M, tot;
    double e;
    replaceleda::mvector<double> dummy;
    std::string line;

    //read N and M
    is >> N; is >> M;
    tot = N*M;
    for(int i = 0; i < tot; ++i){
	is >> e;
	dummy.push_back(e);
    }
    if(dummy.dim() == N * M){
	replaceleda::mmatrix<double> result(N,M,dummy);
	m = result;	
    }
    return is;
}

std::istream& replaceleda::operator>> (std::istream &is, replaceleda::mmatrix<int> &m){
    int N,M, tot;
    int e;
    replaceleda::mvector<int> dummy;
    std::string line;

    //read N and M
    is >> N; is >> M;
    tot = N*M;
    for(int i = 0; i < tot; ++i){
	is >> e;
	dummy.push_back(e);
    }
    if(dummy.dim() == N * M){
	replaceleda::mmatrix<int> result(N,M,dummy);
	m = result;	
    }
    return is;
}



//template<class T>
replaceleda::mmatrix<double> replaceleda::transpose(replaceleda::mmatrix<double> &m){
    return m.trans();
}

replaceleda::mmatrix<int> replaceleda::transpose(replaceleda::mmatrix<int> &m){
    return m.trans();
}

std::ostream& replaceleda::operator<< (std::ostream &os, std::set<replaceleda::node> s){
    os << "{";
    for(std::set<replaceleda::node>::iterator it = s.begin(); it != s.end(); ++it)
	os << (*it)->getIndex() << ", ";
    os << "}" << std::endl;
    return os;
}

//from Vector.hh

//template <class T>
std::ostream& replaceleda::operator<< (std::ostream &os, replaceleda::mvector<double> v){
    os << v.size() << " ";
    for(uint j=0;j<v.size();++j)
	os << v[j] << " ";
    //os << std::endl;
    return os;
} 

std::ostream& replaceleda::operator<< (std::ostream &os, replaceleda::mvector<int> v){
    os << v.size() << " ";
    for(uint j=0;j<v.size();++j)
	os << v[j] << " ";
    //os << std::endl;
    return os;
} 

//template <class T>
std::istream& replaceleda::operator>> (std::istream &is, replaceleda::mvector<double> &v){
    double e;
    uint dim;
    v.clear();

    is >> dim;
    for (uint i = 0; i < dim; ++i){
	is >> e;
	v.push_back(e);
    }
    return is;
}

std::istream& replaceleda::operator>> (std::istream &is, replaceleda::mvector<int> &v){
    int e;
    uint dim;
    v.clear();

    is >> dim;
    for (uint i = 0; i < dim; ++i){
	is >> e;
	v.push_back(e);
    }
    return is;
}


std::ostream& replaceleda::operator<< (std::ostream &os, replaceleda::graph &g){

    node v;
    edge e;
    
    os << "#nodes: " << g.number_of_nodes() << " #edges: " << g.number_of_edges() << std::endl;
    list<node> xxx = g.all_nodes();
    
    forall_nodes(v, g){
	os << v << " " << v->getIndex() << std::endl;
	os << "(" << g.indeg(v) << "," << g.outdeg(v) << "," << g.degree(v) << "):" << std::endl;
	
	forall_out_edges(e, v){
	    node t = target(e);
	    os << v->getIndex() << " --> " << t->getIndex() << std::endl; 
	}
    }

    return os;
}

void replaceleda::printGraph (replaceleda::graph &g, edge_array<double> &weights){

   node v;
    edge e;
    
    std::cerr << "#nodes: " << g.number_of_nodes() << " #edges: " << g.number_of_edges() << std::endl;
    list<node> xxx = g.all_nodes();
    
    forall_nodes(v, g){
	std::cerr << v << " " << v->getIndex() << std::endl;
	std::cerr << "(" << g.indeg(v) << "," << g.outdeg(v) << "," << g.degree(v) << "):" << std::endl;
	
	forall_out_edges(e, v){
	    node t = target(e);
	    std::cerr << v->getIndex() << " -" << weights[e] << "-> " << t->getIndex() << std::endl; 
	}
    }

   

}

//template class replaceleda::mmatrix<double>

