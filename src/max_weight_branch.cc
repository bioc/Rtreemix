/*
  max_weight_branch.c

  Version : 1.1
  Author  : Niko Beerenwinkel
  
  (c) Max-Planck-Institut fuer Informatik, 2004
*/


#include "include/max_weight_branch.h"


// global variables:

static map<node,string> Mut;  // node labels
static map<edge,double> c;    // edge weights

static array< array<node> > cycle_contr;
// cycle_contr[i][k] == v_C <==> the k-th cycle C in B[i] has been contracted to v_C in G[i+1]

static array< array<edge> > cheapest_edge;
// cheapest_edge[i][k] is some cheapest edge on the k-th cycle in B[i]

static map<edge,edge> Phi;
// Phi[e'] == e  <==> e in G[i] has been contracted to e' in G[i+1]

static array< array< map<edge,edge> > > alpha;  
// e = (z, y), z not in C, y in C  <==>  alpha[e] == pred of y in C 



// functions:

int compare_weights(const edge& e1, const edge& e2)
{
  if (c[e1] < c[e2])  return -1;
  if (c[e1] > c[e2])  return  1;
  return 0;
}



double BRANCHING_WEIGHT(list<edge>& B, edge_array<double>& weight)
{
  // total weight of a branching
  
  double total_weight = 0;
  edge e;

  forall(e, B)
    total_weight += weight[e];

  return total_weight;
}



void UNCOVER_BRANCHING(graph& G, list<edge>& B)
{
  // turn G into the branching defined by B (B \subset E(G))

  edge e;

  edge_set E_B(G);  // edge set of B

  forall(e, B)
    E_B.insert(e);

  edge_set T_D(G); //edges to delete

  forall_edges(e, G)
    if (! E_B.member(e))
	T_D.insert(e);
      //G.del_edge(e);

  for(edge_set::iterator it = T_D.begin(); it != T_D.end(); ++it)
      G.del_edge(*it);

}



void DOT(graph& G, map<node,string>& mut, map<edge,double>& cond_prob, char *filestem)
{
  // produces output in dot format (graphviz package)
  // see http://www.research.att.com/sw/tools/graphviz/

  node v;  edge e;

  char filename[128];
  sprintf(filename, "%s.dot", filestem);
  
  std::ofstream dotout(filename);
  dotout << "digraph MWB {" << std::endl << std::endl;
  
  forall_nodes(v, G)
    {
      dotout << "\t \"" << v << "\"";
      dotout << " [ label=\"" << mut[v] << "\", shape=\"plaintext\", height=\"0.3\", fontsize=\"12\", style=\"filled\", fillcolor=\"white\" ];" << std::endl;
    }
  dotout << std::endl;

  forall_edges(e, G)
    {
      node s = G.source(e);
      node t = G.target(e);
      dotout.precision(2);
      dotout << std::showpoint;
      dotout << "\t \"" << s << "\" -> \"" << t << "\"";
      dotout << " [ fontsize=\"10\", label=\"" << cond_prob[e] << "\" ];" << std::endl;
    }

  dotout << "}" << std::endl;
  dotout.close();
}


double branching_weight_intern(list<edge>& B)
{
  // total weight of a branching
  
  double w = 0;

  edge e;
  forall(e, B)
    w += c[e];
  
  return w;
}



list<edge> max_weight_subgraph_indeg_le_1(graph& G)
{
  // find maximum weight subgraph of G with indeg <= 1 for all v in V(G)
  // a list of edges in G of this subgraph is returned
  
  node v;
  list<edge> B;

  forall_nodes(v, G)  // choose heaviest inedge
    {
      list<edge> L = G.in_edges(v);
      if (! L.empty())
	B.append(L.contents(L.max(*compare_weights)));
    }

  return B;
}


array< list<edge> > all_cycles(graph& G, list<edge>& B)
{
  // find all cycles in B

  edge e;

  array< list<edge> > CA;  // array of cycles
  int k = -1;              // cycle index

  // construct subgraph:
  GRAPH<node,edge> G_B;  // the subgraph of G induced by B
  map<node,node> proj;   // proj : V(G) -->> V(G_B)

  forall(e, B)
    {
      node s = G.source(e);
      if (! proj.defined(s))
	{
	  node s_G = G_B.new_node();
	  G_B[s_G] = s;  // needed?
	  proj[s] = s_G;
	}
      node t = G.target(e);
      if (! proj.defined(t))
	{
	  node t_G = G_B.new_node();
	  G_B[t_G] = t;  // needed?
	  proj[t] = t_G;
	}
      G_B[G_B.new_edge(proj[s], proj[t])] = e;
    }
  
  // identify cycles in G_B:
  list<edge> CEL;  // cycle edge list, 
  // will contain one edge (in G_B) from each cycle after Is_Acyclic call

  if (! Is_Acyclic(G_B, CEL))
    {
      edge first_e;
      forall(first_e, CEL)
	{
 	  list<edge> C;  // a cycle is stored as the list of its edges in G
	  
	  // iterate over cycle following inedges:
	  node s = G_B.target(first_e);
	  node first_v = s;  // starting node;
	  edge e = NULL; //edge e = nil;
	  while (s != first_v || e == NULL)//nil)
	    {
	      e = G_B.in_edges(s).head();
	      s = G_B.source(e);  // move on
	      C.push(G_B[e]);  // store corresponding edge in G (!)
	      // the order in C follows the edge directions
	    }

	  k++;  CA.resize(k+1);
 	  CA[k] = C;
 	}
    }
  
  return CA;
}



edge predecessor_in_cycle(node y, list<edge>& C)
{
  // predecessor edge of y in C
  
  edge e;
  forall(e, C)
    {
      if (target(e) == y)
	return e;
    }

  return NULL; //nil;
}


void contract_cycle_edge(edge& e, double alpha_cost, double cheapest_cost, GRAPH<node,edge>& H, node& v_C, map<edge,edge>& edge_contr, node_set& VC)
{
  // note:  e \in E(H)

  node s = H.source(e);  // \in G(H)
  node t = H.target(e);  // \in G(H)

  if (VC.member(s) && (! VC.member(t)))  // e leaves C
    {
      edge e_tilde = H.new_edge(v_C, t);  // new edge
      c[e_tilde] = c[e];
      edge_contr[H[e]] = e_tilde;

      H[e_tilde] = H[e];  // update
      Phi[e_tilde] = H[e];

      H.del_edge(e);
    }

  if ((! VC.member(s)) && VC.member(t))  // e enters C
    {
      edge e_prime = H.new_edge(s, v_C);  // new edge
      c[e_prime] = c[e] - alpha_cost + cheapest_cost;
      edge_contr[H[e]] = e_prime;

      H[e_prime] = H[e];  // update
      Phi[e_prime] = H[e];

      H.del_edge(e);
    }

  if (VC.member(s) && VC.member(t))  // e is in C
    H.del_edge(e);
  
}



string contract_cycle_node(node& v, GRAPH<node,edge>& H, node& v_C, map<node,node>& node_contr, node_set& VC)
{
  string label;

  if (VC.member(v))  // v is on C
    {
      label += Mut[v];
      H.del_node(v);
      node_contr[H[v]] = v_C;
    }
  
  return label;
}



void contract_cycle(int i, array< list<edge> >& CA, int k, GRAPH<node,edge>& H, map<node,node>& node_contr, map<edge,edge>& edge_contr)
{
  node v;  edge e;

  list<edge> C = CA[k];  // cycle to be contracted (indexed by k)
  
  // node set in H indiced by cycle (in G):
  node_set VC(H);  // V(C) \subset H
  forall(e, C)
    {
      VC.insert(H.source(edge_contr[e]));
      VC.insert(H.target(edge_contr[e]));
    }
  
  // some cheapest edge on C
  cheapest_edge.resize(i+1);  cheapest_edge[i].resize(k+1);
  cheapest_edge[i][k] = C.contents(C.min(*compare_weights));
  
  // new node v_C for contracted cycle C
  node v_C = H.new_node();
  cycle_contr.resize(i+1);  cycle_contr[i].resize(k+1);
  cycle_contr[i][k] = v_C;  // in G[i+1]

  // contract edges in H
  list<edge> EL = H.all_edges(); 
  forall(e, EL)
    {
      alpha.resize(i+1);  alpha[i].resize(k+1);
      alpha[i][k][H[e]] = predecessor_in_cycle(target(H[e]), C);  // predecessor edge of target(e) in C
      contract_cycle_edge(e, c[alpha[i][k][H[e]]], c[cheapest_edge[i][k]], H, v_C, edge_contr, VC);
    }
  
  // contract nodes in H
  string new_label("(");
  list<node> NL = H.all_nodes();
  forall(v, NL)
    new_label += contract_cycle_node(v, H, v_C, node_contr, VC);
  Mut[v_C] = new_label + ")";

}



void contract_all_cycles(graph& G, int i, array< list<edge> > CL, GRAPH<node,edge>& H)
{
  // construct the contraction (G, c) --> (H, c)

  node v;  edge e;
  
  // {node, edge}_contr : (G, c) -->> (H, c')
  map<node,node> node_contr;  // node contraction map
  map<edge,edge> edge_contr;  // edge contraction map

  // initialize as identity map
  forall_nodes(v, G)
    {
      node_contr[v] = H.new_node();
      H[node_contr[v]] = v;
      Mut[node_contr[v]] = Mut[v];  // labels
    } 
  forall_edges(e, G)
    {
      edge_contr[e] = H.new_edge(node_contr[G.source(e)], node_contr[G.target(e)]);
      H[edge_contr[e]] = e;
      c[edge_contr[e]] = c[e];  // weights
    }

  // contract cycles in H
  for (int k=CL.low(); k<=CL.high(); k++)
    {
      contract_cycle(i, CL, k, H, node_contr, edge_contr);
      list<edge> empty_list;
    }
}



void reconstruct_branching(int i, GRAPH<node,edge>& H, list<edge>& B_H, list<edge>& B_G, array< list<edge> >& CA)
{
  //  reconstruct branching B_G from the contraction
  //  (G, B_G==B[i-1]) --> (H==G[i], B_H==B[i]) 
  //  CA is the array of cycles in G

  edge e;
  
  B_G.clear();
  
  node_set CCS(H);  // contracted cycles in H
  for (int k=CA.low(); k<=CA.high(); k++)
    CCS.insert(cycle_contr[i-1][k]);
  
  // E(B_G) := E(B_H) \ {e_prime=(z, v_C) | z \in V(B_H), C cycle in B_G}
  forall(e, B_H)
    if (! CCS.member(H.target(e))) 
      B_G.append(H[e]);
 
  // check all cycles in B_G
  for (int k=CA.high(); k>=CA.low(); k--)
    {
      list<edge> C = CA[k];  // cycle edges
      
      // look for e'==(z, v_C) in B_H
      edge e_prime = NULL; //nil;
      forall(e, B_H)
	if (H.target(e) == cycle_contr[i-1][k])  
	  e_prime = e;
      
      if (e_prime != NULL) //nil)  // e' == (z, v_C) in B_H
	{
	  B_G.append(Phi[e_prime]);
	  forall(e, C)
	    if (e != alpha[i-1][k][Phi[e_prime]])
	      B_G.append(e);
	}
      else  // no edge (z, v_C) in B_H
	{
	  forall(e, C)
	    if (e != cheapest_edge[i-1][k])
	      B_G.append(e);
	}
      
    }
}



list<edge> MAX_WEIGHT_BRANCHING(graph& G0, map<node,string>& node_label, edge_array<double>& weight)
{
  // Compute a maximum weight branching of the weighted digraph (G, weight)
  // The list of edges of the branching is returned.
  // The algorithm with all notation is taken from the book
  // Korte, Vygen: Combinatorial Optimization, Springer 2000, pp. 121-125

  node v;  edge e;

  // step 0: data structures:
  array< GRAPH<node,edge> > G(1);  // G[i+1] is subgraph of G[i]
  array< list<edge> > B(1);  // B[i] \subset E(G[i])
  array< array< list<edge> > > CAA(1);  // array of cycle arrays, each cycle is a subset of B[i]
  
  // step 1: initialization.
  CopyGraph(G[0], G0);  
  forall_edges(e, G[0])  // edge weights
    c[e] = weight[G[0][e]];

  forall_nodes(v, G[0])  // node labels
    Mut[v] = node_label[G[0][v]];

  // step 2: construct greedy subgraph B[i].
  B[0] = max_weight_subgraph_indeg_le_1(G[0]);

  // step 3: for all cycles in B[i].
  CAA[0] = all_cycles(G[0], B[0]);

  int i = 0;
  while (CAA[i].size() > 0)
    {
      // step 4: construct (G[i+1], c[i+1]) from (G[i], c[i]).
      G.resize(i+2);  B.resize(i+2);  CAA.resize(i+2);

      contract_all_cycles(G[i], i, CAA[i], G[i+1]);

      i++;
      B[i] = max_weight_subgraph_indeg_le_1(G[i]);  // (step 2)

      CAA[i] = all_cycles(G[i], B[i]);  // (step 3)
    }
  
  // step 5+6: reconstruct the maximum branching B
  while (i > 0)
    {
      reconstruct_branching(i, G[i], B[i], B[i-1], CAA[i-1]);
      i--;
    }

  list<edge> B0;
  forall(e, B[i]){
    B0.append(G[i][e]);
  }

  //clear the global static variables
  Mut.clear();
  c.clear();
  cycle_contr.clear();
  cheapest_edge.clear();
  Phi.clear();
  alpha.clear();

  return B0;
}
