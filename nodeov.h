#ifndef NODEOV_H
#define NODEOV_H

#include "maindata.h"

using namespace std;

/*
                                Classes of descriptors of nodes of Overlap graph
The overlap graph OvGraph of a given pattern H is an oriented graph with the set of nodes
OV(HH) and a set of edges  that is a union of three subsets, left, right and deep edges, that are defined as follows:
  - A left edge links node w to node t iff w = lpred(t);
  - A right edge links node w to node t iff t = rpred(w).
  - A deep edge links node w to node t iff exists a non-empty class H*(w, t),

 where OV(HH) - set of overlap nodes; lpred(t) (rpred(t)) - maximal prefix (suffix) of t from OV(HH);
class H*(w,t) - subset of pattern words h, where lpred(h) = w and rpred(H) = t 

Node w is left (right) deep, iff exists a word h from HH, where lpred(h) = w (rpred(h) = w) 
*/

//Descriptor of node w of OvGraf
class NodeOv {
public:
	//PARAMETERS OF DESCRIPTOR OF w
	int num;			   //number of the node w
	int len;			   //|w|
	bool rdeep;			   //1, if the node is right deep; 0 -otherwise
	bool ldeep;		       //1, if the node is left deep; 0 -otherwise
	
	bool rootchild;			//1, lpred(w)= root; 0 - otherwise
	NodeOv* RParent;		// link to the maximal suffix that is an overlap
	
	int NLChilds;	        // number of left childs
	int NRChilds;		    // Number of right childs	
	int NDLinks;		    // for left deep node number of deep edges leading from the node
	NodeOv* *LChilds;       // List of left links 
    NodeOv* *RChilds;		// List of right links
	NodeOv* *DeepLinks; 	// List of deep links


	///STATIС PARAMETERS (same for all nodes of the graph)
	static NodeOv* Root;			//root of the graph
    static int NumOVNodes;		    // number of nodes in the overlap graph
    static int NumRDNodes;			// number of right deep nodes in the overlap graph
    static int NumLDNodes;			// number of left deep nodes in the overlap graph
    static int NClasses;		    //number of classes

	///METHODS
	NodeOv();						//constructor
	virtual ~NodeOv();				//destructor
	static int CreateGraf(void);	//creating of the overlap graph
	static NodeOv* NewNode();		//creating of a new node of the graph
	static void DeleteGraf(NodeOv* node); //deleting of the graph; here, node - current deleted node
};		


							//AUXILIARY PARAMETERS
//1. PARAMETERS COMMON FOR ALL MODELS

extern NodeOv** Nodes;		//List of nodes of overalp graph
extern int *LLeafPreds;		// list of lpred(h), h in HH
extern int *RLeafPreds;		// list of rpred(h), h in HH

extern int MaxDepth;		//maximal number of nodes on the path leading from root to a terminal node


extern std::vector<int> *NLinkToClass; //during constructing of minimal overlap graph, NLinkToClass(w) is list of classes corresponding to the deep links DeepLinks


//2. PARAMETERS NEEDED FOR HMM AND MARKOV MODELS ONLY
extern vector< vector<int> > CLasses;  //list of overlap classes h*, Classes[h*] is list of numbers of words that are in the h*

extern int* PosNNodeBacks;		//for each node w from OV(HH), PosNNodeBacks is position of Back(w) (w = x.Back(w), x=lpred(w)) in  WNodeBacks
extern std::string WNodeBacks;  //sequence of Back(w), w in OV(HH)

extern int* PosNLeafWords;		//for each word h from HH, PosNLeafWords[h] is position of entry of h in  WLeafWords
extern std::string WLeafWords;  //sequence of entrances of h, w in HH


#endif
