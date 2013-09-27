#ifndef M_TrTree_H
#define M_TrTree_H

#include <vector>
#include <string>
#include "maindata.h"
#include "h_m_node.h"
#include "statedata.h"

//The tree is tree of prefixes of HH of length K and smaller

//A node corresponds to a word in Alp^k (k=1,..,K) 
//Descriptor of node t

class M_TrTree //node t
{
public: 
	int sign;					//symbol on edge leading to t
	int len;					// |t|
	M_TrTree* *Childs;			//list of prefix successors of t, Childs has size |Alp|, if there no node t.a, a in Alp, then Childs[a] = null
	int NStates;				//|AllState(w)|
	int* States;				//list AllState(w)

	static int NumTrNodes;		//number of nodes in the tree
	static M_TrTree* Root;		//root of the tree
	


public:
     M_TrTree();
    virtual ~M_TrTree();

	static void CreateTree(NodeAC* node, M_TrTree* LPnode, int len, string &word);  //Create MTrTree by depth-first traversal of AC trie
};


//Class of K-nodes (nodes of length K) knode
class KM_TrTree:
	public M_TrTree
{
public:

	int NumWLinks;			//Number of links leading from knode to right deep node r of OVGraph, where exists a word h such that 
							//prefix of h is knode and rpred(H)=r
	H_M_Node* *WLinks;		//links leading from knode to right deep node r of OVGraph, where exists a word h such that prefix of h
							//is knode and rpred(h)=r (notation of the set of such words is H~(K,r))
	int* NumWProbs;			// Let there exist a link from knode to r, the number of link is i. Then NumWrobs[i] is number of states
							//q in AllState(r) such that Prob(H~(K,r),q) !=0
	PriorList* *WProbs;		// WProbs[i] is the list of  non-zero probabilities Prob(H~(K,r),q), q in AllState(r) 
	
public:
    KM_TrTree();
    virtual ~KM_TrTree();
};

#endif
