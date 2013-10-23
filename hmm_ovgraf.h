#ifndef HMM_OVGRAF_H
#define HMM_OVGRAF_H
#include "h_m_node.h"
#include "statedata.h"

class hmm_ovgraf
{
public:
	static void CrConsistStatesMatrix(void);											//Set elements of ConsistStMatrix, see h_m_node.h;
	static void States_and_BackProbs(H_M_Node* node, H_M_Node* LPnode, string &back);	//for all w in OV(H) create descriptors of states form AllStates(w) and for all q' in AllStates(w) and q in AllStates(w) compute Prob(q',back,q) 
																						//Input parameters: node - current considered node, LPnode - lpred(node), back - Back(node)
	static void WordProbs(string &word);												//for all w in OV(H) , q' in Q and q in AllState(w) compute Prob(q',word,q)	
	static void FarPartCalc(H_M_Node* rnode);											//compute probabilities of far sets
	static void TransMatrixProduct(int n);												//compute Prob(V^n,q), q in Q
};

#endif
