#ifndef M_OVGRAF_H
#define M_OVGRAF_H_H
#include "m_trtree.h"

class m_ovgraf
{
public:
	static void CrConsistStatesMatrix(void);															//Set elements of ConsistStMatrix, see h_m_node.h;
	static void States_and_BackProbs(H_M_Node* node, H_M_Node* LPnode, string &back, string &word);		//for all w in OV(H) create descriptors of states form AllStates(w) and for all q' in AllStates(w) and q in AllStates(w) compute Prob(q',back,q) 
																										//Input parameters: node - current considered node, LPnode - lpred(node), back - Back(node)
	static void WordProbs(string &word);																//compute word probabilitie
	static void CalEProbOne(M_TrTree* tnode);															//compute Prob(E(n,p,w),q) by depth-first traversal of M_TrTree, w in OV(H), q in AllState(w) 
	static void CalFarProbs(M_TrTree* tnode);															//compute Prob(F(n,p,w),q) by depth-first traversal of M_TrTree, w in OV(H), q in AllState(w) 
};

#endif