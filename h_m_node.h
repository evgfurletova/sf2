
#ifndef H_M_DATA_H
#define H_M_DATA_H

#include "nodeov.h"
#include "statedata.h"
#include "d_hhm_prob.h"
#include "nd_hhm_prob.h"
#include "mmodel_prob.h"

		/*				NOTATIONS
		Len p0 - be desired number of occurences; n0- text length; 
		n - current stage of the algorithm; m - pattern length; Alp - alphabet; HH - pattern; Q - set of HMM states
		п(q',q) - probability from state q' to pass to state q
		*/

////Descriptor of node w of OvGraf in case of  HMM or Markov model 

class H_M_Node:
public NodeOv
{
public:
	///////////////DESCRIPTOR OF A NODE OF OVGRAPH///////////////////////////////////////////////
	int NStates;						//Number of States
	H_M_State* *States;					//list of descriptors of states

	/////////////// SAME PARAMETERS  FOR ALL NODES OF OVGRAPH//////////////////
	
	static int NumAllStates;				//number of all states in HMM
	static double* BnpProbs;				// The 2-dimensional matrix of size Q*p0
									     	// BnpProbs[q][p]= Prob(B(n-m,p),q), q in Q, p = 0,..,p0-1  

	static int** ConsistStMatrix;			//The array of lists. For all each state q  ConsistStMatrix[q] contains states q' (consistent with q) such that п(q',q)>0 
	static int* ConsistStNums;				// The array of size NumAllStates
											//ConsistStNums[q] = (number of states consistent with q)
	
	static double** TransProbMatrix;		//The array of lists.
											//for all q from Q and for all states q' consistent with q 
											//TransProbMatrix[q][q']=п(q',q)

	static double* TransStepProbMatrix;	    //(needed only for HHM)At the stage n of the algorithms work,
											//TransProbMatrix[q][q'] is the probability of transition from q' to q during n steps, i.e Prob(q',V^n,q).
	
	static double* TransStepProbList;		//(needed only for Markov)At the stage n of the algorithms work,
											//TransStepProbList[q][q'] = Prob(V^n,q).
	



////////METHODS////////////////
	H_M_Node();
	~H_M_Node();

	static void ProbCalc(void);					//main function for P-value computation
	static void Preprocessing(void);			//Preprocessing
	static void ClearData(void);				//Clear data structures
};

											//AUXILARY DATA

extern double** HDProbs;						//Let w be processed node, n- current stage; depth - depth of the path leading to w
												//x_0,..x_depth - overlap prefixes of w, x_depth = w
												//HDPobs[k][p] = Prob(D(n-m+|x_k|,p+1,x_k),q), p= 0,...,p0-1, k = 0, ..,depth; q in AllState(w)  

extern double* PrevBnpProbs;					//list to temporarly store data from BnpProbs
extern double* PrevTransStepProbMatrix;			//list to temporarly store data from TransStepProbMatrix
extern double** OrderProbs;						//For Markov model, let t be a prefix of motif words of length len, t1,...,t_len be prefixes of t; 
												//In function, CalEProbOne (see m_ovgraf)  OrderProbs[i][q] = Prob(V^{n-m}.t_i, q), q in AllState(t_i); In function CalFarProbs OrderProbs[i] = Prob(B.t_i,q+p_0 +j);

											
											//FUNCTIONS

extern void SubStr(string &str, int startpos, int nsymb, string &result); //takes substring of string str; 
																		  //Input parameters: str - string, start - starting position, nsymb - nummber of coping symbols, 
																		  //result = str[start, start+nsymb]
extern double TransitionProb(int q1, int q2);							  //computes п(q1,q2)

extern double TermCondProb(int q1, std::string &word,  int wordsize, int q2);  // computes Prob(q1,word[1,wordsize],q2)
																				//Input parameters: q1, q2 -starting and terminal states; word - word, wordsize - length of word

extern double TermProb(std::string word, int wordlen, int q);					//computes Prob(word[1,wordsize],q)
																				//Input parameters: q -starting and terminal states; word - word, wordsize - length of word

extern vector<int> ConsistStates(int q);										//list of states q' such that п(q',q) > 0		



#endif // H_M_DATA_H

