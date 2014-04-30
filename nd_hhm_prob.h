#ifndef ND_HHM_PROB_H
#define ND_HHM_PROB_H


#include <vector>
#include <string>

//Let (q,a,q1) is a probabity beeng in state q to generate a symbol a and go to a state q1
class ND_HHM_Prob
{
public:
	static double HMM_Estimate_Param(int Nseqs, char** seq, double*** HMMProbs, double*** NewHMMProbs);
	static double Find_Best_Model(int NSeqs, char** seq, double*** &BestHMMProbs);
    static double TermCondProb(int q1,const std::string & word, size_t wordsize, size_t from, int q2);	  //calculates Prob_q1(word|q2)
	static double TermProb(std::string word, int q);			  //calculates Prob_q(word
	static void   AllTermProbs(int q0, double* PredProbs, int* PredStates, int npredstates, double* NewProbs, int* NewStates, const std::string & word, size_t wordsize, size_t from);
	static double TransitionProb(int q1, int q2);				  //calculates probability of transition from q1 to q2
	static std::vector<int> ConsistStates(int q);					  //gets all states q' such that exisis a symbol a for that q'=Ф(q,a)

	static void GenRandHMM(double d, double*** HMMProbs);
	static void GenRandHMM1(double d, double*** HMMProbs);
	static void GenRandHMM2(double d, double*** HMMProbs);
	ND_HHM_Prob(void);
	~ND_HHM_Prob(void);

};
#endif
