#ifndef ND_HHM_PROB_H
#define ND_HHM_PROB_H


#include <vector>
#include <string>

//Let (q,a,q1) is a probabity beeng in state q to generate a symbol a and go to a state q1
class ND_HHM_Prob
{
public:

    static double TermCondProb(int q1,const std::string & word, size_t wordsize, size_t from, int q2);	  //calculates Prob_q1(word|q2)
	static double TermProb(std::string word, int q);			  //calculates Prob_q(word

	static double TransitionProb(int q1, int q2);				  //calculates probability of transition from q1 to q2
	static std::vector<int> ConsistStates(int q);					  //gets all states q' such that exisis a symbol a for that q'=Ф(q,a)

	ND_HHM_Prob(void);
	~ND_HHM_Prob(void);

};
#endif
