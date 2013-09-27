#pragma once

#include <vector>
#include <string>

//Let deterministic HMM is given
//Let P(q,a) be a probability function where q - state, a - symbol
//Ô(q,a) - transition function
class D_HHM_Prob {
public:
	static double TermCondProb(int q1, const std::string &word, size_t wordsize, int q2);		// computes Prob(q1,word[1,wordsize],q2)
																								//Input parameters: q1, q2 -starting and terminal states; word - word, wordsize - length of word

	static double TermProb(std::string word, int q);											//calculates Prob(word,q)
	static double TransitionProb(int q1, int q2);												//calculates probability of transition from q1 to q2
	static std::vector<int> ConsistStates(int q);												//gets all states q' such that exisis a symbol a in Alp for that q'=Ô(q,a)
	static void Debug(void);
};
