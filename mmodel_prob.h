#ifndef MMODEL_PROB_H
#define MMODEL_PROB_H
#include <string>
#include <vector>
#include "statedata.h"

class MModel_Prob
{
public:
	static int Power; //Alp^order
	static double* IniProbs; //initial probabilities distribution on the states
public:
	static int NumPower(int num, int step);					//compute num^step
	static int Code(std::string &alpha, int len);           //give substring alpha[1,len]
	static std::string deCode(int k, int cd);				//give word of length k, having code cd (number of the word in prefix order)
	static int PrefixN(int len, std::string &word);			//gives number 
	static int SuffixN(int len, std::string &word, int wordlen);
	static void CalcReachStates(std::string &word, int wordlen, std::vector<int> &states);
	static double TermCondProb(int q1, std::string &word, int wordlen, int q2);
	static double TermProb(std::string &word, int wordlen, int q);
	static std::vector<int> ConsistStates(int q);
	static void SetIniProbs(void);
	static void	TransStepProbCalc(double* &Mass);
	static void CrTrMatrix(void);
	static bool Check_States_Consistence(int q1, std::string &word, int wordlen, int q2);

	MModel_Prob(void);
	~MModel_Prob(void);
};
#endif
