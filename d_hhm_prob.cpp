#include "d_hhm_prob.h"
#include "h_m_node.h"
#include "maindata.h"


//calculates Prob_q1(word|q2)
double D_HHM_Prob::TermCondProb(int q1, const std::string &word, size_t wordsize, int q2){
	if(word.length() == 0){
		if(q1 == q2){
			return 1;
		}
		else{
			return 0;
		}
	}

	int q;
	double p = 1;
	size_t i;
	q = q1;
	size_t s = wordsize;
	for(i = 0; i < s; i++){
		int pos = MainData::AToi(word.at(i));
		p = p*MainData::D_HHMProbs[q][pos];
		q = MainData::D_HHMTrans[q][pos];
		if(q == -1){
			return 0;
		}
	}
	if(q == q2){
		return p;
	}else{
		return 0;
	}
};

//calculates Prob_q(word)
double D_HHM_Prob::TermProb(std::string word, int q){
	size_t s = word.size();
	double p = TermCondProb(0,word,s,q);
	return p;
};


//calculates probability of transition from q1 to q2
double D_HHM_Prob::TransitionProb(int q1, int q2){
	double p = 0;
	int i;
	for(i = 0; i < MainData::AlpSize; i++){
		int state = MainData::D_HHMTrans[q1][i];
		if(state == q2){
			p += MainData::D_HHMProbs[q1][i];
		}
	}
	return p;
};

//gets all states q' such that exisis a symbol a for that q'=Ô(q,a)
vector<int> D_HHM_Prob::ConsistStates(int q){
	vector<int> vec;
	int i,j;
	for(i = 0; i < H_M_Node::NumAllStates; i++){
		for(j = 0; j < MainData::AlpSize; j++){
			if(MainData::D_HHMTrans[i][j] == q){
				vec.push_back(i);
				break;
			}	
		}
	}
	return vec;
};
