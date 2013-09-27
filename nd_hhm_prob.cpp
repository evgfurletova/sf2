#include "nd_hhm_prob.h"
#include "h_m_node.h"
#include "maindata.h"



double ND_HHM_Prob::TermCondProb(int q1,const std::string & word, size_t wordsize, size_t from, int q2){
    if(wordsize-from == 0){
		if(q1 == q2){
			return 1;
		}
		else{
			return 0;
		}
	}
	
	double prob = 0;
	size_t i;
    int pos = MainData::AToi(word.at(from));
	if(wordsize == 1){
		prob = MainData::ND_HHMProbs[q1][q2][pos];
	}
	else{
		size_t s = MainData::ND_HHMTrans[q1][pos].size();
//		word.erase(0,1);
		for(i = 0; i < s; i++){
			int state = MainData::ND_HHMTrans[q1][pos][i];
			double p = MainData::ND_HHMProbs[q1][state][pos];
            prob += TermCondProb(state,word,wordsize,from+1,q2)*p;
		}
	}
	return prob;
};

double ND_HHM_Prob::TermProb(std::string word, int q){
	size_t s = word.size();
    double p = TermCondProb(0,word,s,0,q);
	return p;
};


double ND_HHM_Prob::TransitionProb(int q1, int q2){
	double p = 0;
	int i;
	for(i = 0; i < MainData::AlpSize; i++){
		p += MainData::ND_HHMProbs[q1][q2][i];
	}
	return p;
};

vector<int> ND_HHM_Prob::ConsistStates(int q){
	vector<int> vec;
	int i,j;
	for(i = 0; i < H_M_Node::NumAllStates; i++){
		for(j = 0; j < MainData::AlpSize; j++){
			if(MainData::ND_HHMProbs[i][q][j] != 0){
				vec.push_back(i);
				break;
			}
		}
	}
	return vec;
};

