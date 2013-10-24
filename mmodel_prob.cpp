#include "mmodel_prob.h"
#include "m_trtree.h"
#include "maindata.h"
#include <iostream>


int MModel_Prob::Power = 1; //Alp^order
double* MModel_Prob::IniProbs = nullptr;


MModel_Prob::MModel_Prob(void)
{
}

MModel_Prob::~MModel_Prob(void)
{
} 

////////////////////////////////////////////////
int MModel_Prob::NumPower(int num, int step){
	int i;
	int pow = 1;
	for(i = 0; i < step; i++){
		pow = pow*num;
	}
	return pow;
}

///////////1. functions to work with words//////////////////////

//for a word alpha gives number of the word in the lexicographic order.
int MModel_Prob::Code(std::string &alpha, int len){
	int i;
	int p;
	int cd = 0;

	if(len == 0){
		return -1;
	}
	
	p = 1;
	for(i = 0; i < len -1; i++){
		p = p*MainData::AlpSize;
	}
	
	for(i = 0; i < len; i++){
		cd = cd + p*(MainData::AToi(alpha.at(i)));
		p = p/MainData::AlpSize;
	}
	return cd; 
}

//gives words of length k with number cd in the lexicographic order
std::string MModel_Prob::deCode(int k, int cd){
	std::string alpha;
	if(cd == -1){
	   return alpha;
	}
	int p = 1;
	int s, i;
	for(i = 0; i < k -1; i++){
		p = p*MainData::AlpSize;
	}
	for(i = 0; i < k; i++){
		s = cd/p;
		cd = cd - s*p;
		p = p/MainData::AlpSize;
		alpha += MainData::IToa(s);
	}
	return alpha;
};




int MModel_Prob::PrefixN(int len, std::string &word){
	int cd = 0;
	int pow = NumPower(MainData::AlpSize,len -1);
	int i;
	for(i = 0; i < len; i++){
		cd += MainData::AToi(word.at(i))*pow;
		pow = pow/MainData::AlpSize;
	}
	return cd;
};



//gives suffix of length k
int MModel_Prob::SuffixN(int len, std::string &word, int wordlen){
	int cd = 0;
	int pow = NumPower(MainData::AlpSize,len -1);
	int i;
	for(i = 0; i < len; i++){
		cd += MainData::AToi(word.at(wordlen - len + i))*pow;
		pow = pow/MainData::AlpSize;
	}
	return cd;
};
///////////////////////////////////
 void MModel_Prob::CalcReachStates(std::string &word, int wordlen, vector<int> &states){
	if(wordlen > MainData::order){
		int cd = SuffixN(MainData::order,word,wordlen);
		states.push_back(cd);
	}
	else{
		int i;
		int d = MainData::order - wordlen;
		int cd1,cd2, cd;
		cd1 = 0;
		cd2 = NumPower(MainData::AlpSize,d) - 1;
		cd = Code(word, wordlen);
		int pow = NumPower(MainData::AlpSize,wordlen);
		for( i = cd1; i <= cd2; i++){
			int state = i*pow;
			state += cd;
			states.push_back(state);
		}
	}
	return;
};



//gives probability Prob_q2(word|q1)
double MModel_Prob::TermCondProb(int q1, std::string &word, int wordlen, int q2){
	if(wordlen != 0){
		int i;
		int cd = q1;
		int pow = Power/MainData::AlpSize;
		double p = 1;

		for(i = 0; i < wordlen; i++){
			p = p*MainData::MarkovProbs[MainData::AToi(word.at(i))][cd];
			int cd1 = cd/pow;
			cd = cd - cd1*pow;
			cd = cd*MainData::AlpSize;
			cd += MainData::AToi(word.at(i));
		}
		return p;
	}
	else{
		return 1;
	}
}

//gives Prob_q(word)
string temp;
double MModel_Prob::TermProb(std::string &word, int wordlen, int q){
	if( wordlen!= 0){
		double p = 0;
		if(wordlen < MainData::order){
			p = IniProbs[q];
		}
		else{
			int cd = PrefixN(MainData::order, word);
			if(temp.size() == 0){
				temp.resize(MainData::WordLen);
			}
			SubStr(word,MainData::order,wordlen - MainData::order,temp);

			p = IniProbs[cd]*TermCondProb(cd,temp,wordlen - MainData::order,q);
		}
		return p;
	}
	else{
		return 1;
	}
}


bool MModel_Prob::Check_States_Consistence(int q1, std::string &word, int wordlen, int q2){
	if(wordlen == 0){
		if(q1 == q2) return 1;
		else return 0;
	}
	int sufcd,cd,rest, pow, pow1;
	pow = MModel_Prob::NumPower(MainData::AlpSize, MainData::order - wordlen);
	pow1 = MModel_Prob::Power/pow;
	if(wordlen >= MainData::order){
		sufcd = MModel_Prob::SuffixN(MainData::order, word,wordlen);
		if(sufcd == q2) return 1;
		else return 0;
	}
	else{
		int cd1 = MModel_Prob::Code(word,wordlen);
		rest = q1 % pow;
		cd = rest*pow1 + cd1; 
		if(q2 == cd) return 1;
		else return 0;
	}
};


//gives list of states consistent with q
vector<int> MModel_Prob::ConsistStates(int q){
	int i;
	vector<int> vec;
	int pow = Power/MainData::AlpSize;
	for(i = 0; i < MainData::AlpSize; i++){
		int cd = q/MainData::AlpSize;
		cd += i*pow;
		vec.push_back(cd);
	}
	return vec;
} 

/////////////////////Initial distribution////////////////////////////

int** TrMatrix = NULL;

void MModel_Prob::CrTrMatrix(void){
	int i,j,cd;
	int pow = MModel_Prob::Power/MainData::AlpSize;
	TrMatrix = new int*[MModel_Prob::Power];
	for(i = 0; i < MModel_Prob::Power; i++){
		TrMatrix[i] = new int[MainData::AlpSize];
		for(j = 0; j < MainData::AlpSize; j++){
			cd = i/MainData::AlpSize;
			cd += j*pow;
			TrMatrix[i][j] = cd;
		}
	}
	return;
}



double VecProduct(void){
	int i,j,pos;
	int s;
	double* vec = new double[MModel_Prob::Power];
	int pow = MModel_Prob::Power/MainData::AlpSize;

	for(i = 0; i < MModel_Prob::Power; i++){
		vec[i] = 0;
		s = i/MainData::AlpSize;
		s = i - s*MainData::AlpSize;
		for(j = 0; j < MainData::AlpSize; j++){
			pos = TrMatrix[i][j];
			vec[i] += MModel_Prob::IniProbs[pos]*MainData::MarkovProbs[s][pos];
		}
	}
	double norm = 0;
	double p;
	for(i = 0; i < MModel_Prob::Power; i++){
		p = MModel_Prob::IniProbs[i] - vec[i];
		if(p < 0){
			p = -1*p;
		}
		norm += p;
	}

	for(i = 0; i < MModel_Prob::Power; i++){
		MModel_Prob::IniProbs[i] = vec[i];
	}
	delete[] vec;
	vec = nullptr;
	return norm;
}

void MModel_Prob::SetIniProbs(void){
	int i;
	for( i = 1; i < Power; i++){
		IniProbs[i] = 0;
	}
	IniProbs[0] = 1;

	double norm = 1;
	i = 0;
	
	while((i < 1000)&&(norm > 0.0000000001)){
		norm = VecProduct();
		i++;
	}
	
	/*
	for(i = 0; i < Power; i++){
		delete[] TrMatrix[i];
	}
	delete[] TrMatrix;
	*/
	return;
}

void MModel_Prob::TransStepProbCalc(double* &Mass){
int i,j,pos,s;
	if(MainData::order <= 5){
		double vec[NMaxStates];
	
		for(i = 0; i < MModel_Prob::Power; i++){
			vec[i] = 0;
			s = i/MainData::AlpSize;
			s = i - s*MainData::AlpSize;
			for(j = 0; j < MainData::AlpSize; j++){
				pos = TrMatrix[i][j];
				vec[i] += Mass[pos]*MainData::MarkovProbs[s][pos];
			}
		}

		for(i = 0; i < MModel_Prob::Power; i++){
			Mass[i] = vec[i];
		}
	}
	else{
		double* vec = new double[MModel_Prob::Power];
	
		for(i = 0; i < MModel_Prob::Power; i++){
			vec[i] = 0;
			s = i/MainData::AlpSize;
			s = i - s*MainData::AlpSize;
			for(j = 0; j < MainData::AlpSize; j++){
				pos = TrMatrix[i][j];
				vec[i] += Mass[pos]*MainData::MarkovProbs[s][pos];
			}
		}

		for(i = 0; i < MModel_Prob::Power; i++){
			Mass[i] = vec[i];
		}	
		delete[] vec;
		vec = nullptr;
	}
	return;
}

