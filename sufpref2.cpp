#include "maindata.h"
#include "nodeac.h"
#include "nodeov.h"
#include "nodebern.h"
#include "h_m_node.h"
#include "m_trtree.h"
#include "mmodel_prob.h"
#include <list>

#if !defined(WIN32) && !defined(_WIN32)
// POSIX-functions for time measure
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif


						//NOTATIONS
		//Len p0 - be desired number of occurences; n0- text length; Alp - alphabet
		//m - motif length; Alp - alphabet; HH - motif
		//definition of text sets and sets of states see in our paper xxx;
int main(int argc, char* argv[])
{
	////////////Preliminary actions//////////////
	int i;
	time_t time0 = time(NULL);
	int ExitFlag = 0;
    int Error = MainData::ComLineParse(argc, argv);  //parsing of command line (see maindata.h (cpp))
    if(Error > 0)
        ExitFlag = 1;

	if(ExitFlag == 0){
		Error = MainData::GetInput();				//creating of data structures to storing: 1) alphabet; 2) model parameters; 3) motif (creating of nodes AC trie and prefix links) (see maindata.h (cpp))
		if(Error > 0)								//considering of case p0 == 1 and n0 == m, in the case PValue = Prob(HH) 
			ExitFlag = 1;
	}
	if(ExitFlag == 0){
		Error = MainData::Check_Out_of_Range();   //Check corectness of input parameters (see maindata.h (cpp))
		if(Error > 0){
			ExitFlag = 1;
		}
	}

	//open of out file
	ofstream test;
	test.open("Test.txt", ios::app);
	ofstream out;
	if(MainData::Out_mode ==0){
		out.open(MainData::OutName.c_str(), ios::app);
	}
	else{
		out.open(MainData::OutName.c_str());
	}

	if((MainData::NOccur == 0)&&(ExitFlag == 0)){	//if p0 == 0 then Pvalue = 1 
		MainData::Pvalue =1;
		ExitFlag = 2;
	}
	
	if((MainData::TLen < MainData::WordLen)&&(ExitFlag == 0)){ //if n0 < m  then Pvalue = 0
		MainData::Pvalue = 0;
		ExitFlag = 2;
	}
	if((MainData::TLen == MainData::WordLen)&&(MainData::NOccur > 1)&&(ExitFlag == 0)){ //if n0 == m and p0 > 1 then Pvalue = 0
		MainData::Pvalue = 0;
		ExitFlag = 2;
	};

	if((MainData::AlpSize == 1)&&(ExitFlag == 0)){			//if |Alp| == 1 
		if(MainData::TLen >= MainData::WordLen + MainData::NOccur - 1){
			MainData::Pvalue = 1;
			ExitFlag = 2;
		}
		else{
			MainData::Pvalue = 0;
			ExitFlag = 2;
		}
	}

///////// 

		MainData::ErrorDetect(Error);		//error checking (see maindata.h (cpp))

/////////    
	/////////////Main computations//////////
	if(ExitFlag == 0){
		NodeAC::CreateTrie(); //Create suffix links of nodes of AhoCoracis trie (see nodeac.h, ac_trie.cpp)
		
		if(MainData::Out_mode == 3){
			MainData::ResWords.open(MainData::OutWordsName.c_str()); 
			NodeAC::PrintWords(&MainData::ResWords,NodeAC::ACRoot,"");
		    MainData::ResWords.close();
		}


	   if(MainData::order > 0){ //if Markov model of order K is given
		   MModel_Prob::Power = MModel_Prob::NumPower(MainData::AlpSize, MainData::order);
		   H_M_Node::NumAllStates = MModel_Prob::Power;
			M_TrTree::Root = new M_TrTree;
			std::string word; 
			word.resize(MainData::order);
			M_TrTree::CreateTree(NodeAC::ACRoot,M_TrTree::Root,0,word); // creating of M_TrTree (see m_trtree.h, m_ovgraf.cpp)
	   }
	    NodeOv::CreateGraf(); //Create  Overlap graph (see nodeov.h, ov_graf.cpp)
		
			
		if(MainData::order == 0)
			NodeBern::ProbCalc(); //P-value computation for Bernoulli (see nodebern.h, bprob_graf.cpp)
		else 
			H_M_Node::ProbCalc();  //P-value computation for other models (see h_m_node.h, h_m_ovgraf.cpp)

	
		if(MainData::Out_mode > 1){			// if type of out >1 then printing of report about input and output parameters  (see maindata.h (cpp))
			 MainData::PrintMain(&out,0);
		 }
	}


	//Delete data structures

	if(MainData::mode == 1){
		delete[] MainData::RandPatProbs;
	}
	if((MainData::order > 0)&&(MainData::CrDistribFlag == 1)){
		for(i = 0; i < MainData::AlpSize; i++){
			delete[] MainData::MarkovProbs[i];
		}
		delete[] MainData::MarkovProbs;
	}
	if((MainData::Out_mode <2)&&(ExitFlag != 1)){
		out<<std::setprecision(15)<<MainData::Pvalue<<'\n';
	}
	
	
		int j;
	if((MainData::order == -2)&&(MainData::CrDistribFlag == 1)){
		for(i = 0; i < H_M_Node::NumAllStates; i++){
			for(j = 0; j < H_M_Node::NumAllStates; j++){
				delete[] MainData::ND_HHMProbs[i][j];
			}
			delete[] MainData::ND_HHMProbs[i];
			MainData::ND_HHMTrans[i]->clear();
			delete[] MainData::ND_HHMTrans[i];
		}
		delete[] MainData::ND_HHMProbs;
		delete[] MainData::ND_HHMTrans;
	}
	if((MainData::order == -1)&&(MainData::CrDistribFlag == 1)){
		for(i = 0; i < H_M_Node::NumAllStates; i++){	
			delete[] MainData::D_HHMProbs[i];
			delete[] MainData::D_HHMTrans[i];
		}
		delete[] MainData::D_HHMProbs;
		delete[] MainData::D_HHMTrans;
	}
	
	if(ExitFlag > 0){
		if((MainData::Out_mode > 1)&&(ExitFlag == 2)){
			 MainData::PrintMain(&out,2);
		}
		return 1;
	}

	if((MainData::mode == 2)||(MainData::mode == 3)){
		for(i = 0; i < MainData::WordLen; i++){
			delete[] MainData::PssmMas[i];
		}
		delete[] MainData::PssmMas;
		MainData::PssmMas = NULL;
	}
	if((MainData::mode == 4)&(MainData::ConstPositions != NULL)){
		delete[] MainData::ConstPositions;
	}

	//time0 = time(NULL) - time0;
	//out<<'\t'<<time0<<'\n';
	
	out.close();


#if !defined(WIN32) && !defined(_WIN32)
	// Print resources used
	int who = RUSAGE_SELF;
	struct rusage usage;
	int ret;
	ret = getrusage(who, &usage);
	if (ret==0) {
	    std::cout << "Used resources: \n";
	    struct timeval utime = usage.ru_utime;
	    struct timeval stime = usage.ru_stime;
	    long maxrss = usage.ru_maxrss;
	    long u_millsecs = (utime.tv_sec * 1000000 + utime.tv_usec) / 1000;
	    long s_millsecs = (stime.tv_sec * 1000000 + stime.tv_usec) / 1000;
	    long total_millsecs = u_millsecs + s_millsecs;
	    std::cout << "user time: " << u_millsecs << "ms; system time: " << s_millsecs << "ms" << std::endl;
	    std::cout << "total time: " << total_millsecs << "ms" << std::endl;
	    long pagesize = sysconf(_SC_PAGE_SIZE);
	    long memused = maxrss * pagesize;
	    std::cout << "Used " << maxrss << " KB " << std::endl;
	}
	else {
	    std::cout << "Can't use getrusage!\n";
	}
#endif



	return 0;


}




