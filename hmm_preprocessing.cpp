
#include "h_m_node.h"
#include "nodeov.h"
#include "nd_hhm_prob.h"
#include "d_hhm_prob.h"
#include "hmm_ovgraf.h"
#include "m_ovgraf.h"
#include "m_trtree.h"
#include "statedata.h"
#include "type_malloc.h"
using namespace std;


int H_M_Node::NumAllStates = 0;							//number of all states in HMM
double* H_M_Node::BnpProbs = nullptr;				    // The 2-dimensional matrix of size Q*NOccur
													 	// BnpProbs[q][p-1]= Prob(B(n-m,p),q)  

int** H_M_Node::ConsistStMatrix = nullptr;				//The array of lists. For all each state q  ConsistStMatrix[q] contains states q' (consistent with q) such that ï(q',q)>0 
int* H_M_Node::ConsistStNums = nullptr;					// The array of size NumAllStates
														//ConsistStNums[q] = (number of states consistent with q)
	
double** H_M_Node::TransProbMatrix = nullptr;		    //The array of lists.
														//for all q from Q and for all states q' consistent with q 
														//TransProbMatrix[q][q']=ï(q',q)
double* H_M_Node::TransStepProbMatrix = nullptr;	    //(needed only for HHM)At the stage n of the algorithms work,
														//TransProbMatrix[q][q'] is the probability of transition from q' to q during n steps.

double* H_M_Node::TransStepProbList = nullptr;			//(needed only for Markov)At the stage n of the algorithms work,
														//TransStepProbList[q][q'] = Prob(V^n,q).

					////AUXILLARY MATRECES//////
double** HDProbs = nullptr;								//Let w be processed node, n- current stage; depth - depth of the path leading to w
														//x_0,..x_depth - overlap prefixes of w, x_depth = w
														//HDPobs[k][p] = Prob(D(n-m+|x_k|,p+1,x_k),q), p= 0,...,p0-1, k = 0, ..,depth; q in AllState(w)  

double* PrevBnpProbs = nullptr;							//list to temporarly store data from BnpProbs

double* PrevTransStepProbMatrix = nullptr;				//list to temporarly store data from TransStepProbMatrix


double** OrderProbs = nullptr;							//For Markov model, let t be a prefix of motif words of length len, t1,...,t_len be prefixes of t; 
														//In function, CalEProbOne (see m_ovgraf)  OrderProbs[i][q] = Prob(V^{n-m}.t_i, q), q in AllState(t_i); 
														//In function CalFarProbs OrderProbs[i] = Prob(B.t_i,q+p_0 +j);



// computes Prob(q1,word[1,wordsize],q2)
//Input parameters: q1, q2 -starting and terminal states; word - word, wordsize - length of word

double TermCondProb(int q1, std::string &word,  int wordsize, int q2){
	double p = 0;
	if(MainData::order > 0){
		p = MModel_Prob::TermCondProb(q1,word,wordsize,q2);
	}
	if(MainData::order == -1){
        p = D_HHM_Prob::TermCondProb(q1,word, wordsize, q2);
	}
	if(MainData::order == -2){
        p = ND_HHM_Prob::TermCondProb(q1,word,wordsize,0,q2);
	}
	return p;
} 


//computes Prob(word[1,wordsize],q)
//Input parameters: q -starting and terminal states; word - word, wordsize - length of word

double TermProb(std::string word, int wordlen, int q){
	double p = 0;
	if(MainData::order > 0){
		p = MModel_Prob::TermProb(word,wordlen,q); 
	}
	if(MainData::order == -1){
		p = D_HHM_Prob::TermProb(word,q);
	}
	if(MainData::order == -2){
		p = ND_HHM_Prob::TermProb(word,q);
	}
	return p;
}


//calculates probability of transition from q1 to q2
double TransitionProb(int q1, int q2){
	double p = 0;
	if(MainData::order > 0){
		int s = q2/MainData::AlpSize;
		s = q2 - s*MainData::AlpSize;
		p = MainData::MarkovProbs[s][q1];
	}
	if(MainData::order == -1){
		p = D_HHM_Prob::TransitionProb(q1,q2);
	}
	if(MainData::order == -2){
		p = ND_HHM_Prob::TransitionProb(q1,q2);
	}
	return p;
}



//for a state q gets all states q' such that exisis a symbol 'a' for that Ï(q',a,q)>0
vector<int> ConsistStates(int q){
	vector<int> vec;
	if(MainData::order > 0){
	//	vec = MModel_Prob::ConsistStates(q);
	}
	if(MainData::order == -1){
		vec = D_HHM_Prob::ConsistStates(q);
	}
	if(MainData::order == -2){
		vec = ND_HHM_Prob::ConsistStates(q);
	}
	return vec;
}


//Initialization of parameters depending on model
void ModelInit(void){
	if(MainData::order < 0){							//IF given HMM
		hmm_ovgraf::CrConsistStatesMatrix();			//Create matrix  ConsistStMatrix,  For all each state q  ConsistStMatrix[q] contains states q' 
														//(consistent with q) such that ï(q',q)>0 
	}else{ //IF given Markov model
		MModel_Prob::Power = MModel_Prob::NumPower(MainData::AlpSize,MainData::order);  //For Markov model of order K, MModel_Prob::Power = |Alp|^K
		m_ovgraf::CrConsistStatesMatrix();
		MModel_Prob::CrTrMatrix();						// Create TransProbMatrix; for all q from Q and for all states q' consistent with q 
														//TransProbMatrix[q][q']=ï(q',q)
		if(MainData::MarkovType == 0){					//the model is stationary
			MModel_Prob::IniProbs = Malloc<double>(H_M_Node::NumAllStates); //MModel_Prob::IniProbs - list with stationary probabilities of states
			MModel_Prob::SetIniProbs();					//Initial probabilities are computed as eigen vector of TransProbMatrix
		}
	}
	return;
}

// INITIALIZATION OF DATA IN DESCRIPTORS OF STATES (see StateData.h, paper)

///////////////////////
//takes substring of string str; 
//Input parameters: str - string, start - starting position, nsymb - nummber of coping symbols, 
//result = str[start, start+nsymb]

void SubStr(string &str, int startpos, int nsymb, string &result){
	int pos,i;
	pos = startpos;
	for(i = 0; i < nsymb; i++){
		result[i] = str[pos];
		pos++;
	}
}

///////////////////////


//Computes deep probabilities by depth-first traversal of overlap graph, for each node w and for all right deep nodes r 
//s.t. exist h*(w,r) computes Prob(q',h*(w,r),q), q' in AllState(w), q in ReachState(q',Back(h*(w,r)))

//Input parameters: node - current processing node w; back - Back(h*(w,r)
void DeepProbs(H_M_Node* node, string &back){
	int i,j;

	//Deep probabilities
	if(node->ldeep == 1){	//if w is left deep node
		int k,l; 

		for(i = 0; i < node->NStates; i++){  //for all q in AllState(w)
            typedef PriorList* PriorListPtr;
            node->States[i]->DeepProbs = new PriorListPtr[node->NDLinks](); //initialization of data structures
			node->States[i]->NumDeepProbs = Malloc<int>(node->NDLinks);
			int state = node->States[i]->ID;								//state is q'

			for(j = 0; j < node->NDLinks; j++){								//for r s.t. exists h*(w,r) (deep link)
				int ncl = NLinkToClass[node->num][j];
				int clsize = CLasses[ncl-1].size();
				H_M_Node* DL = static_cast<H_M_Node*>(node->DeepLinks[j]);  //DL is r
					
				vector<PriorList> deepprobs;
				deepprobs.resize(DL->NStates);
				int flag = 1;
				for(l = 0; l <clsize; l++){									//for each h from h*(w,r)
					 int leaf = CLasses[ncl-1][l];							
					 int pos1 = PosNLeafWords[leaf] + node->len;
					 int backlen = MainData::WordLen - node->len;
			         SubStr(WLeafWords, pos1, backlen, back);				//take Back(h)

				     for(k = 0; k <  DL->NStates; k++){						// for each q in ReachState(q',Back(h*(w,r)))
						 int dstate = DL->States[k]->ID;					//q
						 if(MainData::order > 0){
							 flag = MModel_Prob::Check_States_Consistence(state, back, backlen, dstate); //flag == 1, if Prob(q',Back(h),q) >0, 0 - otherwise
						 }
						 if(flag == 1){
							deepprobs[k].pos = k;
							double p = TermCondProb(state, back, backlen, dstate);						//compute Prob(q',Back(h*(w,r)),q)
						    deepprobs[k].prob += p;
						 }
					}
				}


				//Recorde Prob(q',Back(h*(w,r)),q) to <w,q'>.DeepProbs
				for(k = 0; k <  DL->NStates; k++)
					if(deepprobs[k].prob > 0)
						node->States[i]->NumDeepProbs[j] ++;	
                        
				int s = node->States[i]->NumDeepProbs[j];
                node->States[i]->DeepProbs[j] = new PriorList[s]();
				s = 0;
				for(k = 0; k <  DL->NStates; k++){
					if(deepprobs[k].prob > 0){
						node->States[i]->DeepProbs[j][s].pos = deepprobs[k].pos;	
						node->States[i]->DeepProbs[j][s].prob = deepprobs[k].prob;
						s++;
					}
				}
			}
		}/*for(j = 0; j < node->NDLinks; j++)*/
	}

	///Recursion. Processing of left successors of w
   for( i = 0; i < node->NLChilds; i++){
		H_M_Node* LC = static_cast<H_M_Node*>(node->LChilds[i]);
		DeepProbs(LC, back);
	}
}



/////////////////////DEBUG!!!//////////////////////////////////
void DebPrint(H_M_Node* node, H_M_Node* LParent){
	int i,j,k;	
	/*
	if(node->num == 0){
		cout<<"Step "<<step<<"\n-------\n--------\n\n";
	}
	*/
	
	//int row = MainData::WordLen - node->len;

	//HMM_State* States = reinterpret_cast<HM_State*>(node->States);

	cout<<"Num of node: "<<node->num<<'\n';
	cout<<"Num of states: "<<node->NStates<<'\n';
	for(i = 0; i < node->NStates; i++){
		cout<<"Data related to the state "<<node->States[i]->ID<<'\n';
		cout<<"RP pos: " <<node->States[i]->RParentPos<<'\n';

		cout<<"Left Probs \n";
		for(j = 0; j < node->States[i]->NumBackProbs; j++){
			int pos = node->States[i]->BackProbs[j].pos;
			int lpstate = LParent->States[pos]->ID;
			cout<<"Prob("<<lpstate<<",Back,"<<node->States[i]->ID<<") = "<<std::setprecision(15)<<node->States[i]->BackProbs[j].prob<<'\n';
		}
		cout<<"\n\n";

	
		cout<<"Deep Probs \n";
		for(k = 0; k < node->NDLinks; k++){
		//	cout<<"Class: "<<NLinkToClass[node->num][k]<<'\n';
			for(j = 0; j < node->States[i]->NumDeepProbs[k]; j++){
				int pos = node->States[i]->DeepProbs[k][j].pos;
				H_M_Node* rnode = static_cast<H_M_Node*>(node->DeepLinks[k]);
				int	rstate = rnode->States[pos]->ID;
				cout<<"Prob("<<node->States[i]->ID<<",HBack,"<<rstate <<") = "<<std::setprecision(15)<<node->States[i]->DeepProbs[k][j].prob<<'\n';
			}
		cout<<"\n\n";
		}
		
		cout<<"\n\n";
/*		if(MainData::order > 0){
			cout<<"Word Probs \n";
			for(j = 0; j < States[i].NumWordProbs; j++){
				int startstate = States[i].WordProbs[j].pos;
				cout<<"Prob("<<startstate<<",H~(w),"<<node->States[i].ID <<") = "<<std::setprecision(15)<<node->States[i].WordProbs[j].prob<<'\n';
			}
			cout<<"\n\n";
		}
	}

	cout<<"ProbMark \n";
	for(i = 0; i< row; i++){ 
		for(j = 0; j < MainData::NOccur; j++){
			cout<<std::setprecision(15)<<node->ProbMark[j*row + i]<<'\t';
		}
			cout<<'\n';
	}
	cout<<"\n\n";

*/
	}
		cout<<"\n\n";
	for ( i = 0; i < node->NLChilds; i++){
		H_M_Node* G = static_cast<H_M_Node*>(node->LChilds[i]);
		DebPrint(G,node);
	}
	
}





// find position of a state q from  AllState(w) in list of descriptors of states of rpred(w) 
//by bottom-up traversal of overlap graph
void FindRPPos(H_M_Node* node){

	int i,j,k;

	for(i = 0; i < node->NRChilds; i++){
		H_M_Node* RC = static_cast<H_M_Node*>(node->RChilds[i]);
		for(j = 0; j < node->NStates; j++){
			for(k = 0; k < RC->NStates; k++){
				if(node->States[j]->ID == RC->States[k]->ID){
					RC->States[k]->RParentPos = j;
				}
			}
		}

		FindRPPos(RC);
	} 
	return;
}


//Initialization of ProbMark and Firsttemp
void IniProbs(H_M_Node* node){

	
	int i,row;
	
	row = MainData::WordLen - node->len;  

	for(i = 0; i < node->NStates; i++){
		node->States[i]->ProbMark = Malloc<double>(MainData::NOccur*row);
		node->States[i]->FirstTemp = Malloc<double>(MainData::NOccur);
	}	

	for( i = 0; i < node->NLChilds; i++){
		H_M_Node* LC = static_cast<H_M_Node*>(node->LChilds[i]);
		IniProbs(LC);
	}
	
	return;
}


//Set parameters in descriptors of pairs <w,q>, w in OV(HH), q in AllState(w)
void SetStatesData(void){
		H_M_Node* root = static_cast<H_M_Node*>(NodeOv::Root); //root
		string back;
		back.resize(MainData::WordLen);							//Back(w)
		string word;											//word w	
		word.resize(MainData::WordLen);
		if(MainData::order < 0)
			hmm_ovgraf::States_and_BackProbs(root,root,back);	//Create list AllState(w), compute back  probabilities 
		else 
			m_ovgraf::States_and_BackProbs(root,root,back,word);
		


		///free///
		free(PosNNodeBacks);
		WNodeBacks.clear();
		///////
		
		DeepProbs(root,back);									//computation of deep probabilities
		
		///free//////
		CLasses.clear();
		delete[] NLinkToClass;
		///////

		IniProbs(root);											//Initialization of ProbMark and Firsttemp
		///////////
		if(MainData::order < 0)
			hmm_ovgraf::WordProbs(word);						//compute word probabilities	
		else
			m_ovgraf::WordProbs(word);
		////free//////////
		free(Nodes);
		free(RLeafPreds);
		free(PosNLeafWords);
		WLeafWords.clear();
		/////////////////////////////


		FindRPPos(root);										 //find position of a state of node w in list of descriptors of states of rpred(w) 
		
		
}

/////////////INITIALIZATION OF STATIC MATRICES (see h_m_node.h)//////////////////////////////
void IniStaticMas(void){
	int i;

	if(MainData::order > 0){
		H_M_Node::TransStepProbList = Malloc<double>(H_M_Node::NumAllStates);
		for(i = 0; i < H_M_Node::NumAllStates; i++){
			H_M_Node::TransStepProbList[i] = MModel_Prob::IniProbs[i];
		}

		OrderProbs = Malloc<double*>(MainData::order);
		for(i = 0; i < MainData::order; i++){
			int pow = MModel_Prob::NumPower(MainData::AlpSize,MainData::order -i-1);
			OrderProbs[i] = Malloc<double>(pow*MainData::NOccur);
		} 
	}
		
	else{
		H_M_Node::TransStepProbMatrix =  Malloc<double>(H_M_Node::NumAllStates* H_M_Node::NumAllStates);
		PrevTransStepProbMatrix =  Malloc<double>(H_M_Node::NumAllStates* H_M_Node::NumAllStates);
	}
	


	H_M_Node::BnpProbs=  Malloc<double>(H_M_Node::NumAllStates*MainData::NOccur);
	PrevBnpProbs=  Malloc<double>(H_M_Node::NumAllStates*MainData::NOccur);



	HDProbs = Malloc<double*>(MaxDepth);
	for(i = 0; i < MaxDepth; i++)
		HDProbs[i] = Malloc<double>(H_M_Node::NumAllStates* MainData::NOccur);
	

	return;
}






////////////////////////////////////

///////////PREPROCESSING/////////////////////////////
//preprocessing

void H_M_Node::Preprocessing(void){
	
	ModelInit();			// initialisation of model parameters

	SetStatesData();		//initialization of data in descriptors of pairs <w,q>, w in OV(HH), q in AllState(w)
	IniStaticMas();			//Initialisation of static matrices
	
	return;
}

///////////////////////////////////////////////////

//Free data structures
void Del_HMM_Mas(void){
	int i;
	int j;		
	free(H_M_Node::BnpProbs);
	free(PrevBnpProbs);
	
	for(i = 0; i < H_M_Node::NumAllStates; i++){
		free(H_M_Node::ConsistStMatrix[i]);
		free(H_M_Node::TransProbMatrix[i]);
	}
	free(H_M_Node::ConsistStMatrix);
	free(H_M_Node::TransProbMatrix);
	free(H_M_Node::TransStepProbMatrix);
	free(PrevTransStepProbMatrix);
	free(H_M_Node::ConsistStNums);

	H_M_Node::ConsistStMatrix  = nullptr;
	H_M_Node::TransProbMatrix = nullptr;
	H_M_Node::TransStepProbMatrix = nullptr;
	PrevTransStepProbMatrix = nullptr;
	H_M_Node::ConsistStNums = nullptr;

	for(i = 0; i < MaxDepth; i++)
		free(HDProbs[i]);	
	free(HDProbs);

	if(MainData::order > 0){
		free(MModel_Prob::IniProbs);
		free(H_M_Node::TransStepProbList);
		for(i = 0; i < MainData::order; i++)
			free(OrderProbs[i]); 
		free(OrderProbs);
	}
	
	/*
	if(MainData::order == -2){
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
	if(MainData::order == -1){
		for(i = 0; i < H_M_Node::NumAllStates; i++){	
			delete[] MainData::D_HHMProbs[i];
			delete[] MainData::D_HHMTrans[i];
		}
		delete[] MainData::D_HHMProbs;
		delete[] MainData::D_HHMTrans;
	}
	*/
	return;
}


//Detete overlap graph
void H_M_Node::ClearData(void){
	Del_HMM_Mas();
	if(MainData::order > 0)
		delete M_TrTree::Root;
	  NodeOv::DeleteGraf(NodeOv::Root);
	return;
}


////////////////////////////////////////////






