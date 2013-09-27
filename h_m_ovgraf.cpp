
#include "nd_hhm_prob.h"
#include "d_hhm_prob.h"
#include "maindata.h"
#include "h_m_node.h"
#include "m_ovgraf.h"
#include "hmm_ovgraf.h"
#include "m_trtree.h"




////////////////////////////////////////////

						//NOTATIONS
		//Len p0 - be desired number of occurences; n0- text length; 
		//n - current stage of the algorithm; m - pattern length; Alp - alphabet; HH - pattern
		//definition of text sets and sets of states see in our paper xxx; Q is the set af HMM states


//calculate initial probabilities Prob(R(m,p,H(w))) and add it to w.ProbMark, 
//by bottom-up traversal of overlap graph following by right links

//Input parameters: node - processing node w
void Stepm(H_M_Node* node){
	int i,k;
    int pos_in_matrix = MainData::WordLen - node->len - 1;

	if(MainData::order < 0){		//is the model is HMM
        if(node->rdeep == 1){		//is w is right deep node
			for(k = 0; k < node->NStates; k++){	  //for all states q in AllStates(w)
				HMM_State* State = static_cast<HMM_State*>(node->States[k]);
				State->ProbMark[pos_in_matrix] = 0;
				
				if(State->NumWordProbs != 0){	
					if(State->WordProbs[0].pos == 0){
						State->ProbMark[pos_in_matrix] = State->WordProbs[0].prob;  //record Prob(H~(w),q) to <w,q>.ProbMark
					}
				}
			}
		}
	}

	//compute Prob(H(w),q) = Sum_{x in OV(HH), rpred(w)=x}Prob(H(x)) and record it to <r,w>.ProbMark 
	for ( i = 0; i < node->NRChilds; i++){	
		H_M_Node* RC = static_cast<H_M_Node*>(node->RChilds[i]);
		Stepm(RC);		//Run the procedure for right successor x of  w
        int pos1 = MainData::WordLen - RC->len - 1;
		for(k = 0; k < RC->NStates; k++){
			int state = RC->States[k]->RParentPos;
			node->States[state]->ProbMark[pos_in_matrix]  += RC->States[k]->ProbMark[pos1]; 
		}
	}
}


//Updating ProbMark by bottom up traversal of  ROG following suffix links. Computing  Prob(R(n,p,H(w)))
void ROGUpdate(H_M_Node* node, int n){
	int i,j,k;
	int row = MainData::WordLen - node->len;
	int pos = n % row;						//position in w.ProbMark to record Prob(R(n,p,H(w)))

	
	for(i = 0; i < node->NRChilds; i++){
		H_M_Node* RC = static_cast<H_M_Node*>(node->RChilds[i]);
		ROGUpdate(RC,n);				   //processing of right successors of w
		
		for(k = 0; k < RC->NStates; k++){
			int state = RC->States[k]->RParentPos;
			for(j = 0; j < MainData::NOccur; j++){
				node->States[state]->FirstTemp[j] += RC->States[k]->FirstTemp[j]; //Prob(R(n,p,H(w))) = Prob(RE(n,p,H(w))) + Sum_{t,w=rpred(t)}Prob(R(n,p,H(t)))
				RC->States[k]->FirstTemp[j] = 0;
			}	
		}

	}
	//updating of data structures
	for(j = 0; j < MainData::NOccur; j++){
		for(k = 0; k < node->NStates; k++){
			node->States[k]->ProbMark[j*row + pos] = node->States[k]->FirstTemp[j];
			if(node == NodeOv::Root){
				node->States[k]->FirstTemp[j] = 0;
			}
		}
	}



	return;
}








//For a node w compute  Prob(D(n-|Back(w)|,p-1,w),q)
//Input parameters: node - current processing node w; n - stage of the algorithm minus motif length (we will think that n is algorithm stage)
//depth - depth of branch leading from root to the node
void DProbCalc(H_M_Node* node, int n, int depth){
	int i,j;

	int row = MainData::WordLen - node->len;
	int pos = n % row;						//position of Prob(R(n-|Back(w)|,p-1,w),q) in <w,q>.ProbMark (independent from q)
	memset(HDProbs[depth], 0x00, H_M_Node::NumAllStates*MainData::NOccur * sizeof(double));	//structure to store  Prob(D(n-|Back(w)|,p-1,w),q)

	for(i = 0; i < node->NStates; i++){			//for all q in AllState(w)
		for(j = 0; j < MainData::NOccur; j++){	// for all j = 0,...,p0-1
			if(node->rootchild == 1){	//if lpred(w) == root
				HDProbs[depth][i*MainData::NOccur + j] = node->States[i]->ProbMark[j*row + pos]; //Prob(D(n-|Back(w)|,p-1,w),q) = Prob(R(n-|Back(w)|,p-1,w),q)
			}
			else{
				int k;
				double Sum = 0;														//variable to compute Sum_{q' in PriorState(w,q)}Prob(D(n-|Back(w)|-|Back(lpred(w))|,p-1,lpred(w)),q')*Prob(q',Back(w),q)
				for(k = 0; k < node->States[i]->NumBackProbs; k++){					//for all q' in PriorState(w,q)
					int state = node->States[i]->BackProbs[k].pos;
					double p = node->States[i]->BackProbs[k].prob;
					Sum = Sum + HDProbs[depth-1][state*MainData::NOccur + j]*p;		//Sum = Sum + Prob(D(n-|Back(w)|-|Back(lpred(w))|,p-1,lpred(w)),q')*Prob(q',Back(w),q)
				}
						
				HDProbs[depth][i*MainData::NOccur + j] = Sum + node->States[i]->ProbMark[j*row + pos]; //Prob(D(n-|Back(w)|,p-1,w),q) = Sum + Prob(R(n-|Back(w)|,p-1,w),q)
			}
		}
	}

}



//For all right deep nodes r compute  close probabilities Prob(C(n,p,r),q) and add it to <r,q>.FirstTemp[p-1]
//Input parameters: node - current processing node w; n - stage of the algorithm minus motif length (we will think that n is algorithm stage)
//depth - depth of branch leading from root to the node
void ClosePartCalc(H_M_Node* node, int n,int depth){
	int i,j,k,l;
	for(i = 0; i < node->NDLinks; i++){			// for all r such that h*(w,r) is not empty
		H_M_Node* rnode = static_cast<H_M_Node*>(node->DeepLinks[i]); //node r
		for(l = 0; l < node->NStates; l++){		//for all q' in AllState(w)						
			double Sum = 0;
			for(k = 0; k < node->States[l]->NumDeepProbs[i]; k++){	//for all q in ReachState(q',h*(w,r)) 
				int state = node->States[l]->DeepProbs[i][k].pos;  
				double p = node->States[l]->DeepProbs[i][k].prob;
				for(j = 0; j < MainData::NOccur; j++){				//for all j = 0,..,p0
					if(j == 0){
						rnode->States[state]->FirstTemp[0] -= HDProbs[depth][l*MainData::NOccur + j]* p; //<r,q>.FirsTemp[0] -= Prob(D(n-|w|,1,w),q')*Back(q',Back(h*(w,r)),q)
					}
					else{
						rnode->States[state]->FirstTemp[j] += (HDProbs[depth][l*MainData::NOccur + j-1] - HDProbs[depth][l*MainData::NOccur + j])*p; //<r,q>.FirsTemp[j] += (Prob(D(n-|w|,j-1,w),q') - Prob(D(n-|w|,j,w),q')*Back(q',Back(h*(w,r)),q) 
					}
				}
			}
			
		}
	}
}





//Computation of RE-sets probabilities, store them in FirstTemp
//Input parameters: node - current processing node w; n - stage of the algorithm minus motif length (we will think that n is algorithm stage)
//depth - depth of branch leading from root to the node
void ComputeMainInduction(H_M_Node* node, int n, int depth){
	int i;
	if(node->num == 18){
		i = 8;
	}
	//////////////////PROCESSING OF NODE //////////////////////////////////////
	if(node != NodeOv::Root){
		//////// 1. Computation of D-sets probabilities //////////////////////
		DProbCalc(node, n, depth);
		
		//////// 2. Computation of RE-sets probabilities, store them in FirstTemp//////////////////////////	
	
		//////////A. Computation of Prob(C(n,p-1,r))-Prob(C(n,p,r)) for all r such that exist a class h*(w,r)
		if(node->ldeep == 1){ //if w is left deep node
			ClosePartCalc(node, n, depth);
		}
	}
	
	//////////B. If model is HMM, computation of Prob(F(n,p-1,w))-Prob(F(n,p,w))
	if((node->rdeep == 1)&&(MainData::order < 0)){
		hmm_ovgraf::FarPartCalc(node);
	}

	if(node != NodeOv::Root)
		depth++;
	////////////////////////PROCESSING OF CHILDS OF NODE//////////////////////////////////////////////
	
	for( i = 0; i < node->NLChilds; i++){
		H_M_Node* LC = static_cast<H_M_Node*>(node->LChilds[i]);
		ComputeMainInduction(LC, n, depth);
	}
	return;
}





// Compute Prob(B(n-m,p-1))
void CalcBnp(int n){
	int i,j,k,s;
	double sum;
	
	int pos = n % MainData::WordLen;
	H_M_Node* root = static_cast<H_M_Node*>(NodeOv::Root);

	double* temp = PrevBnpProbs; //structure for temporary storing probabilities Prob(B(n-m-1,p),q)
	PrevBnpProbs = H_M_Node::BnpProbs;
	H_M_Node::BnpProbs = temp;
	memset(H_M_Node::BnpProbs, 0x00, MainData::NOccur*H_M_Node::NumAllStates * sizeof(double));

	if(MainData::order > 0)
			s = MainData::AlpSize;

	//internal loop on j = 0,..,p0-1
	for(j = 0; j < MainData::NOccur; j++){ 
		
		//compute Prob(B(n-m,j).V),q)
		for(i = 0; i < H_M_Node::NumAllStates; i++){ //for all q in Q
			if(MainData::order < 0){
				s = H_M_Node::ConsistStNums[i];		//for all q' s.t. ï(q',q)>0
			}
			
			sum = 0; //variable to compute Sum_{q'}Prob(B(n-m-1,j),q')*ï(q',q)
	        for(k = 0; k < s; k++){
				int state = H_M_Node::ConsistStMatrix[i][k];
				sum += PrevBnpProbs[state* MainData::NOccur + j]*H_M_Node::TransProbMatrix[i][k];  //Sum += Prob(B(n-m-1,j),q')*ï(q',q)
			}
				
			H_M_Node::BnpProbs[i*MainData::NOccur + j]= sum;
		}
	

		
		for(i = 0; i < root->NStates; i++){
			int  state = root->States[i]->ID;
			H_M_Node::BnpProbs[state*MainData::NOccur + j] += root->States[i]->ProbMark[j*MainData::WordLen + pos]; // +=Prob(R(n-m,j,q))
		}
	}

	return;
}
////////////////////////////////////////////////////

void DebPrint1(H_M_Node* node, H_M_Node* LParent){
	int i,j,k;	
	/*
	if(node->num == 0){
		cout<<"Step "<<step<<"\n-------\n--------\n\n";
	}
	*/
	
	//int row = MainData::WordLen - node->len;

	//H_M_State* States = reinterpret_cast<H_M_State*>(node->States);

	cout<<"Num of node: "<<node->num<<'\n';
	cout<<"Num of states: "<<node->NStates<<'\n';
	for(i = 0; i < node->NStates; i++){
		cout<<"Data related to the state "<<node->States[i]->ID<<'\n';
		cout<<"RP pos: " <<node->States[i]->RParentPos<<'\n';
	
		if(node->rdeep == 1){
			cout<<"Word Probs \n";
			HMM_State* State = static_cast<HMM_State*>(node->States[i]); 
			for(j = 0; j < State->NumWordProbs; j++){
				int startstate = State->WordProbs[j].pos;
				cout<<"Prob("<<startstate<<",H~(w),"<<State->ID <<") = "<<State->WordProbs[j].prob<<'\n';
			}
		
			cout<<"\n\n";
		}
		int row = MainData::WordLen - node->len;
		cout<<"ProbMark \n";
		for(k = 0; k < MainData::NOccur; k++){
			for(j = 0; j< row; j++){ 
				cout<<std::setprecision(15)<<node->States[i]->ProbMark[k*row + j]<<'\t';
			}
			cout<<'\n';
		}
		cout<<"\n\n";
	}
	cout<<"\n\n\n";


	for ( i = 0; i < node->NLChilds; i++){
		H_M_Node* G = static_cast<H_M_Node*>(node->LChilds[i]);
		DebPrint1(G,node);
	}
}


void DebHMMPrint(H_M_Node* node, int step){
	int i,j,k;	

	if(node->num == 0){
		cout<<"Step "<<step<<"\n-------\n--------\n\n";
	}
	
	
	//int row = MainData::WordLen - node->len;

	//H_M_State* States = reinterpret_cast<H_M_State*>(node->States);

	cout<<"Num of node: "<<node->num<<'\n';
	cout<<"Num of states: "<<node->NStates<<'\n';
	for(i = 0; i < node->NStates; i++){
		cout<<"Data related to the state "<<node->States[i]->ID<<'\n';
	
		int row = MainData::WordLen - node->len;
		cout<<"ProbMark \n";
		for(k = 0; k < MainData::NOccur; k++){
			for(j = 0; j< row; j++){ 
				cout<<std::setprecision(15)<<node->States[i]->ProbMark[k*row + j]<<'\n';
			}
			cout<<'\n';
		}
		cout<<"\n\n";
	}
	cout<<"\n\n\n";


	for ( i = 0; i < node->NLChilds; i++){
		H_M_Node* G = static_cast<H_M_Node*>(node->LChilds[i]);
		DebHMMPrint(G,step);
	}
}








//////////////////////////////////////////////////////////////////


//main algorithm 
void H_M_Node::ProbCalc(void){
	int i;
	H_M_Node* root = static_cast<H_M_Node*>(NodeOv::Root);

	
	/////////////////PRE-PROCESSING/////////////////////
	Preprocessing();							 //initialization of data in descriptors on nodes and states of ovgraph
	

	if(MainData::TLen > MainData::WordLen-1){
		Stepm(root);							//compute initial probabilities Prob(R(m,p,w))
	}
	
	if(MainData::TLen > MainData::WordLen){
		int n = 0;


   //////////////////MAIN LOOP///////////////////////////////
		for(i = MainData::WordLen+1; i<= MainData::TLen; i++){ //i = m+1,..,n0
			if(n >= MainData::WordLen - 1){
				CalcBnp(n);		//computation of B-sets probabilities
			}

			if(MainData::order < 0){
				hmm_ovgraf::TransMatrixProduct(n);	//for HMM, computation of Prob(V^n-m,q), V - alphabet
			}
			else{
				MModel_Prob::TransStepProbCalc(H_M_Node::TransStepProbList); //for Markov mode, computation of Prob(V^n-m,q)
			}	
		
			
			if(MainData::order > 0){  //Let Markov model is given
				m_ovgraf::CalEProbOne(M_TrTree::Root);	 //computation of Prob(E(n,1,H~(r))), r is right deep nodes
				m_ovgraf::CalFarProbs(M_TrTree::Root);	 //computation of far sets probabilities Prob(F(n,1,H~(r)))
			}
		
			//////////////////////
	
			ComputeMainInduction(root, n, 0); //computation of RE-sets probabilities by depth first traversal of graph
		
			ROGUpdate(root,n);				  //computation of R-sets probabilities by bottom up traversal of graph, updating of data in ProbMark

			n++;

		}
		
		//POST-PROCESSING
		for(i = n; i < MainData::TLen; i++){
			CalcBnp(i);						//computation of probabilities Prob(B(n,p0),q), n = n0-m, ...,n0, q in Q
		}
		double p = 0;
		for(i = 0; i < NumAllStates; i++){  //computation of P-value = Sum_q(Prob(B(n,p0),q))
			p += BnpProbs[i*MainData::NOccur + MainData::NOccur - 1]; 
		}
		MainData::Pvalue = p;

	}


	ClearData(); //delete descriptors of nodes of OVGraph
	
	return;
}



