#include "hmm_ovgraf.h"
#include "h_m_node.h"


////////////////////For Preprocessing////////////////////

//Initialization of matrices H_M_Node::ConsistStMatrix and H_M_Node::TransProbMatrix (see h_m_node.h)
void hmm_ovgraf::CrConsistStatesMatrix(void){
	H_M_Node::ConsistStMatrix = Malloc<int*>(H_M_Node::NumAllStates);
	H_M_Node::TransProbMatrix = Malloc<double*>(H_M_Node::NumAllStates);
	H_M_Node::ConsistStNums = Malloc<int>(H_M_Node::NumAllStates);
	int i,j;
	for(i = 0; i < H_M_Node::NumAllStates; i++){
		vector<int> vec = ConsistStates(i);
		int s= (int)vec.size();
		H_M_Node::ConsistStNums[i] = s;
		H_M_Node::ConsistStMatrix[i] = Malloc<int>(s);
		H_M_Node::TransProbMatrix[i] = Malloc<double>(s);
		for(j = 0; j < s; j++){
			H_M_Node::ConsistStMatrix[i][j] = vec[j];
			H_M_Node::TransProbMatrix[i][j]= TransitionProb(vec[j],i);
		}
	}
};

//Design set of states for all nodes of ovgraph, compute back probabilities 
//Input parameters: node - node w, LPnode - lpred(w), back - Back(w)
void hmm_ovgraf::States_and_BackProbs(H_M_Node* node, H_M_Node* LPnode, string &back){
	int i,j;
	if(node->num == 0){ //if w == root
		node->States = Malloc<H_M_State*>(H_M_Node::NumAllStates); //AllStates(root) == Q
		node->NStates = H_M_Node::NumAllStates;

		for(i = 0; i < H_M_Node::NumAllStates; i++){
			if((MainData::order < 0)&&(node->rdeep == 1)) node->States[i] = new HMM_State;
			else node->States[i] = new H_M_State;
			
			node->States[i]->ID = i;
		}
	}
	else{
		vector<H_M_State> States; //w is not root
		States.reserve(H_M_Node::NumAllStates);

		int pos1 = PosNNodeBacks[node->num];
		int backlen = node->len - LPnode->len;
				////////Back Probs/////////
		SubStr(WNodeBacks,pos1,backlen,back);  //take Back(w)


		for(i = 0; i < H_M_Node::NumAllStates; i++){  //for all q in Q
			H_M_State NewStateData; 
			vector<PriorList> backprobs;
			backprobs.reserve(H_M_Node::NumAllStates);
			int flag = 0;

			for(j = 0; j < LPnode->NStates; j++){  //for all q' in AllState(lpred(w))
				int lpstate = LPnode->States[j]->ID;	
				
				double prob = TermCondProb(lpstate,back,backlen, i); //Compute Prob(q',Back(w),q)
				if(prob > 0){										//if Prob(q',Back(w),q) > 0 then q in AllState(w), q' in PriorList(w,q)
					PriorList List;									
					List.pos = j;
					List.prob = prob;
					flag = 1;
					backprobs.push_back(List);
				}
			}

    		if(flag == 1){	//create descriptor of <w,q>
				NewStateData.ID = i;
				NewStateData.NumBackProbs = backprobs.size(); //|PriorList(w,q)|
				NewStateData.BackProbs = Malloc<PriorList>(NewStateData.NumBackProbs);
				int k;
				for(k = 0; k <  NewStateData.NumBackProbs; k++) //for all q' in PriorList(w,q)
					NewStateData.BackProbs[k] = backprobs[k]; //<w,q>.BackProbs[k] = q'; <w,q>.BackProbs[k] = Prob(q',Back(w),q); 

				States.push_back(NewStateData);
			}
		}

		node->NStates = States.size();
		node->States = Malloc<H_M_State*>(node->NStates);
	//updating of data structures
		for(i = 0; i < node->NStates; i++){
			if((MainData::order < 0)&&(node->rdeep == 1)) node->States[i] = new HMM_State;
			else node->States[i] = new H_M_State;
			node->States[i]->ID = States[i].ID;
			node->States[i]->NumBackProbs = States[i].NumBackProbs;
			node->States[i]->BackProbs = States[i].BackProbs;
		}
		//States.clear();
	} /*if(node->num != 0)*/
	  for( i = 0; i < node->NLChilds; i++){
		H_M_Node* LC = static_cast<H_M_Node*>(node->LChilds[i]);
		States_and_BackProbs(LC,node,back);
	}
}



//compute word probabilities Prob(H~(r),q)
void hmm_ovgraf::WordProbs(string &word){
	int i,j,k;
	//Let r be right deep node, q in AllState(r), WordProbs[r][q] = Prob(H~(r),q)

	double***  WordProbs=  new double**[NodeOv::NumOVNodes]();
	for(i = 0; i < NodeOv::NumOVNodes; i++){
		H_M_Node* node = static_cast<H_M_Node*>(Nodes[i]);
		WordProbs[i] = new double*[node->NStates]();
		 for(j = 0; j < node->NStates; j++){
			 WordProbs[i][j] = new double[H_M_Node::NumAllStates]();
		 }
	}
	

	///////////// computation of word probs Prob(q',H~(r),q)/////////
	for(i = 0; i < MainData::NWords; i++){ //for all h in HH
		int pos1 = PosNLeafWords[i];
		int wordlen = MainData::WordLen;
	    SubStr(WLeafWords, pos1, wordlen, word); //take word h
		H_M_Node* RP = static_cast<H_M_Node*>(Nodes[RLeafPreds[i]]); //rpred(h) r
		

		for(j = 0; j <  RP->NStates; j++){		//for all q in AllState(r)

			int rpstate =  RP->States[j]->ID;
			for(k = 0; k <H_M_Node::NumAllStates; k++){ //for all q' in Q
				double prob = TermCondProb(k, word, wordlen, rpstate); //Prob(q',h,q)
				if(prob > 0){
					WordProbs[RP->num][j][k] += prob; 
				}
			
			}
			
		}

	}


	//Updating of data structures
	for(i = 0; i < NodeOv::NumOVNodes; i++){ //for all r in OV(HH)
		H_M_Node* node = static_cast<H_M_Node*>(Nodes[i]);
		if(node->rdeep == 1){  //if r is right deep
			for(j = 0; j < node->NStates; j++){ //for all q in AllState(r)
				HMM_State* State = static_cast<HMM_State*>(node->States[j]);
				int s = 0;
				for(k = 0; k < H_M_Node::NumAllStates; k++){
					if(WordProbs[node->num][j][k] >0) s++;
				}
			
				State->NumWordProbs = s;
				State->WordProbs = Malloc<PriorList>(s);
			
				if(State->NumWordProbs > 0){
					s = 0;
					for(k = 0; k < H_M_Node::NumAllStates; k++){
						if(WordProbs[node->num][j][k] >0){
							State->WordProbs[s].pos = k; 
							State->WordProbs[s].prob  = WordProbs[node->num][j][k]; 
							s++;
						}
					}
				}
			}
		}
	}

	//delete WordProbs;

	for(i = 0; i < NodeOv::NumOVNodes; i++){
		H_M_Node* node = static_cast<H_M_Node*>(Nodes[i]);
		for(j = 0; j < node->NStates; j++){

			delete[] WordProbs[i][j];
		}
		delete[] WordProbs[i];
	}
	delete[] WordProbs;
	return;
}

////////////////////////////////////////////////////////////

//Compute Prob(V^n-m)
void hmm_ovgraf::TransMatrixProduct(int n){
	int i, j;
	if(n == 0){ // if algorithm stage is m 
		for(i = 0; i < H_M_Node::NumAllStates; i++){
			for(j = 0; j < H_M_Node::ConsistStNums[i]; j++){
				int state = H_M_Node::ConsistStMatrix[i][j];
				double p = H_M_Node::TransProbMatrix[i][j];
				H_M_Node::TransStepProbMatrix[state* H_M_Node::NumAllStates + i] = p;
			}
		}	
	}
	else{
		double* temp = PrevTransStepProbMatrix;
		PrevTransStepProbMatrix = H_M_Node::TransStepProbMatrix;
		H_M_Node::TransStepProbMatrix = temp;
		memset(H_M_Node::TransStepProbMatrix, 0x00, H_M_Node::NumAllStates*H_M_Node::NumAllStates * sizeof(double));

		 int k;
		 for( i =0; i < H_M_Node::NumAllStates; i++){
			 for(j = 0; j < H_M_Node::ConsistStNums[i]; j++){
				 int state = H_M_Node::ConsistStMatrix[i][j];
				 double p = H_M_Node::TransProbMatrix[i][j];
				 for(k = 0; k < H_M_Node::NumAllStates; k++){
					H_M_Node::TransStepProbMatrix[k* H_M_Node::NumAllStates + i]  += PrevTransStepProbMatrix[k* H_M_Node::NumAllStates + state]*p;
				 }		 
			 }
		 }
	
	}
	return;
};





//For a node rnode compute  far probabilities Prob(F(n,p,rnode))
//node - current processing right deep node r
void hmm_ovgraf::FarPartCalc(H_M_Node* rnode){
	int i,j,k;

	for(i = 0; i < rnode->NStates; i++){ //for all q in AllState(r)
		HMM_State* State = static_cast<HMM_State*>(rnode->States[i]); 
		for(j = 0; j < MainData::NOccur; j++){
			double fprob = 0; //variable to compute far probabilities
			for(k = 0; k < State->NumWordProbs; k++){ //for all q' in Start(H~(r),q)
				int startstate =State->WordProbs[k].pos; //q'
				double prob = State->WordProbs[k].prob; //Prob(q',H~(r),q)
				if(j == 0){
					fprob +=(H_M_Node::TransStepProbMatrix[startstate] - H_M_Node::BnpProbs[startstate*MainData::NOccur])*prob; //+= (Prob(V^{n-m},q') - Prob(B(n-m,1),q'))*Prob(q',H~(r),q)
				}
				else{
					fprob +=(H_M_Node::BnpProbs[startstate*MainData::NOccur + j-1] - H_M_Node::BnpProbs[startstate*MainData::NOccur + j])*prob; //+= (Prob(B(n-m,j-1),q') - Prob(B(n-m,j),q'))*Prob(q',H~(r),q)
				}
			}
			State->FirstTemp[j] +=fprob;
		}
	}
};
