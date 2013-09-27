
#include "mmodel_prob.h"
#include "m_trtree.h"
#include "maindata.h"
#include "m_ovgraf.h"

//Let Markov model of order K is given

int M_TrTree::NumTrNodes = 0;			//Number of nodes of MTrTree
M_TrTree* M_TrTree::Root = nullptr;		//root of MTrTree

//Design set of states for all nodes of ovgraph, compute back probabilities 
//Input parameters: node - node w, LPnode - lpred(w), back - Back(w), word - word w

void m_ovgraf::States_and_BackProbs(H_M_Node* node, H_M_Node* LPnode, string &back, string &word){
	int i,j;
	if(node->num == 0){ //if w is root
		node->States = Malloc<H_M_State*>(H_M_Node::NumAllStates);  //initialization of w.States
		node->NStates = H_M_Node::NumAllStates;
		for(i = 0; i < H_M_Node::NumAllStates; i++){
			node->States[i] = new H_M_State;
			node->States[i]->ID = i;
		}
	}
	else{

		int pos1 = PosNNodeBacks[node->num];
		int backlen = node->len - LPnode->len;					

		SubStr(WNodeBacks,pos1,backlen,back);						//take Back(w)

		for(i = 0; i < backlen; i++)
			word[LPnode->len + i] = back[i];						//take word w

		vector<int> statesid;
		statesid.reserve(H_M_Node::NumAllStates);
		MModel_Prob::CalcReachStates(word, node->len, statesid);	//statesid contains ids of states from AllState(w)
		
		node->NStates = (int)statesid.size();
		
		node->States = Malloc<H_M_State*>(node->NStates);			
		for(i = 0; i < node->NStates; i++)
			node->States[i] = new H_M_State;

						////////Back Probs/////////;
		int priorsize;								//|PriorState(w,q)|
		if(node->len >= MainData::order)			//if |w| >= K
			priorsize = LPnode->NStates;
		else
			priorsize = MModel_Prob::NumPower(MainData::AlpSize, backlen); 

		int pow = MModel_Prob::Power/priorsize;
		int pow1 = MModel_Prob::NumPower(MainData::AlpSize, LPnode->len);

		int lpstate, pos;
		double prob;
		for(i = 0; i < node->NStates; i++){         // for all  q in AllStates(w)
			
			node->States[i]->NumBackProbs = priorsize;
			node->States[i]->BackProbs = Malloc<PriorList>(priorsize);
			int state = statesid[i];
			int cd = state/priorsize;
			node->States[i]->ID = state;
			
			
			if(node->len >= MainData::order){		//if |w|>K
				for(j = 0; j < LPnode->NStates; j++){  //for all q' in AllState(lpred(w))
					lpstate = LPnode->States[j]->ID;	//q'
				    prob = MModel_Prob::TermCondProb(lpstate,back,backlen, state); //compute Prob(q',Back(w),q)
					node->States[i]->BackProbs[j].pos = j;
					node->States[i]->BackProbs[j].prob = prob;
				}
			}
			else{
				for(j = 0; j < priorsize; j++){
					lpstate = j*pow + cd;
					pos = lpstate/pow1;       //position q' int lpred(w).States
					prob = MModel_Prob::TermCondProb(lpstate,back,backlen, state); //compute Prob(q',Back(w),q)
					node->States[i]->BackProbs[j].pos = pos;
					node->States[i]->BackProbs[j].prob = prob;
				}
			}
		}

	} /*if(node->num != 0)*/
	  for( i = 0; i < node->NLChilds; i++){
		H_M_Node* LC = static_cast<H_M_Node*>(node->LChilds[i]);
		States_and_BackProbs(LC,node,back,word);
	}
}






///////////Functions to calculate states associated to an overlap w from OV(HH)///////////////////


//Initialization of matrices H_M_Node::ConsistStMatrix and H_M_Node::TransProbMatrix (see h_m_node.h)

void m_ovgraf::CrConsistStatesMatrix(void){
	H_M_Node::ConsistStMatrix = new int* [H_M_Node::NumAllStates]();
    typedef double * doublePtr;
    H_M_Node::TransProbMatrix = new doublePtr[H_M_Node::NumAllStates]();
	int i;
	for(i = 0; i < H_M_Node::NumAllStates; i++){
		H_M_Node::ConsistStMatrix[i] = new int[MainData::AlpSize]();
		H_M_Node::TransProbMatrix[i] = new double[MainData::AlpSize]();

		vector<int> vec = MModel_Prob::ConsistStates(i);
		size_t j;
		for(j = 0; j < vec.size(); j++){
			H_M_Node::ConsistStMatrix[i][j] = vec[j];
			double p = TransitionProb(vec[j],i);
			H_M_Node::TransProbMatrix[i][j] = p;
		}
	}
};


KM_TrTree* *KNodes;
////////////////////////////////////////////////////

//Create MTrTree by depth-first traversal of AC trie
//Input parameters: node t - current processed node of AC trie; LPnode - prefix predecessor of t; len - |t|; word - word t
void M_TrTree::CreateTree(NodeAC* node, M_TrTree* LPnode, int len, string &word){
	M_TrTree* tnode;
	std::list<NodeAC*>::iterator i;
	InternAC* node1 = static_cast<InternAC*>(node);
	int j;
	M_TrTree::NumTrNodes ++;

	if(node == NodeAC::ACRoot){	//if t is root of AC trie
		M_TrTree::Root = new M_TrTree;
		tnode = M_TrTree::Root;
        typedef M_TrTree* M_TrTreePtr;
        tnode->Childs = new M_TrTreePtr[MainData::AlpSize]();
		tnode->NStates = MModel_Prob::Power;
		tnode->States = Malloc<int>(MModel_Prob::Power);
		for(j = 0; j < MModel_Prob::Power; j++){			//list of states for root contains all  possible words of length K
			tnode->States[j] = j;
		}
        typedef KM_TrTree* KM_TrTreePtr;
        KNodes = new KM_TrTreePtr[MModel_Prob::Power]();
	}
	else{
		word[len - 1] = MainData::IToa(node->sign);

		if(len <=MainData::order){							// if |t|<=K
			if(len == MainData::order)
                tnode = new KM_TrTree();					//creating of new node
			else
                tnode = new M_TrTree();
			//initialization of parameters in descriptor of t
			tnode->len = len;
			tnode->sign = node->sign;
            typedef M_TrTree* M_TrTreePtr;
            tnode->Childs = new M_TrTreePtr[MainData::AlpSize]();
			LPnode->Childs[tnode->sign] = tnode;

			
			vector<int> statesid;
			int pow = MModel_Prob::NumPower(MainData::AlpSize, MainData::order - len);
			statesid.reserve(pow);
			MModel_Prob::CalcReachStates(word, len, statesid);

			//fill list of states
			if(len < MainData::order){
				tnode->NStates = statesid.size();
				tnode->States = Malloc<int>(tnode->NStates);
				for(j = 0; j < tnode->NStates; j++){
					tnode->States[j] = statesid[j];
				}
			}
			else{	
				int cd = MModel_Prob::Code(word,len);
				tnode->NStates = 1;
				tnode->States = Malloc<int>(1);
				tnode->States[0] = cd;

                KNodes[cd] = static_cast<KM_TrTree*>(tnode);
			}
			
			
		}
	}

	///recursion
	if(len < MainData::order){
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* LC = static_cast<NodeAC*>(*i);
			CreateTree(LC, tnode, len +1, word);
		}
	}
	return;
}




//compute word probabilities
void m_ovgraf::WordProbs(string &word){
	int i,j,k;

	/////////////Word Probs/////////
	int restlen = MainData::WordLen - MainData::order;
	string wordrest; // wordrest = word[K+1,|word|]
	wordrest.resize(restlen);

	vector<vector<int>> RPList;				//for each w in OV(HH), RPList(w) is list h in HH s.t. w=rpred(h)
	RPList.resize(NodeOv::NumOVNodes);
	//Remark. Words in RPList(w) have prefix order 

	vector<int> RPSizes;					//RPSizes(w) is size of RPList(w)
	RPSizes.resize(NodeOv::NumOVNodes);
	
	vector<int> KNumLinks;					//for each word t of length K, KNumLinks(w) is number of right deep nodes r s.t. exists h in HH having prefix t and rpred(h) = t 
	KNumLinks.resize(MModel_Prob::Power);
	
	vector<int> Prefixes;					//for each h in HH, Prefixes(h)  - prefix of length K of h
	Prefixes.resize(MainData::NWords);

	vector<int> Suffixes;					//for each h in HH, Suffixes(h)  - suffix of length K of h
	Suffixes.resize(MainData::NWords);


	////////Create RPList; Suffixes; Prefixes//////////////////


	for(i = 0; i < MainData::NWords; i++){
		H_M_Node* RP = static_cast<H_M_Node*>(Nodes[RLeafPreds[i]]);
		RPSizes[RP->num] ++;
		
		int pos1 = PosNLeafWords[i];
		SubStr(WLeafWords, pos1, MainData::WordLen, word);
		
		int prefcd = MModel_Prob::PrefixN(MainData::order, word);
		int sufcd = MModel_Prob::SuffixN(MainData::order, word, MainData::WordLen);
		Prefixes[i] = prefcd;
		Suffixes[i] = sufcd;
		
		////////////////Compute probabilities Prob(H~(r),q), H~(r) - set of motif words h s.t. t = rpred(h); record the probability in <r,q>.ProbMark ///////////////////////////////////////
		double prob = MModel_Prob::TermProb(word, MainData::WordLen, sufcd);
		int pos_in_matrix = MainData::WordLen - RP->len - 1;

		if(RP->len >= MainData::order){
			RP->States[0]->ProbMark[pos_in_matrix] += prob;
		}else{
			int pow = MModel_Prob::NumPower(MainData::AlpSize,RP->len);
			int rppos = sufcd/pow;
			RP->States[rppos]->ProbMark[pos_in_matrix] += prob;
		}

	}

	
	for(i = 0; i < NodeOv::NumOVNodes; i++){
		RPList[i].reserve(RPSizes[i]);
	}

	RPSizes.clear();

	for(i = 0; i < MainData::NWords; i++){
		H_M_Node* RP = static_cast<H_M_Node*>(Nodes[RLeafPreds[i]]);
		RPList[RP->num].push_back(i);
	}

	////////////////Create KNumLinks////////////////////////////

	bool* Flags = Malloc<bool>(MModel_Prob::Power);

	for(i = 0; i < NodeOv::NumOVNodes; i++){
		int s = RPList[i].size();
		memset(Flags, 0x00, MModel_Prob::Power * sizeof(bool));
		
		for(j = 0; j < s; j++){
			int nword = RPList[i][j];
			int prefcd = Prefixes[nword];

			if(Flags[prefcd] == 0){ 
				KNumLinks[prefcd] ++;
				Flags[prefcd] = 1;
			}
		}
	}
	
	free(Flags);

	///////////////////Initialization of WordProbs data structures (see MTrTree.h)///////////////////////
	for(i = 0; i < MModel_Prob::Power; i++){
		if(KNodes[i] != nullptr){
			KNodes[i]->NumWLinks = 0;
			KNodes[i]->WLinks = Malloc<H_M_Node*>(KNumLinks[i]);
            typedef PriorList* PriorListPtr;
            KNodes[i]->WProbs = new PriorListPtr[KNumLinks[i]]();
			KNodes[i]->NumWProbs = Malloc<int>(KNumLinks[i]);
		} 	
	}

	KNumLinks.clear();

	////////////////Compute Word Probabilities, for each r in OV(HH) computing of  Sum_{h in HH, r = rpred(h)}(Prob(pref_K(h),h[K+1,|h|],suf_K(h)))//////////////////////////////////////////
	//let pref_K(h) (suf_K(h)) be prefix (suffix) of h of length K

	double* Probs = Malloc<double>(MModel_Prob::Power);
	 

	int curprefcd, prefcd, s, pow, sufcd, rppos,nword, numprobs;
	H_M_Node* RP;
	KM_TrTree* knode;
	double prob;

	for(i = 0; i < NodeOv::NumOVNodes; i++){  //for r in OV(HH)
		RP = static_cast<H_M_Node*>(Nodes[i]);
		s = (int)RPList[i].size();
		pow = MModel_Prob::NumPower(MainData::AlpSize, RP->len);
		j = 0;
		while(j < s){      //loop of h in RPList(r)
			nword = RPList[i][j];
			curprefcd = Prefixes[RPList[i][j]];
			prefcd = curprefcd;

			knode =KNodes[prefcd];  //node in MTRTree cooresponding prefixes of h of length K
			knode->NumWLinks ++;
			knode->WLinks[knode->NumWLinks -1] = RP;

			memset(Probs, 0x00, RP->NStates * sizeof(double));

			numprobs  = 0;
			
															
			while((curprefcd == prefcd)&&(j < s)){		//loop on h in RPList(r) having the same prefix t of length K 

				
				int pos1 = PosNLeafWords[nword] + MainData::order;
				SubStr(WLeafWords, pos1, restlen, wordrest);
				sufcd  = Suffixes[nword];
	
				if(RP->len >= MainData::order)	rppos = 0;
				else rppos = sufcd/pow; 
				
				prob = MModel_Prob::TermCondProb(prefcd, wordrest, restlen, sufcd); //computing of Prob(t,h[K+1,|h|],suf_K(h))

				if(Probs[rppos] == 0) numprobs ++;

				Probs[rppos] += prob;
				j++;
				if(j < s){
					nword = RPList[i][j];
					prefcd  = Prefixes[nword];
				}
			}

			//fill word probs to t.WProbs
			knode->NumWProbs[knode->NumWLinks -1] = numprobs;
            knode->WProbs[knode->NumWLinks -1] = new PriorList[numprobs];
			int l = 0;
			for(k = 0; k < RP->NStates; k++){
				if(Probs[k] != 0){
					knode->WProbs[knode->NumWLinks -1][l].pos = k;
					knode->WProbs[knode->NumWLinks -1][l].prob = Probs[k];
					l++;
				}
			}
		}
	}
	free(Probs);
	

}

/////////////////////////////////////////////////////////////////



//for all right deep nodes r in OV(HH) of OVGraf computation of Prob(E(n,1,r),q), q in AllState(w)
//by depth-first traversal of MTrTree

//Input parameters: tnode - current processing node t of MTrTree
void  m_ovgraf::CalEProbOne(M_TrTree* tnode){
	
	int i,j, k,predstate, state;
	

	if(tnode->len != 0){	//if t is not root of MTrTree 
		int pow = MModel_Prob::Power/MainData::AlpSize;
		int pow1 = MModel_Prob::NumPower(MainData::AlpSize, tnode->len -1);
		memset(OrderProbs[tnode->len-1], 0x00, tnode->NStates*sizeof(double));
		for(i = 0; i < tnode->NStates; i++){  //for all states q in AllState(t)
			int cd = tnode->States[i]/MainData::AlpSize;
			for(k = 0; k < MainData::AlpSize; k++){	//for all a in Alp
													// let t =t'.a, q' be a state in AllState(t') s.t. q = q'[2,K].a
				predstate = cd + k*pow;				//position of q' in t'.States
				double mp = MainData::MarkovProbs[tnode->sign][predstate];
				if(tnode->len == 1){	//if |t| ==1
					OrderProbs[0][i] += H_M_Node::TransStepProbList[predstate]*mp;  //Compute Prob(V^{n-m}.a,q) = Prob(V^{n-m},q')*Prob(q',a,q) and add it to  OrderProbs[0][i]
				}
				else{
					int pos = predstate/pow1;
					OrderProbs[tnode->len-1][i] += OrderProbs[tnode->len-2][pos]*mp;	//Compute Prob(V^{n-m}.t,q) = Prob(V^{n-m}.t',q')*Prob(q',a,q) and add it to  OrderProbs[|t|-1][i]
				}
			}
		}
	}
	
	if(tnode->len == MainData::order){													//if |t|==K
		KM_TrTree* knode = static_cast<KM_TrTree*>(tnode);
		predstate = knode->States[0];
		for(i = 0; i < knode->NumWLinks; i++){											//for all r s.t. exists h in HH, having prefix t and rpred(h) = r 
			H_M_Node* rovnode = knode->WLinks[i];
			for(j = 0; j < knode->NumWProbs[i]; j++){									//for all rq in AllState(r), s.t. exists h in HH, having prefix t and rpred(h) = r, where rq is h-reachable  
				state = knode->WProbs[i][j].pos;
				double prob = knode->WProbs[i][j].prob;
				rovnode->States[state]->FirstTemp[0] += OrderProbs[MainData::order -1][0]*prob; //compute Prob(V^{n-m}.H~(r),rq) and add it to <r,q>.Firstemp[0]
			}
		}
	}
	else{//processing of prefix childs of t
		for(i = 0; i < MainData::AlpSize; i++){
			if(tnode->Childs[i] != nullptr)
				CalEProbOne(tnode->Childs[i]);
		}
	}
	return;
}



//for all right deep nodes r in OV(HH) of OVGraf computation of Prob(F(n,p,r),q), q in AllState(w), p = 1,..,p0
//by depth-first traversal of MTrTree

//Input parameters: tnode - current processing node t of MTrTree
void  m_ovgraf::CalFarProbs(M_TrTree* tnode){
	int i,j, k,predstate, state;
	

	if(tnode->len != 0){												//if t is not root of MTrTree 
		int pow = MModel_Prob::Power/MainData::AlpSize;
		int pow1 = MModel_Prob::NumPower(MainData::AlpSize, tnode->len -1);
		memset(OrderProbs[tnode->len-1], 0x00, tnode->NStates*MainData::NOccur*sizeof(double));
		for(i = 0; i < tnode->NStates; i++){								//for all states q in AllState(t)
			int cd = tnode->States[i]/MainData::AlpSize;
			for(k = 0; k < MainData::AlpSize; k++){							//for all a in Alp
																			// let t =t'.a, q' be a state in AllState(t') s.t. q = q'[2,K].a
				predstate = cd + k*pow;									
				double mp = MainData::MarkovProbs[tnode->sign][predstate];	//position of q' in t'.States
				for(j = 0; j < MainData::NOccur; j++){						//for all j = 0, .., p0-1
					if(tnode->len == 1){
						OrderProbs[0][i*MainData::NOccur + j] += H_M_Node::BnpProbs[predstate*MainData::NOccur + j]*mp;		//Compute Prob(B(n-m,j).a,q) = Prob(B(n-m),q')*Prob(q',a,q) and add it to  OrderProbs[0][i*p0+j]
					}
					else{
						int pos = predstate/pow1;
						OrderProbs[tnode->len-1][i*MainData::NOccur + j] += OrderProbs[tnode->len-2][pos*MainData::NOccur + j]*mp;  //Compute Prob(B(n-m,j).t,q) = Prob(B(n-m,j).t',q')*Prob(q',a,q) and add it to  OrderProbs[|t|-1][i*p0+j]
					}
				}
			}
		}
	}
	
	if(tnode->len == MainData::order){						//if |t|==K
		KM_TrTree* knode = static_cast<KM_TrTree*>(tnode);
		predstate = knode->States[0];
		for(i = 0; i < knode->NumWLinks; i++){				//for all r s.t. exists h in HH, having prefix t and rpred(h) = r 
			H_M_Node* rovnode = knode->WLinks[i];
			for(k = 0; k < knode->NumWProbs[i]; k++){		//for all rq in AllState(r), s.t. exists h in HH, having prefix t and rpred(h) = r, where rq is h-reachable  
				state = knode->WProbs[i][k].pos;
				double prob = knode->WProbs[i][k].prob;
				for(j = 0; j < MainData::NOccur; j++){
					//compute  SUmProb(j) = Sum_{h, pref_k(h) = t, rpred(h)=r}(Prob(B(n-m,j+1).h,rq) = Sum_{/../}Prob(B(n-m,j+1).pref_K(h),pref_K(h))*Prob(Pref_K(h),h[k+1,K],suf_k(h))
					if(j == 0){
						rovnode->States[state]->FirstTemp[0] -= OrderProbs[MainData::order -1][0]*prob;   // subtract FProb(0) from <r,q>.Firstemp[0]
					}
					else{
						rovnode->States[state]->FirstTemp[j] += (OrderProbs[MainData::order -1][j -1] - OrderProbs[MainData::order -1][j])*prob; // add (FProb(j-1) - FProb(j)) to <r,q>.Firstemp[j]
					}
				}
			}
		}
	}
	else{//processing of prefix childs of t
		for(i = 0; i < MainData::AlpSize; i++){
			 if(tnode->Childs[i] != nullptr)
				CalFarProbs(tnode->Childs[i]);
		}
	}
	return;
}




