#include "nodeac.h"
#include "nodeov.h"
#include "nodebern.h"
//#include "H_M_Data.h"



////////////AUXILIARY PARAMETERS NEEDED ON THE PREPROCESSING PART OF THE PROGRAM////////////////////////////	
//1. PARAMETERS COMMON FOR ALL MODELS

NodeOv** Nodes = nullptr;			//List of nodes of overalp graph
int* LPreds = nullptr;				// list of lpred(w), w in OV(HH)
int* RPreds = nullptr;				// list of rpred(h), w in OV(HH) 

int *LLeafPreds = nullptr;			// list of lpred(h), h in HH
int *RLeafPreds = nullptr;			// list of rpred(h), h in HH 

int* NumLLeafChilds = nullptr;		//NumLLeafChilds(w) is numbers of h from HH, where w = lpred(h)
int* *LLeafChilds = nullptr;		// NumLLeafChilds(w) is lists of numbers of h from HH, where w = lpred(h)

vector< vector<int> > CLasses;		//list of overlap classes h*, Classes[h*] is list of numbers of words that are in the h*
int* CLassNum = nullptr;			//CLussNum(h) number of overlap class of h, h in HH

//2. PARAMETERS USED DURING GRAPH MINIMIZATION
vector<int> *NLinkToClass;			//NLinkToClass(w) is list of numbers of overlap classes corresponding to the deep links from DeepLinks
//during processing of current word h from HH, let w = lpred(h) and r = rpred(h)
int* NLink = nullptr;				//NLink[r] is number of deep link corresponding to h*(w,r) in w.DeepLinks  
int* NLOG = nullptr;				//NLOG[r] = w
int* NCLass = nullptr;				//NCLass[r] = number of class corresponding to h

//3. PARAMETERS NEDED FOR BERNOULLI MODEL 

double* LeafProbs = nullptr;		//list of LeafProbs[h] = Prob(h), h in HH
double* LeafBacks = nullptr;		//list of LeafBacks[h] =  Prob(Back(h)), h in HH


//4. PARAMETERS NEDED FOR MARKOV MODELS AND HMM
int* PosNNodeBacks;		//for each node w from OV(HH), PosNNodeBacks is position of Back(w) (w = x.Back(w), x=lpred(w)) in  WNodeBacks
std::string WNodeBacks;  //sequence of Back(w), w in OV(HH)
int WNodeBacks_len = 0;  // length of WNodeBacks

int* PosNLeafWords;		//for each word h from HH, PosNLeafWords[h] is position of entry of h in  WLeafWords
std::string WLeafWords;  //sequence of entrances of h, w in HH

int WNodeBacks_size = 0;  //current size of WNodeBacks
int WLeafWords_size = 0;   //current size of WLeafWords
int    MaxDepth = 1;									//maximal number of nodes on the path leading from root to a terminal node




//Computing of  sum  WNodeBacks_len of |Back(w)|, w in OV(HH), during depth-first traversal of AC trie  
//Input parameters:
//node - current considered node
//lplen - length of lpred(node)
//len - length of the node

void Compute_WNodeBacks_len(InternAC* node, int lplen, int len){

	if(node->main == 1){
		WNodeBacks_len += len - lplen;
	}
	
	if(len < MainData::WordLen -1){
		std::list<NodeAC*>::iterator i;
		for(i = node->LChilds.begin(); i != node->LChilds.end(); i++){
			InternAC* LC = static_cast<InternAC*>(*i);
			if(node->main == 1)	
				Compute_WNodeBacks_len(LC,len,len+1);
			else
				Compute_WNodeBacks_len(LC,lplen,len+1);
		}
	}
	return;
}


//Fill lists RPreds and  RLeafPreds with links to rpred(t), r in OV(HH) or t in  HH,
//during depth-first traversal of AC trie  
//Input parameters:
//node - current considered node


void  ComputeRPred(NodeAC* node){
	InternAC* RP = static_cast<InternAC*>(node->RParent);
	if(node->leaf == 0){
		InternAC* node1 = static_cast<InternAC*>(node);

		if(node1->main == 1){
			RPreds[node1->num] = RP->num;
		}
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			ComputeRPred(node2);
		}
	}
	else{
		LeafAC* leaf = static_cast<LeafAC*>(node);
		RLeafPreds[leaf->leafnum -1] = RP->num;
	}
}



//Creating of overlap graph OVGraf by depth-first traversal of AC trie, deleting of  AC trie
//Let ovnode be new creating node
//Input parameters:
//T - node of AC trie corresponding to ovnode 
//LParent - lpred(newnode), before creating node of OVGraf  
//prob - Prob(newnode) according to Bernoulli model, needed for Bern model
//Back - Prob(Back(newnode)) according to Bernoulli model, needed for Bern model
//len - length of newnode
//BackW - Back(newnode), needed for HMM and Markov models
//Word  - word corresponding to ovnode, needed for HMM and Markov models

void CreateLOG(NodeAC* T, NodeOv* LParent, double prob, double Back, int len, string &BackW, string &Word){
	std::list<NodeAC*>::iterator i;
	if(T != NodeAC::ACRoot) len++;
	int backlen = len - LParent->len;
	char symb;
	//Computing probabilities  Prob(Back(T)) and Prob(T) for the case of Bernoulli
	if(T!= NodeAC::ACRoot){
		if(MainData::order == 0){
			prob = prob* MainData::BernProb[T->sign];
			Back = Back*MainData::BernProb[T->sign];
		}else{
		//Computing of Back(T) and word corresponding to T for the case of Markov and HMM
			Word[len - 1] = MainData::IToa(T->sign);
			symb = BackW[backlen - 1];
			BackW[backlen - 1] = MainData::IToa(T->sign);
		}
	}

	//CREATING OF NEW NODE ovnode
	if(T->leaf == 0){ //if T is internal node
		InternAC* node = static_cast<InternAC*>(T);             // node == T
		NodeOv* ovnode;											//
		
		if((node->main != 0) && (node != NodeAC::ACRoot)){		//If T corresponds to an overlap 
		   
		   
			if(Nodes[node->num] == nullptr){
				ovnode = NodeOv::NewNode();  //creating ovnode, if node of OVGraf corresponding to T is not created before
				Nodes[node->num] = ovnode;
			}
			else{
				ovnode = Nodes[node->num];	//else, ovnode is node in Nodes coorresponding to T
			}

			//Set parameters in descriptor of ovnode
			ovnode->num = node->num;
			ovnode->len = len;

			LPreds[ovnode->num] = LParent->num;
			LParent->NLChilds ++;

			int rnum  = RPreds[ovnode->num];
			if(Nodes[rnum] == nullptr){
				Nodes[rnum] = NodeOv::NewNode();
			}
			ovnode->RParent = Nodes[rnum];
			ovnode->RParent->NRChilds ++;
		

			if(MainData::order == 0){ //if model is Bernolli, computing Prob(Back(ovnode))
				NodeBern* ovnode1 = static_cast<NodeBern*>(ovnode);
				ovnode1->BackProb = Back;
				Back = 1;
			}else{ //otherwise, fill PosNNodeBacks and WNodeBacks
				int s = WNodeBacks_size + backlen;
				PosNNodeBacks[node->num] = WNodeBacks_size;
				int pos = 0;
				int k;
				for( k = WNodeBacks_size; k< s; k++){
					WNodeBacks[k] = BackW.at(pos);
					pos++;
				}
				WNodeBacks_size = s;
			}

			//RECURSION, MAINTAINING OF PREFIX SUCCESSORS OF T 
		    std::list<NodeAC*> LChilds = node->LChilds;
			delete node; //delete processed node T of the trie

			for(i = LChilds.begin(); i != LChilds.end(); i++){
				NodeAC* LCnode = *i;
				CreateLOG(LCnode,ovnode,prob,Back,len,BackW, Word);
			}

		} /* if((node->main != 0) && (node != NodeAC::ACRoot))*/
		else{ 
			//processing of internal nodes that are not overlaps, maintaning of prefix successors of T

			std::list<NodeAC*> LChilds = node->LChilds;
			if(node != NodeAC::ACRoot) delete node; //delete processed node T of the trie

			for(i = LChilds.begin(); i != LChilds.end(); i++){
				NodeAC* LCnode = *i;
				CreateLOG(LCnode,LParent,prob,Back,len,BackW, Word);
			}
		} 
	}/* if(T->leaf == 0) */
	else{ //T is leaf, fill lists LLeafPreds, LeafBacks, LeafProbs, PosNLeafWords, WLeafWords
		LeafAC* leaf = static_cast<LeafAC*>(T);
		LLeafPreds[leaf->leafnum-1] = LParent->num;
		if(MainData::order == 0){
			LeafBacks[leaf->leafnum-1] = Back;
			LeafProbs[leaf->leafnum-1] = prob; 
		}
		else{
			PosNLeafWords[leaf->leafnum - 1] = WLeafWords_size;
			int k; 
			for(k = 0; k< MainData::WordLen; k++)
				WLeafWords[WLeafWords_size + k] = Word[k];

			WLeafWords_size += MainData::WordLen;
		}
		delete T; 
	}
	if((T!= NodeAC::ACRoot)&&(MainData::order != 0)){
		BackW[backlen - 1] = symb;
	}
}








//ofstream debff("nwwcdee.txt");
/////////////////////////////////////////////
void DebPrint1(void){
	int i,j;
	for(i = 0; i < NodeOv::NumOVNodes; i++){
		NodeOv* node = Nodes[i];
		//cout<<"/////////////Next node///////////" <<'\n';
		cout<<"Number: "<<i<<'\t'<<"Len: "<<node->len<<'\n';
		cout<<"lpred: "<<LPreds[i]<<'\t';
		cout<<"rpred: "<<RPreds[i]<<'\t';

		NodeBern* node1 = static_cast<NodeBern*>(node);
		cout<<"Back: "<<node1->BackProb<<"\n";
		if(node1->rdeep == 1) cout<<"WordProb:  "<<node1->WordProb<<'\n';
			
		if(node1->ldeep == 1){ 
			int s = node1->NDLinks;
			cout<<"(DLink,DProb): ";
			for(j = 0; j < s; j++){
				cout<<'('<<node1->DeepLinks[j]->num<<','<<node1->DeepProbs[j]<<')'<<'\t';
			}
			cout<<"\n\n";
		}

		cout<<"LChilds: ";
		for(j = 0; j<node->NLChilds; j++){
			cout<<node->LChilds[j]->num <<'\t';
		}
		cout<<'\n';

		cout<<"RChilds: ";
		for(j = 0; j<node->NRChilds; j++){
			cout<<node->RChilds[j]->num <<'\t';
		}
		cout<<'\n';
		for(j = 0; j < NumLLeafChilds[i]; j++){
			cout<<LLeafChilds[i][j]+1 <<'\t';
		}
		cout<<"\n\n\n";
	}
	cout<<"\n\n\n\n\n";
	/*
	cout<<"========CLASSES=========\n\n";
	for(i = 0; i < NodeOv::NClasses; i++){
		cout<<"CLass "<<i<<'\n';
		for(j = 0; j < CLasses[i].size(); j++){
			cout<<CLasses[i][j]<<'\t';
		}
		cout<<"\n\n";
	}

	cout<<"\n\n\n\n\n";
*/
	cout<<"========LEAVES=========\n\n";
	
	for(i = 0; i < MainData::NWords; i++){
		cout<<"Number: "<<i+1<<'\n';
		cout<<"lpred: "<<LLeafPreds[i]<<'\t'<<"rpred: "<<RLeafPreds[i]<<'\n';
		cout<<"Back: "<<LeafBacks[i]<<'\n';
		cout<<"Prob: "<<LeafProbs[i]<<"\n\n\n";
	}
	//debff.close();

}


///////////////////////////////
//For each node w of OvGraf, creates links to left and right successors of w, fill lists LChilds and RChilds;   
//Fill lists  LLeafChilds and computing other parameters
int CrChilds(void){
	int i, lnum, rnum;

	//1. FOR EACH NODE w of OvGraf, CREATES LINKS TO LEFT AND RIGHT SUCCESSORS of w
	//FILL LISTS w.LChilds and w.RChilds
	int* CurNumLChilds = Malloc<int>(NodeOv::NumOVNodes);
	int* CurNumRChilds = Malloc<int>(NodeOv::NumOVNodes);

	// 1.1. initialization of lists w.LChilds and w.RChilds
	for(i = 0; i < NodeOv::NumOVNodes; i++){
		lnum = Nodes[i]->NLChilds;
		rnum = Nodes[i]->NRChilds;
		Nodes[i]->LChilds = Malloc<NodeOv*>(lnum); 
		Nodes[i]->RChilds = Malloc<NodeOv*>(rnum); 
	}
	
	//1.2. fill lists w.LChilds and w.RChilds
	int lpos,rpos;
	for(i = 1; i < NodeOv::NumOVNodes; i++){
		lnum = LPreds[i];
		lpos = CurNumLChilds[lnum];

		Nodes[lnum]->LChilds[lpos] = Nodes[i];

		CurNumLChilds[lnum]++;

        if(lnum == 0) Nodes[i]->rootchild = true;

		rnum = Nodes[i]->RParent->num;
		rpos = CurNumRChilds[rnum];

		Nodes[rnum]->RChilds[rpos] = Nodes[i];

		CurNumRChilds[rnum]++;

	}
	
	free(CurNumRChilds);
	


	//2. FILL LIST LLeafChilds and computing other parameters
	NumLLeafChilds =  Malloc<int>(NodeOv::NumOVNodes);

	int lid, rid;
	for(i = 0; i < MainData::NWords; i++){ //for each w, computing of number of words h in HH, where w= lpred(h), i.e. size of LLeafChilds(w) 
		lid = LLeafPreds[i];
		rid = RLeafPreds[i];
		NumLLeafChilds[lid]++;
		if(MainData::order == 0){ //For bernoulli model, computing of word probabilities w.WordProb for all nodes, see nodebern.h
			NodeBern* RPnode= static_cast<NodeBern*>(Nodes[rid]);
			RPnode->WordProb +=LeafProbs[i];
		}

		//computing of numbers of left and right deep nodes 
		if(Nodes[lid]->ldeep == 0){
			Nodes[lid]->ldeep = 1;
			NodeOv::NumLDNodes ++;
		}
		if(Nodes[rid]->rdeep == 0){
		   Nodes[rid]->rdeep = 1;
		   NodeOv::NumRDNodes ++;
		}
	}

	//fill list LLeafChilds 
    LLeafChilds =new int* [NodeOv::NumOVNodes]();

	for(i = 0; i < NodeOv::NumOVNodes; i++){
	    int lnum = NumLLeafChilds[i]; 
		LLeafChilds[i] = new int[lnum]();
	}
	
	
	memset(CurNumLChilds, 0x00, NodeOv::NumOVNodes* sizeof(int));
	
	for(i = 0; i < MainData::NWords; i++){
		lid = LLeafPreds[i];
		lpos =  CurNumLChilds[lid];
        LLeafChilds[lid][lpos] = i;
		CurNumLChilds[lid]++;
	}
	
	free(CurNumLChilds);

	return 0;
}
///////////////////////////////

//minimization of overlap graph by depth-first traversal using left links

//let w be current processed left deep node and h we word from motif such that w = lpred(h) 
//RECALL! during processing of current word h from HH, let w = lpred(h) and r = rpred(h)
//NLink[rpred(h)] is number of deep link corresponding to h*(w,r) in w.DeepLinks  
//NLOG[rpred(h)] = w
//NCLass[r] = number of class corresponding to h

//Input parameters: node - current processing node w


void MinGraf(NodeOv* node){
    int i;

	if(node->ldeep == 1){ // w is left deep node

		//1. initialazing some auxilary lists 
		vector<NodeOv*> deeplinks;
		vector<double> deepprobs;
		int size= (int) NumLLeafChilds[node->num];
		
		deeplinks.reserve(size);
		NodeBern* bnode;
		if(MainData::order == 0){
		    bnode = static_cast<NodeBern*>(node);
			deepprobs.reserve(size);
		}
		else{
			NLinkToClass[node->num].reserve(size);
			
		}

		//2. Processing of motif words h such that w = lpred(h)
		for(i = 0; i< size; i++){
			int leaf = LLeafChilds[node->num][i]; //number of h
			int rid = RLeafPreds[leaf];
			NodeOv* RP = Nodes[rid];              //rpred(h)
			if(NLOG[RP->num] != node->num){		  //NLOG[rpred(h)] != w, there is no class for h.  
												  //Creating of new class h*(w,rpred(h))
				NodeOv::NClasses++;				  
				deeplinks.push_back(RP);		  //creating new deep link corresponding to  h*(w,rpred(h))
				int nlink = deeplinks.size()-1;
				NLink[RP->num] = nlink;			  //NLink[rpred(h)] = number of the link is list w.DeepLinks
				NLOG[RP->num] = node->num;		  //NLOG[rpred(h)] = w
				
				//fill other parameters
				if(MainData::order == 0){
					deepprobs.push_back(LeafBacks[leaf]);
				}
				else{	
					CLassNum[leaf] = NodeOv::NClasses;
					NCLass[RP->num] = NodeOv::NClasses;
					NLinkToClass[node->num].push_back(NodeOv::NClasses);
					vector<int> clas;
					clas.push_back(leaf);
					CLasses.push_back(clas);
				}
			
			}
			else{						        //NLOG[rpred(h)] == w, there exists class for h. 
				int nlink = NLink[RP->num];

				//fill some parameters
				if(MainData::order == 0){
					deepprobs[nlink] += LeafBacks[leaf];
				}
				else{
					int nclass = NCLass[RP->num];
					CLassNum[leaf] = nclass;
					CLasses[nclass -1].push_back(leaf);
				}
			}			
		}/*for(i = 0; i< size; i++)*/

		//3.Creating lists w.DeepLinks and w.DeepProbs (for Bernoulli)
		size  = (int)deeplinks.size();
		node->NDLinks = size;
		node->DeepLinks =Malloc<NodeOv*>(size);
		if(MainData::order == 0){
			bnode->DeepProbs = Malloc<double>(size);
		}

		for(i = 0; i < size; i++){
			node->DeepLinks[i]= deeplinks[i];
			if(MainData::order == 0){
				bnode->DeepProbs[i] = deepprobs[i];
			}
		}
	}

	for (  i = 0; i < node->NLChilds; i++) {
			MinGraf(node->LChilds[i]);
	}
	return;
}


//////////////////////////////////////////////////

//computing of maximal number of nodes MaxDepth on the path leading from root to a terminal node
void ComputeMaxDepth(NodeOv* node, int depth){
	if(MaxDepth < depth)
		MaxDepth = depth;
	depth++;
	int i;
	for( i = 0; i < node->NLChilds; i++) {
		ComputeMaxDepth(node->LChilds[i],depth);
	}
}
/////////////////////////


// It creates overlap graph
int NodeOv::CreateGraf(void){
// Given: AC-trie, see nodeac.h
// Aim:   Overlap Graph, see nodeov.h
	//		1. INITIALIZATION.
	int i;
	NodeOv::Root = NodeOv::NewNode(); 

	Nodes =Malloc<NodeOv*>(NodeOv::NumOVNodes);		// Node descriptors
	Nodes[0] = NodeOv::Root;						// Root

	LPreds =  Malloc<int>(NodeOv::NumOVNodes);		//left parents of nodes
	RPreds =  Malloc<int>(NodeOv::NumOVNodes);		//left parents of nodes

	LLeafPreds = Malloc<int>(MainData::NWords);		//left parents of leaves
	RLeafPreds = Malloc<int>(MainData::NWords);		//right parents of leaves

	if(MainData::order == 0){
		LeafProbs = Malloc<double>(MainData::NWords);       // Probabilities of words
	    LeafBacks = Malloc<double>(MainData::NWords);	    // Back for leaves
	}
	else{
		PosNNodeBacks = Malloc<int>(NodeOv::NumOVNodes);

		PosNLeafWords =Malloc<int>(MainData::NWords);
		NLinkToClass = new vector<int>[NodeOv::NumOVNodes];
	    WLeafWords.resize(MainData::WordLen*MainData::NWords);

	}

	ComputeRPred(NodeAC::ACRoot);

	InternAC* acroot = static_cast<InternAC*>(NodeAC::ACRoot);
	Compute_WNodeBacks_len(acroot, 0, 0);
	
	string word;
	word.resize(MainData::WordLen);
	string BackW;
	BackW.resize(MainData::WordLen);
	WNodeBacks.resize(WNodeBacks_len);

		//	2. CREATING OF OVERLAP GRAPH
    CreateLOG(NodeAC::ACRoot, NodeOv::Root, 1,1,0,word,BackW);	//create overlaps, compute some data in descriptors of nodes, delete AC trie
	CrChilds();			// Create arrays with left and right childs of overlap graph


	//free and initialize some auxilary data
	free(LPreds);
	free(RPreds);
	free(LLeafPreds);
	
	if(MainData::order == 0) free(LeafProbs);

	CLassNum = Malloc<int>(MainData::NWords);
	NCLass = Malloc<int>(NodeOv::NumOVNodes);
    NLOG = Malloc<int>(NodeOv::NumOVNodes);
    NLink = Malloc<int>(NodeOv::NumOVNodes);
	
	for(i = 0; i < NodeOv::NumOVNodes; i++){
		NLOG[i] = -1;
		NLink[i] = -1;
	}
    ////

    MinGraf(NodeOv::Root);	//overlap graph minimization


	ComputeMaxDepth(NodeOv::Root, 0); //computing of  MaxDepth

		//3. FREE SOME DATA STRUCTURES

	if(MainData::order == 0){
		free(Nodes);
		 free(LeafBacks);
		 free(RLeafPreds);
	}

	free(NumLLeafChilds);
	for(i = 0; i < NodeOv::NumOVNodes; i++){
		delete[] LLeafChilds[i];
	}
	delete[] LLeafChilds;

	free(NLOG);
	free(NLink);
	
	
	return 0; 
}




//////////////////////////////////////////////////

//delete overlap graph
void NodeOv::DeleteGraf(NodeOv* node){
	int i;
	
	for( i = 0; i < node->NLChilds; i++){
		NodeOv::DeleteGraf(node->LChilds[i]);
	}
	delete node;
	return;
};
